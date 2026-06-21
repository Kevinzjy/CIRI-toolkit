//! Scan2: second-pass rescue, candidate validation, and FSJ counting.
//!
//! This stage consumes the Java-compatible BSJ output from Scan1, rebuilds
//! the de-duplicated candidate indexes expected by Java CIRI3, and then revisits
//! the input alignments to rescue additional support while counting FSJ evidence.

use crate::circ_catalog::CircCatalogRecord;
use crate::is_bsj_hg2::{report_scan2_hg_profile, IsBSJHg2};
use crate::misd::misd;
use crate::runtime::{
    emit_perf_line, emit_trace_line, scan2_profile_enabled, should_trace_read, with_trace_hg2_scope,
};
use crate::utils::{
    alignment_short_cs, bam_shard_count, bsj_is_summary_priority, bsj_payload_start,
    cigar_is_full_match, clip_sequence_payload, local_clip_evidence_lines, parse_cigar_ops_basic,
    part_path, reverse_complement, AlignmentRecord,
};
use anyhow::{bail, Result};
use indicatif::{ProgressBar, ProgressStyle};
use memmap2::Mmap;
use noodles::sam::{
    self,
    alignment::{
        record::{
            data::field::{Tag, Value},
            Data as _, Sequence as _,
        },
        Record as _,
    },
};
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::env;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Seek, Write};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc;
use std::sync::OnceLock;
use std::thread;
use std::time::Duration;
use std::time::Instant;

const NON_BSJ_SEGMENT_MAPQ_THRES: i32 = 5;
const NON_BSJ_SEGMENT_MAX_SPAN: i32 = 200000;
const NON_BSJ_OUTWARD_MIN_TERMINAL_CLIP: i32 = 19;
const NON_BSJ_OUTWARD_MIN_PAIR_OFFSET: i32 = 19;
const NON_BSJ_OUTWARD_RETENTION_ENV: &str = "CIRI_SCAN2_NON_BSJ_OUTWARD_RETENTION";
const STRONG_FSJ_ENV: &str = "CIRI_SCAN2_STRONG_FSJ";
const STRONG_FSJ_MIN_ANCHOR: i32 = 19;
const DIRECT_SIDE_FSJ_ANCHOR_ENV: &str = "CIRI_SCAN2_DIRECT_SIDE_FSJ_ANCHOR";
const DIRECT_SIDE_FSJ_RESCUE_ANCHOR_ENV: &str = "CIRI_SCAN2_DIRECT_SIDE_FSJ_RESCUE_ANCHOR";
const DIRECT_SIDE_FSJ_SPLICE_TOLERANCE: i32 = 6;
const FSJ_PAIR_DUMP_ENV: &str = "CIRI_SCAN2_FSJ_PAIR_DUMP";

/// Outward-retention modes for the non-BSJ sidecar prefilter.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum NonBsjOutwardRetention {
    /// Keep the legacy primary-pair-only spill policy.
    PrimaryOnly,
    /// Also allow high-MAPQ supplementary pairs to trigger sidecar retention.
    Supplementary,
    /// Also allow one low-MAPQ mate when the opposite mate has usable MAPQ.
    Relaxed,
    /// Allow high-MAPQ primary/supplementary pairs for stricter precision testing.
    Filtered,
}

/// Core logic for the second scan pass.
pub struct Scan2 {
    /// Minimum mapping quality for a read.
    pub min_mapq_uni: i32,
    /// Minimum linear range size for competition check.
    pub linear_range_size_min: i32,
    /// Index for BSJ candidates (indexed by Site1).
    ///
    /// Java parity note:
    /// - Backed by de-duplicated `siteInfor` records (site1..sum_q), not raw BSJ1 lines.
    /// - Candidate payload layout is `[site1, site2, strand, signal1, signal2, sum_q]`.
    pub index1: HashMap<String, Vec<CandidateBreakpoint>>,
    /// Index for BSJ candidates (indexed by Site2).
    pub index2: HashMap<String, Vec<CandidateBreakpoint>>,
    /// Java parity gate: existing buckets for Site1 index.
    pub site_array1: HashMap<String, HashSet<i32>>,
    /// Java parity gate: existing buckets for Site2 index.
    pub site_array2: HashMap<String, HashSet<i32>>,
    /// Map to store Forward Spliced Junction (FSJ) counts.
    pub fsj_map: HashMap<String, i32>,
    /// Memory limit per thread for Page Cache eviction.
    pub mem_limit: u64,
    /// Read IDs already assigned in Scan1 (Java parity: skip in Scan2).
    pub scan1_ids: HashSet<String>,
    /// Java Scan2 bucket size (`seqLen`) used for directional candidate traversal.
    pub seq_len: i32,
    /// Unique read IDs rescued by Scan2 and appended to the final BSJ file.
    pub rescued_reads: usize,
    /// Total unique BSJ-supporting read IDs in the final merged BSJ file.
    pub final_bsj_reads: usize,
}

/// Temporary sidecar paths produced while Scan2 feeds the segments stage.
///
/// These files are not part of the CIRI3-compatible BSJ output contract. They
/// let the CLI hand large read-level extension evidence directly to the
/// post-Summary segments stage without first merging it into a second huge text
/// file that would immediately be mmap-parsed again.
#[derive(Default)]
pub struct Scan2SegmentArtifacts {
    /// Shard-local non-BSJ topology sidecars retained for direct segments input.
    pub non_bsj_segment_evidence_paths: Vec<String>,
}

/// Returns a stable text payload for one optional `XA:Z` field.
///
/// Scan2 writes `*` for missing alternatives, but otherwise preserves raw BWA
/// XA syntax so the post-Summary segments phase can parse alternatives without
/// a second pass over BAM/SAM.
fn xa_payload(raw: &str) -> &str {
    if raw.is_empty() {
        "*"
    } else {
        raw
    }
}

/// Returns the 1-based mate bucket encoded by a SAM flag.
fn mate_bucket_from_flag(flag: i32) -> i32 {
    if flag & 0x40 != 0 {
        1
    } else if flag & 0x80 != 0 {
        2
    } else {
        0
    }
}

/// Returns whether Scan2 should require strong side anchors for FSJ counting.
///
/// The switch is cached because the check sits on the Scan2 hot path. Leaving
/// it unset preserves Java-compatible span-based FSJ counts.
fn strong_fsj_filter_enabled() -> bool {
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| {
        env::var(STRONG_FSJ_ENV)
            .map(|value| matches!(value.as_str(), "1" | "true" | "TRUE" | "yes" | "on"))
            .unwrap_or(false)
    })
}

/// Returns the optional direct side-spanning FSJ anchor threshold.
///
/// This is an experimental minibwa evaluation mode. When set, Scan2 replaces
/// Java's broad span-based FSJ counting with direct side support from CIGAR
/// `M` blocks and splice-like `N` gaps. Leaving the variable unset keeps the
/// Java-compatible path byte-for-byte reachable for parity checks.
fn direct_side_fsj_anchor() -> Option<i32> {
    static ANCHOR: OnceLock<Option<i32>> = OnceLock::new();
    *ANCHOR.get_or_init(|| {
        env::var(DIRECT_SIDE_FSJ_ANCHOR_ENV)
            .ok()
            .and_then(|value| value.parse::<i32>().ok())
            .filter(|value| *value > 0)
    })
}

/// Returns the optional direct side-spanning FSJ rescue anchor threshold.
///
/// Unlike `CIRI_SCAN2_DIRECT_SIDE_FSJ_ANCHOR`, this mode keeps the
/// Java-compatible CIRI3 FSJ set and only adds extra side-spanning linear
/// evidence. It is the preferred experiment for minibwa because it improves
/// recall without discarding CIRI3's already high-precision default FSJs.
fn direct_side_fsj_rescue_anchor() -> Option<i32> {
    static ANCHOR: OnceLock<Option<i32>> = OnceLock::new();
    *ANCHOR.get_or_init(|| {
        env::var(DIRECT_SIDE_FSJ_RESCUE_ANCHOR_ENV)
            .ok()
            .and_then(|value| value.parse::<i32>().ok())
            .filter(|value| *value > 0)
    })
}

/// Returns the optional read-level FSJ audit path.
///
/// The main `.out` format only contains circ-level FSJ counts. This debug-only
/// sink records the exact read groups that incremented those counts, allowing
/// simulator precision/recall audits without changing the public output
/// contract.
fn fsj_pair_dump_path() -> Option<String> {
    static PATH: OnceLock<Option<String>> = OnceLock::new();
    PATH.get_or_init(|| {
        env::var(FSJ_PAIR_DUMP_ENV)
            .ok()
            .filter(|value| !value.trim().is_empty())
    })
    .clone()
}

/// Returns side-specific FSJ keys with strong linear support.
///
/// The CIRI3-compatible FSJ counter treats a broad alignment span as enough
/// evidence. This enhanced predicate is stricter: each candidate side is counted
/// separately only when it is in the middle of a read chain with enough aligned
/// sequence on both sides. Continuous match blocks and CIGAR `N` splice gaps can
/// provide that evidence; large `D` artifacts cannot.
fn strong_fsj_side_keys(key: &str, aln_pos: i32, cigar: &str) -> Vec<String> {
    let mut fields = key.split('\t');
    let Some(chr) = fields.next() else {
        return Vec::new();
    };
    let Some(site1) = fields.next().and_then(|value| value.parse::<i32>().ok()) else {
        return Vec::new();
    };
    let Some(site2) = fields.next().and_then(|value| value.parse::<i32>().ok()) else {
        return Vec::new();
    };
    let mut out = Vec::with_capacity(2);
    if cigar_has_strong_side_spanning(cigar, aln_pos, site1) {
        out.push(format!("{chr}\t{site1}\t{site2}\tS"));
    }
    if cigar_has_strong_side_spanning(cigar, aln_pos, site2) {
        out.push(format!("{chr}\t{site1}\t{site2}\tE"));
    }
    out
}

/// Tests one candidate side against the read-chain encoded by one CIGAR string.
///
/// A side can be supported either inside one match block or across an `N` splice
/// gap whose flanking match blocks both provide enough query sequence. Deletions
/// do not count as splice support because minibwa can use large `D` operations
/// to absorb junction gaps that BWA would represent with split alignments.
fn cigar_has_strong_side_spanning(cigar: &str, aln_pos: i32, side: i32) -> bool {
    if cigar == "*" || aln_pos <= 0 {
        return false;
    }
    let mut ref_pos = aln_pos;
    let mut prev_match: Option<(i32, i32)> = None;
    let mut pending_skip_after: Option<(i32, i32)> = None;
    let mut count = 0_i32;
    let mut has_count = false;
    for op in cigar.chars() {
        if op.is_ascii_digit() {
            count = match count
                .checked_mul(10)
                .and_then(|value| value.checked_add(op.to_digit(10).unwrap_or(0) as i32))
            {
                Some(value) => value,
                None => return false,
            };
            has_count = true;
            continue;
        }
        if !has_count {
            return false;
        }
        match op {
            'M' | '=' | 'X' => {
                let block_start = ref_pos;
                let block_end = ref_pos + count - 1;
                if side >= block_start
                    && side <= block_end
                    && side - block_start + 1 >= STRONG_FSJ_MIN_ANCHOR
                    && block_end - side + 1 >= STRONG_FSJ_MIN_ANCHOR
                {
                    return true;
                }
                if let Some((skip_start, skip_end)) = pending_skip_after.take() {
                    if let Some((prev_start, prev_end)) = prev_match {
                        let left_anchor = prev_end - prev_start + 1;
                        let right_anchor = block_end - block_start + 1;
                        if side >= skip_start - 6
                            && side <= skip_end + 6
                            && (side - prev_end).abs().min((block_start - side).abs()) <= 6
                            && left_anchor >= STRONG_FSJ_MIN_ANCHOR
                            && right_anchor >= STRONG_FSJ_MIN_ANCHOR
                        {
                            return true;
                        }
                    }
                }
                prev_match = Some((block_start, block_end));
                ref_pos += count;
            }
            'N' => {
                pending_skip_after = prev_match.map(|(_start, end)| (end + 1, ref_pos + count - 1));
                ref_pos += count;
            }
            'D' => {
                pending_skip_after = None;
                ref_pos += count;
            }
            'I' | 'S' | 'H' | 'P' => {
                pending_skip_after = None;
            }
            _ => return false,
        }
        count = 0;
        has_count = false;
    }
    false
}

/// Collapses optional side-specific strong-FSJ keys back to circ-level keys.
fn fsj_count_key(key: &str) -> String {
    let mut fields = key.split('\t');
    match (fields.next(), fields.next(), fields.next(), fields.next()) {
        (Some(chr), Some(start), Some(end), Some(_side)) => format!("{chr}\t{start}\t{end}"),
        _ => key.to_string(),
    }
}

/// Collapses side-specific FSJ audit keys to one read-pair-level count key.
///
/// Strong-FSJ mode can record separate start- and end-side evidence for one
/// read. Summary FSJ counts are still read-pair-level, so side keys must be
/// deduplicated after collapsing back to the circRNA key.
fn collapsed_fsj_count_keys(keys: HashSet<String>) -> HashSet<String> {
    keys.into_iter().map(|key| fsj_count_key(&key)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(flag: i32, pos: i32, mapq: i32, cigar: &'static str) -> AlignmentRecord<'static> {
        AlignmentRecord {
            flag,
            chrom: Cow::Borrowed("chr1"),
            pos,
            mapq,
            cigar: Cow::Borrowed(cigar),
            seq: Cow::Borrowed("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
            xa: Cow::Borrowed(""),
        }
    }

    #[test]
    fn non_bsj_outward_retention_defaults_to_recognition_first_relaxed() {
        assert_eq!(
            Scan2::parse_non_bsj_outward_retention_mode(None),
            NonBsjOutwardRetention::Relaxed
        );
        assert_eq!(
            Scan2::parse_non_bsj_outward_retention_mode(Some("primary")),
            NonBsjOutwardRetention::PrimaryOnly
        );
        assert_eq!(
            Scan2::parse_non_bsj_outward_retention_mode(Some(" primary-only ")),
            NonBsjOutwardRetention::PrimaryOnly
        );
        assert_eq!(
            Scan2::parse_non_bsj_outward_retention_mode(Some("bogus")),
            NonBsjOutwardRetention::PrimaryOnly
        );
        assert_eq!(
            Scan2::parse_non_bsj_outward_retention_mode(Some("filtered")),
            NonBsjOutwardRetention::Filtered
        );
    }

    #[test]
    fn non_bsj_outward_prefilter_accepts_gap_facing_primary_pair() {
        let records = vec![
            aln(97, 94_955_296, 60, "80M70S"),
            aln(145, 94_953_252, 60, "97M53S"),
        ];

        assert!(Scan2::may_support_outward_segments(&records));
    }

    #[test]
    fn non_bsj_outward_prefilter_rejects_linear_fr_pair() {
        let records = vec![aln(65, 100, 60, "150M"), aln(145, 300, 60, "150M")];

        assert!(!Scan2::may_support_outward_segments(&records));
    }

    #[test]
    fn non_bsj_outward_prefilter_keeps_mapq_gate() {
        let records = vec![aln(97, 1_000, 60, "100M50S"), aln(145, 500, 4, "100M50S")];

        assert!(!Scan2::may_support_outward_segments(&records));
    }

    #[test]
    fn non_bsj_outward_prefilter_accepts_same_span_terminal_clips() {
        let records = vec![aln(65, 100, 60, "131M19S"), aln(145, 100, 60, "131M19S")];

        assert!(Scan2::may_support_outward_segments(&records));
    }

    #[test]
    fn strong_fsj_side_anchor_accepts_match_block_with_two_anchors() {
        assert_eq!(
            strong_fsj_side_keys("chr1\t150\t300", 100, "101M"),
            vec!["chr1\t150\t300\tS".to_string()]
        );
    }

    #[test]
    fn strong_fsj_side_anchor_accepts_spliced_read_chain() {
        assert_eq!(
            strong_fsj_side_keys("chr1\t141\t300", 100, "40M100N40M"),
            vec!["chr1\t141\t300\tS".to_string()]
        );
    }

    #[test]
    fn strong_fsj_side_anchor_rejects_large_deletion_span() {
        assert!(strong_fsj_side_keys("chr1\t141\t200", 100, "40M117D80M").is_empty());
    }

    #[test]
    fn strong_fsj_side_anchor_rejects_short_flanking_anchor() {
        assert!(strong_fsj_side_keys("chr1\t150\t300", 140, "60M").is_empty());
    }

    #[test]
    fn strong_fsj_side_anchor_counts_start_and_end_separately() {
        assert_eq!(
            strong_fsj_side_keys("chr1\t150\t180", 100, "101M"),
            vec![
                "chr1\t150\t180\tS".to_string(),
                "chr1\t150\t180\tE".to_string()
            ]
        );
    }

    #[test]
    fn collapsed_fsj_count_keys_deduplicates_side_specific_keys() {
        let mut keys = HashSet::new();
        keys.insert("chr1\t150\t180\tS".to_string());
        keys.insert("chr1\t150\t180\tE".to_string());

        let collapsed = collapsed_fsj_count_keys(keys);

        assert_eq!(collapsed.len(), 1);
        assert!(collapsed.contains("chr1\t150\t180"));
    }

    #[test]
    fn non_bsj_outward_supplementary_retention_keeps_high_mapq_split_pair() {
        let records = vec![
            aln(97, 32_740_922, 1, "31S99M149D20M"),
            aln(2145, 32_740_654, 34, "32M118S"),
            aln(145, 32_716_883, 60, "91M59S"),
        ];

        assert!(!Scan2::may_support_outward_segments(&records));
        assert!(Scan2::may_support_outward_segments_with_retention(
            &records,
            NonBsjOutwardRetention::Supplementary
        ));
    }

    #[test]
    fn non_bsj_outward_relaxed_retention_keeps_one_low_mapq_mate() {
        let records = vec![
            aln(97, 29_650_146, 60, "150M"),
            aln(145, 29_649_893, 1, "10S115M134D25M"),
        ];

        assert!(!Scan2::may_support_outward_segments(&records));
        assert!(!Scan2::may_support_outward_segments_with_retention(
            &records,
            NonBsjOutwardRetention::Supplementary
        ));
        assert!(Scan2::may_support_outward_segments_with_retention(
            &records,
            NonBsjOutwardRetention::Relaxed
        ));
    }
}

/// Extracts a raw BWA `XA:Z` payload from one SAM record.
///
/// XA never feeds Scan2 rescue or FSJ counting. It is carried only as sidecar
/// evidence for later circ-context segment chain selection.
fn raw_xa_from_sam_record(record: &sam::Record) -> Result<String> {
    let tag = Tag::new(b'X', b'A');
    let data = record.data();
    let Some(value) = data.get(&tag).transpose()? else {
        return Ok(String::new());
    };
    let Value::String(raw) = value else {
        return Ok(String::new());
    };
    let Ok(raw) = std::str::from_utf8(raw.as_ref()) else {
        return Ok(String::new());
    };
    Ok(raw.to_string())
}

/// Extracts a raw BWA `XA:Z` payload from one BAM record.
///
/// Malformed optional tags are treated as absent because XA is a post-Summary
/// representation hint, not required evidence for CIRI3 parity decisions.
fn raw_xa_from_bam_record(record: &noodles::bam::Record) -> Result<String> {
    let tag = Tag::new(b'X', b'A');
    let data = record.data();
    let Some(value) = data.get(&tag).transpose()? else {
        return Ok(String::new());
    };
    let Value::String(raw) = value else {
        return Ok(String::new());
    };
    let Ok(raw) = std::str::from_utf8(raw.as_ref()) else {
        return Ok(String::new());
    };
    Ok(raw.to_string())
}

/// Returns a fixed-width reference window, padding out-of-range bases with `N`.
///
/// External BED catalogs can legally point close to chromosome ends. Scan2 needs
/// two signal bases in the same payload slots as Scan1-derived candidates, so
/// padding keeps the payload well-formed while the downstream sequence validator
/// still decides whether any read truly supports the boundary.
fn reference_window_or_n(seq: &str, start: i32, end: i32) -> String {
    let bytes = seq.as_bytes();
    let mut out = String::with_capacity((end - start).max(0) as usize);
    for pos in start..end {
        if pos >= 0 {
            let idx = pos as usize;
            if idx < bytes.len() {
                out.push(bytes[idx] as char);
                continue;
            }
        }
        out.push('N');
    }
    out
}

#[derive(Clone)]
pub struct CandidateBreakpoint {
    /// Genomic coordinate used for bucket/range lookup.
    pub site: i32,
    /// Insertion order after Java-style de-duplication.
    /// Used as a stable tie-breaker when two candidates share the same `site`.
    pub order: usize,
    /// Candidate payload in Java `siteInfor` layout:
    /// `[site1, site2, strand, signal1, signal2, sum_q]`.
    pub data: Vec<String>,
}

type OwnedAlignmentRecord = AlignmentRecord<'static>;
type OwnedStandMap = HashMap<i32, (char, Cow<'static, str>)>;

struct SamOwnedScan2Group {
    read_id: String,
    alignments: Vec<OwnedAlignmentRecord>,
    all_alignments: Vec<OwnedAlignmentRecord>,
    stand_map: OwnedStandMap,
}

#[derive(Clone, Copy)]
struct NonBsjMsid {
    kind: i32,
    clip1: i32,
    clip2: i32,
    ref_len: i32,
}

impl NonBsjMsid {
    fn invalid(ref_len: i32) -> Self {
        Self {
            kind: 0,
            clip1: 0,
            clip2: 0,
            ref_len,
        }
    }
}

/// Advises the kernel that a processed `mmap` slice can be evicted from cache.
///
/// The offset is aligned manually because `madvise` requires page alignment on
/// Linux. Without that alignment the call can fail silently from the pipeline's
/// point of view and the full BAM/SAM page cache footprint keeps growing.
fn advise_dontneed(mmap: &Mmap, offset: usize, len: usize) {
    if len == 0 {
        return;
    }
    let page_size = 4096usize;
    let aligned_offset = (offset / page_size) * page_size;
    let aligned_len = ((offset + len + page_size - 1) / page_size) * page_size - aligned_offset;
    unsafe {
        // SAFETY: `aligned_offset` and `aligned_len` are page-aligned within the
        // live mmap, and `madvise` only receives a non-mutating cache hint.
        let ptr = mmap.as_ptr().add(aligned_offset);
        libc::madvise(ptr as *mut libc::c_void, aligned_len, libc::MADV_DONTNEED);
    }
}

/// Optional Scan2 profiler used when CLI `--perf` or legacy
/// `CIRI_PROFILE_SCAN2=1` is enabled.
///
/// The counters stay intentionally coarse so the profiler can be left inside the
/// candidate hot path without materially perturbing release timings.
#[derive(Default)]
pub(crate) struct Scan2Profile {
    wall_total_ns: AtomicU64,
    shard_total_ns: AtomicU64,
    group_process_ns: AtomicU64,
    validator_call_ns: AtomicU64,
    display_ns: AtomicU64,
    segment_sidecar_ns: AtomicU64,
    non_bsj_sidecar_ns: AtomicU64,
    write_ns: AtomicU64,
    merge_ns: AtomicU64,
    records: AtomicU64,
    groups: AtomicU64,
    candidate_checks: AtomicU64,
    candidate_hits: AtomicU64,
    segment_sidecar_rows: AtomicU64,
    segment_sidecar_bytes: AtomicU64,
    non_bsj_sidecar_rows: AtomicU64,
    non_bsj_sidecar_bytes: AtomicU64,
}

/// Mate-level Scan1 claims used to gate Scan2 display-only rescue rows.
///
/// The display path is intentionally outside Summary counting, but it runs for
/// every Scan2 read group. Storing claims as a compact per-read bitmask avoids
/// allocating `read_id\tmate` strings in that hot path and lets fully claimed
/// groups skip display reconstruction before any per-mate alignment indexing.
struct Scan2DisplayClaims {
    by_read: HashMap<String, u8>,
}

impl Scan2DisplayClaims {
    /// Creates an empty display claim index.
    fn new() -> Self {
        Self {
            by_read: HashMap::new(),
        }
    }

    /// Records that Scan1 already owns one displayed mate for `read_id`.
    fn insert(&mut self, read_id: &str, mate_label: &str) {
        let bit = match mate_label {
            "R1" => 0b01,
            "R2" => 0b10,
            _ => 0,
        };
        if bit != 0 {
            *self.by_read.entry(read_id.to_string()).or_insert(0) |= bit;
        }
    }

    /// Tests whether the display path should skip this mate.
    fn contains_mate(&self, read_id: &str, mate_label: &str) -> bool {
        let bit = match mate_label {
            "R1" => 0b01,
            "R2" => 0b10,
            _ => 0,
        };
        bit != 0
            && self
                .by_read
                .get(read_id)
                .is_some_and(|mask| mask & bit != 0)
    }

    /// Tests whether both mates are already represented by Scan1 display rows.
    fn contains_both_mates(&self, read_id: &str) -> bool {
        self.by_read
            .get(read_id)
            .is_some_and(|mask| mask & 0b11 == 0b11)
    }
}

/// Borrowed identity for a main Scan2 hit that display rescue would rediscover.
///
/// The display path emits rows without the priority column, but the candidate
/// identity before validation is the same `(mate, cigar, chr, sites, signals)`
/// payload. Keeping this as borrowed fields lets Scan2 suppress exact display
/// duplicates before the expensive HG2 validator call without allocating a
/// formatted key for every display candidate.
#[derive(Clone, Copy)]
struct Scan2DisplaySkipKey<'a> {
    mate: &'a str,
    cigar: &'a str,
    chrom: &'a str,
    site1: &'a str,
    site2: &'a str,
    signal1: &'a str,
    signal2: &'a str,
    sum_q: &'a str,
}

impl<'a> Scan2DisplaySkipKey<'a> {
    /// Tests whether a display candidate would duplicate this main Scan2 hit.
    fn matches(&self, mate: &str, cigar: &str, chrom: &str, cand: &CandidateBreakpoint) -> bool {
        self.mate == mate
            && self.cigar == cigar
            && self.chrom == chrom
            && self.site1 == cand.data[0]
            && self.site2 == cand.data[1]
            && self.signal1 == cand.data[2]
            && self.signal2 == cand.data[3]
            && self.sum_q == cand.data[4]
    }
}

impl Scan2Profile {
    /// Checks whether release profiling is enabled for the current Scan2 run.
    fn enabled_from_env() -> bool {
        scan2_profile_enabled()
    }

    /// Emits the aggregated Scan2 timing summary.
    fn report(&self) {
        let wall_total_ns = self.wall_total_ns.load(Ordering::Relaxed);
        let shard_total_ns = self.shard_total_ns.load(Ordering::Relaxed);
        let group_process_ns = self.group_process_ns.load(Ordering::Relaxed);
        let validator_call_ns = self.validator_call_ns.load(Ordering::Relaxed);
        let display_ns = self.display_ns.load(Ordering::Relaxed);
        let segment_sidecar_ns = self.segment_sidecar_ns.load(Ordering::Relaxed);
        let non_bsj_sidecar_ns = self.non_bsj_sidecar_ns.load(Ordering::Relaxed);
        let write_ns = self.write_ns.load(Ordering::Relaxed);
        let merge_ns = self.merge_ns.load(Ordering::Relaxed);
        let records = self.records.load(Ordering::Relaxed);
        let groups = self.groups.load(Ordering::Relaxed);
        let candidate_checks = self.candidate_checks.load(Ordering::Relaxed);
        let candidate_hits = self.candidate_hits.load(Ordering::Relaxed);
        let segment_sidecar_rows = self.segment_sidecar_rows.load(Ordering::Relaxed);
        let segment_sidecar_bytes = self.segment_sidecar_bytes.load(Ordering::Relaxed);
        let non_bsj_sidecar_rows = self.non_bsj_sidecar_rows.load(Ordering::Relaxed);
        let non_bsj_sidecar_bytes = self.non_bsj_sidecar_bytes.load(Ordering::Relaxed);
        let accounted_ns = group_process_ns
            .saturating_add(display_ns)
            .saturating_add(segment_sidecar_ns)
            .saturating_add(non_bsj_sidecar_ns)
            .saturating_add(write_ns);
        let other_shard_ns = shard_total_ns.saturating_sub(accounted_ns);
        let pct = |part: u64, whole: u64| -> f64 {
            if whole == 0 {
                0.0
            } else {
                part as f64 * 100.0 / whole as f64
            }
        };
        emit_perf_line(&format!(
            "[PROFILE_SCAN2] wall_ms={:.3} shard_work_ms={:.3} merge_ms={:.3} records={} groups={} candidate_checks={} candidate_hits={}",
            wall_total_ns as f64 / 1_000_000.0,
            shard_total_ns as f64 / 1_000_000.0,
            merge_ns as f64 / 1_000_000.0,
            records,
            groups,
            candidate_checks,
            candidate_hits,
        ));
        emit_perf_line(&format!(
            "[PROFILE_SCAN2] shard_breakdown_ms group_process={:.3} ({:.1}%) validator={:.3} ({:.1}% of group) display={:.3} ({:.1}%) segment_sidecar={:.3} ({:.1}%) non_bsj_sidecar={:.3} ({:.1}%) write={:.3} ({:.1}%) other={:.3} ({:.1}%)",
            group_process_ns as f64 / 1_000_000.0,
            pct(group_process_ns, shard_total_ns),
            validator_call_ns as f64 / 1_000_000.0,
            pct(validator_call_ns, group_process_ns),
            display_ns as f64 / 1_000_000.0,
            pct(display_ns, shard_total_ns),
            segment_sidecar_ns as f64 / 1_000_000.0,
            pct(segment_sidecar_ns, shard_total_ns),
            non_bsj_sidecar_ns as f64 / 1_000_000.0,
            pct(non_bsj_sidecar_ns, shard_total_ns),
            write_ns as f64 / 1_000_000.0,
            pct(write_ns, shard_total_ns),
            other_shard_ns as f64 / 1_000_000.0,
            pct(other_shard_ns, shard_total_ns),
        ));
        emit_perf_line(&format!(
            "[PROFILE_SCAN2_SIDECAR] bsj_rows={} bsj_bytes={} non_bsj_rows={} non_bsj_bytes={}",
            segment_sidecar_rows,
            segment_sidecar_bytes,
            non_bsj_sidecar_rows,
            non_bsj_sidecar_bytes,
        ));
    }
}

impl Scan2 {
    /// Mirrors Java `String.hashCode` so Scan2 can rebuild Java-like `HashSet`
    /// iteration order when constructing the candidate index.
    #[inline]
    fn java_string_hash(s: &str) -> i32 {
        // Java String.hashCode(): h = 31*h + ch (32-bit overflow).
        let mut h: i32 = 0;
        for ch in s.chars() {
            h = h.wrapping_mul(31).wrapping_add(ch as i32);
        }
        h
    }

    /// Applies the same spread step Java `HashMap` uses before bucket selection.
    #[inline]
    fn java_hash_spread(h: i32) -> u32 {
        let x = h as u32;
        x ^ (x >> 16)
    }

    /// Computes the minimum Java-style `HashSet` backing capacity for `n` items.
    #[inline]
    fn java_hashset_capacity(n: usize) -> usize {
        // Java HashMap defaults: initial capacity 16, load factor 0.75.
        // Capacity is power-of-two and grows so that size <= cap * 0.75.
        let mut cap = 16usize;
        while n > ((cap as f64) * 0.75f64) as usize {
            cap = cap.saturating_mul(2);
        }
        cap.max(1)
    }

    /// Creates a Scan2 runner with CIRI3-compatible thresholds and bucket size.
    pub fn new(min_mapq_uni: i32, linear_range_size_min: i32, seq_len: i32) -> Self {
        Self {
            min_mapq_uni,
            linear_range_size_min,
            index1: HashMap::new(),
            index2: HashMap::new(),
            site_array1: HashMap::new(),
            site_array2: HashMap::new(),
            fsj_map: HashMap::new(),
            mem_limit: 2 * 1024 * 1024 * 1024,
            scan1_ids: HashSet::new(),
            seq_len,
            rescued_reads: 0,
            final_bsj_reads: 0,
        }
    }

    /// Binary-search lower bound used by the sorted candidate buckets.
    fn lower_bound(list: &[CandidateBreakpoint], key: i32) -> usize {
        let mut lo = 0usize;
        let mut hi = list.len();
        while lo < hi {
            let mid = (lo + hi) / 2;
            if list[mid].site < key {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    }

    /// Returns the candidate slice covering one Java bucket.
    ///
    /// The index vectors are sorted once during `build_index`, so hot-path bucket
    /// scans can stay branch-light and avoid re-filtering the whole chromosome list.
    fn bucket_range(list: &[CandidateBreakpoint], bucket: i32, bucket_size: i32) -> (usize, usize) {
        if bucket < 0 || bucket_size <= 0 {
            return (0, 0);
        }
        let start = bucket.saturating_mul(bucket_size);
        let end_exclusive = (bucket + 1).saturating_mul(bucket_size);
        let l = Self::lower_bound(list, start);
        let r = Self::lower_bound(list, end_exclusive);
        (l, r)
    }

    /// Overrides the per-thread page-cache eviction threshold.
    pub fn set_mem_limit(&mut self, limit: u64) {
        self.mem_limit = limit;
    }

    /// Drops Scan2-only indexing state once rescue is complete.
    ///
    /// Summary only needs the merged `fsj_map`, so the large candidate indexes and
    /// Scan1 read-id filter can be released explicitly before the final stage to
    /// keep peak RSS lower on whole-genome inputs.
    pub fn release_working_set(&mut self) {
        self.index1.clear();
        self.index1.shrink_to_fit();
        self.index2.clear();
        self.index2.shrink_to_fit();
        self.site_array1.clear();
        self.site_array1.shrink_to_fit();
        self.site_array2.clear();
        self.site_array2.shrink_to_fit();
        self.scan1_ids.clear();
        self.scan1_ids.shrink_to_fit();
    }

    /// Returns the temporary path used for one shard-local FSJ spill file.
    #[inline]
    fn shard_fsj_path(output_fsj: &str, shard_idx: usize) -> String {
        part_path(output_fsj, shard_idx)
    }

    /// Returns the shard-local read-level FSJ audit path.
    fn shard_fsj_pair_dump_path(output_fsj: &str, shard_idx: usize) -> String {
        format!("{}.pairs", Self::shard_fsj_path(output_fsj, shard_idx))
    }

    /// Writes one shard-local FSJ map to disk so the main thread can merge it
    /// later without retaining every shard map in memory at once.
    fn write_fsj_shard(path: &str, local_fsj: &HashMap<String, i32>) -> Result<()> {
        let mut writer = BufWriter::with_capacity(256 * 1024, File::create(path)?);
        for (key, count) in local_fsj {
            writeln!(writer, "{}\t{}", key, count)?;
        }
        writer.flush()?;
        Ok(())
    }

    /// Loads Scan1 BSJ1 output and builds the de-duplicated Scan2 candidate index.
    ///
    /// This function intentionally mirrors Java's unique-site indexing semantics:
    /// the traversal order of candidates must match Java, because Scan2 returns on
    /// the first valid non-`2` tag it encounters.
    ///
    /// Two details here are easy to "simplify" incorrectly:
    /// - Scan2 de-duplicates by circ site payload before building the index, so
    ///   repeated BSJ1 rows for the same circ do not create extra candidates.
    /// - The insertion order is preserved through `order`, then reused after the
    ///   final per-bucket sort, because Java's `HashSet -> ArrayList` path still
    ///   leaves a stable first-hit order once the same inputs are replayed.
    pub fn build_index(&mut self, bsj1_file: &str) -> Result<()> {
        // Java parity: Scan2 index is built from unique circ sites, not all BSJ1 rows.
        // Equivalent Java flow:
        //   BSJ1 -> scan1IdMap + chrCircSiteMap(HashSet) -> siteArray/siteMap.
        let file = File::open(bsj1_file)?;
        let reader = BufReader::new(file);
        let mut chr_circ_site_seen: HashMap<String, HashSet<String>> = HashMap::new();
        let mut chr_circ_site_insertion: HashMap<String, Vec<String>> = HashMap::new();
        for line_res in reader.lines() {
            let line = line_res?;
            let p: Vec<&str> = line.split('\t').collect();
            if p.len() < 10 || !bsj_is_summary_priority(&p) {
                continue;
            }
            let payload_start = bsj_payload_start(&p);
            if p.len() <= payload_start + 8 {
                continue;
            }
            self.scan1_ids.insert(p[0].to_string());
            let chr = p[payload_start + 2].to_string();
            let site_infor = p[payload_start + 3..].join("\t");
            let seen = chr_circ_site_seen
                .entry(chr.clone())
                .or_insert_with(HashSet::new);
            if seen.insert(site_infor.clone()) {
                chr_circ_site_insertion
                    .entry(chr)
                    .or_insert_with(Vec::new)
                    .push(site_infor);
            }
        }

        let mut order_counter: usize = 0;
        for (chr, insertion_list) in chr_circ_site_insertion {
            // Java parity: emulate HashSet iteration order:
            // bucket index ascending, and insertion order within each bucket.
            let cap = Self::java_hashset_capacity(insertion_list.len());
            let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); cap];
            for (idx, s) in insertion_list.iter().enumerate() {
                let spread = Self::java_hash_spread(Self::java_string_hash(s));
                let bucket_idx = (spread as usize) & (cap - 1);
                buckets[bucket_idx].push(idx);
            }

            for bucket in buckets {
                for idx in bucket {
                    let arr: Vec<String> = insertion_list[idx]
                        .split('\t')
                        .map(|s| s.to_string())
                        .collect();
                    if arr.len() < 6 {
                        continue;
                    }
                    let site1 = arr[0].parse().unwrap_or(0);
                    let site2 = arr[1].parse().unwrap_or(0);
                    self.fsj_map
                        .entry(format!("{}\t{}\t{}", chr, arr[0], arr[1]))
                        .or_insert(0);
                    self.index1
                        .entry(chr.clone())
                        .or_insert_with(Vec::new)
                        .push(CandidateBreakpoint {
                            site: site1,
                            order: order_counter,
                            data: arr.clone(),
                        });
                    self.index2
                        .entry(chr.clone())
                        .or_insert_with(Vec::new)
                        .push(CandidateBreakpoint {
                            site: site2,
                            order: order_counter,
                            data: arr,
                        });
                    let bucket_size = self.seq_len.max(1);
                    self.site_array1
                        .entry(chr.clone())
                        .or_insert_with(HashSet::new)
                        .insert(site1 / bucket_size);
                    self.site_array2
                        .entry(chr.clone())
                        .or_insert_with(HashSet::new)
                        .insert(site2 / bucket_size);
                    order_counter += 1;
                }
            }
        }
        for list in self.index1.values_mut() {
            list.sort_by_key(|x| (x.site, x.order));
        }
        for list in self.index2.values_mut() {
            list.sort_by_key(|x| (x.site, x.order));
        }
        Ok(())
    }

    /// Adds cohort-level external BSJ sites to the Scan2 candidate indexes.
    ///
    /// `--circ` candidates come from BED6 catalog coordinates rather than from
    /// Scan1 read evidence. They are therefore appended after the Java-parity
    /// Scan1 candidates: default runs stay unchanged, while a second-pass run can
    /// test additional cohort-supported BSJ boundaries against this sample's own
    /// reads. The left/right signal fields are reconstructed from the reference
    /// so the existing `is_bsj_hg2` validator can remain the single source of
    /// sequence-match decisions.
    pub fn add_external_circ_candidates(
        &mut self,
        candidates: &[CircCatalogRecord],
        reference: &HashMap<String, String>,
    ) -> Result<usize> {
        let mut seen = HashSet::new();
        let mut next_order = 0usize;
        for (chrom, list) in &self.index1 {
            for candidate in list {
                next_order = next_order.max(candidate.order + 1);
                if candidate.data.len() >= 3 {
                    seen.insert(format!(
                        "{}\t{}\t{}\t{}",
                        chrom, candidate.data[0], candidate.data[1], candidate.data[2]
                    ));
                }
            }
        }

        let mut added = 0usize;
        let bucket_size = self.seq_len.max(1);
        for candidate in candidates {
            if candidate.start < 1 || candidate.end < candidate.start {
                bail!(
                    "invalid external circ candidate {}:{}-{}",
                    candidate.chrom,
                    candidate.start,
                    candidate.end
                );
            }
            let Some(chr_seq) = reference.get(&candidate.chrom) else {
                bail!(
                    "external circ candidate {}:{}-{} uses chromosome absent from reference",
                    candidate.chrom,
                    candidate.start,
                    candidate.end
                );
            };
            if candidate.end as usize > chr_seq.len() {
                bail!(
                    "external circ candidate {}:{}-{} exceeds reference length {}",
                    candidate.chrom,
                    candidate.start,
                    candidate.end,
                    chr_seq.len()
                );
            }
            let dedup_key = format!(
                "{}\t{}\t{}\t{}",
                candidate.chrom, candidate.start, candidate.end, candidate.strand
            );
            if !seen.insert(dedup_key) {
                continue;
            }
            let signal_left =
                reference_window_or_n(chr_seq, candidate.start - 3, candidate.start - 1);
            let signal_right = reference_window_or_n(chr_seq, candidate.end, candidate.end + 2);
            let payload = vec![
                candidate.start.to_string(),
                candidate.end.to_string(),
                candidate.strand.clone(),
                signal_left,
                signal_right,
                "1".to_string(),
            ];
            self.fsj_map
                .entry(format!(
                    "{}\t{}\t{}",
                    candidate.chrom, candidate.start, candidate.end
                ))
                .or_insert(0);
            self.index1
                .entry(candidate.chrom.clone())
                .or_insert_with(Vec::new)
                .push(CandidateBreakpoint {
                    site: candidate.start,
                    order: next_order,
                    data: payload.clone(),
                });
            self.index2
                .entry(candidate.chrom.clone())
                .or_insert_with(Vec::new)
                .push(CandidateBreakpoint {
                    site: candidate.end,
                    order: next_order,
                    data: payload,
                });
            self.site_array1
                .entry(candidate.chrom.clone())
                .or_insert_with(HashSet::new)
                .insert(candidate.start / bucket_size);
            self.site_array2
                .entry(candidate.chrom.clone())
                .or_insert_with(HashSet::new)
                .insert(candidate.end / bucket_size);
            next_order += 1;
            added += 1;
        }

        for list in self.index1.values_mut() {
            list.sort_by_key(|x| (x.site, x.order));
        }
        for list in self.index2.values_mut() {
            list.sort_by_key(|x| (x.site, x.order));
        }
        Ok(added)
    }

    /// Builds a mate-level display index from the post-Summary Scan1 display rows.
    ///
    /// This path is intentionally isolated from the parity index above. It keeps
    /// the user-facing `.bsj` display self-consistent without feeding any of the
    /// extra mate-level rows back into `.out`.
    fn build_display_index(&mut self, bsj1_display_file: &str) -> Result<Scan2DisplayClaims> {
        self.index1.clear();
        self.index2.clear();
        self.site_array1.clear();
        self.site_array2.clear();
        let file = File::open(bsj1_display_file)?;
        let reader = BufReader::new(file);
        let mut claims = Scan2DisplayClaims::new();
        let mut chr_circ_site_seen: HashMap<String, HashSet<String>> = HashMap::new();
        let mut chr_circ_site_insertion: HashMap<String, Vec<String>> = HashMap::new();
        for line_res in reader.lines() {
            let line = line_res?;
            let p: Vec<&str> = line.split('\t').collect();
            if p.len() < 10 {
                continue;
            }
            let payload_start = bsj_payload_start(&p);
            let mate_label = if payload_start == 3 { p[1] } else { "NA" };
            if p.len() <= payload_start + 8 {
                continue;
            }
            claims.insert(p[0], mate_label);
            let chr = p[payload_start + 2].to_string();
            let site_infor = p[payload_start + 3..].join("\t");
            let seen = chr_circ_site_seen
                .entry(chr.clone())
                .or_insert_with(HashSet::new);
            if seen.insert(site_infor.clone()) {
                chr_circ_site_insertion
                    .entry(chr)
                    .or_insert_with(Vec::new)
                    .push(site_infor);
            }
        }

        let mut order_counter: usize = 0;
        for (chr, sites) in chr_circ_site_insertion {
            let mut list1 = Vec::new();
            let mut list2 = Vec::new();
            let mut set1 = HashSet::new();
            let mut set2 = HashSet::new();
            for site in sites {
                let data: Vec<String> = site.split('\t').map(|s| s.to_string()).collect();
                if data.len() != 6 {
                    continue;
                }
                let site1 = data[0].parse::<i32>().unwrap_or(0);
                let site2 = data[1].parse::<i32>().unwrap_or(0);
                let num1 = site1 / self.seq_len.max(1);
                let num2 = site2 / self.seq_len.max(1);
                set1.insert(num1);
                set2.insert(num2);
                list1.push(CandidateBreakpoint {
                    site: site1,
                    order: order_counter,
                    data: data.clone(),
                });
                list2.push(CandidateBreakpoint {
                    site: site2,
                    order: order_counter,
                    data,
                });
                order_counter += 1;
            }
            list1.sort_by_key(|c| (c.site, c.order));
            list2.sort_by_key(|c| (c.site, c.order));
            self.index1.insert(chr.clone(), list1);
            self.index2.insert(chr.clone(), list2);
            self.site_array1.insert(chr.clone(), set1);
            self.site_array2.insert(chr, set2);
        }
        Ok(claims)
    }

    /// Collects FSJ keys overlapped by a linear alignment span.
    ///
    /// This mirrors Java `GetFSJClass.getFSJ(...)`: only the two buckets touched
    /// by `[start_tem, end_tem]` are scanned, with reverse traversal on `num1`
    /// and forward traversal on `num2`. The broader lower-bound scan is cheaper
    /// to write but not Java-compatible, and it over-counts FSJs on chr1 while
    /// leaving BSJ rescue unchanged.
    ///
    /// `style` is kept in Java's original encoding (`0`, `1`, `-1`, `10`) on
    /// purpose. Those magic-looking values control which bucket edge is treated
    /// as inclusive for full-match, MS, SM, and SMS alignments, and normalizing
    /// them into a more abstract enum made earlier parity checks harder to audit.
    fn collect_fsj_keys_in_range(
        &self,
        chr: &str,
        start_tem: i32,
        end_tem: i32,
        style: i32,
        out: &mut HashSet<String>,
    ) {
        let bucket_size = self.seq_len.max(1);
        let num1 = start_tem / bucket_size;
        let num2 = end_tem / bucket_size;
        if let Some(list1) = self.index1.get(chr) {
            if self.site_array1.get(chr).is_some_and(|s| s.contains(&num1)) {
                let (l1, r1) = Self::bucket_range(list1, num1, bucket_size);
                for idx in (l1..r1).rev() {
                    let cand = &list1[idx];
                    if cand.site >= start_tem {
                        if cand.site <= end_tem || style == 10 || style == 1 {
                            out.insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                        }
                    } else {
                        break;
                    }
                }
            }
            if num2 != num1 && self.site_array1.get(chr).is_some_and(|s| s.contains(&num2)) {
                let (l2, r2) = Self::bucket_range(list1, num2, bucket_size);
                for idx in l2..r2 {
                    let cand = &list1[idx];
                    if cand.site <= end_tem {
                        if cand.site >= start_tem || style == 10 || style == -1 {
                            out.insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                        }
                    } else {
                        break;
                    }
                }
            }
        }
        if let Some(list2) = self.index2.get(chr) {
            if self.site_array2.get(chr).is_some_and(|s| s.contains(&num1)) {
                let (l1, r1) = Self::bucket_range(list2, num1, bucket_size);
                for idx in (l1..r1).rev() {
                    let cand = &list2[idx];
                    if cand.site >= start_tem {
                        if cand.site <= end_tem || style == 10 || style == 1 {
                            out.insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                        }
                    } else {
                        break;
                    }
                }
            }
            if num2 != num1 && self.site_array2.get(chr).is_some_and(|s| s.contains(&num2)) {
                let (l2, r2) = Self::bucket_range(list2, num2, bucket_size);
                for idx in l2..r2 {
                    let cand = &list2[idx];
                    if cand.site <= end_tem {
                        if cand.site >= start_tem || style == 10 || style == -1 {
                            out.insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                        }
                    } else {
                        break;
                    }
                }
            }
        }
    }

    /// Collects FSJ keys directly from side-spanning linear alignment evidence.
    ///
    /// This minibwa-oriented experimental mode is intentionally separate from
    /// Java-compatible bucket spans: `M` blocks contribute only the interval
    /// where both flanks have at least `anchor` aligned bases, and `N` gaps
    /// contribute only small windows around splice boundaries with adequate
    /// flanking match blocks. `D` gaps are ignored because large deletion spans
    /// were a recurrent false-FSJ source in minibwa.
    fn collect_direct_side_fsj_keys(
        &self,
        chr: &str,
        aln_pos: i32,
        cigar: &str,
        anchor: i32,
        out: &mut HashSet<String>,
    ) {
        if cigar == "*" || aln_pos <= 0 {
            return;
        }
        let mut ref_pos = aln_pos;
        let mut prev_match: Option<(i32, i32)> = None;
        let mut pending_skip_after: Option<(i32, i32, i32, i32)> = None;
        let mut count = 0_i32;
        let mut has_count = false;
        for op in cigar.chars() {
            if op.is_ascii_digit() {
                count = match count
                    .checked_mul(10)
                    .and_then(|value| value.checked_add(op.to_digit(10).unwrap_or(0) as i32))
                {
                    Some(value) => value,
                    None => return,
                };
                has_count = true;
                continue;
            }
            if !has_count {
                return;
            }
            match op {
                'M' | '=' | 'X' => {
                    let block_start = ref_pos;
                    let block_end = ref_pos + count - 1;
                    let direct_start = block_start + anchor - 1;
                    let direct_end = block_end - anchor + 1;
                    if direct_start <= direct_end {
                        self.collect_fsj_keys_in_range(chr, direct_start, direct_end, 0, out);
                    }
                    if let Some((_skip_start, _skip_end, left_start, left_end)) =
                        pending_skip_after.take()
                    {
                        if left_end - left_start + 1 >= anchor
                            && block_end - block_start + 1 >= anchor
                        {
                            self.collect_fsj_keys_in_range(
                                chr,
                                left_end - DIRECT_SIDE_FSJ_SPLICE_TOLERANCE,
                                left_end + DIRECT_SIDE_FSJ_SPLICE_TOLERANCE,
                                0,
                                out,
                            );
                            self.collect_fsj_keys_in_range(
                                chr,
                                block_start - DIRECT_SIDE_FSJ_SPLICE_TOLERANCE,
                                block_start + DIRECT_SIDE_FSJ_SPLICE_TOLERANCE,
                                0,
                                out,
                            );
                        }
                    }
                    prev_match = Some((block_start, block_end));
                    ref_pos += count;
                }
                'N' => {
                    pending_skip_after =
                        prev_match.map(|(start, end)| (end + 1, ref_pos + count - 1, start, end));
                    ref_pos += count;
                }
                'D' => {
                    pending_skip_after = None;
                    ref_pos += count;
                }
                'I' | 'S' | 'H' | 'P' => {
                    pending_skip_after = None;
                }
                _ => return,
            }
            count = 0;
            has_count = false;
        }
    }

    /// Collects FSJ keys with optional strong side-anchor filtering.
    ///
    /// The default CIRI3-compatible path keeps Java's span-based FSJ counting.
    /// `CIRI_SCAN2_STRONG_FSJ=1` is an evaluation mode for minibwa: a candidate
    /// circRNA is counted only when the current linear alignment has one BSJ side
    /// inside a match block with enough matched bases on both sides. Large `D`
    /// spans therefore no longer masquerade as strong forward evidence.
    fn collect_fsj_keys_for_alignment(
        &self,
        chr: &str,
        start_tem: i32,
        end_tem: i32,
        style: i32,
        aln_pos: i32,
        cigar: &str,
        out: &mut HashSet<String>,
    ) {
        if let Some(anchor) = direct_side_fsj_anchor() {
            self.collect_direct_side_fsj_keys(chr, aln_pos, cigar, anchor, out);
            return;
        }
        if !strong_fsj_filter_enabled() {
            self.collect_fsj_keys_in_range(chr, start_tem, end_tem, style, out);
            if let Some(anchor) = direct_side_fsj_rescue_anchor() {
                self.collect_direct_side_fsj_keys(chr, aln_pos, cigar, anchor, out);
            }
            return;
        }
        let mut candidates = HashSet::new();
        self.collect_fsj_keys_in_range(chr, start_tem, end_tem, style, &mut candidates);
        for key in candidates {
            for side_key in strong_fsj_side_keys(&key, aln_pos, cigar) {
                out.insert(side_key);
            }
        }
        if let Some(anchor) = direct_side_fsj_rescue_anchor() {
            self.collect_direct_side_fsj_keys(chr, aln_pos, cigar, anchor, out);
        }
    }

    /// Inserts one validator-derived FSJ competitor under the strong-FSJ policy.
    ///
    /// Scan2 adds `tag==0` candidates as linear competitors. In enhanced
    /// strong-FSJ mode, those candidates must pass the same side-anchor test as
    /// span-derived FSJ keys; otherwise BSJ-like split reads can inflate the
    /// forward denominator without a strong linear explanation.
    fn insert_fsj_key_for_alignment(
        chr: &str,
        site1: &str,
        site2: &str,
        aln_pos: i32,
        cigar: &str,
        out: &mut HashSet<String>,
    ) {
        if direct_side_fsj_anchor().is_some() {
            return;
        }
        let key = format!("{chr}\t{site1}\t{site2}");
        if !strong_fsj_filter_enabled() {
            out.insert(key);
        } else {
            for side_key in strong_fsj_side_keys(&key, aln_pos, cigar) {
                out.insert(side_key);
            }
        }
    }

    /// Dispatches to the SAM or BAM Scan2 implementation.
    pub fn run(
        &mut self,
        sam_file: &str,
        output_bsj2: &str,
        output_fsj: &str,
        chr_tcga_map: &HashMap<String, String>,
    ) -> Result<()> {
        self.run_with_display(sam_file, output_bsj2, output_fsj, None, None, chr_tcga_map)
    }

    /// Runs Scan2 and optionally emits mate-level `priority=0` rescue rows into
    /// the same `.bsj2` stream during the scan.
    pub fn run_with_display(
        &mut self,
        sam_file: &str,
        output_bsj2: &str,
        output_fsj: &str,
        scan1_display_path: Option<&str>,
        display_output_path: Option<&str>,
        chr_tcga_map: &HashMap<String, String>,
    ) -> Result<()> {
        self.run_with_display_and_segments(
            sam_file,
            output_bsj2,
            output_fsj,
            scan1_display_path,
            display_output_path,
            None,
            None,
            chr_tcga_map,
        )
        .map(|_| ())
    }

    /// Runs Scan2 while also writing a segments-evidence sidecar for rescued BSJ reads.
    ///
    /// The extra sidecar is independent from `.bsj2` and FSJ counting. It records
    /// mapper alignment blocks for read groups with Scan2 BSJ evidence so the
    /// post-Summary segments stage can reconstruct confirmed BSJ rows without
    /// scanning the input alignment file again.
    pub fn run_with_display_and_segments(
        &mut self,
        sam_file: &str,
        output_bsj2: &str,
        output_fsj: &str,
        scan1_display_path: Option<&str>,
        display_output_path: Option<&str>,
        segments_output_path: Option<&str>,
        non_bsj_segments_output_path: Option<&str>,
        chr_tcga_map: &HashMap<String, String>,
    ) -> Result<Scan2SegmentArtifacts> {
        use crate::sam_bam::{detect_format, InputFormat};
        let display_scan2 = if let Some(scan1_path) = scan1_display_path {
            let mut helper =
                Scan2::new(self.min_mapq_uni, self.linear_range_size_min, self.seq_len);
            let claims = helper.build_display_index(scan1_path)?;
            Some((helper, claims))
        } else {
            None
        };
        let format = detect_format(sam_file)?;
        match format {
            InputFormat::Sam => self.run_sam_with_display(
                sam_file,
                output_bsj2,
                output_fsj,
                display_output_path,
                segments_output_path,
                non_bsj_segments_output_path,
                display_scan2.as_ref(),
                chr_tcga_map,
            ),
            InputFormat::Bam => self.run_bam_with_display(
                sam_file,
                output_bsj2,
                output_fsj,
                display_output_path,
                segments_output_path,
                non_bsj_segments_output_path,
                display_scan2.as_ref(),
                chr_tcga_map,
            ),
        }
    }

    /// Runs Scan2 on SAM input.
    ///
    /// The SAM path preserves queryname groups before batch-parallel processing,
    /// matching the BAM path's read-group semantics while avoiding shard-boundary
    /// FSJ duplication.
    pub fn run_sam(
        &mut self,
        sam_file: &str,
        output_bsj2: &str,
        output_fsj: &str,
        chr_tcga_map: &HashMap<String, String>,
    ) -> Result<()> {
        self.run_sam_with_display(
            sam_file,
            output_bsj2,
            output_fsj,
            None,
            None,
            None,
            None,
            chr_tcga_map,
        )
        .map(|_| ())
    }

    fn run_sam_with_display(
        &mut self,
        sam_file: &str,
        output_bsj2: &str,
        output_fsj: &str,
        display_output_path: Option<&str>,
        segments_output_path: Option<&str>,
        non_bsj_segments_output_path: Option<&str>,
        display_scan2: Option<&(Scan2, Scan2DisplayClaims)>,
        chr_tcga_map: &HashMap<String, String>,
    ) -> Result<Scan2SegmentArtifacts> {
        let run_started = Instant::now();
        let profile = if Scan2Profile::enabled_from_env() {
            Some(Scan2Profile::default())
        } else {
            None
        };
        let profile_ref = profile.as_ref();
        let file_size = std::fs::metadata(sam_file)?.len();

        let pb = ProgressBar::new(file_size as u64);
        pb.set_style(ProgressStyle::default_bar().template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {percent:>3}% ({eta}) {msg}")?.progress_chars("#>-"));
        pb.set_message("");
        let shard_out = part_path(output_bsj2, 0);
        let display_shard_out = display_output_path.map(|path| part_path(path, 0));
        let segments_shard_out = segments_output_path.map(|path| part_path(path, 0));
        let non_bsj_segments_shard_out =
            non_bsj_segments_output_path.map(|path| part_path(path, 0));
        let fsj_out = Self::shard_fsj_path(output_fsj, 0);
        self.process_sam_file_to_file(
            sam_file,
            chr_tcga_map,
            &pb,
            &shard_out,
            &fsj_out,
            display_scan2,
            display_shard_out.as_deref(),
            segments_shard_out.as_deref(),
            non_bsj_segments_shard_out.as_deref(),
            profile_ref,
        )?;

        pb.set_style(ProgressStyle::default_bar().template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {msg}")?.progress_chars("#>-"));
        pb.finish_with_message("Completed");
        let merge_started = Instant::now();
        self.merge_shards_and_fsj(output_bsj2, output_fsj, 1)?;
        if let Some(profile) = profile_ref {
            profile
                .merge_ns
                .fetch_add(merge_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile
                .wall_total_ns
                .fetch_add(run_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile.report();
            report_scan2_hg_profile();
        }
        if let Some(path) = display_output_path {
            self.merge_display_shards(path, 1)?;
        }
        if let Some(path) = segments_output_path {
            self.merge_segment_shards(path, 1)?;
        }
        let non_bsj_segment_evidence_paths = non_bsj_segments_output_path
            .map(|path| vec![part_path(path, 0)])
            .unwrap_or_default();
        Ok(Scan2SegmentArtifacts {
            non_bsj_segment_evidence_paths,
        })
    }

    /// Merges per-shard Scan2 BSJ2 outputs and FSJ spill files.
    ///
    /// `rescued_reads` and `final_bsj_reads` are stage-level accounting only.
    /// The user-facing "final BSJ reads" summary is recomputed later from the
    /// clustered `.out`, because Summary can still merge circ families without
    /// changing the raw `.bsj1/.bsj2` membership.
    fn merge_shards_and_fsj(
        &mut self,
        output_bsj2: &str,
        output_fsj: &str,
        num_threads: usize,
    ) -> Result<()> {
        let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(output_bsj2)?);
        let mut pair_dump_writer = if let Some(path) = fsj_pair_dump_path() {
            Some(BufWriter::with_capacity(1024 * 1024, File::create(path)?))
        } else {
            None
        };
        let mut rescued_ids = HashSet::new();
        for i in 0..num_threads {
            let shard_path = part_path(output_bsj2, i);
            if let Ok(shard_file) = File::open(&shard_path) {
                let mut shard_reader = BufReader::new(shard_file);
                let mut line = String::new();
                while shard_reader.read_line(&mut line)? != 0 {
                    let trimmed = line.trim_end();
                    if !trimmed.is_empty() {
                        let parts: Vec<&str> = trimmed.split('\t').collect();
                        if !parts.is_empty() && bsj_is_summary_priority(&parts) {
                            rescued_ids.insert(parts[0].to_string());
                        }
                        writeln!(writer, "{}", trimmed)?;
                    }
                    line.clear();
                }
            }
            let _ = std::fs::remove_file(shard_path);
        }
        for i in 0..num_threads {
            let fsj_path = Self::shard_fsj_path(output_fsj, i);
            if let Ok(fsj_file) = File::open(&fsj_path) {
                let mut fsj_reader = BufReader::new(fsj_file);
                let mut line = String::new();
                while fsj_reader.read_line(&mut line)? != 0 {
                    let trimmed = line.trim_end();
                    if !trimmed.is_empty() {
                        let mut parts = trimmed.rsplitn(2, '\t');
                        let count = parts
                            .next()
                            .and_then(|s| s.parse::<i32>().ok())
                            .unwrap_or(0);
                        if let Some(key) = parts.next() {
                            *self.fsj_map.entry(key.to_string()).or_insert(0) += count;
                        }
                    }
                    line.clear();
                }
            }
            let _ = std::fs::remove_file(fsj_path);

            let pair_path = Self::shard_fsj_pair_dump_path(output_fsj, i);
            if let Some(pair_dump_writer) = pair_dump_writer.as_mut() {
                if let Ok(pair_file) = File::open(&pair_path) {
                    let mut pair_reader = BufReader::new(pair_file);
                    let mut line = String::new();
                    while pair_reader.read_line(&mut line)? != 0 {
                        pair_dump_writer.write_all(line.as_bytes())?;
                        line.clear();
                    }
                }
            }
            let _ = std::fs::remove_file(pair_path);
        }
        writer.flush()?;
        if let Some(mut pair_dump_writer) = pair_dump_writer {
            pair_dump_writer.flush()?;
        }
        self.rescued_reads = rescued_ids.len();
        self.final_bsj_reads = self.scan1_ids.len() + self.rescued_reads;
        Ok(())
    }

    /// Runs Scan2 on BAM input using BGZF-aware sharding.
    pub fn run_bam(
        &mut self,
        bam_file: &str,
        output_bsj2: &str,
        output_fsj: &str,
        chr_tcga_map: &HashMap<String, String>,
    ) -> Result<()> {
        self.run_bam_with_display(
            bam_file,
            output_bsj2,
            output_fsj,
            None,
            None,
            None,
            None,
            chr_tcga_map,
        )
        .map(|_| ())
    }

    fn run_bam_with_display(
        &mut self,
        bam_file: &str,
        output_bsj2: &str,
        output_fsj: &str,
        display_output_path: Option<&str>,
        segments_output_path: Option<&str>,
        non_bsj_segments_output_path: Option<&str>,
        display_scan2: Option<&(Scan2, Scan2DisplayClaims)>,
        chr_tcga_map: &HashMap<String, String>,
    ) -> Result<Scan2SegmentArtifacts> {
        use noodles::bam;
        let run_started = Instant::now();
        let profile = if Scan2Profile::enabled_from_env() {
            Some(Scan2Profile::default())
        } else {
            None
        };
        let profile_ref = profile.as_ref();
        let file = File::open(bam_file)?;
        let file_size = std::fs::metadata(bam_file)?.len();
        // SAFETY: the file descriptor remains open while the mapping is created;
        // the returned `Mmap` owns the mapping for the rest of this scope.
        let mmap = unsafe { Mmap::map(&file)? };
        let num_threads = bam_shard_count(mmap.len(), rayon::current_num_threads());
        let shard_size = mmap.len() / num_threads;
        let mut reader = bam::io::Reader::new(file);
        let header = reader.read_header()?;

        let pb = ProgressBar::new(file_size);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {percent:>3}% ({eta}) {msg}")?
                .progress_chars("#>-"),
        );
        pb.set_message("");
        pb.enable_steady_tick(Duration::from_millis(120));

        unsafe {
            // SAFETY: the pointer and length come from the live mapping, and
            // `madvise` only changes the kernel's read-ahead/cache policy.
            libc::madvise(
                mmap.as_ptr() as *mut libc::c_void,
                mmap.len(),
                libc::MADV_SEQUENTIAL,
            );
        }

        let header_ref = &header;
        (0..num_threads).into_par_iter().try_for_each(|i| {
            let start = i * shard_size;
            let end = if i == num_threads - 1 {
                mmap.len()
            } else {
                (i + 1) * shard_size
            };
            let shard_out = part_path(output_bsj2, i);
            let display_shard_out = display_output_path.map(|path| part_path(path, i));
            let segments_shard_out = segments_output_path.map(|path| part_path(path, i));
            let non_bsj_segments_shard_out =
                non_bsj_segments_output_path.map(|path| part_path(path, i));
            let fsj_out = Self::shard_fsj_path(output_fsj, i);
            self.process_bam_shard_to_file(
                &mmap,
                start,
                end,
                header_ref,
                chr_tcga_map,
                &pb,
                &shard_out,
                &fsj_out,
                display_scan2,
                display_shard_out.as_deref(),
                segments_shard_out.as_deref(),
                non_bsj_segments_shard_out.as_deref(),
                profile_ref,
            )
        })?;

        pb.set_position(file_size);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {msg}")?
                .progress_chars("#>-"),
        );
        pb.finish_with_message("Completed");
        let merge_started = Instant::now();
        self.merge_shards_and_fsj(output_bsj2, output_fsj, num_threads)?;
        if let Some(profile) = profile_ref {
            profile
                .merge_ns
                .fetch_add(merge_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile
                .wall_total_ns
                .fetch_add(run_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile.report();
            report_scan2_hg_profile();
        }
        if let Some(path) = display_output_path {
            self.merge_display_shards(path, num_threads)?;
        }
        if let Some(path) = segments_output_path {
            self.merge_segment_shards(path, num_threads)?;
        }
        let non_bsj_segment_evidence_paths = non_bsj_segments_output_path
            .map(|path| {
                (0..num_threads)
                    .map(|idx| part_path(path, idx))
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        Ok(Scan2SegmentArtifacts {
            non_bsj_segment_evidence_paths,
        })
    }

    /// Merges display shard files in shard order into one final temporary file.
    fn merge_display_shards(&self, output_path: &str, num_threads: usize) -> Result<()> {
        let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(output_path)?);
        for i in 0..num_threads {
            let shard_path = part_path(output_path, i);
            if let Ok(shard_file) = File::open(&shard_path) {
                let mut reader = BufReader::new(shard_file);
                let mut line = String::new();
                loop {
                    line.clear();
                    if reader.read_line(&mut line)? == 0 {
                        break;
                    }
                    writer.write_all(line.as_bytes())?;
                }
            }
            let _ = std::fs::remove_file(&shard_path);
        }
        writer.flush()?;
        Ok(())
    }

    /// Merges shard-local Scan2 segment-evidence sidecars in stable shard order.
    ///
    /// The merged file is not consumed by Summary. It is only a post-Summary
    /// reconstruction input, so keeping it separate from `.bsj2` preserves the
    /// CIRI3-compatible rescue/counting contract.
    fn merge_segment_shards(&self, output_path: &str, num_threads: usize) -> Result<()> {
        let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(output_path)?);
        for i in 0..num_threads {
            let shard_path = part_path(output_path, i);
            if let Ok(shard_file) = File::open(&shard_path) {
                let mut reader = BufReader::new(shard_file);
                let mut line = String::new();
                while reader.read_line(&mut line)? != 0 {
                    writer.write_all(line.as_bytes())?;
                    line.clear();
                }
            }
            let _ = std::fs::remove_file(&shard_path);
        }
        writer.flush()?;
        Ok(())
    }

    /// Formats read-group alignments for the Scan2 `<prefix>.segments2` sidecar.
    ///
    /// These rows intentionally mirror Scan1's sidecar protocol: raw mapper
    /// blocks are captured only for read groups with Scan2 BSJ evidence, and the
    /// post-Summary segments stage later applies circ confirmation, mate labels,
    /// and chain-level CIGAR reconstruction.
    fn segment_evidence_lines<'a>(
        read_id: &str,
        stage: &str,
        alignments: &[AlignmentRecord<'a>],
        reference: &HashMap<String, String>,
    ) -> Vec<String> {
        let mut rows = Vec::new();
        for aln in alignments {
            if aln.chrom.as_ref() == "*" || aln.cigar.as_ref() == "*" {
                continue;
            }
            let mate = if aln.flag & 0x40 != 0 { "R1" } else { "R2" };
            let clips = clip_sequence_payload(aln.cigar.as_ref(), aln.seq.as_ref());
            let cs = alignment_short_cs(
                aln.cigar.as_ref(),
                aln.seq.as_ref(),
                aln.chrom.as_ref(),
                aln.pos,
                reference,
            );
            rows.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                read_id,
                stage,
                mate,
                aln.flag,
                aln.chrom,
                aln.pos,
                aln.mapq,
                aln.cigar,
                aln.seq.len(),
                clips,
                cs,
                xa_payload(aln.xa.as_ref())
            ));
        }
        rows
    }

    /// Formats Scan2 read-level non-BSJ topology candidates.
    ///
    /// This sidecar is deliberately independent of the BSJ candidate index and
    /// Summary circRNA spans. It only uses read-level split/outward geometry to
    /// decide whether a group is worth re-evaluating after Summary. Soft clips
    /// are stored as compact `L:/R:` payloads so local clip placement remains
    /// available without repeating full read sequences across every alignment.
    /// The short `N` stage and compact alignment payload also carries
    /// short-form cs so sequence-aware `.segments.bam` can be reconstructed
    /// without rescanning FASTQ/BAM.
    fn non_bsj_segment_evidence_lines<'a>(
        read_id: &str,
        alignments: &[AlignmentRecord<'a>],
        reference: &HashMap<String, String>,
    ) -> Vec<String> {
        if !Self::may_support_non_bsj_segments(alignments) {
            return Vec::new();
        }
        let mut records = Vec::new();
        for aln in alignments {
            if aln.chrom.as_ref() == "*" || aln.cigar.as_ref() == "*" {
                continue;
            }
            let clips = clip_sequence_payload(aln.cigar.as_ref(), aln.seq.as_ref());
            let cs = alignment_short_cs(
                aln.cigar.as_ref(),
                aln.seq.as_ref(),
                aln.chrom.as_ref(),
                aln.pos,
                reference,
            );
            records.push(format!(
                "{}|{}|{}|{}|{}|{}|{}",
                aln.flag, aln.chrom, aln.pos, aln.mapq, aln.cigar, clips, cs
            ));
        }
        if records.is_empty() {
            Vec::new()
        } else {
            vec![format!("{}\t{}\t{}", read_id, "N", records.join(";"))]
        }
    }

    /// Cheap read-level filter for supplemental backward/outward candidates.
    fn may_support_non_bsj_segments<'a>(alignments: &[AlignmentRecord<'a>]) -> bool {
        Self::may_support_outward_segments(alignments)
            || Self::may_support_outward_segments_with_retention(
                alignments,
                Self::non_bsj_outward_retention_mode(),
            )
            || Self::may_support_backward_segments(alignments)
    }

    /// Reads the non-BSJ outward retention mode from the environment.
    ///
    /// The default is the recognition-first relaxed mode because simulator
    /// evaluation showed that it restores minibwa outward recall to the BWA
    /// control range while preserving BSJ and isoform accuracy. `primary` and
    /// `filtered` remain explicit overrides for speed- or precision-first
    /// audits.
    fn non_bsj_outward_retention_mode() -> NonBsjOutwardRetention {
        let value = env::var(NON_BSJ_OUTWARD_RETENTION_ENV).ok();
        Self::parse_non_bsj_outward_retention_mode(value.as_deref())
    }

    /// Parses a configured non-BSJ outward retention mode.
    ///
    /// Keeping parsing separate from the environment read makes the
    /// recognition-first default testable despite the process-wide environment.
    /// Only an unset environment variable opts into that default; unknown
    /// explicit values fall back to the legacy primary-only policy so typoed
    /// speed or precision overrides do not silently enable extra retention.
    fn parse_non_bsj_outward_retention_mode(value: Option<&str>) -> NonBsjOutwardRetention {
        let Some(value) = value else {
            return NonBsjOutwardRetention::Relaxed;
        };
        let value = value.trim();
        if value.eq_ignore_ascii_case("supplementary") {
            NonBsjOutwardRetention::Supplementary
        } else if value.eq_ignore_ascii_case("relaxed") {
            NonBsjOutwardRetention::Relaxed
        } else if value.eq_ignore_ascii_case("filtered") {
            NonBsjOutwardRetention::Filtered
        } else if value.eq_ignore_ascii_case("primary")
            || value.eq_ignore_ascii_case("primary_only")
            || value.eq_ignore_ascii_case("primary-only")
        {
            NonBsjOutwardRetention::PrimaryOnly
        } else {
            NonBsjOutwardRetention::PrimaryOnly
        }
    }

    /// Tests 3' outward pair geometry before writing non-BSJ segment evidence.
    ///
    /// This must stay aligned with the stricter post-Summary outward detector in
    /// `ciri_as`: Scan2 only decides whether the full read group is worth
    /// spilling to the sidecar, while the later segment stage still applies the
    /// linear-mate negative filter and builds only `type=outward` rows. Keeping
    /// the same coordinate rule here prevents gap-facing outward pairs from
    /// being dropped before that safer final check can run.
    fn may_support_outward_segments<'a>(alignments: &[AlignmentRecord<'a>]) -> bool {
        let mut r1: Option<&AlignmentRecord<'_>> = None;
        let mut r2: Option<&AlignmentRecord<'_>> = None;
        for aln in alignments {
            if aln.flag & 0x4 != 0 || aln.flag & 0x100 != 0 || aln.flag & 0x800 != 0 {
                continue;
            }
            match aln.flag & 0x40 != 0 {
                true if r1.is_none() => r1 = Some(aln),
                true => return false,
                false if r2.is_none() => r2 = Some(aln),
                false => return false,
            }
        }
        let (Some(r1), Some(r2)) = (r1, r2) else {
            return false;
        };
        if r1.chrom != r2.chrom
            || r1.chrom.as_ref() == "*"
            || r1.mapq < NON_BSJ_SEGMENT_MAPQ_THRES
            || r2.mapq < NON_BSJ_SEGMENT_MAPQ_THRES
        {
            return false;
        }
        let Some(r1_span) = Self::alignment_ref_span(r1) else {
            return false;
        };
        let Some(r2_span) = Self::alignment_ref_span(r2) else {
            return false;
        };
        if !Self::has_3p_outward_pair_geometry(r1, r1_span, r2, r2_span) {
            return false;
        }
        let span_start = r1_span.0.min(r2_span.0);
        let span_end = r1_span.1.max(r2_span.1);
        span_end - span_start + 1 <= NON_BSJ_SEGMENT_MAX_SPAN
    }

    /// Tests supplementary-aware outward geometry for opt-in sidecar retention.
    ///
    /// This function never runs in the default path. It is intentionally more
    /// permissive about supplementary records because Scan2 only decides whether
    /// to spill the read group; the later segments stage still applies full
    /// linear negative filters, circ containment and chain ranking before
    /// emitting `type=outward`.
    fn may_support_outward_segments_with_retention<'a>(
        alignments: &[AlignmentRecord<'a>],
        mode: NonBsjOutwardRetention,
    ) -> bool {
        if mode == NonBsjOutwardRetention::PrimaryOnly {
            return false;
        }
        for (left_idx, left) in alignments.iter().enumerate() {
            if !Self::is_non_bsj_outward_retention_record(left) {
                continue;
            }
            let Some(left_span) = Self::alignment_ref_span(left) else {
                continue;
            };
            for right in alignments.iter().skip(left_idx + 1) {
                if !Self::is_non_bsj_outward_retention_record(right)
                    || left.chrom != right.chrom
                    || mate_bucket_from_flag(left.flag) == mate_bucket_from_flag(right.flag)
                {
                    continue;
                }
                if !Self::outward_retention_mapq_pass(left, right, mode) {
                    continue;
                }
                let Some(right_span) = Self::alignment_ref_span(right) else {
                    continue;
                };
                if !Self::has_3p_outward_pair_geometry(left, left_span, right, right_span) {
                    continue;
                }
                let span_start = left_span.0.min(right_span.0);
                let span_end = left_span.1.max(right_span.1);
                if span_end - span_start + 1 <= NON_BSJ_SEGMENT_MAX_SPAN {
                    return true;
                }
            }
        }
        false
    }

    /// Returns whether one alignment can participate in opt-in outward retention.
    fn is_non_bsj_outward_retention_record<'a>(aln: &AlignmentRecord<'a>) -> bool {
        aln.flag & 0x4 == 0
            && aln.flag & 0x100 == 0
            && aln.chrom.as_ref() != "*"
            && aln.cigar.as_ref() != "*"
            && mate_bucket_from_flag(aln.flag) != 0
    }

    /// Applies MAPQ thresholds for supplementary-aware outward retention modes.
    fn outward_retention_mapq_pass<'a>(
        left: &AlignmentRecord<'a>,
        right: &AlignmentRecord<'a>,
        mode: NonBsjOutwardRetention,
    ) -> bool {
        match mode {
            NonBsjOutwardRetention::PrimaryOnly => false,
            NonBsjOutwardRetention::Supplementary => {
                left.mapq >= NON_BSJ_SEGMENT_MAPQ_THRES && right.mapq >= NON_BSJ_SEGMENT_MAPQ_THRES
            }
            NonBsjOutwardRetention::Filtered => {
                left.mapq >= NON_BSJ_SEGMENT_MAPQ_THRES && right.mapq >= NON_BSJ_SEGMENT_MAPQ_THRES
            }
            NonBsjOutwardRetention::Relaxed => {
                left.mapq.max(right.mapq) >= NON_BSJ_SEGMENT_MAPQ_THRES
            }
        }
    }

    /// Mirrors the post-Summary outward geometry used by CIRI-AS segments.
    ///
    /// A non-identical pair is accepted when the reverse-strand mate lies left
    /// of the forward-strand mate by both start and end coordinates. Exact
    /// same-span pairs need paired 3' terminal clips so fully overlapping
    /// artifacts are not spilled into the non-BSJ sidecar.
    fn has_3p_outward_pair_geometry<'a>(
        r1: &AlignmentRecord<'a>,
        r1_span: (i32, i32),
        r2: &AlignmentRecord<'a>,
        r2_span: (i32, i32),
    ) -> bool {
        let r1_reverse = r1.flag & 0x10 != 0;
        let r2_reverse = r2.flag & 0x10 != 0;
        if r1_reverse == r2_reverse {
            return false;
        }
        let (reverse_span, forward_span) = if r1_reverse {
            (r1_span, r2_span)
        } else {
            (r2_span, r1_span)
        };
        let (reverse_record, forward_record) = if r1_reverse { (r1, r2) } else { (r2, r1) };
        if reverse_span == forward_span {
            return Self::has_paired_outward_terminal_clip(reverse_record, forward_record);
        }
        forward_span.0 - reverse_span.0 >= NON_BSJ_OUTWARD_MIN_PAIR_OFFSET
            && forward_span.1 - reverse_span.1 >= NON_BSJ_OUTWARD_MIN_PAIR_OFFSET
    }

    /// Returns whether same-span mate pairs carry outward-specific clip tails.
    ///
    /// CIGAR terminal clips are in query order. Subtracting the opposite mate's
    /// 5' clip keeps ordinary mate-overlap clipping from masquerading as paired
    /// 3' outward evidence.
    fn has_paired_outward_terminal_clip<'a>(
        reverse_record: &AlignmentRecord<'a>,
        forward_record: &AlignmentRecord<'a>,
    ) -> bool {
        let Some((reverse_5p, reverse_3p)) =
            Self::terminal_query_clips(reverse_record.cigar.as_ref())
        else {
            return false;
        };
        let Some((forward_5p, forward_3p)) =
            Self::terminal_query_clips(forward_record.cigar.as_ref())
        else {
            return false;
        };
        reverse_3p >= NON_BSJ_OUTWARD_MIN_TERMINAL_CLIP
            && forward_3p >= NON_BSJ_OUTWARD_MIN_TERMINAL_CLIP
            && reverse_3p - forward_5p >= NON_BSJ_OUTWARD_MIN_TERMINAL_CLIP
            && forward_3p - reverse_5p >= NON_BSJ_OUTWARD_MIN_TERMINAL_CLIP
    }

    /// Returns 5' and 3' terminal query clip lengths from a SAM CIGAR.
    fn terminal_query_clips(cigar: &str) -> Option<(i32, i32)> {
        let ops = parse_cigar_ops_basic(cigar)?;
        let five_prime = ops
            .first()
            .copied()
            .filter(|(_, op)| matches!(op, 'S' | 'H'))
            .map_or(0, |(len, _)| len);
        let three_prime = ops
            .last()
            .copied()
            .filter(|(_, op)| matches!(op, 'S' | 'H'))
            .map_or(0, |(len, _)| len);
        Some((five_prime, three_prime))
    }

    /// Tests whether split alignments can form a CIRI-AS-style backward candidate.
    ///
    /// The filter intentionally mirrors the cheap, circ-independent part of
    /// `mapping_check2`: same chromosome and strand are necessary but not enough;
    /// the paired CIGAR classifiers must also satisfy a backward junction
    /// geometry. This keeps Scan2 from writing ordinary supplementary-heavy
    /// groups into the non-BSJ sidecar while still allowing reads with no known
    /// BSJ or circRNA assignment.
    fn may_support_backward_segments<'a>(alignments: &[AlignmentRecord<'a>]) -> bool {
        for i in 0..alignments.len().saturating_sub(1) {
            let left = &alignments[i];
            if left.chrom.as_ref() == "*" || left.cigar.as_ref() == "*" {
                continue;
            }
            for right in &alignments[i + 1..] {
                if left.chrom != right.chrom
                    || right.cigar.as_ref() == "*"
                    || (left.flag & 0x10 != 0) != (right.flag & 0x10 != 0)
                {
                    continue;
                }
                let read_len = left.seq.len().max(right.seq.len()) as i32;
                if read_len <= 0 {
                    continue;
                }
                let left_msid = Self::non_bsj_msid(left.cigar.as_ref(), read_len);
                let right_msid = Self::non_bsj_msid(right.cigar.as_ref(), read_len);
                if Self::non_bsj_backward_pair_possible(
                    left, right, left_msid, right_msid, read_len,
                ) || Self::non_bsj_backward_pair_possible(
                    right, left, right_msid, left_msid, read_len,
                ) {
                    return true;
                }
            }
        }
        false
    }

    /// Tests the CIRI-AS backward candidate geometry without circ-boundary checks.
    fn non_bsj_backward_pair_possible<'a>(
        left: &AlignmentRecord<'a>,
        right: &AlignmentRecord<'a>,
        left_msid: NonBsjMsid,
        right_msid: NonBsjMsid,
        read_len: i32,
    ) -> bool {
        let product = left_msid.kind * right_msid.kind;
        if product == -1 {
            let cir_scale = left_msid.kind * (left.pos + left_msid.clip2)
                + right_msid.kind * (right.pos + right_msid.clip2);
            return (left_msid.clip1 - right_msid.clip1).abs() <= 6 && cir_scale < 0;
        }
        if product.abs() != 10 {
            return false;
        }
        let (mx, my, rx, ry) = if left_msid.kind <= right_msid.kind {
            (left_msid, right_msid, left, right)
        } else {
            (right_msid, left_msid, right, left)
        };
        if mx.kind == -1 {
            let cir_scale = ry.pos + my.ref_len - 1 - rx.pos;
            (read_len - my.clip2 - mx.clip1).abs() <= 6 && cir_scale < 0
        } else {
            let cir_scale = rx.pos + mx.ref_len - 1 - ry.pos;
            (mx.clip1 - my.clip1).abs() <= 6 && cir_scale < 0
        }
    }

    /// Classifies CIGARs using the subset of CIRI-AS `MSID` needed for prefiltering.
    fn non_bsj_msid(cigar: &str, read_len: i32) -> NonBsjMsid {
        if cigar == "*" || cigar.is_empty() {
            return NonBsjMsid::invalid(-1);
        }
        let mut counts = Vec::new();
        let mut styles = Vec::new();
        let mut count = 0i32;
        let mut has_count = false;
        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                let Some(next) = count
                    .checked_mul(10)
                    .and_then(|value| value.checked_add(ch.to_digit(10).unwrap_or(0) as i32))
                else {
                    return NonBsjMsid::invalid(-2);
                };
                count = next;
                has_count = true;
            } else {
                if !has_count {
                    return NonBsjMsid::invalid(-2);
                }
                counts.push(count);
                styles.push(if ch == 'H' { 'S' } else { ch });
                count = 0;
                has_count = false;
            }
        }
        if has_count || styles.is_empty() {
            return NonBsjMsid::invalid(-2);
        }

        if counts.len() == 1 {
            return if styles[0] == 'M' && counts[0] == read_len {
                NonBsjMsid {
                    kind: 0,
                    clip1: 0,
                    clip2: 0,
                    ref_len: read_len,
                }
            } else {
                NonBsjMsid::invalid(-1)
            };
        }

        match counts.len() {
            2 => match (styles[0], styles[1]) {
                ('M', 'S') => NonBsjMsid {
                    kind: 1,
                    clip1: counts[0],
                    clip2: counts[0] - 1,
                    ref_len: counts[0],
                },
                ('S', 'M') => NonBsjMsid {
                    kind: -1,
                    clip1: counts[0],
                    clip2: 0,
                    ref_len: counts[1],
                },
                _ => NonBsjMsid::invalid(-2),
            },
            3 => match (styles[0], styles[1], styles[2]) {
                ('S', _, 'S') => NonBsjMsid {
                    kind: 10,
                    clip1: counts[0],
                    clip2: counts[2],
                    ref_len: counts[1],
                },
                ('M', 'D', 'M') => NonBsjMsid {
                    kind: 0,
                    clip1: 0,
                    clip2: 0,
                    ref_len: read_len + counts[1],
                },
                ('M', 'I', 'M') => NonBsjMsid {
                    kind: 0,
                    clip1: 0,
                    clip2: 0,
                    ref_len: read_len - counts[1],
                },
                _ => NonBsjMsid::invalid(-2),
            },
            _ if styles[0] == 'M' && *styles.last().unwrap() == 'S' => {
                let (m_sum, d_sum) =
                    Self::non_bsj_sum_md(&styles[..styles.len() - 1], &counts[..counts.len() - 1]);
                NonBsjMsid {
                    kind: 1,
                    clip1: read_len - counts[counts.len() - 1],
                    clip2: m_sum + d_sum - 1,
                    ref_len: m_sum + d_sum,
                }
            }
            _ if styles[0] == 'S' && *styles.last().unwrap() == 'M' => {
                let (m_sum, d_sum) = Self::non_bsj_sum_md(&styles[1..], &counts[1..]);
                NonBsjMsid {
                    kind: -1,
                    clip1: counts[0],
                    clip2: 0,
                    ref_len: m_sum + d_sum,
                }
            }
            _ if styles[0] == 'M' && *styles.last().unwrap() == 'M' => {
                let (m_sum, d_sum) = Self::non_bsj_sum_md(&styles, &counts);
                NonBsjMsid {
                    kind: 0,
                    clip1: 0,
                    clip2: 0,
                    ref_len: m_sum + d_sum,
                }
            }
            _ if styles[0] == 'S' && *styles.last().unwrap() == 'S' => {
                let (m_sum, d_sum) = Self::non_bsj_sum_md(
                    &styles[1..styles.len() - 1],
                    &counts[1..counts.len() - 1],
                );
                NonBsjMsid {
                    kind: 10,
                    clip1: counts[0],
                    clip2: counts[counts.len() - 1],
                    ref_len: m_sum + d_sum,
                }
            }
            _ => NonBsjMsid::invalid(-2),
        }
    }

    /// Sums reference match and deletion lengths for the local MSID classifier.
    fn non_bsj_sum_md(styles: &[char], counts: &[i32]) -> (i32, i32) {
        let mut match_sum = 0;
        let mut deletion_sum = 0;
        for (style, count) in styles.iter().zip(counts.iter()) {
            match style {
                'M' | '=' | 'X' => match_sum += count,
                'D' => deletion_sum += count,
                _ => {}
            }
        }
        (match_sum, deletion_sum)
    }

    /// Returns the genomic span covered by reference-consuming CIGAR operators.
    fn alignment_ref_span(aln: &AlignmentRecord<'_>) -> Option<(i32, i32)> {
        let mut ref_len = 0i32;
        let mut count = 0i32;
        let mut has_count = false;
        for op in aln.cigar.as_ref().chars() {
            if op.is_ascii_digit() {
                count = count
                    .checked_mul(10)?
                    .checked_add(op.to_digit(10)? as i32)?;
                has_count = true;
                continue;
            }
            if !has_count {
                return None;
            }
            if matches!(op, 'M' | 'D' | 'N' | '=' | 'X') {
                ref_len = ref_len.checked_add(count)?;
            }
            count = 0;
            has_count = false;
        }
        if has_count || ref_len <= 0 {
            return None;
        }
        Some((aln.pos, aln.pos + ref_len - 1))
    }

    /// Formats mapper rows plus validator-accepted local clip pseudo rows.
    fn segment_evidence_lines_with_local<'a>(
        read_id: &str,
        stage: &str,
        local_stage: &str,
        alignments: &[AlignmentRecord<'a>],
        bsj_lines: &[String],
        reference: &HashMap<String, String>,
    ) -> Vec<String> {
        let mut rows = Self::segment_evidence_lines(read_id, stage, alignments, reference);
        let refs: Vec<&AlignmentRecord<'_>> = alignments.iter().collect();
        rows.extend(local_clip_evidence_lines(
            read_id,
            local_stage,
            &refs,
            bsj_lines,
            reference,
            10,
        ));
        rows
    }

    /// Builds a `(read_id, mate_label, legacy_payload)` key from a priority row.
    ///
    /// The key lets the expanded Scan2 writer avoid emitting the same mate-level
    /// rescue twice when the display path rediscovers the Java first-hit row.
    fn scan2_priority_key(line: &str) -> Option<String> {
        let mut parts = line.splitn(4, '\t');
        let read_id = parts.next()?;
        let mate_label = parts.next()?;
        let _priority = parts.next()?;
        let payload = parts.next()?;
        Some(format!("{read_id}\t{mate_label}\t{payload}"))
    }

    /// Adds a `priority` field to a display row of the form
    /// `read_id, mate_label, legacy_payload...`.
    fn scan2_display_key_and_priority_line(line: &str, priority: &str) -> Option<(String, String)> {
        let mut parts = line.splitn(3, '\t');
        let read_id = parts.next()?;
        let mate_label = parts.next()?;
        let payload = parts.next()?;
        Some((
            format!("{read_id}\t{mate_label}\t{payload}"),
            format!("{read_id}\t{mate_label}\t{priority}\t{payload}"),
        ))
    }

    /// Builds a borrowed display-duplicate key from one priority Scan2 row.
    ///
    /// Display rescue is only allowed to skip candidates that would be filtered
    /// later as exact duplicates of main Scan2 output. The validator-derived tag
    /// fields are deliberately excluded because the pre-validator candidate
    /// identity is enough to prove the display branch is rediscovering the same
    /// first-hit row.
    fn scan2_display_skip_key(line: &str) -> Option<Scan2DisplaySkipKey<'_>> {
        let parts: Vec<&str> = line.split('\t').collect();
        let payload_start = bsj_payload_start(&parts);
        if payload_start != 3 || parts.len() < payload_start + 9 {
            return None;
        }
        Some(Scan2DisplaySkipKey {
            mate: parts[1],
            cigar: parts[payload_start],
            chrom: parts[payload_start + 2],
            site1: parts[payload_start + 3],
            site2: parts[payload_start + 4],
            signal1: parts[payload_start + 5],
            signal2: parts[payload_start + 6],
            sum_q: parts[payload_start + 7],
        })
    }

    /// Returns whether a display candidate is an exact duplicate of a main hit.
    fn display_candidate_is_main_duplicate(
        skip_keys: &[Scan2DisplaySkipKey<'_>],
        mate: &str,
        cigar: &str,
        chrom: &str,
        cand: &CandidateBreakpoint,
    ) -> bool {
        skip_keys
            .iter()
            .any(|key| key.matches(mate, cigar, chrom, cand))
    }

    /// Streams a SAM file into one shard output while batching grouped reads for
    /// parallel Scan2 rescue evaluation.
    fn process_sam_file_to_file(
        &self,
        sam_file: &str,
        chr_tcga_map: &HashMap<String, String>,
        pb: &ProgressBar,
        out_path: &str,
        fsj_path: &str,
        display_scan2: Option<&(Scan2, Scan2DisplayClaims)>,
        display_out_path: Option<&str>,
        segments_out_path: Option<&str>,
        non_bsj_segments_out_path: Option<&str>,
        profile: Option<&Scan2Profile>,
    ) -> Result<()> {
        let mut writer = BufWriter::with_capacity(256 * 1024, File::create(out_path)?);
        let mut display_writer = if let Some(path) = display_out_path {
            Some(BufWriter::with_capacity(256 * 1024, File::create(path)?))
        } else {
            None
        };
        let mut segments_writer = if let Some(path) = segments_out_path {
            Some(BufWriter::with_capacity(256 * 1024, File::create(path)?))
        } else {
            None
        };
        let mut non_bsj_segments_writer = if let Some(path) = non_bsj_segments_out_path {
            Some(BufWriter::with_capacity(256 * 1024, File::create(path)?))
        } else {
            None
        };
        let mut fsj_pair_writer = if fsj_pair_dump_path().is_some() {
            Some(BufWriter::with_capacity(
                256 * 1024,
                File::create(format!("{fsj_path}.pairs"))?,
            ))
        } else {
            None
        };
        let batch_size = (rayon::current_num_threads().max(1) * 256).max(1024);
        let (tx, rx) =
            mpsc::sync_channel::<Vec<SamOwnedScan2Group>>(rayon::current_num_threads().max(2));
        let pb_clone = pb.clone();

        let local_fsj = thread::scope(|scope| -> Result<HashMap<String, i32>> {
            let producer = scope.spawn(|| {
                self.stream_sam_group_batches_for_scan2(
                    sam_file, &pb_clone, tx, batch_size, profile,
                )
            });
            let mut merged_fsj = HashMap::new();
            for mut batch in rx {
                self.process_sam_group_batch_for_scan2(
                    &mut writer,
                    display_writer.as_mut(),
                    segments_writer.as_mut(),
                    non_bsj_segments_writer.as_mut(),
                    fsj_pair_writer.as_mut(),
                    &mut batch,
                    &mut merged_fsj,
                    chr_tcga_map,
                    display_scan2,
                    profile,
                )?;
            }
            producer
                .join()
                .map_err(|_| anyhow::anyhow!("SAM Scan2 group producer thread panicked"))??;
            Ok(merged_fsj)
        })?;

        writer.flush()?;
        if let Some(writer) = display_writer.as_mut() {
            writer.flush()?;
        }
        if let Some(writer) = segments_writer.as_mut() {
            writer.flush()?;
        }
        if let Some(writer) = non_bsj_segments_writer.as_mut() {
            writer.flush()?;
        }
        if let Some(writer) = fsj_pair_writer.as_mut() {
            writer.flush()?;
        }
        Self::write_fsj_shard(fsj_path, &local_fsj)
    }

    /// Sequentially parses SAM records into grouped Scan2 inputs.
    fn stream_sam_group_batches_for_scan2(
        &self,
        sam_file: &str,
        pb: &ProgressBar,
        tx: mpsc::SyncSender<Vec<SamOwnedScan2Group>>,
        batch_size: usize,
        profile: Option<&Scan2Profile>,
    ) -> Result<()> {
        let sam_reader = BufReader::with_capacity(1024 * 1024, File::open(sam_file)?);
        let mut reader = sam::io::Reader::new(sam_reader);
        let header = reader.read_header()?;
        let mut record = sam::Record::default();
        let mut current_id: Vec<u8> = Vec::new();
        let mut alignments: Vec<OwnedAlignmentRecord> = Vec::with_capacity(16);
        let mut all_alignments: Vec<OwnedAlignmentRecord> = Vec::with_capacity(16);
        let mut stand_map: OwnedStandMap = HashMap::with_capacity(4);
        let mut one_read_key: i32 = -1;
        let mut batch = Vec::with_capacity(batch_size);
        let mut records_since_progress = 0usize;
        let mut last_progress_pos = reader.get_mut().stream_position()?;
        let mut cigar_buf = String::with_capacity(64);
        let mut seq_buf = String::with_capacity(256);

        while reader.read_record(&mut record)? != 0 {
            let read_id = record
                .name()
                .ok_or_else(|| anyhow::anyhow!("Missing read name"))?;
            if let Some(profile) = profile {
                profile.records.fetch_add(1, Ordering::Relaxed);
            }
            if read_id != current_id.as_slice() {
                if !current_id.is_empty() {
                    self.push_sam_owned_group_for_scan2(
                        &mut batch,
                        &current_id,
                        &mut alignments,
                        &mut all_alignments,
                        &mut stand_map,
                    );
                    if batch.len() >= batch_size {
                        tx.send(std::mem::take(&mut batch))
                            .map_err(|_| anyhow::anyhow!("SAM Scan2 consumer dropped"))?;
                        batch = Vec::with_capacity(batch_size);
                    }
                }
                current_id.clear();
                current_id.extend_from_slice(read_id);
                one_read_key = -1;
            }

            let flag = i32::from(u16::from(record.flags()?));
            let chrom = match record.reference_sequence(&header) {
                Some(Ok((name, _))) => String::from_utf8_lossy(name).to_string(),
                _ => "*".to_string(),
            };
            let start_pos = record
                .alignment_start()
                .transpose()?
                .map(|p| p.get() as i32)
                .unwrap_or(0);
            let mapq = record
                .mapping_quality()
                .transpose()?
                .map(u8::from)
                .unwrap_or(0) as i32;
            cigar_buf.clear();
            for result in record.cigar().iter() {
                let op = result?;
                use noodles::sam::alignment::record::cigar::op::Kind;
                let op_char = match op.kind() {
                    Kind::Match => 'M',
                    Kind::Insertion => 'I',
                    Kind::Deletion => 'D',
                    Kind::Skip => 'N',
                    Kind::SoftClip => 'S',
                    Kind::HardClip => 'H',
                    Kind::Pad => 'P',
                    Kind::SequenceMatch => '=',
                    Kind::SequenceMismatch => 'X',
                };
                let _ = write!(&mut cigar_buf, "{}{}", op.len(), op_char);
            }
            seq_buf.clear();
            for b in record.sequence().iter() {
                seq_buf.push(char::from(b));
            }
            let seq = seq_buf.clone();
            let s_idx = if flag & 0x40 != 0 { 1 } else { 0 };
            let st_c = if flag & 0x10 != 0 { '1' } else { '0' };
            if s_idx != one_read_key {
                one_read_key = s_idx;
                alignments.retain(|a| {
                    let idx = if a.flag & 0x40 != 0 { 1 } else { 0 };
                    idx != s_idx
                });
                if !seq.is_empty() && seq != "*" {
                    stand_map.insert(s_idx, (st_c, Cow::Owned(seq.clone())));
                }
            }
            let alignment = AlignmentRecord {
                flag,
                chrom: Cow::Owned(chrom),
                pos: start_pos,
                mapq,
                cigar: Cow::Owned(cigar_buf.clone()),
                seq: Cow::Owned(seq),
                xa: Cow::Owned(raw_xa_from_sam_record(&record)?),
            };
            all_alignments.push(alignment.clone());
            alignments.push(alignment);

            records_since_progress += 1;
            if records_since_progress >= 4096 {
                let pos = reader.get_mut().stream_position()?;
                if pos > last_progress_pos {
                    pb.inc(pos - last_progress_pos);
                    last_progress_pos = pos;
                }
                records_since_progress = 0;
            }
        }

        if !current_id.is_empty() {
            self.push_sam_owned_group_for_scan2(
                &mut batch,
                &current_id,
                &mut alignments,
                &mut all_alignments,
                &mut stand_map,
            );
        }
        if !batch.is_empty() {
            tx.send(batch)
                .map_err(|_| anyhow::anyhow!("SAM Scan2 consumer dropped"))?;
        }

        let final_pos = reader.get_mut().stream_position()?;
        if final_pos > last_progress_pos {
            pb.inc(final_pos - last_progress_pos);
        }
        Ok(())
    }

    /// Moves the currently accumulated SAM Scan2 group into the worker batch.
    fn push_sam_owned_group_for_scan2(
        &self,
        batch: &mut Vec<SamOwnedScan2Group>,
        current_id: &[u8],
        alignments: &mut Vec<OwnedAlignmentRecord>,
        all_alignments: &mut Vec<OwnedAlignmentRecord>,
        stand_map: &mut OwnedStandMap,
    ) {
        batch.push(SamOwnedScan2Group {
            read_id: String::from_utf8_lossy(current_id).into_owned(),
            alignments: std::mem::take(alignments),
            all_alignments: std::mem::take(all_alignments),
            stand_map: std::mem::take(stand_map),
        });
        *alignments = Vec::with_capacity(16);
        *all_alignments = Vec::with_capacity(16);
        *stand_map = HashMap::with_capacity(4);
    }

    /// Evaluates one batch of grouped SAM reads in parallel and writes Scan2 output.
    fn process_sam_group_batch_for_scan2(
        &self,
        writer: &mut BufWriter<File>,
        mut display_writer: Option<&mut BufWriter<File>>,
        mut segments_writer: Option<&mut BufWriter<File>>,
        mut non_bsj_segments_writer: Option<&mut BufWriter<File>>,
        mut fsj_pair_writer: Option<&mut BufWriter<File>>,
        batch: &mut Vec<SamOwnedScan2Group>,
        merged_fsj: &mut HashMap<String, i32>,
        chr_tcga_map: &HashMap<String, String>,
        display_scan2: Option<&(Scan2, Scan2DisplayClaims)>,
        profile: Option<&Scan2Profile>,
    ) -> Result<()> {
        let groups = std::mem::take(batch);
        let collect_fsj_pair_rows = fsj_pair_writer.is_some();
        let suppress_display_duplicates = display_writer.is_none();
        let results: Vec<(
            Vec<String>,
            HashMap<String, i32>,
            Vec<String>,
            Vec<String>,
            Vec<String>,
            Vec<String>,
        )> = groups
            .into_par_iter()
            .map(|owned| {
                let mut local_lines = Vec::new();
                let mut local_fsj = HashMap::new();
                let mut fsj_pair_rows = Vec::new();
                let fsj_pair_dump = collect_fsj_pair_rows.then_some(&mut fsj_pair_rows);
                let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
                let _ = self.process_group_view(
                    &owned.read_id,
                    &owned.alignments,
                    &owned.stand_map,
                    &mut local_lines,
                    &mut local_fsj,
                    fsj_pair_dump,
                    chr_tcga_map,
                    &mut validator,
                    profile,
                );
                let display_lines = if let Some((display_helper, scan1_claims)) = display_scan2 {
                    let mut display_validator =
                        IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
                    let display_skip_keys: Vec<_> = if suppress_display_duplicates {
                        local_lines
                            .iter()
                            .filter_map(|line| Self::scan2_display_skip_key(line))
                            .collect()
                    } else {
                        Vec::new()
                    };
                    display_helper
                        .process_group_view_display(
                            &owned.read_id,
                            &owned.alignments,
                            &owned.stand_map,
                            scan1_claims,
                            &display_skip_keys,
                            chr_tcga_map,
                            &mut display_validator,
                        )
                        .unwrap_or_default()
                } else {
                    Vec::new()
                };
                let has_bsj = !local_lines.is_empty() || !display_lines.is_empty();
                let evidence = if has_bsj {
                    let mut bsj_lines = local_lines.clone();
                    bsj_lines.extend(display_lines.iter().cloned());
                    Self::segment_evidence_lines_with_local(
                        &owned.read_id,
                        "scan2",
                        "scan2_local",
                        &owned.alignments,
                        &bsj_lines,
                        chr_tcga_map,
                    )
                } else {
                    Vec::new()
                };
                let scan1_claimed = display_scan2.as_ref().is_some_and(|(_, scan1_claims)| {
                    scan1_claims.by_read.contains_key(owned.read_id.as_str())
                });
                let non_bsj_evidence =
                    if local_lines.is_empty() && display_lines.is_empty() && !scan1_claimed {
                        Self::non_bsj_segment_evidence_lines(
                            &owned.read_id,
                            &owned.all_alignments,
                            chr_tcga_map,
                        )
                    } else {
                        Vec::new()
                    };
                (
                    local_lines,
                    local_fsj,
                    display_lines,
                    evidence,
                    non_bsj_evidence,
                    fsj_pair_rows,
                )
            })
            .collect();

        for (lines, fsj_map, display_lines, evidence, non_bsj_evidence, fsj_pair_rows) in results {
            let main_keys: HashSet<String> = lines
                .iter()
                .filter_map(|line| Self::scan2_priority_key(line))
                .collect();
            for line in lines {
                writeln!(writer, "{}", line)?;
            }
            if let Some(display_writer) = display_writer.as_deref_mut() {
                for line in display_lines {
                    writeln!(display_writer, "{}", line)?;
                }
            } else {
                for line in display_lines {
                    if let Some((key, priority_line)) =
                        Self::scan2_display_key_and_priority_line(&line, "0")
                    {
                        if !main_keys.contains(&key) {
                            writeln!(writer, "{}", priority_line)?;
                        }
                    }
                }
            }
            for (key, count) in fsj_map {
                *merged_fsj.entry(key).or_insert(0) += count;
            }
            if let Some(fsj_pair_writer) = fsj_pair_writer.as_deref_mut() {
                for row in fsj_pair_rows {
                    writeln!(fsj_pair_writer, "{}", row)?;
                }
            }
            if let Some(segments_writer) = segments_writer.as_deref_mut() {
                let mut rows = 0_u64;
                let mut bytes = 0_u64;
                for line in evidence {
                    rows += 1;
                    bytes += line.len() as u64 + 1;
                    writeln!(segments_writer, "{}", line)?;
                }
                if let Some(profile) = profile {
                    profile
                        .segment_sidecar_rows
                        .fetch_add(rows, Ordering::Relaxed);
                    profile
                        .segment_sidecar_bytes
                        .fetch_add(bytes, Ordering::Relaxed);
                }
            }
            if let Some(non_bsj_segments_writer) = non_bsj_segments_writer.as_deref_mut() {
                let mut rows = 0_u64;
                let mut bytes = 0_u64;
                for line in non_bsj_evidence {
                    rows += 1;
                    bytes += line.len() as u64 + 1;
                    writeln!(non_bsj_segments_writer, "{}", line)?;
                }
                if let Some(profile) = profile {
                    profile
                        .non_bsj_sidecar_rows
                        .fetch_add(rows, Ordering::Relaxed);
                    profile
                        .non_bsj_sidecar_bytes
                        .fetch_add(bytes, Ordering::Relaxed);
                }
            }
        }
        Ok(())
    }

    /// Processes one BAM shard and emits rescued BSJ2 lines plus one local FSJ spill file.
    ///
    /// The small buffer reuses inside this function are intentional low-risk
    /// optimizations: they reduce per-record allocation churn without changing the
    /// order or content of any candidate checks.
    fn process_bam_shard_to_file(
        &self,
        mmap: &Mmap,
        start: usize,
        end: usize,
        header: &sam::Header,
        chr_tcga_map: &HashMap<String, String>,
        pb: &ProgressBar,
        out_path: &str,
        fsj_path: &str,
        display_scan2: Option<&(Scan2, Scan2DisplayClaims)>,
        display_out_path: Option<&str>,
        segments_out_path: Option<&str>,
        non_bsj_segments_out_path: Option<&str>,
        profile: Option<&Scan2Profile>,
    ) -> Result<()> {
        use noodles::bam;
        let shard_started = profile.map(|_| Instant::now());
        let mut writer = BufWriter::with_capacity(256 * 1024, File::create(out_path)?);
        let mut display_writer = if let Some(path) = display_out_path {
            Some(BufWriter::with_capacity(256 * 1024, File::create(path)?))
        } else {
            None
        };
        let mut segments_writer = if let Some(path) = segments_out_path {
            Some(BufWriter::with_capacity(256 * 1024, File::create(path)?))
        } else {
            None
        };
        let mut non_bsj_segments_writer = if let Some(path) = non_bsj_segments_out_path {
            Some(BufWriter::with_capacity(256 * 1024, File::create(path)?))
        } else {
            None
        };
        let mut fsj_pair_writer = if fsj_pair_dump_path().is_some() {
            Some(BufWriter::with_capacity(
                256 * 1024,
                File::create(format!("{fsj_path}.pairs"))?,
            ))
        } else {
            None
        };
        let mut local_fsj = HashMap::new();
        let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
        let mut display_validator = display_scan2
            .as_ref()
            .map(|_| IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni));
        let block_start = if start == 0 {
            0
        } else {
            let mut found = None;
            for i in start..mmap.len().saturating_sub(3) {
                if &mmap[i..i + 4] == b"\x1f\x8b\x08\x04" {
                    found = Some(i);
                    break;
                }
            }
            found.unwrap_or(mmap.len())
        };
        let pos = if start == 0 {
            0
        } else {
            // Decode from the previous BGZF block so we can distinguish a true
            // new read group from one that started in the earlier shard.
            let mut found = None;
            let mut i = block_start
                .saturating_sub(1)
                .min(mmap.len().saturating_sub(4));
            loop {
                if i + 3 < mmap.len() && &mmap[i..i + 4] == b"\x1f\x8b\x08\x04" {
                    found = Some(i);
                    break;
                }
                if i == 0 {
                    break;
                }
                i -= 1;
            }
            found.unwrap_or(0)
        };
        if pos >= mmap.len() {
            writer.flush()?;
            if let Some(display_writer) = display_writer.as_mut() {
                display_writer.flush()?;
            }
            if let Some(segments_writer) = segments_writer.as_mut() {
                segments_writer.flush()?;
            }
            if let Some(non_bsj_segments_writer) = non_bsj_segments_writer.as_mut() {
                non_bsj_segments_writer.flush()?;
            }
            if let Some(fsj_pair_writer) = fsj_pair_writer.as_mut() {
                fsj_pair_writer.flush()?;
            }
            return Self::write_fsj_shard(fsj_path, &local_fsj);
        }

        let mut reader = bam::io::Reader::new(&mmap[pos..]);
        if start == 0 {
            let _ = reader.read_header()?;
        }
        let mut record = bam::Record::default();
        let mut current_id: Vec<u8> = Vec::new();
        let mut alignments: Vec<AlignmentRecord> = Vec::with_capacity(16);
        let mut all_alignments: Vec<AlignmentRecord> = Vec::with_capacity(16);
        let mut stand_map: HashMap<i32, (char, Cow<str>)> = HashMap::with_capacity(4);
        let mut one_read_key: i32 = -1;
        let mut crossed_start = start == 0;
        let mut leading_partial_id: Option<Vec<u8>> = None;
        let mut last_compressed_pos = block_start.saturating_sub(pos);
        let mut last_evicted_pos = pos;
        // Keep page-cache eviction on a small sliding window. Large thresholds let
        // the whole mapped BAM drift into RSS before any MADV_DONTNEED call fires.
        let eviction_threshold =
            (self.mem_limit as usize / 8).clamp(8 * 1024 * 1024, 64 * 1024 * 1024);
        let mut cigar_buf = String::with_capacity(64);
        let mut seq_buf = String::with_capacity(256);
        let mut res_batch = Vec::new();
        let mut fsj_pair_rows = Vec::new();

        while reader.read_record(&mut record)? != 0 {
            if let Some(profile) = profile {
                profile.records.fetch_add(1, Ordering::Relaxed);
            }
            let curr_c_pos = reader.get_ref().virtual_position().compressed() as usize;
            pb.inc(curr_c_pos.saturating_sub(last_compressed_pos) as u64);
            last_compressed_pos = curr_c_pos;

            let evicted_rel = last_evicted_pos.saturating_sub(pos);
            if curr_c_pos.saturating_sub(evicted_rel) > eviction_threshold {
                advise_dontneed(
                    mmap,
                    last_evicted_pos,
                    curr_c_pos.saturating_sub(evicted_rel),
                );
                last_evicted_pos = pos + curr_c_pos;
            }

            let read_id = record
                .name()
                .ok_or_else(|| anyhow::anyhow!("Missing read name"))?;
            let abs_c_pos = pos + curr_c_pos;
            if !crossed_start {
                if abs_c_pos < block_start {
                    leading_partial_id = Some(read_id.to_vec());
                    continue;
                }
                crossed_start = true;
            }
            if let Some(partial_id) = &leading_partial_id {
                if read_id == partial_id.as_slice() {
                    // Let the earlier shard own the cross-boundary read group
                    // completely; otherwise Scan2 can silently drop the first
                    // full read group visible in this shard.
                    continue;
                }
                leading_partial_id = None;
            }

            if read_id != current_id.as_slice() {
                if !current_id.is_empty() {
                    let id_str = String::from_utf8_lossy(&current_id);
                    res_batch.clear();
                    fsj_pair_rows.clear();
                    let fsj_pair_dump = if fsj_pair_writer.is_some() {
                        Some(&mut fsj_pair_rows)
                    } else {
                        None
                    };
                    self.process_group_view(
                        &id_str,
                        &alignments,
                        &stand_map,
                        &mut res_batch,
                        &mut local_fsj,
                        fsj_pair_dump,
                        chr_tcga_map,
                        &mut validator,
                        profile,
                    )?;
                    if let Some(fsj_pair_writer) = fsj_pair_writer.as_mut() {
                        for row in &fsj_pair_rows {
                            writeln!(fsj_pair_writer, "{}", row)?;
                        }
                    }
                    let write_started = profile.map(|_| Instant::now());
                    let main_keys: HashSet<String> = res_batch
                        .iter()
                        .filter_map(|line| Self::scan2_priority_key(line))
                        .collect();
                    let display_skip_keys: Vec<_> = if display_writer.is_none() {
                        res_batch
                            .iter()
                            .filter_map(|line| Self::scan2_display_skip_key(line))
                            .collect()
                    } else {
                        Vec::new()
                    };
                    for line in &res_batch {
                        writeln!(writer, "{}", line)?;
                    }
                    if let (Some(profile), Some(write_started)) = (profile, write_started) {
                        profile.write_ns.fetch_add(
                            write_started.elapsed().as_nanos() as u64,
                            Ordering::Relaxed,
                        );
                    }
                    let display_started = profile.map(|_| Instant::now());
                    let display_lines =
                        if let (Some((display_helper, scan1_claims)), Some(display_validator)) =
                            (display_scan2, display_validator.as_mut())
                        {
                            display_helper.process_group_view_display(
                                &id_str,
                                &alignments,
                                &stand_map,
                                scan1_claims,
                                &display_skip_keys,
                                chr_tcga_map,
                                display_validator,
                            )?
                        } else {
                            Vec::new()
                        };
                    if let (Some(profile), Some(display_started)) = (profile, display_started) {
                        profile.display_ns.fetch_add(
                            display_started.elapsed().as_nanos() as u64,
                            Ordering::Relaxed,
                        );
                    }
                    if !display_lines.is_empty() {
                        let write_started = profile.map(|_| Instant::now());
                        if let Some(display_writer) = display_writer.as_mut() {
                            for line in &display_lines {
                                writeln!(display_writer, "{}", line)?;
                            }
                        } else {
                            for line in &display_lines {
                                if let Some((key, priority_line)) =
                                    Self::scan2_display_key_and_priority_line(line, "0")
                                {
                                    if !main_keys.contains(&key) {
                                        writeln!(writer, "{}", priority_line)?;
                                    }
                                }
                            }
                        }
                        if let (Some(profile), Some(write_started)) = (profile, write_started) {
                            profile.write_ns.fetch_add(
                                write_started.elapsed().as_nanos() as u64,
                                Ordering::Relaxed,
                            );
                        }
                    }
                    let scan1_claimed = display_scan2.as_ref().is_some_and(|(_, scan1_claims)| {
                        scan1_claims.by_read.contains_key(id_str.as_ref())
                    });
                    let has_bsj =
                        !res_batch.is_empty() || !display_lines.is_empty() || scan1_claimed;
                    if let Some(segments_writer) = segments_writer.as_mut() {
                        let segment_started = profile.map(|_| Instant::now());
                        if has_bsj {
                            let mut bsj_lines = res_batch.clone();
                            bsj_lines.extend(display_lines.iter().cloned());
                            let mut rows = 0_u64;
                            let mut bytes = 0_u64;
                            for line in Self::segment_evidence_lines_with_local(
                                &id_str,
                                "scan2",
                                "scan2_local",
                                &alignments,
                                &bsj_lines,
                                chr_tcga_map,
                            ) {
                                rows += 1;
                                bytes += line.len() as u64 + 1;
                                writeln!(segments_writer, "{}", line)?;
                            }
                            if let Some(profile) = profile {
                                profile
                                    .segment_sidecar_rows
                                    .fetch_add(rows, Ordering::Relaxed);
                                profile
                                    .segment_sidecar_bytes
                                    .fetch_add(bytes, Ordering::Relaxed);
                            }
                        }
                        if let (Some(profile), Some(segment_started)) = (profile, segment_started) {
                            profile.segment_sidecar_ns.fetch_add(
                                segment_started.elapsed().as_nanos() as u64,
                                Ordering::Relaxed,
                            );
                        }
                    }
                    if !has_bsj {
                        if let Some(non_bsj_segments_writer) = non_bsj_segments_writer.as_mut() {
                            let non_bsj_started = profile.map(|_| Instant::now());
                            let mut rows = 0_u64;
                            let mut bytes = 0_u64;
                            for line in Self::non_bsj_segment_evidence_lines(
                                &id_str,
                                &all_alignments,
                                chr_tcga_map,
                            ) {
                                rows += 1;
                                bytes += line.len() as u64 + 1;
                                writeln!(non_bsj_segments_writer, "{}", line)?;
                            }
                            if let (Some(profile), Some(non_bsj_started)) =
                                (profile, non_bsj_started)
                            {
                                profile.non_bsj_sidecar_ns.fetch_add(
                                    non_bsj_started.elapsed().as_nanos() as u64,
                                    Ordering::Relaxed,
                                );
                                profile
                                    .non_bsj_sidecar_rows
                                    .fetch_add(rows, Ordering::Relaxed);
                                profile
                                    .non_bsj_sidecar_bytes
                                    .fetch_add(bytes, Ordering::Relaxed);
                            }
                        }
                    }
                    if abs_c_pos > end {
                        current_id.clear();
                        break;
                    }
                }
                current_id.clear();
                current_id.extend_from_slice(read_id);
                alignments.clear();
                all_alignments.clear();
                stand_map.clear();
                one_read_key = -1;
            }

            let flag = i32::from(u16::from(record.flags()));
            let chrom = match record.reference_sequence(header) {
                Some(Ok((name, _))) => String::from_utf8_lossy(name).to_string(),
                _ => "*".to_string(),
            };
            let start_pos = record
                .alignment_start()
                .transpose()?
                .map(|p| p.get() as i32)
                .unwrap_or(0);
            let mapq = record.mapping_quality().map(u8::from).unwrap_or(0) as i32;
            cigar_buf.clear();
            for result in record.cigar().iter() {
                let op = result?;
                use noodles::sam::alignment::record::cigar::op::Kind;
                let op_char = match op.kind() {
                    Kind::Match => 'M',
                    Kind::Insertion => 'I',
                    Kind::Deletion => 'D',
                    Kind::Skip => 'N',
                    Kind::SoftClip => 'S',
                    Kind::HardClip => 'H',
                    Kind::Pad => 'P',
                    Kind::SequenceMatch => '=',
                    Kind::SequenceMismatch => 'X',
                };
                let _ = write!(&mut cigar_buf, "{}{}", op.len(), op_char);
            }
            seq_buf.clear();
            for b in record.sequence().iter() {
                seq_buf.push(char::from(b));
            }
            let seq = seq_buf.clone();
            let s_idx = if flag & 0x40 != 0 { 1 } else { 0 };
            let st_c = if flag & 0x10 != 0 { '1' } else { '0' };
            if s_idx != one_read_key {
                one_read_key = s_idx;
                // Java parity: BAM Scan2 rebuilds the per-mate alignment list each
                // time the iterator switches between R1 and R2 within the same
                // read group. Keeping every earlier alignment for that mate looks
                // harmless, but it changes which candidate payloads contribute to
                // `temFSJId` on complex supplementary-heavy reads and was one of
                // the reasons full hg38 FSJ counting drifted while chr1 stayed
                // clean. This overwrite is therefore intentional and verified.
                alignments.retain(|a| {
                    let idx = if a.flag & 0x40 != 0 { 1 } else { 0 };
                    idx != s_idx
                });
                // Java parity: `standMap` is overwritten on mate switches with the
                // current record's sequence; it is not a "longest-sequence wins"
                // cache in Scan2. This only affects Scan2 candidate validation,
                // not Scan1 representative-sequence handling. The distinction is
                // easy to miss because Scan1 does keep the longest representative
                // sequence, but carrying that policy into Scan2 perturbs FSJ-only
                // behavior on supplementary-heavy BAM families.
                if !seq.is_empty() && seq != "*" {
                    stand_map.insert(s_idx, (st_c, Cow::Owned(seq.clone())));
                }
            }
            let alignment = AlignmentRecord {
                flag,
                chrom: Cow::Owned(chrom),
                pos: start_pos,
                mapq,
                cigar: Cow::Owned(cigar_buf.clone()),
                seq: Cow::Owned(seq),
                xa: Cow::Owned(raw_xa_from_bam_record(&record)?),
            };
            all_alignments.push(alignment.clone());
            alignments.push(alignment);
        }
        if !current_id.is_empty() {
            let id_str = String::from_utf8_lossy(&current_id);
            res_batch.clear();
            fsj_pair_rows.clear();
            let fsj_pair_dump = if fsj_pair_writer.is_some() {
                Some(&mut fsj_pair_rows)
            } else {
                None
            };
            self.process_group_view(
                &id_str,
                &alignments,
                &stand_map,
                &mut res_batch,
                &mut local_fsj,
                fsj_pair_dump,
                chr_tcga_map,
                &mut validator,
                profile,
            )?;
            if let Some(fsj_pair_writer) = fsj_pair_writer.as_mut() {
                for row in &fsj_pair_rows {
                    writeln!(fsj_pair_writer, "{}", row)?;
                }
            }
            let write_started = profile.map(|_| Instant::now());
            let main_keys: HashSet<String> = res_batch
                .iter()
                .filter_map(|line| Self::scan2_priority_key(line))
                .collect();
            let display_skip_keys: Vec<_> = if display_writer.is_none() {
                res_batch
                    .iter()
                    .filter_map(|line| Self::scan2_display_skip_key(line))
                    .collect()
            } else {
                Vec::new()
            };
            for line in &res_batch {
                writeln!(writer, "{}", line)?;
            }
            if let (Some(profile), Some(write_started)) = (profile, write_started) {
                profile
                    .write_ns
                    .fetch_add(write_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            }
            let display_started = profile.map(|_| Instant::now());
            let display_lines =
                if let (Some((display_helper, scan1_claims)), Some(display_validator)) =
                    (display_scan2, display_validator.as_mut())
                {
                    display_helper.process_group_view_display(
                        &id_str,
                        &alignments,
                        &stand_map,
                        scan1_claims,
                        &display_skip_keys,
                        chr_tcga_map,
                        display_validator,
                    )?
                } else {
                    Vec::new()
                };
            if let (Some(profile), Some(display_started)) = (profile, display_started) {
                profile.display_ns.fetch_add(
                    display_started.elapsed().as_nanos() as u64,
                    Ordering::Relaxed,
                );
            }
            if !display_lines.is_empty() {
                let write_started = profile.map(|_| Instant::now());
                if let Some(display_writer) = display_writer.as_mut() {
                    for line in &display_lines {
                        writeln!(display_writer, "{}", line)?;
                    }
                } else {
                    for line in &display_lines {
                        if let Some((key, priority_line)) =
                            Self::scan2_display_key_and_priority_line(line, "0")
                        {
                            if !main_keys.contains(&key) {
                                writeln!(writer, "{}", priority_line)?;
                            }
                        }
                    }
                }
                if let (Some(profile), Some(write_started)) = (profile, write_started) {
                    profile
                        .write_ns
                        .fetch_add(write_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
                }
            }
            let scan1_claimed = display_scan2.as_ref().is_some_and(|(_, scan1_claims)| {
                scan1_claims.by_read.contains_key(id_str.as_ref())
            });
            let has_bsj = !res_batch.is_empty() || !display_lines.is_empty() || scan1_claimed;
            if let Some(segments_writer) = segments_writer.as_mut() {
                let segment_started = profile.map(|_| Instant::now());
                if has_bsj {
                    let mut bsj_lines = res_batch.clone();
                    bsj_lines.extend(display_lines.iter().cloned());
                    let mut rows = 0_u64;
                    let mut bytes = 0_u64;
                    for line in Self::segment_evidence_lines_with_local(
                        &id_str,
                        "scan2",
                        "scan2_local",
                        &alignments,
                        &bsj_lines,
                        chr_tcga_map,
                    ) {
                        rows += 1;
                        bytes += line.len() as u64 + 1;
                        writeln!(segments_writer, "{}", line)?;
                    }
                    if let Some(profile) = profile {
                        profile
                            .segment_sidecar_rows
                            .fetch_add(rows, Ordering::Relaxed);
                        profile
                            .segment_sidecar_bytes
                            .fetch_add(bytes, Ordering::Relaxed);
                    }
                }
                if let (Some(profile), Some(segment_started)) = (profile, segment_started) {
                    profile.segment_sidecar_ns.fetch_add(
                        segment_started.elapsed().as_nanos() as u64,
                        Ordering::Relaxed,
                    );
                }
            }
            if !has_bsj {
                if let Some(non_bsj_segments_writer) = non_bsj_segments_writer.as_mut() {
                    let non_bsj_started = profile.map(|_| Instant::now());
                    let mut rows = 0_u64;
                    let mut bytes = 0_u64;
                    for line in
                        Self::non_bsj_segment_evidence_lines(&id_str, &all_alignments, chr_tcga_map)
                    {
                        rows += 1;
                        bytes += line.len() as u64 + 1;
                        writeln!(non_bsj_segments_writer, "{}", line)?;
                    }
                    if let (Some(profile), Some(non_bsj_started)) = (profile, non_bsj_started) {
                        profile.non_bsj_sidecar_ns.fetch_add(
                            non_bsj_started.elapsed().as_nanos() as u64,
                            Ordering::Relaxed,
                        );
                        profile
                            .non_bsj_sidecar_rows
                            .fetch_add(rows, Ordering::Relaxed);
                        profile
                            .non_bsj_sidecar_bytes
                            .fetch_add(bytes, Ordering::Relaxed);
                    }
                }
            }
        }
        let evicted_rel = last_evicted_pos.saturating_sub(pos);
        if last_compressed_pos > evicted_rel {
            advise_dontneed(mmap, last_evicted_pos, last_compressed_pos - evicted_rel);
        }
        writer.flush()?;
        if let Some(display_writer) = display_writer.as_mut() {
            display_writer.flush()?;
        }
        if let Some(segments_writer) = segments_writer.as_mut() {
            segments_writer.flush()?;
        }
        if let Some(non_bsj_segments_writer) = non_bsj_segments_writer.as_mut() {
            non_bsj_segments_writer.flush()?;
        }
        if let Some(fsj_pair_writer) = fsj_pair_writer.as_mut() {
            fsj_pair_writer.flush()?;
        }
        if let (Some(profile), Some(shard_started)) = (profile, shard_started) {
            profile
                .shard_total_ns
                .fetch_add(shard_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
        Self::write_fsj_shard(fsj_path, &local_fsj)
    }

    /// Rescues BSJ candidates and counts overlapping FSJs for one read group.
    ///
    /// This is the Scan2 hot path. The preserved optimizations are intentionally
    /// modest:
    /// - mate-segment lookup is hoisted once per segment side
    /// - oriented read/mate strings are computed once per alignment, not once per
    ///   candidate
    /// - bucket traversal order remains Java-compatible
    pub(crate) fn process_group_view<'a>(
        &self,
        id: &str,
        alignments: &[AlignmentRecord<'a>],
        stand_map: &HashMap<i32, (char, Cow<'a, str>)>,
        results: &mut Vec<String>,
        local_fsj: &mut HashMap<String, i32>,
        mut fsj_pair_dump: Option<&mut Vec<String>>,
        chr_tcga_map: &HashMap<String, String>,
        is_bsj_hg2: &mut IsBSJHg2,
        profile: Option<&Scan2Profile>,
    ) -> Result<()> {
        let started = profile.map(|_| Instant::now());
        let trace_read = should_trace_read(id);
        let trace_all_candidates =
            trace_read && std::env::var("CIRI_TRACE_ALL_CANDS").ok().as_deref() == Some("1");
        if self.scan1_ids.contains(id) {
            return Ok(());
        }
        let mut segments: HashMap<i32, Vec<&AlignmentRecord<'a>>> = HashMap::new();
        for aln in alignments {
            segments
                .entry(if aln.flag & 0x40 != 0 { 1 } else { 0 })
                .or_default()
                .push(aln);
        }
        // Java counts FSJ support per read group, not per alignment record. The
        // temporary set intentionally de-duplicates all linear evidence gathered
        // from one read before the local shard counter is incremented.
        let mut tem_fsj_keys = HashSet::new();
        for seg_idx in [0_i32, 1_i32] {
            let Some(seg_alns) = segments.get(&seg_idx) else {
                continue;
            };
            let (read_strand, seq) = match stand_map.get(&seg_idx) {
                Some(&(st, ref s)) => (st, s.as_ref()),
                None => continue,
            };
            let slen = seq.len() as i32;
            // Mate alignment lookup is invariant for all candidates generated from
            // this read side, so keep the reference once and reuse it below.
            let mate_seg_alns = segments.get(&(1 - seg_idx));
            for aln in seg_alns {
                let chr = aln.chrom.as_ref();
                if !self.index1.contains_key(chr) {
                    continue;
                }
                let cigar_ref: Cow<'_, str> = if aln.cigar.contains('H') {
                    Cow::Owned(aln.cigar.replace('H', "S"))
                } else {
                    Cow::Borrowed(aln.cigar.as_ref())
                };
                if cigar_ref == "*" {
                    continue;
                }
                let c = misd(cigar_ref.as_ref(), slen);
                if cigar_is_full_match(cigar_ref.as_ref(), slen) {
                    let start_tem = aln.pos + 6;
                    let end_tem = aln.pos + slen - 7;
                    self.collect_fsj_keys_for_alignment(
                        chr,
                        start_tem,
                        end_tem,
                        0,
                        aln.pos,
                        cigar_ref.as_ref(),
                        &mut tem_fsj_keys,
                    );
                    continue;
                }
                if c[0] == -1 || c[0] == 10 {
                    if let Some(list) = self.index1.get(chr) {
                        let curr_strand = if aln.flag & 0x10 != 0 { '1' } else { '0' };
                        // Orientation depends on the alignment strand of the
                        // current record, so the normalized strings are hoisted to
                        // alignment scope but not further.
                        let str_e_owned = if aln.flag & 0x10 != 0 {
                            if read_strand == '0' {
                                reverse_complement(seq)
                            } else {
                                seq.to_string()
                            }
                        } else if read_strand == '0' {
                            seq.to_string()
                        } else {
                            reverse_complement(seq)
                        };
                        let p_str_owned =
                            if let Some(&(p_strand, ref p_seq)) = stand_map.get(&(1 - seg_idx)) {
                                if p_strand != curr_strand {
                                    p_seq.to_string()
                                } else {
                                    reverse_complement(p_seq)
                                }
                            } else {
                                String::new()
                            };
                        let new_num_site = aln.pos;
                        // Java parity: fixed seqLen from constructor.
                        let bucket_size = self.seq_len.max(1);
                        let num1 = (new_num_site - 6) / bucket_size;
                        let num2 = (new_num_site + 6) / bucket_size;

                        let mut eval_candidate = |cand: &CandidateBreakpoint| -> Result<bool> {
                            let e_idx = c[1] + (cand.site - new_num_site);
                            if e_idx <= 0 || e_idx > slen {
                                return Ok(false);
                            }
                            let str_f = &str_e_owned[0..e_idx as usize];
                            let str3 = if c[0] == 10 {
                                let si = slen - c[2];
                                if si >= 0 && si <= slen {
                                    &str_e_owned[si as usize..]
                                } else {
                                    ""
                                }
                            } else {
                                "*"
                            };
                            let site1 = cand.data[0].parse::<i32>().unwrap_or(0);
                            let site2 = cand.data[1].parse::<i32>().unwrap_or(0);
                            let mut s2_ok = if mate_seg_alns.is_some() { 0 } else { 1 };
                            if let Some(mate_seg_alns) = mate_seg_alns {
                                for mate_aln in mate_seg_alns {
                                    if mate_aln.chrom.as_ref() == chr
                                        && mate_aln.mapq >= self.min_mapq_uni
                                    {
                                        let c_ano = misd(&mate_aln.cigar, slen);
                                        let mate_strand =
                                            if mate_aln.flag & 0x10 != 0 { '1' } else { '0' };
                                        if mate_strand != curr_strand
                                            && mate_aln.pos >= site1 - 6
                                            && mate_aln.pos + c_ano[3] - 1 <= site2 + 6
                                        {
                                            s2_ok = 1;
                                        } else {
                                            s2_ok = -1;
                                        }
                                        break;
                                    }
                                }
                            }
                            let circ_c = vec![
                                (if aln.flag & 0x10 != 0 { "1" } else { "0" }).to_string(),
                                chr.to_string(),
                                "sm".to_string(),
                                cand.data[0].clone(),
                                cand.data[1].clone(),
                                str_f.to_string(),
                                p_str_owned.clone(),
                                str3.to_string(),
                                s2_ok.to_string(),
                                cand.data[2].clone(),
                                cand.data[3].clone(),
                                cand.data[4].clone(),
                                aln.mapq.to_string(),
                            ];
                            let validator_started = profile.map(|_| Instant::now());
                            let tag = with_trace_hg2_scope(trace_read, || {
                                is_bsj_hg2.is_bsj_hg2(&circ_c, chr_tcga_map.get(chr).unwrap())
                            });
                            if let Some(profile) = profile {
                                profile.candidate_checks.fetch_add(1, Ordering::Relaxed);
                                if let Some(validator_started) = validator_started {
                                    profile.validator_call_ns.fetch_add(
                                        validator_started.elapsed().as_nanos() as u64,
                                        Ordering::Relaxed,
                                    );
                                }
                            }
                            if trace_read {
                                emit_trace_line(&format!(
                                    "[TRACE_SCAN2_CAND] id={} type=sm seg={} aln_pos={} chr={} site1={} site2={} cand_site={} cigar={} mapq={} s2_ok={} str_len={} pair_len={} tag={}",
                                    id,
                                    seg_idx,
                                    aln.pos,
                                    chr,
                                    cand.data[0],
                                    cand.data[1],
                                    cand.site,
                                    cigar_ref.as_ref(),
                                    aln.mapq,
                                    circ_c[8],
                                    circ_c[5].len(),
                                    circ_c[6].len(),
                                    tag
                                ));
                            }
                            if tag == "0" {
                                // Java treats validator tag `0` as "this linear
                                // candidate is not a rescuable BSJ, but it still
                                // supports the circ as an FSJ competitor". These
                                // insertions are one of the main places where
                                // BSJ parity can already be perfect while FSJ
                                // counts still drift.
                                Self::insert_fsj_key_for_alignment(
                                    chr,
                                    &cand.data[0],
                                    &cand.data[1],
                                    aln.pos,
                                    cigar_ref.as_ref(),
                                    &mut tem_fsj_keys,
                                );
                                Ok(false)
                            } else if tag != "2" {
                                if let Some(profile) = profile {
                                    profile.candidate_hits.fetch_add(1, Ordering::Relaxed);
                                }
                                let tag_body = &tag[0..tag.len() - 1];
                                let mate_label = if seg_idx == 1 { "R1" } else { "R2" };
                                results.push(format!(
                                    "{}\t{}\t1\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                    id,
                                    mate_label,
                                    cigar_ref.as_ref(),
                                    tag_body,
                                    chr,
                                    cand.data[0],
                                    cand.data[1],
                                    cand.data[2],
                                    cand.data[3],
                                    cand.data[4],
                                    tag.chars().next_back().unwrap()
                                ));
                                Ok(true)
                            } else {
                                Ok(false)
                            }
                        };

                        // Java parity: num1 bucket reverse order.
                        if num1 > 0 && self.site_array1.get(chr).is_some_and(|s| s.contains(&num1))
                        {
                            let (l1, r1) = Self::bucket_range(list, num1, bucket_size);
                            for idx in (l1..r1).rev() {
                                let cand = &list[idx];
                                let bias = cand.site - new_num_site;
                                if bias >= -6 {
                                    if bias <= 6 && eval_candidate(cand)? && !trace_all_candidates {
                                        return Ok(());
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                        // Java parity: num2 bucket forward order.
                        if num2 != num1
                            && self.site_array1.get(chr).is_some_and(|s| s.contains(&num2))
                        {
                            let (l2, r2) = Self::bucket_range(list, num2, bucket_size);
                            for idx in l2..r2 {
                                let cand = &list[idx];
                                let bias = cand.site - new_num_site;
                                if bias <= 6 {
                                    if eval_candidate(cand)? && !trace_all_candidates {
                                        return Ok(());
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                }
                if c[0] == 1 || c[0] == 10 {
                    let new_site = aln.pos + c[3] - 1;
                    if let Some(list) = self.index2.get(chr) {
                        let curr_strand = if aln.flag & 0x10 != 0 { '1' } else { '0' };
                        let str_e_owned = if aln.flag & 0x10 != 0 {
                            if read_strand == '0' {
                                reverse_complement(seq)
                            } else {
                                seq.to_string()
                            }
                        } else if read_strand == '0' {
                            seq.to_string()
                        } else {
                            reverse_complement(seq)
                        };
                        let p_str_owned =
                            if let Some(&(p_strand, ref p_seq)) = stand_map.get(&(1 - seg_idx)) {
                                if p_strand != curr_strand {
                                    p_seq.to_string()
                                } else {
                                    reverse_complement(p_seq)
                                }
                            } else {
                                String::new()
                            };
                        // Java parity: fixed seqLen from constructor.
                        let bucket_size = self.seq_len.max(1);
                        let num1 = (new_site - 6) / bucket_size;
                        let num2 = (new_site + 6) / bucket_size;

                        let mut eval_candidate = |cand: &CandidateBreakpoint| -> Result<bool> {
                            let s_idx = if c[0] == 10 {
                                // Java parity: SMS right-anchor slicing uses
                                // `seqLen - tail_soft_clip + bias`, not
                                // `start_soft_clip + mapped_len + bias`.
                                // The latter drifts whenever the middle CIGAR
                                // contains insertions or deletions.
                                slen - c[2] + (cand.site - new_site)
                            } else {
                                c[1] + (cand.site - new_site)
                            };
                            if s_idx < 0 || s_idx >= slen {
                                return Ok(false);
                            }
                            let str_f = &str_e_owned[s_idx as usize..];
                            let str3 = if c[0] == 10 {
                                &str_e_owned[0..c[1] as usize]
                            } else {
                                "*"
                            };
                            let site1 = cand.data[0].parse::<i32>().unwrap_or(0);
                            let site2 = cand.data[1].parse::<i32>().unwrap_or(0);
                            let mut s2_ok = if mate_seg_alns.is_some() { 0 } else { 1 };
                            if let Some(mate_seg_alns) = mate_seg_alns {
                                for mate_aln in mate_seg_alns {
                                    if mate_aln.chrom.as_ref() == chr
                                        && mate_aln.mapq >= self.min_mapq_uni
                                    {
                                        let c_ano = misd(&mate_aln.cigar, slen);
                                        let mate_strand =
                                            if mate_aln.flag & 0x10 != 0 { '1' } else { '0' };
                                        if mate_strand != curr_strand
                                            && mate_aln.pos >= site1 - 6
                                            && mate_aln.pos + c_ano[3] - 1 <= site2 + 6
                                        {
                                            s2_ok = 1;
                                        } else {
                                            s2_ok = -1;
                                        }
                                        break;
                                    }
                                }
                            }
                            let circ_c = vec![
                                (if aln.flag & 0x10 != 0 { "1" } else { "0" }).to_string(),
                                chr.to_string(),
                                "ms".to_string(),
                                cand.data[0].clone(),
                                cand.data[1].clone(),
                                str_f.to_string(),
                                p_str_owned.clone(),
                                str3.to_string(),
                                s2_ok.to_string(),
                                cand.data[2].clone(),
                                cand.data[3].clone(),
                                cand.data[4].clone(),
                                aln.mapq.to_string(),
                            ];
                            let validator_started = profile.map(|_| Instant::now());
                            let tag = with_trace_hg2_scope(trace_read, || {
                                is_bsj_hg2.is_bsj_hg2(&circ_c, chr_tcga_map.get(chr).unwrap())
                            });
                            if let Some(profile) = profile {
                                profile.candidate_checks.fetch_add(1, Ordering::Relaxed);
                                if let Some(validator_started) = validator_started {
                                    profile.validator_call_ns.fetch_add(
                                        validator_started.elapsed().as_nanos() as u64,
                                        Ordering::Relaxed,
                                    );
                                }
                            }
                            if trace_read {
                                emit_trace_line(&format!(
                                    "[TRACE_SCAN2_CAND] id={} type=ms seg={} aln_pos={} chr={} site1={} site2={} cand_site={} cigar={} mapq={} s2_ok={} str_len={} pair_len={} tag={}",
                                    id,
                                    seg_idx,
                                    aln.pos,
                                    chr,
                                    cand.data[0],
                                    cand.data[1],
                                    cand.site,
                                    cigar_ref.as_ref(),
                                    aln.mapq,
                                    circ_c[8],
                                    circ_c[5].len(),
                                    circ_c[6].len(),
                                    tag
                                ));
                            }
                            if tag == "0" {
                                Self::insert_fsj_key_for_alignment(
                                    chr,
                                    &cand.data[0],
                                    &cand.data[1],
                                    aln.pos,
                                    cigar_ref.as_ref(),
                                    &mut tem_fsj_keys,
                                );
                                Ok(false)
                            } else if tag != "2" {
                                if let Some(profile) = profile {
                                    profile.candidate_hits.fetch_add(1, Ordering::Relaxed);
                                }
                                let tag_body = &tag[0..tag.len() - 1];
                                let mate_label = if seg_idx == 1 { "R1" } else { "R2" };
                                results.push(format!(
                                    "{}\t{}\t1\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                    id,
                                    mate_label,
                                    cigar_ref.as_ref(),
                                    tag_body,
                                    chr,
                                    cand.data[0],
                                    cand.data[1],
                                    cand.data[2],
                                    cand.data[3],
                                    cand.data[4],
                                    tag.chars().next_back().unwrap()
                                ));
                                Ok(true)
                            } else {
                                Ok(false)
                            }
                        };

                        // Java parity: num1 bucket reverse order.
                        if self.site_array2.get(chr).is_some_and(|s| s.contains(&num1)) {
                            let (l1, r1) = Self::bucket_range(list, num1, bucket_size);
                            for idx in (l1..r1).rev() {
                                let cand = &list[idx];
                                let bias = cand.site - new_site;
                                if bias >= -6 {
                                    if bias <= 6 && eval_candidate(cand)? && !trace_all_candidates {
                                        return Ok(());
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                        // Java parity: num2 bucket forward order.
                        if num2 != num1
                            && self.site_array2.get(chr).is_some_and(|s| s.contains(&num2))
                        {
                            let (l2, r2) = Self::bucket_range(list, num2, bucket_size);
                            for idx in l2..r2 {
                                let cand = &list[idx];
                                let bias = cand.site - new_site;
                                if bias <= 6 {
                                    if eval_candidate(cand)? && !trace_all_candidates {
                                        return Ok(());
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                }
                // Besides explicit tag==0 additions above, every linear segment
                // also contributes bucket-gated FSJ overlaps from its covered
                // genomic span, exactly like Java `GetFSJClass.getFSJ(...)`.
                let start_tem = aln.pos + 6;
                let end_tem = aln.pos + c[3] - 7;
                self.collect_fsj_keys_for_alignment(
                    chr,
                    start_tem,
                    end_tem,
                    c[0],
                    aln.pos,
                    cigar_ref.as_ref(),
                    &mut tem_fsj_keys,
                );
            }
        }
        for count_key in collapsed_fsj_count_keys(tem_fsj_keys) {
            if let Some(dump) = fsj_pair_dump.as_deref_mut() {
                dump.push(format!("{count_key}\t{id}"));
            }
            *local_fsj.entry(count_key).or_insert(0) += 1;
        }
        if let (Some(profile), Some(started)) = (profile, started) {
            profile
                .group_process_ns
                .fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile.groups.fetch_add(1, Ordering::Relaxed);
        }
        Ok(())
    }

    /// Enumerates mate-level Scan2 rescue hits for the post-Summary `.bsj`.
    ///
    /// This helper keeps the expanded mate view out of the parity path. It uses
    /// the display-only claim set from Scan1 and returns at most one first-hit
    /// rescue per mate.
    fn process_group_view_display<'a>(
        &self,
        id: &str,
        alignments: &[AlignmentRecord<'a>],
        stand_map: &HashMap<i32, (char, Cow<'a, str>)>,
        scan1_claims: &Scan2DisplayClaims,
        main_skip_keys: &[Scan2DisplaySkipKey<'_>],
        chr_tcga_map: &HashMap<String, String>,
        is_bsj_hg2: &mut IsBSJHg2,
    ) -> Result<Vec<String>> {
        let trace_read = should_trace_read(id);
        let trace_all_candidates =
            trace_read && std::env::var("CIRI_TRACE_ALL_CANDS").ok().as_deref() == Some("1");
        if scan1_claims.contains_both_mates(id) {
            return Ok(Vec::new());
        }
        let mut results = Vec::new();
        let mut segments: HashMap<i32, Vec<&AlignmentRecord<'a>>> = HashMap::new();
        for aln in alignments {
            segments
                .entry(if aln.flag & 0x40 != 0 { 1 } else { 0 })
                .or_default()
                .push(aln);
        }

        'mate: for seg_idx in [1_i32, 0_i32] {
            let mate_label = if seg_idx == 1 { "R1" } else { "R2" };
            if scan1_claims.contains_mate(id, mate_label) {
                continue;
            }
            let Some(seg_alns) = segments.get(&seg_idx) else {
                continue;
            };
            let (read_strand, seq) = match stand_map.get(&seg_idx) {
                Some(&(st, ref s)) => (st, s.as_ref()),
                None => continue,
            };
            let slen = seq.len() as i32;
            let mate_seg_alns = segments.get(&(1 - seg_idx));
            for aln in seg_alns {
                let chr = aln.chrom.as_ref();
                if !self.index1.contains_key(chr) {
                    continue;
                }
                let cigar_ref: Cow<'_, str> = if aln.cigar.contains('H') {
                    Cow::Owned(aln.cigar.replace('H', "S"))
                } else {
                    Cow::Borrowed(aln.cigar.as_ref())
                };
                if cigar_ref == "*" {
                    continue;
                }
                let c = misd(cigar_ref.as_ref(), slen);
                if c[0] == -1 || c[0] == 10 {
                    if let Some(list) = self.index1.get(chr) {
                        let curr_strand = if aln.flag & 0x10 != 0 { '1' } else { '0' };
                        let str_e_owned = if aln.flag & 0x10 != 0 {
                            if read_strand == '0' {
                                reverse_complement(seq)
                            } else {
                                seq.to_string()
                            }
                        } else if read_strand == '0' {
                            seq.to_string()
                        } else {
                            reverse_complement(seq)
                        };
                        let p_str_owned =
                            if let Some(&(p_strand, ref p_seq)) = stand_map.get(&(1 - seg_idx)) {
                                if p_strand != curr_strand {
                                    p_seq.to_string()
                                } else {
                                    reverse_complement(p_seq)
                                }
                            } else {
                                String::new()
                            };
                        let new_num_site = aln.pos;
                        let bucket_size = self.seq_len.max(1);
                        let num1 = (new_num_site - 6) / bucket_size;
                        let num2 = (new_num_site + 6) / bucket_size;

                        let mut eval_candidate = |cand: &CandidateBreakpoint| -> Result<bool> {
                            let e_idx = c[1] + (cand.site - new_num_site);
                            if e_idx <= 0 || e_idx > slen {
                                return Ok(false);
                            }
                            let str_f = &str_e_owned[0..e_idx as usize];
                            let str3 = if c[0] == 10 {
                                let si = slen - c[2];
                                if si >= 0 && si <= slen {
                                    &str_e_owned[si as usize..]
                                } else {
                                    ""
                                }
                            } else {
                                "*"
                            };
                            let site1 = cand.data[0].parse::<i32>().unwrap_or(0);
                            let site2 = cand.data[1].parse::<i32>().unwrap_or(0);
                            let mut s2_ok = if mate_seg_alns.is_some() { 0 } else { 1 };
                            if let Some(mate_seg_alns) = mate_seg_alns {
                                for mate_aln in mate_seg_alns {
                                    if mate_aln.chrom.as_ref() == chr
                                        && mate_aln.mapq >= self.min_mapq_uni
                                    {
                                        let c_ano = misd(&mate_aln.cigar, slen);
                                        let mate_strand =
                                            if mate_aln.flag & 0x10 != 0 { '1' } else { '0' };
                                        if mate_strand != curr_strand
                                            && mate_aln.pos >= site1 - 6
                                            && mate_aln.pos + c_ano[3] - 1 <= site2 + 6
                                        {
                                            s2_ok = 1;
                                        } else {
                                            s2_ok = -1;
                                        }
                                        break;
                                    }
                                }
                            }
                            if Self::display_candidate_is_main_duplicate(
                                main_skip_keys,
                                mate_label,
                                cigar_ref.as_ref(),
                                chr,
                                cand,
                            ) {
                                if trace_read {
                                    emit_trace_line(&format!(
                                        "[TRACE_SCAN2_DISPLAY] id={} type=sm seg={} aln_pos={} chr={} site1={} site2={} cand_site={} cigar={} skipped=main_duplicate",
                                        id,
                                        seg_idx,
                                        aln.pos,
                                        chr,
                                        cand.data[0],
                                        cand.data[1],
                                        cand.site,
                                        cigar_ref.as_ref()
                                    ));
                                }
                                return Ok(true);
                            }
                            let circ_c = vec![
                                (if aln.flag & 0x10 != 0 { "1" } else { "0" }).to_string(),
                                chr.to_string(),
                                "sm".to_string(),
                                cand.data[0].clone(),
                                cand.data[1].clone(),
                                str_f.to_string(),
                                p_str_owned.clone(),
                                str3.to_string(),
                                s2_ok.to_string(),
                                cand.data[2].clone(),
                                cand.data[3].clone(),
                                cand.data[4].clone(),
                                aln.mapq.to_string(),
                            ];
                            let tag = with_trace_hg2_scope(trace_read, || {
                                is_bsj_hg2.is_bsj_hg2(&circ_c, chr_tcga_map.get(chr).unwrap())
                            });
                            if trace_read {
                                emit_trace_line(&format!(
                                    "[TRACE_SCAN2_DISPLAY] id={} type=sm seg={} aln_pos={} chr={} site1={} site2={} cand_site={} cigar={} mapq={} tag={}",
                                    id,
                                    seg_idx,
                                    aln.pos,
                                    chr,
                                    cand.data[0],
                                    cand.data[1],
                                    cand.site,
                                    cigar_ref.as_ref(),
                                    aln.mapq,
                                    tag
                                ));
                            }
                            if tag != "0" && tag != "2" {
                                let tag_body = &tag[0..tag.len() - 1];
                                results.push(format!(
                                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                    id,
                                    mate_label,
                                    cigar_ref.as_ref(),
                                    tag_body,
                                    chr,
                                    cand.data[0],
                                    cand.data[1],
                                    cand.data[2],
                                    cand.data[3],
                                    cand.data[4],
                                    tag.chars().next_back().unwrap()
                                ));
                                Ok(true)
                            } else {
                                Ok(false)
                            }
                        };

                        if num1 > 0 && self.site_array1.get(chr).is_some_and(|s| s.contains(&num1))
                        {
                            let (l1, r1) = Self::bucket_range(list, num1, bucket_size);
                            for idx in (l1..r1).rev() {
                                let cand = &list[idx];
                                let bias = cand.site - new_num_site;
                                if bias >= -6 {
                                    if bias <= 6 && eval_candidate(cand)? && !trace_all_candidates {
                                        continue 'mate;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                        if num2 != num1
                            && self.site_array1.get(chr).is_some_and(|s| s.contains(&num2))
                        {
                            let (l2, r2) = Self::bucket_range(list, num2, bucket_size);
                            for idx in l2..r2 {
                                let cand = &list[idx];
                                let bias = cand.site - new_num_site;
                                if bias <= 6 {
                                    if eval_candidate(cand)? && !trace_all_candidates {
                                        continue 'mate;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                }
                if c[0] == 1 || c[0] == 10 {
                    let new_site = aln.pos + c[3] - 1;
                    if let Some(list) = self.index2.get(chr) {
                        let curr_strand = if aln.flag & 0x10 != 0 { '1' } else { '0' };
                        let str_e_owned = if aln.flag & 0x10 != 0 {
                            if read_strand == '0' {
                                reverse_complement(seq)
                            } else {
                                seq.to_string()
                            }
                        } else if read_strand == '0' {
                            seq.to_string()
                        } else {
                            reverse_complement(seq)
                        };
                        let p_str_owned =
                            if let Some(&(p_strand, ref p_seq)) = stand_map.get(&(1 - seg_idx)) {
                                if p_strand != curr_strand {
                                    p_seq.to_string()
                                } else {
                                    reverse_complement(p_seq)
                                }
                            } else {
                                String::new()
                            };
                        let bucket_size = self.seq_len.max(1);
                        let num1 = (new_site - 6) / bucket_size;
                        let num2 = (new_site + 6) / bucket_size;

                        let mut eval_candidate = |cand: &CandidateBreakpoint| -> Result<bool> {
                            let s_idx = if c[0] == 10 {
                                slen - c[2] + (cand.site - new_site)
                            } else {
                                c[1] + (cand.site - new_site)
                            };
                            if s_idx < 0 || s_idx >= slen {
                                return Ok(false);
                            }
                            let str_f = &str_e_owned[s_idx as usize..];
                            let str3 = if c[0] == 10 {
                                &str_e_owned[0..c[1] as usize]
                            } else {
                                "*"
                            };
                            let site1 = cand.data[0].parse::<i32>().unwrap_or(0);
                            let site2 = cand.data[1].parse::<i32>().unwrap_or(0);
                            let mut s2_ok = if mate_seg_alns.is_some() { 0 } else { 1 };
                            if let Some(mate_seg_alns) = mate_seg_alns {
                                for mate_aln in mate_seg_alns {
                                    if mate_aln.chrom.as_ref() == chr
                                        && mate_aln.mapq >= self.min_mapq_uni
                                    {
                                        let c_ano = misd(&mate_aln.cigar, slen);
                                        let mate_strand =
                                            if mate_aln.flag & 0x10 != 0 { '1' } else { '0' };
                                        if mate_strand != curr_strand
                                            && mate_aln.pos >= site1 - 6
                                            && mate_aln.pos + c_ano[3] - 1 <= site2 + 6
                                        {
                                            s2_ok = 1;
                                        } else {
                                            s2_ok = -1;
                                        }
                                        break;
                                    }
                                }
                            }
                            if Self::display_candidate_is_main_duplicate(
                                main_skip_keys,
                                mate_label,
                                cigar_ref.as_ref(),
                                chr,
                                cand,
                            ) {
                                if trace_read {
                                    emit_trace_line(&format!(
                                        "[TRACE_SCAN2_DISPLAY] id={} type=ms seg={} aln_pos={} chr={} site1={} site2={} cand_site={} cigar={} skipped=main_duplicate",
                                        id,
                                        seg_idx,
                                        aln.pos,
                                        chr,
                                        cand.data[0],
                                        cand.data[1],
                                        cand.site,
                                        cigar_ref.as_ref()
                                    ));
                                }
                                return Ok(true);
                            }
                            let circ_c = vec![
                                (if aln.flag & 0x10 != 0 { "1" } else { "0" }).to_string(),
                                chr.to_string(),
                                "ms".to_string(),
                                cand.data[0].clone(),
                                cand.data[1].clone(),
                                str_f.to_string(),
                                p_str_owned.clone(),
                                str3.to_string(),
                                s2_ok.to_string(),
                                cand.data[2].clone(),
                                cand.data[3].clone(),
                                cand.data[4].clone(),
                                aln.mapq.to_string(),
                            ];
                            let tag = with_trace_hg2_scope(trace_read, || {
                                is_bsj_hg2.is_bsj_hg2(&circ_c, chr_tcga_map.get(chr).unwrap())
                            });
                            if trace_read {
                                emit_trace_line(&format!(
                                    "[TRACE_SCAN2_DISPLAY] id={} type=ms seg={} aln_pos={} chr={} site1={} site2={} cand_site={} cigar={} mapq={} tag={}",
                                    id,
                                    seg_idx,
                                    aln.pos,
                                    chr,
                                    cand.data[0],
                                    cand.data[1],
                                    cand.site,
                                    cigar_ref.as_ref(),
                                    aln.mapq,
                                    tag
                                ));
                            }
                            if tag != "0" && tag != "2" {
                                let tag_body = &tag[0..tag.len() - 1];
                                results.push(format!(
                                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                    id,
                                    mate_label,
                                    cigar_ref.as_ref(),
                                    tag_body,
                                    chr,
                                    cand.data[0],
                                    cand.data[1],
                                    cand.data[2],
                                    cand.data[3],
                                    cand.data[4],
                                    tag.chars().next_back().unwrap()
                                ));
                                Ok(true)
                            } else {
                                Ok(false)
                            }
                        };

                        if self.site_array2.get(chr).is_some_and(|s| s.contains(&num1)) {
                            let (l1, r1) = Self::bucket_range(list, num1, bucket_size);
                            for idx in (l1..r1).rev() {
                                let cand = &list[idx];
                                let bias = cand.site - new_site;
                                if bias >= -6 {
                                    if bias <= 6 && eval_candidate(cand)? && !trace_all_candidates {
                                        continue 'mate;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                        if num2 != num1
                            && self.site_array2.get(chr).is_some_and(|s| s.contains(&num2))
                        {
                            let (l2, r2) = Self::bucket_range(list, num2, bucket_size);
                            for idx in l2..r2 {
                                let cand = &list[idx];
                                let bias = cand.site - new_site;
                                if bias <= 6 {
                                    if eval_candidate(cand)? && !trace_all_candidates {
                                        continue 'mate;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        Ok(results)
    }
}
