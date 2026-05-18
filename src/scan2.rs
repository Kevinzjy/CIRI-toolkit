//! Scan2: second-pass rescue, candidate validation, and FSJ counting.
//!
//! This stage consumes the Java-compatible BSJ output from Scan1, rebuilds
//! the de-duplicated candidate indexes expected by Java CIRI3, and then revisits
//! the input alignments to rescue additional support while counting FSJ evidence.

use crate::is_bsj_hg2::{report_scan2_hg_profile, IsBSJHg2};
use crate::misd::misd;
use crate::runtime::{
    emit_perf_line, emit_trace_line, scan2_profile_enabled, should_trace_read, with_trace_hg2_scope,
};
use crate::utils::{
    bam_shard_count, bsj_is_summary_priority, bsj_payload_start, clip_sequence_payload,
    local_clip_evidence_lines, part_path, reverse_complement, AlignmentRecord,
};
use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use memmap2::Mmap;
use noodles::sam::{
    self,
    alignment::{record::Sequence as _, Record as _},
};
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Seek, Write};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc;
use std::thread;
use std::time::Duration;
use std::time::Instant;

const NON_BSJ_SEGMENT_MAPQ_THRES: i32 = 5;
const NON_BSJ_SEGMENT_MAX_SPAN: i32 = 200000;

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
        }
        writer.flush()?;
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
    ) -> Vec<String> {
        let mut rows = Vec::new();
        for aln in alignments {
            if aln.chrom.as_ref() == "*" || aln.cigar.as_ref() == "*" {
                continue;
            }
            let mate = if aln.flag & 0x40 != 0 { "R1" } else { "R2" };
            let clips = clip_sequence_payload(aln.cigar.as_ref(), aln.seq.as_ref());
            rows.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                read_id,
                stage,
                mate,
                aln.flag,
                aln.chrom,
                aln.pos,
                aln.mapq,
                aln.cigar,
                aln.seq.len(),
                clips
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
    /// The short `N` stage and six-field alignment payload keep the large
    /// whole-genome sidecar cheaper to write and parse.
    fn non_bsj_segment_evidence_lines<'a>(
        read_id: &str,
        alignments: &[AlignmentRecord<'a>],
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
            records.push(format!(
                "{}|{}|{}|{}|{}|{}",
                aln.flag, aln.chrom, aln.pos, aln.mapq, aln.cigar, clips
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
            || Self::may_support_backward_segments(alignments)
    }

    /// Tests the 5' overlap plus 3' outward pair geometry without circ gating.
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
        reverse_span.0 < forward_span.0
            && forward_span.0 <= reverse_span.1
            && reverse_span.1 < forward_span.1
            && forward_span.1 - reverse_span.0 + 1 <= NON_BSJ_SEGMENT_MAX_SPAN
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
        let mut rows = Self::segment_evidence_lines(read_id, stage, alignments);
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
        batch: &mut Vec<SamOwnedScan2Group>,
        merged_fsj: &mut HashMap<String, i32>,
        chr_tcga_map: &HashMap<String, String>,
        display_scan2: Option<&(Scan2, Scan2DisplayClaims)>,
        profile: Option<&Scan2Profile>,
    ) -> Result<()> {
        let groups = std::mem::take(batch);
        let results: Vec<(
            Vec<String>,
            HashMap<String, i32>,
            Vec<String>,
            Vec<String>,
            Vec<String>,
        )> = groups
            .into_par_iter()
            .map(|owned| {
                let mut local_lines = Vec::new();
                let mut local_fsj = HashMap::new();
                let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
                let _ = self.process_group_view(
                    &owned.read_id,
                    &owned.alignments,
                    &owned.stand_map,
                    &mut local_lines,
                    &mut local_fsj,
                    chr_tcga_map,
                    &mut validator,
                    profile,
                );
                let display_lines = if let Some((display_helper, scan1_claims)) = display_scan2 {
                    let mut display_validator =
                        IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
                    display_helper
                        .process_group_view_display(
                            &owned.read_id,
                            &owned.alignments,
                            &owned.stand_map,
                            scan1_claims,
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
                        Self::non_bsj_segment_evidence_lines(&owned.read_id, &owned.all_alignments)
                    } else {
                        Vec::new()
                    };
                (
                    local_lines,
                    local_fsj,
                    display_lines,
                    evidence,
                    non_bsj_evidence,
                )
            })
            .collect();

        for (lines, fsj_map, display_lines, evidence, non_bsj_evidence) in results {
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
                    self.process_group_view(
                        &id_str,
                        &alignments,
                        &stand_map,
                        &mut res_batch,
                        &mut local_fsj,
                        chr_tcga_map,
                        &mut validator,
                        profile,
                    )?;
                    let write_started = profile.map(|_| Instant::now());
                    let main_keys: HashSet<String> = res_batch
                        .iter()
                        .filter_map(|line| Self::scan2_priority_key(line))
                        .collect();
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
                            for line in
                                Self::non_bsj_segment_evidence_lines(&id_str, &all_alignments)
                            {
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
            };
            all_alignments.push(alignment.clone());
            alignments.push(alignment);
        }
        if !current_id.is_empty() {
            let id_str = String::from_utf8_lossy(&current_id);
            res_batch.clear();
            self.process_group_view(
                &id_str,
                &alignments,
                &stand_map,
                &mut res_batch,
                &mut local_fsj,
                chr_tcga_map,
                &mut validator,
                profile,
            )?;
            let write_started = profile.map(|_| Instant::now());
            let main_keys: HashSet<String> = res_batch
                .iter()
                .filter_map(|line| Self::scan2_priority_key(line))
                .collect();
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
                    for line in Self::non_bsj_segment_evidence_lines(&id_str, &all_alignments) {
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
                if cigar_ref == format!("{}M", slen) {
                    let start_tem = aln.pos + 6;
                    let end_tem = aln.pos + slen - 7;
                    self.collect_fsj_keys_in_range(chr, start_tem, end_tem, 0, &mut tem_fsj_keys);
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
                                tem_fsj_keys
                                    .insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
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
                                tem_fsj_keys
                                    .insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
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
                self.collect_fsj_keys_in_range(chr, start_tem, end_tem, c[0], &mut tem_fsj_keys);
            }
        }
        for k in tem_fsj_keys {
            *local_fsj.entry(k).or_insert(0) += 1;
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
