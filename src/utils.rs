//! Utility module: Common bioinformatics and string processing functions.
//!
//! This module provides the `AlignmentRecord` abstraction and helper functions
//! like reverse complementation and memory unit parsing.

use std::borrow::Cow;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ClipPlacement {
    pub(crate) pos: i32,
    pub(crate) clip_offset: usize,
    pub(crate) len: usize,
}

/// Smallest compressed BAM shard size worth parallelizing.
///
/// BGZF streams are block-compressed, and overly small shards make the
/// "search forward to the next block header" fallback ambiguous on tiny BAMs.
/// That can land non-zero shards inside incomplete compressed members and yield
/// decoder errors such as `failed to fill whole buffer`. Capping the shard count
/// by a conservative compressed-byte minimum keeps tiny regression BAMs stable
/// while leaving large production BAMs fully parallel.
pub const MIN_BAM_SHARD_BYTES: usize = 4 * 1024 * 1024;

/// Unified representation of one alignment record used by Scan1 and Scan2.
///
/// `Cow` keeps the type flexible: BAM/SAM parsers can borrow transient slices or
/// promote data to owned strings when groups must outlive the parser buffer.
#[derive(Debug, Clone)]
pub struct AlignmentRecord<'a> {
    pub flag: i32,
    pub chrom: Cow<'a, str>,
    pub pos: i32,
    pub mapq: i32,
    pub cigar: Cow<'a, str>,
    pub seq: Cow<'a, str>,
}

/// Returns the reverse complement of a DNA sequence.
///
/// This helper preserves non-ACGT characters as-is because the reference and read
/// sequences may contain `N`, and parity code expects those positions to survive.
pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'a' => 't',
            't' => 'a',
            'c' => 'g',
            'g' => 'c',
            _ => c,
        })
        .collect()
}

/// Parses a SAM-like CIGAR into `(length, op)` pairs for sidecar formatting.
///
/// This helper intentionally performs only lightweight syntax parsing and keeps
/// operation semantics with the caller. It is used by Scan1/Scan2 sidecar
/// writers, where pulling in CIRI-AS chain reconstruction would create an
/// unnecessary dependency from the parity-sensitive scan stages.
pub fn parse_cigar_ops_basic(cigar: &str) -> Option<Vec<(i32, char)>> {
    let mut ops = Vec::new();
    let mut number = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            number.push(ch);
            continue;
        }
        let len = number.parse::<i32>().ok()?;
        number.clear();
        ops.push((len, ch));
    }
    if number.is_empty() {
        Some(ops)
    } else {
        None
    }
}

/// Encodes only recoverable soft-clipped bases for segment-evidence sidecars.
///
/// Storing full read sequences would inflate the BSJ sidecar substantially. The
/// final segments stage only needs clipped subsequences, so the compact
/// `L:<seq>,R:<seq>` payload keeps I/O bounded while retaining the local evidence
/// validated by Scan1/Scan2.
pub fn clip_sequence_payload(cigar: &str, seq: &str) -> String {
    if seq.is_empty() || seq == "*" || !cigar.contains('S') {
        return "*".to_string();
    }
    let Some(ops) = parse_cigar_ops_basic(cigar) else {
        return "*".to_string();
    };
    let mut fields = Vec::new();
    if let Some((len, 'S')) = ops.first().copied() {
        if let Some(part) = seq.get(..len as usize) {
            fields.push(format!("L:{}", part));
        }
    }
    if let Some((len, 'S')) = ops.last().copied() {
        if let Some(part) = seq.get(seq.len().saturating_sub(len as usize)..) {
            fields.push(format!("R:{}", part));
        }
    }
    if fields.is_empty() {
        "*".to_string()
    } else {
        fields.join(",")
    }
}

/// Adds pseudo-alignment rows for local clip evidence accepted with BSJ calls.
///
/// This helper is deliberately sidecar-only: Scan1/Scan2 first make their normal
/// CIRI3-compatible BSJ decision, then use the confirmed circ span to place
/// clipped subsequences inside that circ interval. Full-clip exact matches are
/// preferred; if the full clip does not place, the longest prefix/suffix partial
/// exact match is retained when it still satisfies the minimum length. The
/// resulting rows feed post-Summary segment assembly but never Summary itself.
pub fn local_clip_evidence_lines<'a>(
    read_id: &str,
    stage: &str,
    alignments: &[&AlignmentRecord<'a>],
    bsj_lines: &[String],
    reference: &HashMap<String, String>,
    min_clip_len: usize,
) -> Vec<String> {
    let mut circ_spans = Vec::new();
    let mut seen_circ = HashSet::new();
    for line in bsj_lines {
        if let Some((chr, start, end)) = parse_bsj_line_span(line) {
            if seen_circ.insert(format!("{chr}\t{start}\t{end}")) {
                circ_spans.push((chr, start, end));
            }
        }
    }
    if circ_spans.is_empty() {
        return Vec::new();
    }

    let mut rows = Vec::new();
    let mut seen_rows = HashSet::new();
    for (chr, circ_start, circ_end) in circ_spans {
        let Some(chr_seq) = reference.get(&chr) else {
            continue;
        };
        if circ_start <= 0 || circ_end < circ_start || circ_end as usize > chr_seq.len() {
            continue;
        }
        let circ_seq = chr_seq[(circ_start - 1) as usize..circ_end as usize].to_uppercase();
        for aln in alignments {
            if aln.chrom.as_ref() != chr || aln.seq.is_empty() || aln.seq.as_ref() == "*" {
                continue;
            }
            let payload = clip_sequence_payload(aln.cigar.as_ref(), aln.seq.as_ref());
            if payload == "*" {
                continue;
            }
            for (side, clip_seq) in parse_clip_payload(&payload) {
                if clip_seq.len() < min_clip_len || clip_seq.contains('N') {
                    continue;
                }
                // SAM/BAM stores the sequence for this alignment record in the
                // same orientation used by its CIGAR. Reversing clips again for
                // `0x10` records loses valid circ-side local blocks from reverse
                // strand BSJ reads.
                let genomic_query = clip_seq.to_uppercase();
                if genomic_query.contains('N') {
                    continue;
                }
                let mut candidates = clip_match_placements(
                    &circ_seq,
                    &genomic_query,
                    circ_start,
                    circ_end,
                    min_clip_len,
                );
                candidates.truncate(4);
                for placement in candidates {
                    let read_len = aln.seq.len() as i32;
                    let Some(cigar) = clip_placement_cigar(
                        side,
                        clip_seq.len(),
                        placement.clip_offset,
                        placement.len,
                        read_len,
                    ) else {
                        continue;
                    };
                    let pos = placement.pos;
                    let flag = aln.flag | 0x800;
                    let key = format!("{flag}\t{chr}\t{pos}\t{cigar}");
                    if !seen_rows.insert(key) {
                        continue;
                    }
                    let mate = if flag & 0x40 != 0 { "R1" } else { "R2" };
                    rows.push(format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t*",
                        read_id,
                        stage,
                        mate,
                        flag,
                        chr,
                        pos,
                        aln.mapq.saturating_sub(1),
                        cigar,
                        read_len
                    ));
                }
            }
        }
    }
    rows
}

/// Builds a pseudo-CIGAR for one full or partial local clip placement.
///
/// `clip_offset` is zero-based inside the clipped subsequence. Left clips occupy
/// read positions starting at one; right clips occupy the terminal read suffix.
/// Preserving both surrounding soft clips lets the post-Summary chain builder
/// know where a partial local block sits within the mate read.
pub(crate) fn clip_placement_cigar(
    side: char,
    clip_len: usize,
    clip_offset: usize,
    match_len: usize,
    read_len: i32,
) -> Option<String> {
    if match_len == 0 || clip_offset + match_len > clip_len || read_len <= 0 {
        return None;
    }
    let read_len = read_len as usize;
    let read_start = match side {
        'L' => clip_offset + 1,
        'R' => read_len.checked_sub(clip_len)? + clip_offset + 1,
        _ => return None,
    };
    let read_end = read_start + match_len - 1;
    if read_end > read_len {
        return None;
    }
    Some(spliced_single_match_cigar(
        read_start - 1,
        match_len,
        read_len - read_end,
    ))
}

/// Formats one local pseudo-alignment with optional terminal soft clips.
fn spliced_single_match_cigar(leading_s: usize, match_len: usize, trailing_s: usize) -> String {
    let mut cigar = String::new();
    if leading_s > 0 {
        cigar.push_str(&format!("{leading_s}S"));
    }
    cigar.push_str(&format!("{match_len}M"));
    if trailing_s > 0 {
        cigar.push_str(&format!("{trailing_s}S"));
    }
    cigar
}

/// Extracts the circ span from old, display, or priority-annotated BSJ rows.
fn parse_bsj_line_span(line: &str) -> Option<(String, i32, i32)> {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.is_empty() {
        return None;
    }
    let payload_start = if parts.len() >= 3
        && is_bsj_mate_label(parts[1])
        && (parts[2] == "0" || parts[2] == "1")
    {
        3
    } else if parts.len() >= 2 && is_bsj_mate_label(parts[1]) {
        2
    } else {
        1
    };
    if parts.len() < payload_start + 5 {
        return None;
    }
    let chr = parts[payload_start + 2].to_string();
    let start = parts[payload_start + 3].parse::<i32>().ok()?;
    let end = parts[payload_start + 4].parse::<i32>().ok()?;
    Some((chr, start, end))
}

/// Parses compact sidecar clip payload fields.
pub(crate) fn parse_clip_payload(payload: &str) -> Vec<(char, &str)> {
    payload
        .split(',')
        .filter_map(|field| {
            let (side, seq) = field.split_once(':')?;
            let side = side.chars().next()?;
            matches!(side, 'L' | 'R').then_some((side, seq))
        })
        .collect()
}

/// Finds full or longest prefix/suffix partial clip placements.
pub(crate) fn clip_match_placements(
    circ_seq: &str,
    genomic_query: &str,
    circ_start: i32,
    circ_end: i32,
    min_clip_len: usize,
) -> Vec<ClipPlacement> {
    let query_len = genomic_query.len();
    let full: Vec<ClipPlacement> =
        exact_clip_match_positions(circ_seq, genomic_query, circ_start, circ_end)
            .into_iter()
            .map(|pos| ClipPlacement {
                pos,
                clip_offset: 0,
                len: query_len,
            })
            .collect();
    if !full.is_empty() {
        return full;
    }
    if query_len <= min_clip_len {
        return Vec::new();
    }

    for len in (min_clip_len..query_len).rev() {
        let mut placements = Vec::new();
        let mut seen = HashSet::new();
        for clip_offset in [0, query_len - len] {
            if !seen.insert(clip_offset) {
                continue;
            }
            let sub_query = &genomic_query[clip_offset..clip_offset + len];
            placements.extend(
                exact_clip_match_positions(circ_seq, sub_query, circ_start, circ_end)
                    .into_iter()
                    .map(|pos| ClipPlacement {
                        pos,
                        clip_offset,
                        len,
                    }),
            );
        }
        if !placements.is_empty() {
            placements.sort_by_key(|placement| {
                let end = placement.pos + placement.len as i32 - 1;
                (
                    (placement.pos - circ_start)
                        .abs()
                        .min((end - circ_end).abs()),
                    placement.pos,
                    placement.clip_offset,
                )
            });
            return placements;
        }
    }
    Vec::new()
}

/// Finds exact local-clip placements and ranks boundary-near hits first.
pub(crate) fn exact_clip_match_positions(
    circ_seq: &str,
    genomic_query: &str,
    circ_start: i32,
    circ_end: i32,
) -> Vec<i32> {
    let mut positions = Vec::new();
    let mut offset = 0usize;
    while let Some(found) = circ_seq[offset..].find(genomic_query) {
        let absolute = offset + found;
        let pos = circ_start + absolute as i32;
        positions.push(pos);
        offset = absolute + 1;
    }
    let len = genomic_query.len() as i32;
    positions.sort_by_key(|pos| {
        let end = *pos + len - 1;
        ((*pos - circ_start).abs().min((end - circ_end).abs()), *pos)
    });
    positions
}

/// Parses memory strings like `2G`, `512M`, or `1024K` into bytes.
///
/// The fallback defaults intentionally match the historical CLI behavior rather
/// than failing hard on malformed input.
pub fn parse_mem_str(mem_str: &str) -> u64 {
    let s = mem_str.to_uppercase();
    if s.ends_with('G') {
        s[..s.len() - 1].parse::<u64>().unwrap_or(2) * 1024 * 1024 * 1024
    } else if s.ends_with('M') {
        s[..s.len() - 1].parse::<u64>().unwrap_or(2048) * 1024 * 1024
    } else if s.ends_with('K') {
        s[..s.len() - 1].parse::<u64>().unwrap_or(2097152) * 1024
    } else {
        s.parse::<u64>().unwrap_or(2) * 1024 * 1024 * 1024 // Default 2G
    }
}

/// Returns the final Summary output path for one CLI `-o/--out` prefix.
///
/// The CLI now treats `-o` strictly as a prefix so all outputs follow one
/// predictable scheme: `<prefix>.out`, `<prefix>.bsj1`, `<prefix>.bsj`,
/// and `<prefix>.fsj`.
///
/// Keeping the naming centralized here avoids another round of drift between
/// CLI help, logging, temporary-file cleanup, and test fixtures.
pub fn result_path_for_output(output_arg: &str) -> String {
    format!("{}.out", output_arg)
}

/// Returns the sidecar log path for one CLI prefix.
///
/// This file stores the high-level stage summaries emitted by the main pipeline
/// so long runs can be reviewed after the terminal session ends.
pub fn log_path_for_output(output_arg: &str) -> String {
    format!("{}.log", output_arg)
}

/// Returns the targeted debug trace path for one CLI prefix.
///
/// `--debug` writes verbose per-read processing traces here so they do not get
/// mixed into the always-on stage summary log.
pub fn debug_path_for_output(output_arg: &str) -> String {
    format!("{}.debug.log", output_arg)
}

/// Returns the profiling report path for one CLI prefix.
///
/// `--perf` is intentionally a boolean switch, so profiling output always lands
/// in a predictable sibling file without forcing the user to invent another
/// path on the command line.
pub fn perf_path_for_output(output_arg: &str) -> String {
    format!("{}.perf.log", output_arg)
}

/// Returns the merged Scan1 BSJ output path for one CLI prefix.
pub fn bsj1_path_for_output(output_arg: &str) -> String {
    format!("{}.bsj1", output_arg)
}

/// Returns the Scan1 segments-evidence sidecar path for one CLI prefix.
///
/// This file is intentionally separate from `.bsj1`: Summary consumes `.bsj1`
/// for CIRI3 parity, while the segments sidecar keeps read-level alignment
/// blocks for post-Summary internal-structure reconstruction.
pub fn segments1_path_for_output(output_arg: &str) -> String {
    format!("{}.segments1", output_arg)
}

/// Returns the final combined BSJ output path for one CLI prefix.
pub fn bsj_path_for_output(output_arg: &str) -> String {
    format!("{}.bsj", output_arg)
}

/// Returns whether a field is the mate label used by the expanded BSJ protocol.
///
/// The parser helpers below deliberately accept both the historical CIRI3-style
/// BSJ rows and the new mate-level rows. That compatibility keeps `Summary` and
/// `Scan2` focused on their parity logic instead of scattering column-offset
/// checks across hot code paths.
pub fn is_bsj_mate_label(value: &str) -> bool {
    value == "R1" || value == "R2"
}

/// Returns the first field index of the legacy BSJ payload.
///
/// Old rows start with `read_id, cigar, tag...`; new rows start with
/// `read_id, mate_label, priority, cigar, tag...`.
pub fn bsj_payload_start(parts: &[&str]) -> usize {
    if parts.len() >= 3 && is_bsj_mate_label(parts[1]) && (parts[2] == "0" || parts[2] == "1") {
        3
    } else {
        1
    }
}

/// Returns whether a BSJ row should participate in the Java-compatible summary.
///
/// Historical rows have no explicit priority and are therefore treated as
/// summary-priority rows. New mate-level rows must use `priority=1`.
pub fn bsj_is_summary_priority(parts: &[&str]) -> bool {
    let payload_start = bsj_payload_start(parts);
    payload_start == 1 || parts.get(2).is_some_and(|priority| *priority == "1")
}

/// Returns the basename used for Scan2 shard-local BSJ spill files.
pub fn bsj2_path_for_output(output_arg: &str) -> String {
    format!("{}.bsj2", output_arg)
}

/// Returns the Scan2 segments-evidence sidecar path for one CLI prefix.
///
/// This mirrors [`segments1_path_for_output`] for Scan2 rescue reads and lets
/// the final `<prefix>.segments` stage avoid re-scanning BAM/SAM for confirmed
/// BSJ read-level blocks.
pub fn segments2_path_for_output(output_arg: &str) -> String {
    format!("{}.segments2", output_arg)
}

/// Returns the basename used for FSJ shard-local spill files.
///
/// The pipeline no longer keeps a final `.fsj` artifact, but the `<prefix>.fsj`
/// stem is still used to name temporary spill files such as
/// `<prefix>.fsj.part_0001.tmp`.
pub fn fsj_path_for_output(output_arg: &str) -> String {
    format!("{}.fsj", output_arg)
}

/// Returns the temporary spill-file path for one shard-local intermediate.
///
/// Stable zero-padded numbering keeps shard files easy to scan by eye and
/// avoids mixing two historical naming schemes (`.shard_N` vs `fsj_shard_N`).
pub fn part_path(final_path: &str, part_idx: usize) -> String {
    format!("{}.part_{:04}.tmp", final_path, part_idx + 1)
}

/// Chooses a safe BAM shard count for BGZF-parallel Scan1/Scan2 processing.
///
/// This intentionally falls back to fewer shards for tiny BAMs. The goal is not
/// throughput on tiny fixtures, but avoiding shard starts that are too dense to
/// reliably find a distinct next BGZF member.
pub fn bam_shard_count(file_len: usize, requested_threads: usize) -> usize {
    let requested = requested_threads.max(1);
    let by_size = (file_len / MIN_BAM_SHARD_BYTES).max(1);
    requested.min(by_size.max(1))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATGC"), "GCAT");
        assert_eq!(reverse_complement("AAttGGcc"), "ggCCaaTT");
        assert_eq!(reverse_complement("N"), "N");
    }

    #[test]
    fn test_parse_mem_str() {
        assert_eq!(parse_mem_str("1G"), 1024 * 1024 * 1024);
        assert_eq!(parse_mem_str("512M"), 512 * 1024 * 1024);
        assert_eq!(parse_mem_str("2"), 2 * 1024 * 1024 * 1024);
    }

    #[test]
    fn test_bam_shard_count() {
        assert_eq!(bam_shard_count(32 * 1024, 16), 1);
        assert_eq!(bam_shard_count(4 * 1024 * 1024, 16), 1);
        assert_eq!(bam_shard_count(8 * 1024 * 1024, 16), 2);
        assert_eq!(bam_shard_count(128 * 1024 * 1024, 4), 4);
    }

    #[test]
    fn reverse_strand_clip_matches_record_sequence_orientation() {
        let mut chr = vec![b'N'; 300];
        chr[179..195].copy_from_slice(b"TTTTCTAACCTTGTGA");
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), String::from_utf8(chr).unwrap());
        let seq = format!("{}{}", "TTTTCTAACCTTGTGA", "A".repeat(134));
        let alignment = AlignmentRecord {
            flag: 0x40 | 0x10,
            chrom: Cow::Borrowed("chr1"),
            pos: 100,
            mapq: 60,
            cigar: Cow::Borrowed("16S134M"),
            seq: Cow::Owned(seq),
        };
        let bsj_lines = vec!["read1\tx\ty\tchr1\t100\t199".to_string()];

        let rows = local_clip_evidence_lines(
            "read1",
            "scan1_local",
            &[&alignment],
            &bsj_lines,
            &reference,
            10,
        );

        assert_eq!(
            rows,
            vec!["read1\tscan1_local\tR1\t2128\tchr1\t180\t59\t16M134S\t150\t*"]
        );
    }

    #[test]
    fn local_clip_records_longest_prefix_partial_match() {
        let mut chr = vec![b'N'; 300];
        chr[179..194].copy_from_slice(b"AGCCATCTGTGAGGG");
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), String::from_utf8(chr).unwrap());
        let seq = format!("{}{}", "AGCCATCTGTGAGGGA", "A".repeat(134));
        let alignment = AlignmentRecord {
            flag: 0x40,
            chrom: Cow::Borrowed("chr1"),
            pos: 195,
            mapq: 60,
            cigar: Cow::Borrowed("16S134M"),
            seq: Cow::Owned(seq),
        };
        let bsj_lines = vec!["read1\tx\ty\tchr1\t100\t250".to_string()];

        let rows = local_clip_evidence_lines(
            "read1",
            "scan1_local",
            &[&alignment],
            &bsj_lines,
            &reference,
            10,
        );

        assert_eq!(
            rows,
            vec!["read1\tscan1_local\tR1\t2112\tchr1\t180\t59\t15M135S\t150\t*"]
        );
    }

    #[test]
    fn local_clip_records_longest_suffix_partial_match() {
        let mut chr = vec![b'N'; 300];
        chr[199..214].copy_from_slice(b"CACATGTGGACTAAA");
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), String::from_utf8(chr).unwrap());
        let seq = format!("{}{}", "TTTTTTTTTTTTCACATGTGGACTAAA", "A".repeat(123));
        let alignment = AlignmentRecord {
            flag: 0x80,
            chrom: Cow::Borrowed("chr1"),
            pos: 215,
            mapq: 60,
            cigar: Cow::Borrowed("27S123M"),
            seq: Cow::Owned(seq),
        };
        let bsj_lines = vec!["read1\tx\ty\tchr1\t100\t250".to_string()];

        let rows = local_clip_evidence_lines(
            "read1",
            "scan1_local",
            &[&alignment],
            &bsj_lines,
            &reference,
            10,
        );

        assert_eq!(
            rows,
            vec!["read1\tscan1_local\tR2\t2176\tchr1\t200\t59\t12S15M123S\t150\t*"]
        );
    }

    #[test]
    fn test_output_paths_for_prefix() {
        assert_eq!(result_path_for_output("sample.ciri"), "sample.ciri.out");
        assert_eq!(log_path_for_output("sample.ciri"), "sample.ciri.log");
        assert_eq!(bsj1_path_for_output("sample.ciri"), "sample.ciri.bsj1");
        assert_eq!(bsj_path_for_output("sample.ciri"), "sample.ciri.bsj");
        assert_eq!(bsj2_path_for_output("sample.ciri"), "sample.ciri.bsj2");
        assert_eq!(fsj_path_for_output("sample.ciri"), "sample.ciri.fsj");
        assert_eq!(
            part_path("sample.ciri.bsj", 0),
            "sample.ciri.bsj.part_0001.tmp"
        );
    }
}
