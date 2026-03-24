//! Utility module: Common bioinformatics and string processing functions.
//!
//! This module provides the `AlignmentRecord` abstraction and helper functions
//! like reverse complementation and memory unit parsing.

use std::borrow::Cow;

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
}
