//! Utility module: Common bioinformatics and string processing functions.
//!
//! This module provides the `AlignmentRecord` abstraction and helper functions
//! like reverse complementation and memory unit parsing.

use std::borrow::Cow;

/// Unified representation of a sequence alignment record.
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
pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C',
            'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c',
            _ => c,
        })
        .collect()
}

/// Parses memory strings like "2G", "512M", "1024K" into bytes.
pub fn parse_mem_str(mem_str: &str) -> u64 {
    let s = mem_str.to_uppercase();
    if s.ends_with('G') {
        s[..s.len()-1].parse::<u64>().unwrap_or(2) * 1024 * 1024 * 1024
    } else if s.ends_with('M') {
        s[..s.len()-1].parse::<u64>().unwrap_or(2048) * 1024 * 1024
    } else if s.ends_with('K') {
        s[..s.len()-1].parse::<u64>().unwrap_or(2097152) * 1024
    } else {
        s.parse::<u64>().unwrap_or(2) * 1024 * 1024 * 1024 // Default 2G
    }
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
}
