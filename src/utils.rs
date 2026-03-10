//! Utility module: Common bioinformatics and string processing functions.
//!
//! This module provides the `AlignmentRecord` abstraction and helper functions
//! like reverse complementation.

use std::borrow::Cow;

/// Unified representation of a sequence alignment record.
///
/// This struct uses `Cow` (Copy-On-Write) to support zero-copy views for SAM (from Mmap)
/// and owned strings for BAM (from decompression) without code duplication.
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
/// # Arguments
/// * `seq` - The input DNA sequence string.
///
/// # Returns
/// A new `String` representing the reverse complement.
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATGC"), "GCAT");
        assert_eq!(reverse_complement("AAttGGcc"), "ggCCaaTT");
        assert_eq!(reverse_complement("N"), "N");
    }
}
