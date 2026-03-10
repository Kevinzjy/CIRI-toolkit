use std::borrow::Cow;

#[derive(Debug, Clone)]
pub struct AlignmentRecord<'a> {
    pub flag: i32,
    pub chrom: Cow<'a, str>,
    pub pos: i32,
    pub mapq: i32,
    pub cigar: Cow<'a, str>,
    pub seq: Cow<'a, str>,
}

pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            _ => c, // 保持原样，如 'N'
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATCG"), "CGAT");
        assert_eq!(reverse_complement("atcg"), "CGAT");
        assert_eq!(reverse_complement("ATCNG"), "CNGAT");
    }
}
