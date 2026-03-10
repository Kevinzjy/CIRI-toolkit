//! Port of CIRI3 Misd.java
//!
//! This module parses CIGAR strings using an efficient manual state machine
//! instead of regex to determine the alignment type (SM, MS, SMS) and mapped length.

/// Parses a CIGAR string and returns classification metadata.
///
/// # Arguments
/// * `old_cigar` - The CIGAR string from a SAM/BAM record.
/// * `seq_length` - The total length of the read sequence.
///
/// # Returns
/// An array of 4 integers:
/// 0: Type indicator (-1 for SM, 1 for MS, 10 for SMS, 0 for pure M).
/// 1: Soft-clip length at the start (for SM/SMS) or match length (for MS).
/// 2: Relative position index.
/// 3: Total mapped length (sum of M and D operations).
pub fn misd(old_cigar: &str, seq_length: i32) -> [i32; 4] {
    let mut cigar_ite = [0, 0, 0, 0];
    
    // Standardize: Replace H (Hard clip) with S (Soft clip) conceptually for length.
    let mut counts = Vec::with_capacity(4);
    let mut ops = Vec::with_capacity(4);
    
    let mut current_val = 0;
    for c in old_cigar.chars() {
        if c.is_digit(10) {
            current_val = current_val * 10 + c.to_digit(10).unwrap() as i32;
        } else {
            // Only accept valid CIGAR operations
            if matches!(c, 'M' | 'I' | 'D' | 'N' | 'S' | 'H' | 'P' | '=' | 'X') {
                let op = if c == 'H' { 'S' } else { c };
                counts.push(current_val);
                ops.push(op);
                current_val = 0;
            } else {
                // Invalid character in CIGAR string (like '*')
                cigar_ite[3] = -2;
                return cigar_ite;
            }
        }
    }

    // Safety check for empty CIGAR strings
    if ops.is_empty() {
        cigar_ite[3] = -2;
        return cigar_ite;
    }

    if ops.len() == 1 {
        if ops[0] == 'M' { cigar_ite[3] = seq_length; }
        else { cigar_ite[3] = -1; }
    } else if ops.len() == 2 {
        if ops[0] == 'M' && ops[1] == 'S' {
            cigar_ite[0] = 1; cigar_ite[1] = counts[0]; cigar_ite[2] = counts[0] - 1; cigar_ite[3] = counts[0];
        } else if ops[0] == 'S' && ops[1] == 'M' {
            cigar_ite[0] = -1; cigar_ite[1] = counts[0]; cigar_ite[2] = 0; cigar_ite[3] = counts[1];
        } else { cigar_ite[3] = -2; }
    } else if ops.len() == 3 {
        if ops[0] == 'S' && ops[1] == 'M' && ops[2] == 'S' {
            cigar_ite[0] = 10; cigar_ite[1] = counts[0]; cigar_ite[2] = counts[2]; cigar_ite[3] = counts[1];
        } else if ops[0] == 'M' && ops[1] == 'D' && ops[2] == 'M' {
            cigar_ite[3] = seq_length + counts[1];
        } else if ops[0] == 'M' && ops[1] == 'I' && ops[2] == 'M' {
            cigar_ite[3] = seq_length - counts[1];
        } else { cigar_ite[3] = -2; }
    } else if ops[0] == 'M' && *ops.last().unwrap() == 'S' {
        let (mut m_sum, mut d_sum) = (0, 0);
        for i in 0..ops.len() {
            if ops[i] == 'M' { m_sum += counts[i]; }
            else if ops[i] == 'D' { d_sum += counts[i]; }
        }
        cigar_ite[0] = 1; cigar_ite[1] = seq_length - counts.last().unwrap();
        cigar_ite[2] = m_sum + d_sum - 1; cigar_ite[3] = m_sum + d_sum;
    } else if ops[0] == 'S' && *ops.last().unwrap() == 'M' {
        let (mut m_sum, mut d_sum) = (0, 0);
        for i in 1..ops.len() {
            if ops[i] == 'M' { m_sum += counts[i]; }
            else if ops[i] == 'D' { d_sum += counts[i]; }
        }
        cigar_ite[0] = -1; cigar_ite[1] = counts[0]; cigar_ite[2] = 0; cigar_ite[3] = m_sum + d_sum;
    } else if ops[0] == 'M' && *ops.last().unwrap() == 'M' {
        let (mut m_sum, mut d_sum) = (0, 0);
        for i in 0..ops.len() {
            if ops[i] == 'M' { m_sum += counts[i]; }
            else if ops[i] == 'D' { d_sum += counts[i]; }
        }
        cigar_ite[3] = m_sum + d_sum;
    } else if ops[0] == 'S' && *ops.last().unwrap() == 'S' {
        let (mut m_sum, mut d_sum) = (0, 0);
        for i in 1..ops.len()-1 {
            if ops[i] == 'M' { m_sum += counts[i]; }
            else if ops[i] == 'D' { d_sum += counts[i]; }
        }
        cigar_ite[0] = 10; cigar_ite[1] = counts[0]; cigar_ite[2] = *counts.last().unwrap(); cigar_ite[3] = m_sum + d_sum;
    } else { cigar_ite[3] = -2; }
    
    cigar_ite
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_misd_match() {
        let res = misd("100M", 100);
        assert_eq!(res, [0, 0, 0, 100]);
    }

    #[test]
    fn test_misd_invalid() {
        let res = misd("*", 100);
        assert_eq!(res[3], -2);
    }

    #[test]
    fn test_misd_soft_clip_start() {
        let res = misd("20S80M", 100);
        assert_eq!(res, [-1, 20, 0, 80]);
    }

    #[test]
    fn test_misd_soft_clip_end() {
        let res = misd("80M20S", 100);
        assert_eq!(res, [1, 80, 79, 80]);
    }

    #[test]
    fn test_misd_insertion() {
        let res = misd("40M2I58M", 100);
        assert_eq!(res, [0, 0, 0, 98]);
    }

    #[test]
    fn test_misd_deletion() {
        let res = misd("40M2D58M", 98);
        assert_eq!(res, [0, 0, 0, 100]);
    }

    #[test]
    fn test_misd_complex() {
        let res = misd("10S30M2I40M20S", 102);
        assert_eq!(res, [10, 10, 20, 70]);
    }
}
