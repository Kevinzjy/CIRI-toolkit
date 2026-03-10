/// Port of CIRI3 Misd.java
/// Parses CIGAR strings using efficient manual parsing instead of Regex.

pub fn misd(old_cigar: &str, seq_length: i32) -> [i32; 4] {
    let mut cigar_ite = [0, 0, 0, 0];
    
    // Standardize: Replace H with S
    let mut counts = Vec::with_capacity(4);
    let mut ops = Vec::with_capacity(4);
    
    let mut current_val = 0;
    for c in old_cigar.chars() {
        if c.is_digit(10) {
            current_val = current_val * 10 + c.to_digit(10).unwrap() as i32;
        } else {
            let op = if c == 'H' { 'S' } else { c };
            counts.push(current_val);
            ops.push(op);
            current_val = 0;
        }
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
