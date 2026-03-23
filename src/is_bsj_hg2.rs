//! Validator module: Reproduces the core validation logic of Java CIRI3.
//!
//! This module includes comprehensive sequence validation, canonical splice signal 
//! identification, and linear competition checks to distinguish BSJs from linear splicing noise.

use crate::index_compare::IndexCompare;
use std::collections::HashMap;
use std::fmt::Write as FmtWrite;

/// Smith-Waterman local alignment implementation with traceback.
pub struct SmithWaterman {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_penalty: i32,
    pub score: i32,
    /// Number of steps in the traceback, equivalent to `alignment[1].length()` in Java.
    pub aligned_len: i32,
    seq1: String,
    seq2: String,
}

impl SmithWaterman {
    /// Initializes a new SmithWaterman aligner with given scoring parameters.
    pub fn new(m: i32, mis: i32, gap: i32) -> Self {
        Self { 
            match_score: m, 
            mismatch_penalty: mis, 
            gap_penalty: gap, 
            score: 0, 
            aligned_len: 0,
            seq1: String::new(), 
            seq2: String::new() 
        }
    }

    /// Sets the sequences to be aligned.
    pub fn set_seq(&mut self, s1: &str, s2: &str) {
        self.seq1 = s1.to_uppercase();
        self.seq2 = s2.to_uppercase();
    }

    /// Performs local alignment and traceback to calculate the optimal score and aligned length.
    pub fn align(&mut self) {
        // Java parity (`smith` package): table rows = seq2, cols = seq1.
        let cols = self.seq1.len();
        let rows = self.seq2.len();
        if cols == 0 || rows == 0 {
            self.score = 0;
            self.aligned_len = 0;
            return;
        }

        let s1 = self.seq1.as_bytes();
        let s2 = self.seq2.as_bytes();
        let mut score_table = vec![vec![0i32; cols + 1]; rows + 1];
        let mut prev: Vec<Vec<Option<(usize, usize)>>> = vec![vec![None; cols + 1]; rows + 1];
        let mut high_row = 0usize;
        let mut high_col = 0usize;

        for row in 1..=rows {
            for col in 1..=cols {
                let row_space_score = score_table[row - 1][col] + self.gap_penalty;
                let col_space_score = score_table[row][col - 1] + self.gap_penalty;
                let mut match_or_mismatch_score = score_table[row - 1][col - 1];
                if s2[row - 1] == s1[col - 1] {
                    match_or_mismatch_score += self.match_score;
                } else {
                    match_or_mismatch_score += self.mismatch_penalty;
                }

                let mut cell_score = 0i32;
                let mut cell_prev: Option<(usize, usize)> = None;
                if row_space_score >= col_space_score {
                    if match_or_mismatch_score >= row_space_score {
                        if match_or_mismatch_score > 0 {
                            cell_score = match_or_mismatch_score;
                            cell_prev = Some((row - 1, col - 1));
                        }
                    } else if row_space_score > 0 {
                        cell_score = row_space_score;
                        cell_prev = Some((row - 1, col));
                    }
                } else if match_or_mismatch_score >= col_space_score {
                    if match_or_mismatch_score > 0 {
                        cell_score = match_or_mismatch_score;
                        cell_prev = Some((row - 1, col - 1));
                    }
                } else if col_space_score > 0 {
                    cell_score = col_space_score;
                    cell_prev = Some((row, col - 1));
                }

                score_table[row][col] = cell_score;
                prev[row][col] = cell_prev;
                if cell_score > score_table[high_row][high_col] {
                    high_row = row;
                    high_col = col;
                }
            }
        }

        // Java traceback stops at score == 0.
        let mut row = high_row;
        let mut col = high_col;
        let mut align1: Vec<u8> = Vec::new();
        let mut align2: Vec<u8> = Vec::new();
        while score_table[row][col] != 0 {
            let (pr, pc) = match prev[row][col] {
                Some(p) => p,
                None => break,
            };
            if row - pr == 1 {
                align2.push(s2[row - 1]);
            } else {
                align2.push(b'-');
            }
            if col - pc == 1 {
                align1.push(s1[col - 1]);
            } else {
                align1.push(b'-');
            }
            row = pr;
            col = pc;
        }
        self.aligned_len = align2.len() as i32;

        // Java getAlignmentScore recomputes score from traceback alignments.
        let mut total = 0i32;
        for i in 0..align1.len() {
            let c1 = align1[i];
            let c2 = align2[i];
            if c1 == b'-' || c2 == b'-' {
                total += self.gap_penalty;
            } else if c1 == c2 {
                total += self.match_score;
            } else {
                total += self.mismatch_penalty;
            }
        }
        self.score = total;
    }
}

/// Main validator struct for BSJ candidate verification.
pub struct IsBSJHg2 {
    pub linear_range_size_min: i32,
    pub min_mapq_uni: i32,
    pub initial_size1: i32,
    pub aligner: SmithWaterman,
    window_unit: [i32; 5],
}

/// Helper function to mirror Java's `String.substring` behavior (end-exclusive, safe bounds).
pub fn java_substring(s: &str, start: i32, end: i32) -> &str {
    let len = s.len() as i32;
    assert!(start >= 0, "java_substring start<0: start={}, end={}, len={}", start, end, len);
    assert!(end >= 0, "java_substring end<0: start={}, end={}, len={}", start, end, len);
    assert!(start <= end, "java_substring start>end: start={}, end={}, len={}", start, end, len);
    assert!(end <= len, "java_substring end>len: start={}, end={}, len={}", start, end, len);
    &s[start as usize..end as usize]
}

impl IsBSJHg2 {
    /// Creates a new `IsBSJHg2` validator.
    pub fn new(linear_range_size_min: i32, min_mapq_uni: i32) -> Self {
        Self {
            linear_range_size_min,
            min_mapq_uni,
            initial_size1: 7,
            aligner: SmithWaterman::new(1, -1, -1),
            window_unit: [9, 7, 5, 4, 3],
        }
    }

    /// Internal distance check for window-based sequence mapping.
    fn distance_loci_stats(&self, locus_count: usize, locus_sum: i32, locus2_count: usize, locus2_sum: i32, window_step: i32) -> i32 {
        if locus_count < 2 && locus2_count < 2 { return 0; }
        if locus_sum <= window_step * locus_count as i32 && locus_sum * 20 < locus2_sum { 1 } else { 0 }
    }

    /// Linear competition check (Phase 1.1).
    pub fn is_in_circ_rna_1_1(&self, len_str: i32, str_val: &str, circ_range_seq: &str, linear_range: &str) -> i32 {
        for &step in &self.window_unit {
            if len_str < step * 2 { continue; }
            let window_size = step * 2;
            let mut trial = (len_str - window_size) / step;
            let mut locus_count = 0usize;
            let mut locus2_count = 0usize;
            let mut locus_sum = 0i32;
            let mut locus2_sum = 0i32;
            let mut prev_locus: Option<i32> = None;
            let mut prev_locus2: Option<i32> = None;
            let mut miss_count = [0, 0, 0]; 
            let mut miss_count2_total = 0;
            for j in 0..=trial {
                let s_idx = len_str - j * step - window_size;
                let e_idx = len_str - j * step;
                let seq = &str_val[s_idx as usize..e_idx as usize];
                if let Some(pos) = circ_range_seq.rfind(seq) {
                    let pos = pos as i32;
                    if let Some(prev) = prev_locus { locus_sum += (pos - prev).abs(); }
                    prev_locus = Some(pos);
                    locus_count += 1;
                    miss_count[0] = 0;
                }
                else { miss_count[1] += 1; miss_count[0] += 1; if miss_count[0] > miss_count[2] { miss_count[2] = miss_count[0]; } }
                if let Some(pos2) = linear_range.rfind(seq) {
                    let pos2 = pos2 as i32;
                    if let Some(prev) = prev_locus2 { locus2_sum += (pos2 - prev).abs(); }
                    prev_locus2 = Some(pos2);
                    locus2_count += 1;
                } else { miss_count2_total += 1; }
            }
            if len_str % step != 0 {
                trial += 1;
                let seq = &str_val[0..window_size as usize];
                if let Some(pos) = circ_range_seq.rfind(seq) {
                    let pos = pos as i32;
                    if let Some(prev) = prev_locus { locus_sum += (pos - prev).abs(); }
                    prev_locus = Some(pos);
                    locus_count += 1;
                    miss_count[0] = 0;
                }
                else { miss_count[1] += 1; miss_count[0] += 1; if miss_count[0] > miss_count[2] { miss_count[2] = miss_count[0]; } }
                if let Some(pos2) = linear_range.rfind(seq) {
                    let pos2 = pos2 as i32;
                    if let Some(prev) = prev_locus2 { locus2_sum += (pos2 - prev).abs(); }
                    prev_locus2 = Some(pos2);
                    locus2_count += 1;
                } else { miss_count2_total += 1; }
            }
            if miss_count2_total == 0 && miss_count[1] == 0 { if self.distance_loci_stats(locus_count, locus_sum, locus2_count, locus2_sum, step) == 1 { return 1; } else { return 0; } }
            else if miss_count2_total <= miss_count[1] { if locus2_count != 0 { return 0; } }
            else if miss_count[2] > 5 || miss_count[1] * 2 > trial { continue; } else { return 1; }
        }
        0
    }

    /// Linear competition check (Phase 1.2).
    pub fn is_in_circ_rna_1_2(&self, len_str: i32, str_val: &str, circ_range_seq: &str, linear_range: &str) -> i32 {
        for &step in &self.window_unit {
            if len_str < step * 2 { continue; }
            let window_size = step * 2;
            let mut trial = (len_str - window_size) / step;
            let mut locus_count = 0usize;
            let mut locus2_count = 0usize;
            let mut locus_sum = 0i32;
            let mut locus2_sum = 0i32;
            let mut prev_locus: Option<i32> = None;
            let mut prev_locus2: Option<i32> = None;
            let mut miss_count = [0, 0, 0];
            let mut miss_count2_total = 0;
            for j in 0..=trial {
                let s_idx = j * step;
                let e_idx = j * step + window_size;
                let seq = &str_val[s_idx as usize..e_idx as usize];
                if let Some(pos) = circ_range_seq.find(seq) {
                    let pos = pos as i32;
                    if let Some(prev) = prev_locus { locus_sum += (pos - prev).abs(); }
                    prev_locus = Some(pos);
                    locus_count += 1;
                    miss_count[0] = 0;
                }
                else { miss_count[1] += 1; miss_count[0] += 1; if miss_count[0] > miss_count[2] { miss_count[2] = miss_count[0]; } }
                if let Some(pos2) = linear_range.find(seq) {
                    let pos2 = pos2 as i32;
                    if let Some(prev) = prev_locus2 { locus2_sum += (pos2 - prev).abs(); }
                    prev_locus2 = Some(pos2);
                    locus2_count += 1;
                } else { miss_count2_total += 1; }
            }
            if len_str % step != 0 {
                trial += 1;
                let seq = &str_val[(len_str - window_size) as usize..len_str as usize];
                if let Some(pos) = circ_range_seq.find(seq) {
                    let pos = pos as i32;
                    if let Some(prev) = prev_locus { locus_sum += (pos - prev).abs(); }
                    prev_locus = Some(pos);
                    locus_count += 1;
                    miss_count[0] = 0;
                }
                else { miss_count[1] += 1; miss_count[0] += 1; if miss_count[0] > miss_count[2] { miss_count[2] = miss_count[0]; } }
                if let Some(pos2) = linear_range.find(seq) {
                    let pos2 = pos2 as i32;
                    if let Some(prev) = prev_locus2 { locus2_sum += (pos2 - prev).abs(); }
                    prev_locus2 = Some(pos2);
                    locus2_count += 1;
                } else { miss_count2_total += 1; }
            }
            if miss_count2_total == 0 && miss_count[1] == 0 { if self.distance_loci_stats(locus_count, locus_sum, locus2_count, locus2_sum, step) == 1 { return 1; } else { return 0; } }
            else if miss_count2_total <= miss_count[1] { if locus2_count != 0 { return 0; } }
            else if miss_count[2] > 5 || miss_count[1] * 2 > trial { continue; } else { return 1; }
        }
        0
    }

    /// Validates sequence mapping for unmapped segments (Phase 2).
    pub fn is_in_circ_rna_2(&self, unmap_seq: &str, circ_range_seq: &str) -> i32 {
        let seq_len = unmap_seq.len() as i32;
        let window_step = 5; let window_size = 10;
        let mut miss_count3 = [0, 0, 0];
        let trial = (seq_len - window_size) / window_step;
        for j in 0..=trial {
            let seq = if seq_len <= 10 { unmap_seq } else { &unmap_seq[(j * window_step) as usize..(j * window_step + window_size) as usize] };
            if circ_range_seq.contains(seq) { miss_count3[0] = 0; }
            else { miss_count3[1] += 1; miss_count3[0] += 1; if miss_count3[0] > miss_count3[2] { miss_count3[2] = miss_count3[0]; } }
        }
        if miss_count3[2] > 5 || (miss_count3[1] - 1) * 2 > trial { return 0; }
        1
    }

    /// Validates sequence mapping for mate reads (Phase 3).
    pub fn is_in_circ_rna_3(&self, ano_read: &str, pre_judge: &str, circ_range_seq: &str, pem_null_range_seq: &str) -> i32 {
        let ano_read_len = ano_read.len() as i32;
        let window_step = 5; let window_size = 10;
        let trial = (ano_read_len - window_size) / window_step;
        let mut miss_count = [0, 0, 0]; let mut miss_count2_total = 0;
        let mut locus_count = 0usize; let mut locus2_count = 0usize;
        let mut locus_sum = 0i32; let mut locus2_sum = 0i32;
        let mut prev_locus: Option<i32> = None; let mut prev_locus2: Option<i32> = None;
        for j in 0..=trial {
            let seq = &ano_read[(j * window_step) as usize..(window_size + j * window_step) as usize];
            if let Some(pos) = circ_range_seq.find(seq) {
                let pos = pos as i32;
                if let Some(prev) = prev_locus { locus_sum += (pos - prev).abs(); }
                prev_locus = Some(pos);
                locus_count += 1;
                miss_count[0] = 0;
            }
            else { miss_count[1] += 1; miss_count[0] += 1; if miss_count[0] > miss_count[2] { miss_count[2] = miss_count[0]; } }
            if !pem_null_range_seq.is_empty() {
                if let Some(pos2) = pem_null_range_seq.find(seq) {
                    let pos2 = pos2 as i32;
                    if let Some(prev) = prev_locus2 { locus2_sum += (pos2 - prev).abs(); }
                    prev_locus2 = Some(pos2);
                    locus2_count += 1;
                } else { miss_count2_total += 1; }
            }
        }
        if !pem_null_range_seq.is_empty() {
            if miss_count2_total == 0 && miss_count[1] == 0 { if self.distance_loci_stats(locus_count, locus_sum, locus2_count, locus2_sum, window_step) == 1 { return 1; } else { return -2; } }
            else if miss_count2_total <= miss_count[1] { if locus2_count != 0 { return -2; } else { return -1; } }
            else if miss_count[1] * 4 > trial * 3 && pre_judge == "0" { return -2; }
            else if miss_count[2] > 5 || miss_count[1] * 2 > trial { return -1; } else { return 1; }
        } else {
            if miss_count[1] * 4 > trial * 3 && pre_judge == "0" { return -2; }
            else if miss_count[2] > 5 || miss_count[1] * 2 > trial { return -1; }
        }
        1
    }

    /// High-level BSJ identification from Scan 1 metadata.
    pub fn is_bsj_hg1(&mut self, circ_line_arr: &mut [String], chr_taga: &str, sum_q: i32, mitochondrion: &str, sp_label: bool, chr_exon_start_map: &HashMap<String, String>, chr_exon_end_map: &HashMap<String, String>) -> Option<String> {
        let site1 = circ_line_arr[9].parse::<i32>().unwrap_or(0);
        let site2 = circ_line_arr[10].parse::<i32>().unwrap_or(0);
        let end_adjt1 = circ_line_arr[11].parse::<i32>().unwrap_or(0);
        let end_adjt2 = circ_line_arr[12].parse::<i32>().unwrap_or(0);
        let total_adjustment = end_adjt1 + end_adjt2;
        let chr_taga_len = chr_taga.len() as i32;

        let (end_string1, end_string2, tmp_site1, tmp_site2, adjt_bp);
        if end_adjt2 >= 0 {
            tmp_site1 = site1 - end_adjt1 - 1; tmp_site2 = site2 - end_adjt1 - 1; adjt_bp = 2 + total_adjustment;
            if site1 - end_adjt1 - 4 >= 0 {
                end_string1 = java_substring(chr_taga, site1 - end_adjt1 - 4, end_adjt2 + site1).to_string();
            } else {
                let n_pad = (0 - (site1 - end_adjt1 - 4)) as usize;
                end_string1 = format!("{}{}", "N".repeat(n_pad), java_substring(chr_taga, 0, end_adjt2 + site1));
            }
            end_string2 = java_substring(chr_taga, site2 - end_adjt1 - 1, 3 + end_adjt2 + site2).to_string();
        } else {
            tmp_site1 = site1 + end_adjt2 - 1; tmp_site2 = site2 + end_adjt2 - 1; adjt_bp = 2 - total_adjustment;
            if site1 + end_adjt2 - 4 >= 0 {
                end_string1 = java_substring(chr_taga, site1 + end_adjt2 - 4, site1 - end_adjt1).to_string();
            } else {
                let n_pad = (0 - (site1 + end_adjt2 - 4)) as usize;
                end_string1 = format!("{}{}", "N".repeat(n_pad), java_substring(chr_taga, 0, site1 - end_adjt1));
            }
            end_string2 = java_substring(chr_taga, site2 + end_adjt2 - 1, 3 + site2 - end_adjt1).to_string();
        }

        let mut index_strand_map = if circ_line_arr[1] == mitochondrion || sp_label { IndexCompare::index_compare_chrm(&end_string1, &end_string2) } else { IndexCompare::index_compare(&end_string1, &end_string2) };
        let chr = circ_line_arr[1].as_str();
        let mut start_key = String::with_capacity(chr.len() + 16);
        let mut end_key = String::with_capacity(chr.len() + 16);
        for i in 0..=adjt_bp {
            start_key.clear();
            start_key.push_str(chr);
            start_key.push('\t');
            let _ = write!(&mut start_key, "{}", tmp_site1 + i);
            end_key.clear();
            end_key.push_str(chr);
            end_key.push('\t');
            let _ = write!(&mut end_key, "{}", tmp_site2 + i);
            if !index_strand_map.contains_key(&i) {
                if let (Some(gene_start), Some(gene_end)) = (chr_exon_start_map.get(&start_key), chr_exon_end_map.get(&end_key)) {
                    if gene_start == gene_end {
                        let mut parts = gene_start.split('\t');
                        let _ = parts.next();
                        let strand = parts.next().unwrap_or("");
                        index_strand_map.insert(i, format!("{}\t{}\t{}\t{}", i, strand, java_substring(chr_taga, tmp_site1 + i - 3, tmp_site1 + i - 1), java_substring(chr_taga, tmp_site2 + i, tmp_site2 + i + 2)));
                    }
                }
            }
        }

        if !index_strand_map.is_empty() {
            for (&shift, sig) in &index_strand_map {
                let mut shift_parts = sig.split('\t');
                let _shift_idx = shift_parts.next().unwrap_or("");
                let shift_strand = shift_parts.next().unwrap_or("");
                let shift_left = shift_parts.next().unwrap_or("");
                let shift_right = shift_parts.next().unwrap_or("");
                let diff_adjt = if end_adjt2 >= 0 { shift - 1 - end_adjt1 } else { shift - 1 + total_adjustment - end_adjt1 };
                let site1_new = site1 + diff_adjt; let site2_new = site2 + diff_adjt;
                let mut str_new = ["".to_string(), "".to_string()];
                if diff_adjt >= 0 {
                    let str_adj = java_substring(&circ_line_arr[2], 0, diff_adjt);
                    str_new[1] = String::with_capacity(circ_line_arr[3].len() + str_adj.len());
                    str_new[1].push_str(&circ_line_arr[3]);
                    str_new[1].push_str(str_adj);
                    str_new[0] = java_substring(&circ_line_arr[2], diff_adjt, circ_line_arr[2].len() as i32).to_string();
                } else {
                    let str_adj = java_substring(&circ_line_arr[3], circ_line_arr[3].len() as i32 + diff_adjt, circ_line_arr[3].len() as i32);
                    str_new[0] = String::with_capacity(str_adj.len() + circ_line_arr[2].len());
                    str_new[0].push_str(str_adj);
                    str_new[0].push_str(&circ_line_arr[2]);
                    str_new[1] = java_substring(&circ_line_arr[3], 0, circ_line_arr[3].len() as i32 + diff_adjt).to_string();
                }
                {
                    let mut merged = String::with_capacity(shift_left.len() + str_new[0].len());
                    merged.push_str(shift_left);
                    merged.push_str(&str_new[0]);
                    str_new[0] = merged;
                }
                str_new[1].push_str(shift_right);
                let initial_seq1 = java_substring(&str_new[0], 0, self.initial_size1);
                let initial_seq2 = java_substring(&str_new[1], str_new[1].len() as i32 - self.initial_size1, str_new[1].len() as i32);
                let circ_range_seq = if site1_new - 3 < 0 && site2_new + 2 > chr_taga_len { &chr_taga[..] }
                else if site1_new - 3 < 0 { java_substring(chr_taga, 0, site2_new + 2) }
                else if site2_new + 2 > chr_taga_len { java_substring(chr_taga, site1_new - 3, chr_taga_len) }
                else { java_substring(chr_taga, site1_new - 3, site2_new + 2) };
                if circ_range_seq.starts_with(initial_seq1) &&
                   circ_range_seq.ends_with(initial_seq2) {
                    for i in 0..=1 {
                        if circ_line_arr[6 + i] != "1" {
                            let linear_range;
                            if i == 1 {
                                if site2_new - site1_new + 5 >= self.linear_range_size_min { if 2 * site1_new >= site2_new + 6 { linear_range = java_substring(chr_taga, 2 * site1_new - site2_new - 6, site1_new - 1); } else { linear_range = java_substring(chr_taga, 0, site1_new - 1); } }
                                else { if site1_new >= self.linear_range_size_min + 1 { linear_range = java_substring(chr_taga, site1_new - self.linear_range_size_min - 1, site1_new - 1); } else { linear_range = java_substring(chr_taga, 0, site1_new - 1); } }
                            } else {
                                if site2_new - site1_new + 5 >= self.linear_range_size_min { if 2 * site2_new - site1_new + 5 > chr_taga_len { linear_range = java_substring(chr_taga, site2_new, chr_taga_len); } else { linear_range = java_substring(chr_taga, site2_new, 2 * site2_new - site1_new + 5); } }
                                else { if site2_new + self.linear_range_size_min > chr_taga_len { linear_range = java_substring(chr_taga, site2_new, chr_taga_len); } else { linear_range = java_substring(chr_taga, site2_new, site2_new + self.linear_range_size_min); } }
                            }
                            circ_line_arr[6 + i] = self.is_in_circ_rna_1_2(str_new[i].len() as i32, &str_new[i], circ_range_seq, linear_range).to_string();
                        }
                    }
                    if circ_line_arr[6] == "1" && circ_line_arr[7] == "1" {
                        if circ_line_arr[4] != "*" && self.is_in_circ_rna_2(&circ_line_arr[4], circ_range_seq) == 0 { return None; }
                        let mut tag = 1;
                        if circ_line_arr[5].len() > 5 {
                            // Java parity: preserve original branch order/conditions.
                            let pem_null = if circ_line_arr[0] == "1" && site1_new - site1_new + 5 >= self.linear_range_size_min {
                                if 2 * site1_new >= site2_new + 6 {
                                    java_substring(chr_taga, 2 * site1_new - site2_new - 6, site1_new - 1)
                                } else {
                                    java_substring(chr_taga, 0, site1_new - 1)
                                }
                            } else if circ_line_arr[0] == "1" {
                                if site1_new >= self.linear_range_size_min + 1 {
                                    java_substring(chr_taga, site1_new - self.linear_range_size_min - 1, site1_new - 1)
                                } else {
                                    java_substring(chr_taga, 0, site1_new - 1)
                                }
                            } else if circ_line_arr[0] == "0" && site2_new - site1_new + 5 > self.linear_range_size_min {
                                if 2 * site2_new - site1_new + 5 > chr_taga_len {
                                    java_substring(chr_taga, site2_new, chr_taga_len)
                                } else {
                                    java_substring(chr_taga, site2_new, 2 * site2_new - site1_new + 5)
                                }
                            } else {
                                if site2_new + self.linear_range_size_min > chr_taga_len {
                                    java_substring(chr_taga, site2_new, chr_taga_len)
                                } else {
                                    java_substring(chr_taga, site2_new, site2_new + self.linear_range_size_min)
                                }
                            };
                            tag = self.is_in_circ_rna_3(&circ_line_arr[5], &circ_line_arr[8], circ_range_seq, pem_null);
                        }
                        return Some(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", tag, circ_line_arr[1], site1_new, site2_new, shift_strand, shift_left, shift_right, sum_q));
                    }
                }
            }
        }
        None
    }

    /// Low-level BSJ identification from rescue metadata (Scan 2).
    pub fn is_bsj_hg2(&mut self, circ_line_arr: &[String], chr_taga: &str) -> String {
        let mut judge_tag = "3".to_string();
        let start_site = circ_line_arr[3].parse::<i32>().unwrap_or(0);
        let end_site = circ_line_arr[4].parse::<i32>().unwrap_or(0);
        let quant = circ_line_arr[12].parse::<i32>().unwrap_or(0);
        let chr_taga_len = chr_taga.len() as i32;
        let circ_range_seq = if start_site - 3 < 0 && end_site + 2 > chr_taga_len { &chr_taga[..] }
        else if start_site - 3 < 0 { java_substring(chr_taga, 0, end_site + 2) }
        else if end_site + 2 > chr_taga_len { java_substring(chr_taga, start_site - 3, chr_taga_len) }
        else { java_substring(chr_taga, start_site - 3, end_site + 2) };
        let circ_range_len = circ_range_seq.len() as i32;
        let str_full; let mut pem_null_range_seq = "";
        if circ_line_arr[2] == "sm" {
            str_full = format!("{}{}", circ_line_arr[5], circ_line_arr[11]);
            let len_str = str_full.len() as i32;
            if len_str < 7 { if java_substring(circ_range_seq, circ_range_len - len_str, circ_range_len) != str_full { return "0".to_string(); } else { return "2".to_string(); } }
            else {
                let s1_part1 = java_substring(circ_range_seq, circ_range_len - 4, circ_range_len - 2);
                let s1_part2 = java_substring(&str_full, len_str - 4, len_str - 2);
                let s2_part1 = java_substring(circ_range_seq, circ_range_len - 7, circ_range_len - 4);
                let s2_part2 = java_substring(&str_full, len_str - 7, len_str - 4);
                if s1_part1 == s1_part2 || s2_part1 == s2_part2 {
                    let mut label = true; let linear_range;
                    if end_site - start_site + 5 >= self.linear_range_size_min { if 2 * start_site >= end_site + 6 { linear_range = java_substring(chr_taga, 2 * start_site - end_site - 6, start_site - 1); } else { linear_range = java_substring(chr_taga, 0, start_site - 1); } }
                    else { if start_site - 1 >= self.linear_range_size_min { linear_range = java_substring(chr_taga, start_site - self.linear_range_size_min - 1, start_site - 1); } else { linear_range = java_substring(chr_taga, 0, start_site - 1); } }
                    if circ_range_len < len_str { return "2".to_string(); }
                    self.aligner.set_seq(java_substring(circ_range_seq, circ_range_len - len_str, circ_range_len), &str_full); self.aligner.align();
                    if self.aligner.score >= (len_str - (len_str - 2) / 10 * 2) {
                        if quant >= self.min_mapq_uni {
                            let initial_seq = java_substring(&str_full, len_str - self.initial_size1, len_str);
                            if java_substring(circ_range_seq, circ_range_len - self.initial_size1, circ_range_len) == initial_seq {
                                if self.is_in_circ_rna_1_1(len_str, &str_full, circ_range_seq, linear_range) == 1 { judge_tag = "2".to_string(); } else { judge_tag = "3".to_string(); }
                            }
                        }
                        if circ_line_arr[0] == "1" { pem_null_range_seq = linear_range; } label = false;
                    }
                    if label {
                        let initial_seq = java_substring(&str_full, len_str - self.initial_size1, len_str);
                        if java_substring(circ_range_seq, circ_range_len - self.initial_size1, circ_range_len) == initial_seq {
                            if self.is_in_circ_rna_1_1(len_str, &str_full, circ_range_seq, linear_range) == 1 {
                                if quant >= self.min_mapq_uni {
                                    judge_tag = "2".to_string();
                                    if circ_line_arr[0] == "1" {
                                        pem_null_range_seq = linear_range;
                                    }
                                } else {
                                    return "2".to_string();
                                }
                            } else {
                                return "0".to_string();
                            }
                        } else { return "0".to_string(); }
                    }
                } else { return "0".to_string(); }
            }
        } else {
            str_full = format!("{}{}", circ_line_arr[10], circ_line_arr[5]);
            let len_str = str_full.len() as i32;
            if len_str < 7 { if java_substring(circ_range_seq, 0, len_str) != str_full { return "0".to_string(); } else { return "2".to_string(); } }
            else {
                let s1_part1 = java_substring(circ_range_seq, 2, 4);
                let s1_part2 = java_substring(&str_full, 2, 4);
                let s2_part1 = java_substring(circ_range_seq, 4, 7);
                let s2_part2 = java_substring(&str_full, 4, 7);
                if s1_part1 == s1_part2 || s2_part1 == s2_part2 {
                    let mut label = true; let linear_range;
                    if end_site - start_site + 5 >= self.linear_range_size_min { if 2 * end_site - start_site + 5 > chr_taga_len { linear_range = java_substring(chr_taga, end_site, chr_taga_len); } else { linear_range = java_substring(chr_taga, end_site, 2 * end_site - start_site + 5); } }
                    else { if end_site + self.linear_range_size_min > chr_taga_len { linear_range = java_substring(chr_taga, end_site, chr_taga_len); } else { linear_range = java_substring(chr_taga, end_site, end_site + self.linear_range_size_min); } }
                    if circ_range_len < len_str { return "2".to_string(); }
                    self.aligner.set_seq(java_substring(circ_range_seq, 0, len_str), &str_full); self.aligner.align();
                    if self.aligner.score >= (len_str - (len_str - 2) / 10 * 2) {
                        let initial_seq = java_substring(&str_full, 0, self.initial_size1);
                        if java_substring(circ_range_seq, 0, self.initial_size1) == initial_seq {
                            if self.is_in_circ_rna_1_2(len_str, &str_full, circ_range_seq, linear_range) == 1 { if quant >= self.min_mapq_uni { judge_tag = "2".to_string(); } } else { judge_tag = "3".to_string(); }
                        }
                        if circ_line_arr[0] == "0" { pem_null_range_seq = linear_range; } label = false;
                    }
                    if label {
                        let initial_seq = java_substring(&str_full, 0, self.initial_size1);
                        if java_substring(circ_range_seq, 0, self.initial_size1) == initial_seq {
                            if self.is_in_circ_rna_1_2(len_str, &str_full, circ_range_seq, linear_range) == 1 {
                                if quant >= self.min_mapq_uni {
                                    judge_tag = "2".to_string();
                                    if circ_line_arr[0] == "0" {
                                        pem_null_range_seq = linear_range;
                                    }
                                } else {
                                    return "2".to_string();
                                }
                            } else {
                                return "0".to_string();
                            }
                        } else { return "0".to_string(); }
                    }
                } else { return "0".to_string(); }
            }
        }
        if circ_line_arr[7] != "*" && self.is_in_circ_rna_2(&circ_line_arr[7], circ_range_seq) == 0 { return "0".to_string(); }
        if circ_line_arr[6].len() > 5 {
            let res = self.is_in_circ_rna_3(&circ_line_arr[6], &circ_line_arr[8], circ_range_seq, pem_null_range_seq);
            return format!("{}{}", res, judge_tag);
        }
        format!("1{}", judge_tag)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_java_substring() {
        assert_eq!(java_substring("ABCDE", 1, 3), "BC");
        assert_eq!(java_substring("ABCDE", 0, 5), "ABCDE");
    }

    #[test]
    fn test_smith_waterman_exact() {
        let mut sw = SmithWaterman::new(1, -1, -1);
        sw.set_seq("ATGC", "ATGC");
        sw.align();
        assert_eq!(sw.score, 4);
        assert_eq!(sw.aligned_len, 4);
    }

    #[test]
    fn test_smith_waterman_mismatch() {
        let mut sw = SmithWaterman::new(1, -1, -1);
        sw.set_seq("ATGC", "ATGG");
        sw.align();
        assert_eq!(sw.score, 3);
    }

    #[test]
    fn test_smith_waterman_short() {
        let mut sw = SmithWaterman::new(1, -1, -1);
        sw.set_seq("GCAT", "GC");
        sw.align();
        assert_eq!(sw.score, 2);
    }
}
