/// Validator module: Reproduces the core validation logic of Java CIRI3.
/// This module includes comprehensive sequence validation and linear competition checks.

use crate::index_compare::IndexCompare;
use std::collections::HashMap;

/// Smith-Waterman local alignment implementation with traceback.
pub struct SmithWaterman {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_penalty: i32,
    pub score: i32,
    /// Number of steps in the traceback, equivalent to alignment[1].length() in Java.
    pub aligned_len: i32,
    seq1: String,
    seq2: String,
}

impl SmithWaterman {
    /// Initializes a new SmithWaterman aligner.
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

    pub fn set_seq(&mut self, s1: &str, s2: &str) {
        self.seq1 = s1.to_uppercase();
        self.seq2 = s2.to_uppercase();
    }

    /// Performs local alignment and traceback to calculate score and aligned length.
    pub fn align(&mut self) {
        let n = self.seq1.len();
        let m = self.seq2.len();
        if n == 0 || m == 0 { self.score = 0; self.aligned_len = 0; return; }
        
        let mut dp = vec![vec![0; m + 1]; n + 1];
        let mut max_score = 0;
        let mut max_i = 0;
        let mut max_j = 0;
        let s1_bytes = self.seq1.as_bytes();
        let s2_bytes = self.seq2.as_bytes();
        
        for i in 1..=n {
            for j in 1..=m {
                let score = if s1_bytes[i - 1] == s2_bytes[j - 1] { self.match_score } else { self.mismatch_penalty };
                dp[i][j] = (dp[i - 1][j - 1] + score)
                    .max(dp[i - 1][j] + self.gap_penalty)
                    .max(dp[i][j - 1] + self.gap_penalty)
                    .max(0);
                    
                if dp[i][j] >= max_score { 
                    max_score = dp[i][j];
                    max_i = i;
                    max_j = j;
                }
            }
        }
        self.score = max_score;

        // Traceback to find aligned length (parity with Summary.java alignment[1].length())
        let mut curr_i = max_i;
        let mut curr_j = max_j;
        let mut steps = 0;
        while curr_i > 0 && curr_j > 0 && dp[curr_i][curr_j] > 0 {
            let score = if s1_bytes[curr_i - 1] == s2_bytes[curr_j - 1] { self.match_score } else { self.mismatch_penalty };
            if dp[curr_i][curr_j] == dp[curr_i - 1][curr_j - 1] + score {
                curr_i -= 1; curr_j -= 1;
            } else if dp[curr_i][curr_j] == dp[curr_i - 1][curr_j] + self.gap_penalty {
                curr_i -= 1;
            } else {
                curr_j -= 1;
            }
            steps += 1;
        }
        self.aligned_len = steps;
    }
}

pub struct IsBSJHg2 {
    pub linear_range_size_min: i32,
    pub min_mapq_uni: i32,
    pub initial_size1: i32,
    pub aligner: SmithWaterman,
    window_unit: [i32; 5],
}

pub fn java_substring(s: &str, start: i32, end: i32) -> &str {
    let len = s.len() as i32;
    let s_idx = start.max(0).min(len);
    let e_idx = end.max(0).min(len);
    if s_idx >= e_idx { return ""; }
    &s[s_idx as usize..e_idx as usize]
}

impl IsBSJHg2 {
    pub fn new(linear_range_size_min: i32, min_mapq_uni: i32) -> Self {
        Self {
            linear_range_size_min,
            min_mapq_uni,
            initial_size1: 7,
            aligner: SmithWaterman::new(1, -1, -1),
            window_unit: [9, 7, 5, 4, 3],
        }
    }

    fn distance_loci(&self, locus_list: &[i32], locus2_list: &[i32], window_step: i32) -> i32 {
        if locus_list.len() < 2 && locus2_list.len() < 2 { return 0; }
        let mut locus_sum = 0;
        let mut locus2_sum = 0;
        for i in 1..locus_list.len() { locus_sum += (locus_list[i] - locus_list[i-1]).abs(); }
        for i in 1..locus2_list.len() { locus2_sum += (locus2_list[i] - locus2_list[i-1]).abs(); }
        if locus_sum <= window_step * locus_list.len() as i32 && locus_sum * 20 < locus2_sum { 1 } else { 0 }
    }

    pub fn is_in_circ_rna_1_1(&self, len_str: i32, str_val: &str, circ_range_seq: &str, linear_range: &str) -> i32 {
        for &step in &self.window_unit {
            if len_str < step * 2 { continue; }
            let mut locus_list = Vec::new();
            let mut locus2_list = Vec::new();
            let window_size = step * 2;
            let mut trial = (len_str - window_size) / step;
            let mut miss_count = [0, 0, 0]; 
            let mut miss_count2_total = 0;
            for j in 0..=trial {
                let s_idx = len_str - j * step - window_size;
                let e_idx = len_str - j * step;
                let seq = &str_val[s_idx as usize..e_idx as usize];
                if let Some(pos) = circ_range_seq.rfind(seq) { locus_list.push(pos as i32); miss_count[0] = 0; }
                else { miss_count[1] += 1; miss_count[0] += 1; if miss_count[0] > miss_count[2] { miss_count[2] = miss_count[0]; } }
                if let Some(pos2) = linear_range.rfind(seq) { locus2_list.push(pos2 as i32); } else { miss_count2_total += 1; }
            }
            if len_str % step != 0 {
                trial += 1;
                let seq = &str_val[0..window_size as usize];
                if let Some(pos) = circ_range_seq.rfind(seq) { locus_list.push(pos as i32); miss_count[0] = 0; }
                else { miss_count[1] += 1; miss_count[0] += 1; if miss_count[0] > miss_count[2] { miss_count[2] = miss_count[0]; } }
                if let Some(pos2) = linear_range.rfind(seq) { locus2_list.push(pos2 as i32); } else { miss_count2_total += 1; }
            }
            if miss_count2_total == 0 && miss_count[1] == 0 { if self.distance_loci(&locus_list, &locus2_list, step) == 1 { return 1; } else { return 0; } }
            else if miss_count2_total <= miss_count[1] { if !locus2_list.is_empty() { return 0; } }
            else if miss_count[2] > 5 || miss_count[1] * 2 > trial { continue; } else { return 1; }
        }
        0
    }

    pub fn is_in_circ_rna_1_2(&self, len_str: i32, str_val: &str, circ_range_seq: &str, linear_range: &str) -> i32 {
        for &step in &self.window_unit {
            if len_str < step * 2 { continue; }
            let mut locus_list = Vec::new();
            let mut locus2_list = Vec::new();
            let window_size = step * 2;
            let mut trial = (len_str - window_size) / step;
            let mut miss_count = [0, 0, 0];
            let mut miss_count2_total = 0;
            for j in 0..=trial {
                let s_idx = j * step;
                let e_idx = j * step + window_size;
                let seq = &str_val[s_idx as usize..e_idx as usize];
                if let Some(pos) = circ_range_seq.find(seq) { locus_list.push(pos as i32); miss_count[0] = 0; }
                else { miss_count[1] += 1; miss_count[0] += 1; if miss_count[0] > miss_count[2] { miss_count[2] = miss_count[0]; } }
                if let Some(pos2) = linear_range.find(seq) { locus2_list.push(pos2 as i32); } else { miss_count2_total += 1; }
            }
            if len_str % step != 0 {
                trial += 1;
                let seq = &str_val[(len_str - window_size) as usize..len_str as usize];
                if let Some(pos) = circ_range_seq.find(seq) { locus_list.push(pos as i32); miss_count[0] = 0; }
                else { miss_count[1] += 1; miss_count[0] += 1; if miss_count[0] > miss_count[2] { miss_count[2] = miss_count[0]; } }
                if let Some(pos2) = linear_range.find(seq) { locus2_list.push(pos2 as i32); } else { miss_count2_total += 1; }
            }
            if miss_count2_total == 0 && miss_count[1] == 0 { if self.distance_loci(&locus_list, &locus2_list, step) == 1 { return 1; } else { return 0; } }
            else if miss_count2_total <= miss_count[1] { if !locus2_list.is_empty() { return 0; } }
            else if miss_count[2] > 5 || miss_count[1] * 2 > trial { continue; } else { return 1; }
        }
        0
    }

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

    pub fn is_in_circ_rna_3(&self, ano_read: &str, pre_judge: &str, circ_range_seq: &str, pem_null_range_seq: &str) -> i32 {
        let ano_read_len = ano_read.len() as i32;
        let window_step = 5; let window_size = 10;
        let trial = (ano_read_len - window_size) / window_step;
        let mut miss_count = [0, 0, 0]; let mut miss_count2_total = 0;
        let mut locus_list = Vec::new(); let mut locus2_list = Vec::new();
        for j in 0..=trial {
            let seq = &ano_read[(j * window_step) as usize..(window_size + j * window_step) as usize];
            if let Some(pos) = circ_range_seq.find(seq) { locus_list.push(pos as i32); miss_count[0] = 0; }
            else { miss_count[1] += 1; miss_count[0] += 1; if miss_count[0] > miss_count[2] { miss_count[2] = miss_count[0]; } }
            if !pem_null_range_seq.is_empty() {
                if let Some(pos2) = pem_null_range_seq.find(seq) { locus2_list.push(pos2 as i32); } else { miss_count2_total += 1; }
            }
        }
        if !pem_null_range_seq.is_empty() {
            if miss_count2_total == 0 && miss_count[1] == 0 { if self.distance_loci(&locus_list, &locus2_list, window_step) == 1 { return 1; } else { return -2; } }
            else if miss_count2_total <= miss_count[1] { if !locus2_list.is_empty() { return -2; } else { return -1; } }
            else if miss_count[1] * 4 > trial * 3 && pre_judge == "0" { return -2; }
            else if miss_count[2] > 5 || miss_count[1] * 2 > trial { return -1; } else { return 1; }
        } else {
            if miss_count[1] * 4 > trial * 3 && pre_judge == "0" { return -2; }
            else if miss_count[2] > 5 || miss_count[1] * 2 > trial { return -1; }
        }
        1
    }

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
            end_string1 = java_substring(chr_taga, site1 - end_adjt1 - 4, end_adjt2 + site1).to_string();
            end_string2 = java_substring(chr_taga, site2 - end_adjt1 - 1, 3 + end_adjt2 + site2).to_string();
        } else {
            tmp_site1 = site1 + end_adjt2 - 1; tmp_site2 = site2 + end_adjt2 - 1; adjt_bp = 2 - total_adjustment;
            end_string1 = java_substring(chr_taga, site1 + end_adjt2 - 4, site1 - end_adjt1).to_string();
            end_string2 = java_substring(chr_taga, site2 + end_adjt2 - 1, 3 + site2 - end_adjt1).to_string();
        }

        let mut index_strand_map = if circ_line_arr[1] == mitochondrion || sp_label { IndexCompare::index_compare_chrm(&end_string1, &end_string2) } else { IndexCompare::index_compare(&end_string1, &end_string2) };
        for i in 0..=adjt_bp {
            let start_key = format!("{}\t{}", circ_line_arr[1], tmp_site1 + i);
            let end_key = format!("{}\t{}", circ_line_arr[1], tmp_site2 + i);
            if !index_strand_map.contains_key(&i) && chr_exon_start_map.contains_key(&start_key) && chr_exon_end_map.contains_key(&end_key) {
                let gene_stand = chr_exon_start_map.get(&start_key).unwrap();
                if gene_stand == chr_exon_end_map.get(&end_key).unwrap() {
                    let parts: Vec<&str> = gene_stand.split('\t').collect();
                    index_strand_map.insert(i, format!("{}\t{}\t{}\t{}", i, parts[1], java_substring(chr_taga, tmp_site1 + i - 3, tmp_site1 + i - 1), java_substring(chr_taga, tmp_site2 + i, tmp_site2 + i + 2)));
                }
            }
        }

        if !index_strand_map.is_empty() {
            for (&shift, sig) in &index_strand_map {
                let shift_arr: Vec<&str> = sig.split('\t').collect();
                let diff_adjt = if end_adjt2 >= 0 { shift - 1 - end_adjt1 } else { shift - 1 + total_adjustment - end_adjt1 };
                let site1_new = site1 + diff_adjt; let site2_new = site2 + diff_adjt;
                let mut str_new = ["".to_string(), "".to_string()];
                if diff_adjt >= 0 {
                    let str_adj = java_substring(&circ_line_arr[2], 0, diff_adjt);
                    str_new[1] = format!("{}{}", circ_line_arr[3], str_adj);
                    str_new[0] = java_substring(&circ_line_arr[2], diff_adjt, circ_line_arr[2].len() as i32).to_string();
                } else {
                    let str_adj = java_substring(&circ_line_arr[3], circ_line_arr[3].len() as i32 + diff_adjt, circ_line_arr[3].len() as i32);
                    str_new[0] = format!("{}{}", str_adj, circ_line_arr[2]);
                    str_new[1] = java_substring(&circ_line_arr[3], 0, circ_line_arr[3].len() as i32 + diff_adjt).to_string();
                }
                str_new[0] = format!("{}{}", shift_arr[2], str_new[0]); str_new[1] = format!("{}{}", str_new[1], shift_arr[3]);
                let initial_seq1 = java_substring(&str_new[0], 0, self.initial_size1);
                let initial_seq2 = java_substring(&str_new[1], str_new[1].len() as i32 - self.initial_size1, str_new[1].len() as i32);
                let circ_range_seq = if site1_new - 3 < 0 && site2_new + 2 > chr_taga_len { &chr_taga[..] }
                else if site1_new - 3 < 0 { java_substring(chr_taga, 0, site2_new + 2) }
                else if site2_new + 2 > chr_taga_len { java_substring(chr_taga, site1_new - 3, chr_taga_len) }
                else { java_substring(chr_taga, site1_new - 3, site2_new + 2) };
                let circ_range_len = circ_range_seq.len() as i32;

                if java_substring(circ_range_seq, 0, initial_seq1.len() as i32) == initial_seq1 &&
                   java_substring(circ_range_seq, circ_range_len - initial_seq2.len() as i32, circ_range_len) == initial_seq2 {
                    let mut j_ok = [circ_line_arr[6].clone(), circ_line_arr[7].clone()];
                    for i in 0..=1 {
                        if j_ok[i] != "1" {
                            let linear_range;
                            if i == 1 {
                                if site2_new - site1_new + 5 >= self.linear_range_size_min { if 2 * site1_new >= site2_new + 6 { linear_range = java_substring(chr_taga, 2 * site1_new - site2_new - 6, site1_new - 1); } else { linear_range = java_substring(chr_taga, 0, site1_new - 1); } }
                                else { if site1_new >= self.linear_range_size_min + 1 { linear_range = java_substring(chr_taga, site1_new - self.linear_range_size_min - 1, site1_new - 1); } else { linear_range = java_substring(chr_taga, 0, site1_new - 1); } }
                            } else {
                                if site2_new - site1_new + 5 >= self.linear_range_size_min { if 2 * site2_new - site1_new + 5 > chr_taga_len { linear_range = java_substring(chr_taga, site2_new, chr_taga_len); } else { linear_range = java_substring(chr_taga, site2_new, 2 * site2_new - site1_new + 5); } }
                                else { if site2_new + self.linear_range_size_min > chr_taga_len { linear_range = java_substring(chr_taga, site2_new, chr_taga_len); } else { linear_range = java_substring(chr_taga, site2_new, site2_new + self.linear_range_size_min); } }
                            }
                            j_ok[i] = self.is_in_circ_rna_1_2(str_new[i].len() as i32, &str_new[i], circ_range_seq, linear_range).to_string();
                        }
                    }
                    if j_ok[0] == "1" && j_ok[1] == "1" {
                        if circ_line_arr[4] != "*" && self.is_in_circ_rna_2(&circ_line_arr[4], circ_range_seq) == 0 { return None; }
                        let mut tag = 1;
                        if circ_line_arr[5].len() > 5 {
                            let mut pem_null = "";
                            if circ_line_arr[0] == "1" { if 2*site1_new >= site2_new+6 { pem_null = java_substring(chr_taga, 2*site1_new-site2_new-6, site1_new-1); } else { pem_null = java_substring(chr_taga, 0, site1_new-1); } }
                            else { if 2*site2_new-site1_new+5 > chr_taga_len { pem_null = java_substring(chr_taga, site2_new, chr_taga_len); } else { pem_null = java_substring(chr_taga, site2_new, 2*site2_new-site1_new+5); } }
                            tag = self.is_in_circ_rna_3(&circ_line_arr[5], &circ_line_arr[8], circ_range_seq, pem_null);
                        }
                        let mut final_s1 = site1_new; let mut final_s2 = site2_new;
                        if final_s1 > final_s2 { std::mem::swap(&mut final_s1, &mut final_s2); }
                        return Some(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", tag, circ_line_arr[1], final_s1, final_s2, shift_arr[1], shift_arr[2], shift_arr[3], sum_q));
                    }
                }
            }
        }
        None
    }

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
                                if self.is_in_circ_rna_1_1(len_str, &str_full, circ_range_seq, linear_range) == 1 { judge_tag = "1".to_string(); } else { judge_tag = "2".to_string(); }
                            }
                        }
                        if circ_line_arr[0] == "1" { pem_null_range_seq = linear_range; } label = false;
                    }
                    if label {
                        let initial_seq = java_substring(&str_full, len_str - self.initial_size1, len_str);
                        if java_substring(circ_range_seq, circ_range_len - self.initial_size1, circ_range_len) == initial_seq {
                            if self.is_in_circ_rna_1_1(len_str, &str_full, circ_range_seq, linear_range) == 1 { if quant >= self.min_mapq_uni { judge_tag = "1".to_string(); } } else { return "0".to_string(); }
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
                            if self.is_in_circ_rna_1_2(len_str, &str_full, circ_range_seq, linear_range) == 1 { if quant >= self.min_mapq_uni { judge_tag = "1".to_string(); } } else { judge_tag = "2".to_string(); }
                        }
                        if circ_line_arr[0] == "0" { pem_null_range_seq = linear_range; } label = false;
                    }
                    if label {
                        let initial_seq = java_substring(&str_full, 0, self.initial_size1);
                        if java_substring(circ_range_seq, 0, self.initial_size1) == initial_seq {
                            if self.is_in_circ_rna_1_2(len_str, &str_full, circ_range_seq, linear_range) == 1 { if quant >= self.min_mapq_uni { judge_tag = "1".to_string(); } } else { return "0".to_string(); }
                        } else { return "0".to_string(); }
                    }
                } else { return "0".to_string(); }
            }
        }
        if circ_line_arr[7] != "*" && self.is_in_circ_rna_2(&circ_line_arr[7], circ_range_seq) == 0 { return "0".to_string(); }
        if circ_line_arr[6].len() > 5 { let res = self.is_in_circ_rna_3(&circ_line_arr[6], &circ_line_arr[8], circ_range_seq, pem_null_range_seq); return format!("{}{}", res, judge_tag); }
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
