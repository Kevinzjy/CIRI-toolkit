//! Summary module aligned to Java `Summary.java` logic.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use anyhow::Result;

use crate::annotation::Annotation;
use crate::is_bsj_hg2::SmithWaterman;

pub struct Summary {
    pub stringency: i32,
}

struct CircSortItem {
    start_site: i32,
    line: String,
}

impl Summary {
    #[inline]
    /// Mirrors Java `String.hashCode()` for deterministic bucket ordering parity.
    fn java_string_hash(s: &str) -> i32 {
        let mut h: i32 = 0;
        for ch in s.chars() {
            h = h.wrapping_mul(31).wrapping_add(ch as i32);
        }
        h
    }

    #[inline]
    /// Mirrors Java HashMap hash spreading (`h ^ (h >>> 16)`).
    fn java_hash_spread(h: i32) -> u32 {
        let x = h as u32;
        x ^ (x >> 16)
    }

    #[inline]
    /// Mirrors Java HashSet capacity growth (initial 16, load factor 0.75).
    fn java_hashset_capacity(n: usize) -> usize {
        let mut cap = 16usize;
        while n > ((cap as f64) * 0.75f64) as usize {
            cap = cap.saturating_mul(2);
        }
        cap.max(1)
    }

    pub fn new(stringency: i32) -> Self {
        Self { stringency }
    }

    fn cigar_len_between_ms(cigar: &str) -> i32 {
        let ops: Vec<char> = cigar.chars().filter(|c| c.is_ascii_alphabetic()).collect();
        let nums: Vec<i32> = {
            let mut out = Vec::new();
            let mut num = 0i32;
            let mut seen_digit = false;
            for ch in cigar.chars() {
                if ch.is_ascii_digit() {
                    num = num * 10 + (ch as i32 - '0' as i32);
                    seen_digit = true;
                } else if ch.is_ascii_alphabetic() {
                    if seen_digit {
                        out.push(num);
                    }
                    num = 0;
                    seen_digit = false;
                }
            }
            out
        };
        if ops.is_empty() || nums.is_empty() || ops.len() != nums.len() {
            return 0;
        }
        let first_m = ops.iter().position(|&c| c == 'M');
        let last_m = ops.iter().rposition(|&c| c == 'M');
        if first_m.is_none() || last_m.is_none() {
            return 0;
        }
        let first = first_m.unwrap();
        let last = last_m.unwrap();
        let mut sum = 0i32;
        for i in first..=last {
            if ops[i] == 'M' || ops[i] == 'D' {
                sum += nums[i];
            }
        }
        sum
    }

    fn update_cigar_counts(cigar: &str, counts: &mut [i32; 3]) {
        let t: String = cigar.chars().filter(|x| x.is_ascii_alphabetic()).collect();
        if t == "M" {
            return;
        } else if t.find('M') > Some(0) && t.rfind('M') < Some(t.len() - 1) {
            counts[2] += 1;
        } else if t.find('M') == Some(0) {
            counts[1] += 1;
        } else {
            counts[0] += 1;
        }
    }

    fn merge_same_start(
        circ_map: &mut HashMap<String, HashSet<String>>,
        chr_tcga_map: &HashMap<String, String>,
        circ_start_insertion: &HashMap<String, Vec<String>>,
    ) {
        let mut aligner = SmithWaterman::new(1, -1, -3);
        for (chr_start, insertion_list) in circ_start_insertion {
            if insertion_list.len() <= 1 {
                continue;
            }
            // Java parity: circStartMap value is HashSet<String>. Emulate HashSet iteration
            // order using Java hash bucket order and insertion order within bucket.
            let cap = Self::java_hashset_capacity(insertion_list.len());
            let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); cap];
            for (idx, s) in insertion_list.iter().enumerate() {
                let spread = Self::java_hash_spread(Self::java_string_hash(s));
                let bucket_idx = (spread as usize) & (cap - 1);
                buckets[bucket_idx].push(idx);
            }
            let mut circ_names: Vec<String> = Vec::with_capacity(insertion_list.len());
            for bucket in buckets {
                for idx in bucket {
                    let k = &insertion_list[idx];
                    if circ_map.contains_key(k) {
                        circ_names.push(k.clone());
                    }
                }
            }
            if circ_names.len() <= 1 {
                continue;
            }
            let chr = chr_start.split('\t').next().unwrap_or_default();
            let chr_seq = match chr_tcga_map.get(chr) {
                Some(v) => v,
                None => continue,
            };
            let mut score = vec![0i32; circ_names.len()];
            let mut end_site = vec![0i32; circ_names.len()];
            let mut cigars = vec![0i32; circ_names.len()];

            for (idx, name) in circ_names.iter().enumerate() {
                if let Some(lines) = circ_map.get(name) {
                    let mut len = 0i32;
                    for line in lines {
                        let p: Vec<&str> = line.split('\t').collect();
                        if p.len() < 10 || !p[1].contains(';') {
                            continue;
                        }
                        end_site[idx] = p[5].parse().unwrap_or(0);
                        if p[2] == "1" && p[9] == "1" {
                            score[idx] += 1;
                        }
                        let c_arr: Vec<&str> = p[1].split(';').collect();
                        if c_arr.len() >= 2 {
                            let l = Self::cigar_len_between_ms(c_arr[1]);
                            if l > len {
                                len = l;
                            }
                        }
                    }
                    cigars[idx] = len;
                }
            }

            let up: Vec<usize> = score.iter().enumerate().filter_map(|(i, &v)| if v > 0 { Some(i) } else { None }).collect();
            let mut zero: Vec<usize> = score.iter().enumerate().filter_map(|(i, &v)| if v == 0 { Some(i) } else { None }).collect();

            for up_idx in up {
                let mut remove_list = Vec::new();
                for &z_idx in &zero {
                    let len2 = cigars[z_idx];
                    if len2 <= 0 {
                        continue;
                    }
                    let e1 = end_site[up_idx];
                    let e2 = end_site[z_idx];
                    if e1 <= len2 || e2 <= len2 {
                        continue;
                    }
                    let s1 = (e1 - len2) as usize;
                    let s2 = (e2 - len2) as usize;
                    let e1u = e1 as usize;
                    let e2u = e2 as usize;
                    if e1u > chr_seq.len() || e2u > chr_seq.len() || s1 >= e1u || s2 >= e2u {
                        continue;
                    }
                    let seq1 = &chr_seq[s1..e1u];
                    let seq2 = &chr_seq[s2..e2u];
                    aligner.set_seq(seq1, seq2);
                    aligner.align();
                    let pass_score = aligner.score >= seq1.len() as i32 - 2 - ((seq1.len() as i32 - 1) / 10) * 2;
                    let pass_len = aligner.aligned_len > seq1.len() as i32 - 2;
                    if pass_score && pass_len {
                        let up_name = &circ_names[up_idx];
                        let z_name = &circ_names[z_idx];
                        if let (Some(set_up), Some(set_z)) = (circ_map.get(up_name).cloned(), circ_map.get(z_name).cloned()) {
                            let mut merged = set_up;
                            merged.extend(set_z);
                            circ_map.insert(up_name.clone(), merged);
                            circ_map.remove(z_name);
                            remove_list.push(z_idx);
                        }
                    }
                }
                zero.retain(|v| !remove_list.contains(v));
            }
        }
    }

    fn merge_same_end(
        circ_map: &mut HashMap<String, HashSet<String>>,
        chr_tcga_map: &HashMap<String, String>,
    ) {
        let mut circ_end_map: HashMap<String, HashSet<String>> = HashMap::new();
        for k in circ_map.keys() {
            let p: Vec<&str> = k.split('\t').collect();
            if p.len() != 3 {
                continue;
            }
            let key = format!("{}\t{}", p[0], p[2]);
            circ_end_map.entry(key).or_default().insert(k.clone());
        }

        let mut aligner = SmithWaterman::new(1, -1, -3);
        for (chr_end, group_set) in circ_end_map {
            if group_set.len() <= 1 {
                continue;
            }
            let mut circ_names: Vec<String> = group_set.into_iter().collect();
            circ_names.retain(|k| circ_map.contains_key(k));
            if circ_names.len() <= 1 {
                continue;
            }
            let chr = chr_end.split('\t').next().unwrap_or_default();
            let chr_seq = match chr_tcga_map.get(chr) {
                Some(v) => v,
                None => continue,
            };
            let mut score = vec![0i32; circ_names.len()];
            let mut start_site = vec![0i32; circ_names.len()];
            let mut cigars = vec![0i32; circ_names.len()];

            for (idx, name) in circ_names.iter().enumerate() {
                if let Some(lines) = circ_map.get(name) {
                    let mut len = 0i32;
                    for line in lines {
                        let p: Vec<&str> = line.split('\t').collect();
                        if p.len() < 10 || !p[1].contains(';') {
                            continue;
                        }
                        start_site[idx] = p[4].parse().unwrap_or(0);
                        if p[2] == "1" && p[9] == "1" {
                            score[idx] += 1;
                        }
                        let c_arr: Vec<&str> = p[1].split(';').collect();
                        if !c_arr.is_empty() {
                            let l = Self::cigar_len_between_ms(c_arr[0]);
                            if l > len {
                                len = l;
                            }
                        }
                    }
                    cigars[idx] = len;
                }
            }

            let up: Vec<usize> = score.iter().enumerate().filter_map(|(i, &v)| if v > 0 { Some(i) } else { None }).collect();
            let mut zero: Vec<usize> = score.iter().enumerate().filter_map(|(i, &v)| if v == 0 { Some(i) } else { None }).collect();

            for up_idx in up {
                let mut remove_list = Vec::new();
                for &z_idx in &zero {
                    let len2 = cigars[z_idx];
                    if len2 <= 0 {
                        continue;
                    }
                    let s1 = start_site[up_idx];
                    let s2 = start_site[z_idx];
                    if s1 < 0 || s2 < 0 {
                        continue;
                    }
                    let e1 = s1 + len2;
                    let e2 = s2 + len2;
                    if e1 as usize > chr_seq.len() || e2 as usize > chr_seq.len() {
                        continue;
                    }
                    let seq1 = &chr_seq[s1 as usize..e1 as usize];
                    let seq2 = &chr_seq[s2 as usize..e2 as usize];
                    aligner.set_seq(seq1, seq2);
                    aligner.align();
                    let pass_score = aligner.score >= seq1.len() as i32 - 2 - ((seq1.len() as i32 - 1) / 10) * 2;
                    let pass_len = aligner.aligned_len > seq1.len() as i32 - 2;
                    if pass_score && pass_len {
                        let up_name = &circ_names[up_idx];
                        let z_name = &circ_names[z_idx];
                        if let (Some(set_up), Some(set_z)) = (circ_map.get(up_name).cloned(), circ_map.get(z_name).cloned()) {
                            let mut merged = set_up;
                            merged.extend(set_z);
                            circ_map.insert(up_name.clone(), merged);
                            circ_map.remove(z_name);
                            remove_list.push(z_idx);
                        }
                    }
                }
                zero.retain(|v| !remove_list.contains(v));
            }
        }
    }

    pub fn run(
        &mut self,
        bsj1_file: &str,
        out_prefix: &str,
        fsj_map: &HashMap<String, i32>,
        chr_tcga_map: &HashMap<String, String>,
        annotation: &Annotation,
    ) -> Result<()> {
        let file = File::open(bsj1_file)?;
        let reader = BufReader::new(file);
        let mut circ_map: HashMap<String, HashSet<String>> = HashMap::new();
        let mut circ_start_seen: HashMap<String, HashSet<String>> = HashMap::new();
        let mut circ_start_insertion: HashMap<String, Vec<String>> = HashMap::new();
        for line_res in reader.lines() {
            let line = line_res?;
            let p: Vec<&str> = line.split('\t').collect();
            if p.len() < 8 {
                continue;
            }
            let key = format!("{}\t{}\t{}", p[3], p[4], p[5]);
            circ_map.entry(key.clone()).or_default().insert(line.clone());
            let start_key = format!("{}\t{}", p[3], p[4]);
            let seen = circ_start_seen.entry(start_key.clone()).or_default();
            if seen.insert(key.clone()) {
                circ_start_insertion.entry(start_key).or_default().push(key);
            }
        }

        Self::merge_same_start(&mut circ_map, chr_tcga_map, &circ_start_insertion);
        Self::merge_same_end(&mut circ_map, chr_tcga_map);

        let mut final_results: BTreeMap<String, Vec<CircSortItem>> = BTreeMap::new();
        for (chr_start_end, lines) in &circ_map {
            let p_key: Vec<&str> = chr_start_end.split('\t').collect();
            if p_key.len() != 3 {
                continue;
            }

            let mut circ_id_set3: HashSet<String> = HashSet::new();
            let mut cigar_set3: HashSet<String> = HashSet::new();
            let mut false_cigar_set3: HashSet<String> = HashSet::new();
            let mut circ_id_set: HashSet<String> = HashSet::new();
            let mut cigar_set: HashSet<String> = HashSet::new();
            let mut false_cigar_set: HashSet<String> = HashSet::new();
            let mut non_reads = 0i32;
            let mut fp_reads = 0i32;
            let mut tag = 0i32;
            let mut strand = " ".to_string();

            for line in lines {
                let p: Vec<&str> = line.split('\t').collect();
                if p.len() < 10 {
                    continue;
                }
                strand = p[6].to_string();
                if p[2] == "1" {
                    if p[9] == "3" {
                        cigar_set3.insert(p[1].to_string());
                        circ_id_set3.insert(p[0].to_string());
                    } else {
                        for c in p[1].split(';') {
                            cigar_set.insert(c.to_string());
                        }
                        if p[9] == "1" {
                            tag += 1;
                        }
                        circ_id_set.insert(p[0].to_string());
                    }
                } else if p[2] == "-1" {
                    if p[9] == "3" {
                        false_cigar_set3.insert(p[1].to_string());
                    } else {
                        for c in p[1].split(';') {
                            false_cigar_set.insert(c.to_string());
                        }
                        non_reads += 1;
                    }
                } else if p[2] == "-2" {
                    if p[9] == "3" {
                        false_cigar_set3.insert(p[1].to_string());
                    } else {
                        for c in p[1].split(';') {
                            false_cigar_set.insert(c.to_string());
                        }
                        fp_reads += 1;
                    }
                }
            }

            cigar_set3.extend(cigar_set.iter().cloned());
            circ_id_set3.extend(circ_id_set.iter().cloned());
            false_cigar_set3.extend(false_cigar_set.iter().cloned());

            let mut cigar_count3 = [0, 0, 0];
            for c in &cigar_set3 {
                Self::update_cigar_counts(c, &mut cigar_count3);
            }

            let tp_reads = circ_id_set.len() as i32;
            let tp_reads3 = circ_id_set3.len() as i32;

            let passed = if self.stringency == 2 {
                ((tp_reads > 19 * fp_reads || fp_reads <= 1) && tp_reads > non_reads + fp_reads && cigar_set.len() >= 3 && tp_reads >= 2)
                    || (tag > 0 && false_cigar_set3.is_empty() && cigar_set3.len() >= 3 && tp_reads3 >= 2)
            } else if self.stringency == 1 {
                ((tp_reads > 19 * fp_reads || false_cigar_set.len() <= 2) && tp_reads > non_reads + fp_reads && tp_reads >= 2)
                    || (tag > 0 && false_cigar_set3.is_empty() && tp_reads3 >= 2)
            } else {
                ((tp_reads > 19 * fp_reads || false_cigar_set.len() <= 2) && tp_reads > non_reads + fp_reads)
                    || (tag > 0 && false_cigar_set3.is_empty() && tp_reads3 >= 2)
            };

            if !passed {
                continue;
            }

            let fsj = *fsj_map.get(chr_start_end).unwrap_or(&0);
            let ratio = if tp_reads3 * 2 + fsj > 0 {
                (tp_reads3 * 2) as f64 / (tp_reads3 * 2 + fsj) as f64
            } else {
                0.0
            };
            let start = p_key[1].parse::<i32>().unwrap_or(0);
            let mut circ_type = "intergenic_region".to_string();
            let mut gene_id = "NA".to_string();
            if let (Some(v1), Some(v2)) = (
                annotation.chr_exon_start_map.get(&format!("{}\t{}", p_key[0], p_key[1])),
                annotation.chr_exon_end_map.get(&format!("{}\t{}", p_key[0], p_key[2])),
            ) {
                if let (Some((g1, _)), Some((g2, _))) = (v1.split_once('\t'), v2.split_once('\t')) {
                    if g1 == g2 {
                        circ_type = "exon".to_string();
                        gene_id = g1.to_string();
                    }
                }
            }
            let mut ids: Vec<String> = circ_id_set3.into_iter().collect();
            ids.sort();
            let line = format!(
                "{}:{}|{}\t{}\t{}\t{}\t{}\t{}_{}_{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}",
                p_key[0],
                p_key[1],
                p_key[2],
                p_key[0],
                p_key[1],
                p_key[2],
                tp_reads3,
                cigar_count3[0],
                cigar_count3[1],
                cigar_count3[2],
                fsj,
                ratio,
                circ_type,
                gene_id,
                strand,
                ids.join(","),
                tag
            );
            final_results.entry(p_key[0].to_string()).or_default().push(CircSortItem { start_site: start, line });
        }

        let out_file = File::create(format!("{}.result", out_prefix))?;
        let mut writer = BufWriter::with_capacity(1024 * 1024, out_file);
        writeln!(writer, "circRNA_ID\tchr\tcircRNA_start\tcircRNA_end\t#junction_reads\tSM_MS_SMS\t#non_junction_reads\tjunction_reads_ratio\tcircRNA_type\tgene_id\tstrand\tjunction_reads_ID\tScore")?;
        for items in final_results.values_mut() {
            items.sort_by_key(|x| x.start_site);
            for item in items {
                writeln!(writer, "{}", item.line)?;
            }
        }
        writer.flush()?;
        Ok(())
    }
}
