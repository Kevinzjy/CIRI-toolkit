//! Summary module aligned to Java `Summary.java` logic.
//!
//! This pass merges nearby circRNA sites, applies Java-compatible stringency
//! filters, and writes the final circRNA report.

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use anyhow::Result;

use crate::annotation::Annotation;
use crate::is_bsj_hg2::SmithWaterman;
use crate::utils::{bsj_is_summary_priority, bsj_payload_start};

/// Final clustering and reporting stage.
///
/// `Summary` deliberately keeps Java-shaped merge and filtering rules, because
/// even small changes here alter the final circRNA set despite identical Scan1
/// and Scan2 intermediates.
pub struct Summary {
    pub stringency: i32,
    /// Number of circRNA rows written in the most recent run.
    pub circ_count: usize,
    /// Number of unique BSJ-supporting read IDs retained in the most recent final report.
    pub final_bsj_reads: usize,
}

/// One sortable final-result row for a chromosome.
///
/// Results are grouped by chromosome and then sorted by start coordinate to
/// preserve the stable output order expected by existing comparisons.
struct CircSortItem {
    start_site: i32,
    line: String,
}

impl Summary {
    /// Compares chromosome names using a "main chromosomes first" order.
    ///
    /// The old `BTreeMap` output order was pure lexicographic, which placed
    /// scaffolds such as `GL*` and `KI*` before or between the expected main
    /// chromosome blocks on hg38 outputs. For downstream inspection we want all
    /// `chr*` contigs first, with canonical chromosomes ordered naturally
    /// (`chr1..chr22, chrX, chrY, chrM/chrMT`), then the remaining `chr*`
    /// contigs, and only after that non-`chr` scaffolds.
    fn compare_chr_names(a: &str, b: &str) -> Ordering {
        fn chr_rank(chr: &str) -> (u8, u8, i32, &str) {
            if let Some(rest) = chr.strip_prefix("chr") {
                if let Ok(n) = rest.parse::<i32>() {
                    return (0, 0, n, "");
                }
                return match rest {
                    "X" => (0, 1, 23, ""),
                    "Y" => (0, 1, 24, ""),
                    "M" | "MT" => (0, 1, 25, ""),
                    _ => (0, 2, 0, rest),
                };
            }
            (1, 0, 0, chr)
        }

        let ka = chr_rank(a);
        let kb = chr_rank(b);
        ka.0.cmp(&kb.0)
            .then_with(|| ka.1.cmp(&kb.1))
            .then_with(|| ka.2.cmp(&kb.2))
            .then_with(|| ka.3.cmp(kb.3))
    }

    /// Mirrors Java Summary's fallback annotation pass for non-exact exon
    /// boundary matches.
    ///
    /// Java first checks whether both circRNA ends hit the same exon boundary
    /// pair exactly. If not, it falls back to scanning gene spans on the same
    /// chromosome and labels the circRNA as `exon`, `intron`, or
    /// `intergenic_region` based on whether both ends fall inside exon intervals.
    fn annotate_circ(annotation: &Annotation, chr: &str, start: i32, end: i32) -> (String, String) {
        if let (Some(v1), Some(v2)) = (
            annotation
                .chr_exon_start_map
                .get(&format!("{}\t{}", chr, start)),
            annotation
                .chr_exon_end_map
                .get(&format!("{}\t{}", chr, end)),
        ) {
            if let (Some((g1, _)), Some((g2, _))) = (v1.split_once('\t'), v2.split_once('\t')) {
                if g1 == g2 {
                    return ("exon".to_string(), g1.to_string());
                }
            }
        }

        if let Some(genes) = annotation.chr_gene_map.get(chr) {
            for gene in genes {
                if start < gene.start {
                    break;
                }
                if end > gene.end {
                    continue;
                }
                if let Some(exons) = annotation.gene_exon_map.get(&gene.gene_id) {
                    let start_in_exon = exons
                        .iter()
                        .any(|&(exon_start, exon_end)| exon_start <= start && exon_end >= start);
                    let end_in_exon = exons
                        .iter()
                        .any(|&(exon_start, exon_end)| exon_start <= end && exon_end >= end);
                    if start_in_exon && end_in_exon {
                        return ("exon".to_string(), gene.gene_id.clone());
                    }
                    return ("intron".to_string(), gene.gene_id.clone());
                }
            }
        }

        ("intergenic_region".to_string(), "NA".to_string())
    }

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

    /// Creates the Summary stage with a Java-compatible stringency level.
    pub fn new(stringency: i32) -> Self {
        Self {
            stringency,
            circ_count: 0,
            final_bsj_reads: 0,
        }
    }

    /// Extracts the total genomic span covered between the outermost `M`
    /// operations of a split CIGAR fragment.
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

    /// Updates the `[SM, MS, SMS]`-style CIGAR pattern counters used in the final
    /// report.
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

    /// Merges circRNA groups that share the same start site when sequence support
    /// indicates they are equivalent.
    ///
    /// The Java-style HashSet bucket reconstruction is intentional: the merge
    /// winner depends on traversal order, so parity requires deterministic
    /// emulation of Java iteration.
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

            let up: Vec<usize> = score
                .iter()
                .enumerate()
                .filter_map(|(i, &v)| if v > 0 { Some(i) } else { None })
                .collect();
            let mut zero: Vec<usize> = score
                .iter()
                .enumerate()
                .filter_map(|(i, &v)| if v == 0 { Some(i) } else { None })
                .collect();

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
                    let pass_score =
                        aligner.score >= seq1.len() as i32 - 2 - ((seq1.len() as i32 - 1) / 10) * 2;
                    let pass_len = aligner.aligned_len > seq1.len() as i32 - 2;
                    if pass_score && pass_len {
                        let up_name = &circ_names[up_idx];
                        let z_name = &circ_names[z_idx];
                        if let (Some(set_up), Some(set_z)) = (
                            circ_map.get(up_name).cloned(),
                            circ_map.get(z_name).cloned(),
                        ) {
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

    /// Merges circRNA groups that share the same end site under the same sequence
    /// similarity rule used for `merge_same_start`.
    fn merge_same_end(
        circ_map: &mut HashMap<String, HashSet<String>>,
        chr_tcga_map: &HashMap<String, String>,
        circ_insertion: &[String],
        circ_map_capacity: usize,
    ) {
        let mut circ_end_seen: HashMap<String, HashSet<String>> = HashMap::new();
        let mut circ_end_insertion: HashMap<String, Vec<String>> = HashMap::new();
        let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); circ_map_capacity.max(1)];
        for (idx, key) in circ_insertion.iter().enumerate() {
            let spread = Self::java_hash_spread(Self::java_string_hash(key));
            let bucket_idx = (spread as usize) & (circ_map_capacity.max(1) - 1);
            buckets[bucket_idx].push(idx);
        }
        for bucket in buckets {
            for idx in bucket {
                let k = &circ_insertion[idx];
                if !circ_map.contains_key(k) {
                    continue;
                }
                let p: Vec<&str> = k.split('\t').collect();
                if p.len() != 3 {
                    continue;
                }
                let key = format!("{}\t{}", p[0], p[2]);
                let seen = circ_end_seen.entry(key.clone()).or_default();
                if seen.insert(k.clone()) {
                    circ_end_insertion.entry(key).or_default().push(k.clone());
                }
            }
        }

        let mut aligner = SmithWaterman::new(1, -1, -3);
        for (chr_end, insertion_list) in circ_end_insertion {
            if insertion_list.len() <= 1 {
                continue;
            }
            // Java parity: circEndMap stores HashSet<String>. Reconstruct the
            // same bucket-order iteration that Java's HashSet would expose for
            // this same-end group, because the first positive-score circ that
            // encounters a zero-score circ wins the merge.
            let cap = Self::java_hashset_capacity(insertion_list.len());
            let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); cap];
            for (idx, s) in insertion_list.iter().enumerate() {
                let spread = Self::java_hash_spread(Self::java_string_hash(s));
                let bucket_idx = (spread as usize) & (cap - 1);
                buckets[bucket_idx].push(idx);
            }
            let circ_names: Vec<String> = buckets
                .into_iter()
                .flat_map(|bucket| bucket.into_iter())
                .map(|idx| insertion_list[idx].clone())
                .filter(|k| circ_map.contains_key(k))
                .collect();
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

            let up: Vec<usize> = score
                .iter()
                .enumerate()
                .filter_map(|(i, &v)| if v > 0 { Some(i) } else { None })
                .collect();
            let mut zero: Vec<usize> = score
                .iter()
                .enumerate()
                .filter_map(|(i, &v)| if v == 0 { Some(i) } else { None })
                .collect();

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
                    let pass_score =
                        aligner.score >= seq1.len() as i32 - 2 - ((seq1.len() as i32 - 1) / 10) * 2;
                    let pass_len = aligner.aligned_len > seq1.len() as i32 - 2;
                    if pass_score && pass_len {
                        let up_name = &circ_names[up_idx];
                        let z_name = &circ_names[z_idx];
                        if let (Some(set_up), Some(set_z)) = (
                            circ_map.get(up_name).cloned(),
                            circ_map.get(z_name).cloned(),
                        ) {
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

    /// Runs the final Summary stage and writes the requested report file.
    ///
    /// Inputs are assumed to already be parity-aligned with Java at the BSJ1/FSJ
    /// level; this stage only performs Java-compatible merging, stringency
    /// filtering, annotation labeling, and output ordering.
    pub fn run(
        &mut self,
        bsj1_file: &str,
        result_path: &str,
        fsj_map: &HashMap<String, i32>,
        chr_tcga_map: &HashMap<String, String>,
        annotation: &Annotation,
    ) -> Result<()> {
        self.run_from_bsj_files(&[bsj1_file], result_path, fsj_map, chr_tcga_map, annotation)
    }

    /// Runs Summary over one or more BSJ protocol files in the supplied order.
    ///
    /// Mate-level BSJ rows carry `mate_label` and `priority`; only `priority=1`
    /// is converted back to the legacy payload consumed by the Java-compatible
    /// merge and stringency logic.
    pub fn run_from_bsj_files(
        &mut self,
        bsj_files: &[&str],
        result_path: &str,
        fsj_map: &HashMap<String, i32>,
        chr_tcga_map: &HashMap<String, String>,
        annotation: &Annotation,
    ) -> Result<()> {
        self.circ_count = 0;
        self.final_bsj_reads = 0;
        let mut circ_map: HashMap<String, HashSet<String>> = HashMap::new();
        let mut circ_start_seen: HashMap<String, HashSet<String>> = HashMap::new();
        let mut circ_start_insertion: HashMap<String, Vec<String>> = HashMap::new();
        let mut circ_seen: HashSet<String> = HashSet::new();
        let mut circ_insertion: Vec<String> = Vec::new();
        for bsj_file in bsj_files {
            let file = File::open(bsj_file)?;
            let reader = BufReader::new(file);
            for line_res in reader.lines() {
                let line = line_res?;
                let p: Vec<&str> = line.split('\t').collect();
                if p.len() < 8 || !bsj_is_summary_priority(&p) {
                    continue;
                }
                let payload_start = bsj_payload_start(&p);
                if p.len() <= payload_start + 8 {
                    continue;
                }
                let legacy_line = if payload_start == 1 {
                    line.clone()
                } else {
                    format!("{}\t{}", p[0], p[payload_start..].join("\t"))
                };
                let legacy: Vec<&str> = legacy_line.split('\t').collect();
                let key = format!("{}\t{}\t{}", legacy[3], legacy[4], legacy[5]);
                circ_map
                    .entry(key.clone())
                    .or_default()
                    .insert(legacy_line.clone());
                if circ_seen.insert(key.clone()) {
                    circ_insertion.push(key.clone());
                }
                let start_key = format!("{}\t{}", legacy[3], legacy[4]);
                let seen = circ_start_seen.entry(start_key.clone()).or_default();
                if seen.insert(key.clone()) {
                    circ_start_insertion.entry(start_key).or_default().push(key);
                }
            }
        }

        Self::merge_same_start(&mut circ_map, chr_tcga_map, &circ_start_insertion);
        let circ_map_capacity = Self::java_hashset_capacity(circ_insertion.len());
        Self::merge_same_end(
            &mut circ_map,
            chr_tcga_map,
            &circ_insertion,
            circ_map_capacity,
        );

        let mut final_results: HashMap<String, Vec<CircSortItem>> = HashMap::new();
        let mut final_read_ids: HashSet<String> = HashSet::new();
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
                ((tp_reads > 19 * fp_reads || fp_reads <= 1)
                    && tp_reads > non_reads + fp_reads
                    && cigar_set.len() >= 3
                    && tp_reads >= 2)
                    || (tag > 0
                        && false_cigar_set3.is_empty()
                        && cigar_set3.len() >= 3
                        && tp_reads3 >= 2)
            } else if self.stringency == 1 {
                ((tp_reads > 19 * fp_reads || false_cigar_set.len() <= 2)
                    && tp_reads > non_reads + fp_reads
                    && tp_reads >= 2)
                    || (tag > 0 && false_cigar_set3.is_empty() && tp_reads3 >= 2)
            } else {
                ((tp_reads > 19 * fp_reads || false_cigar_set.len() <= 2)
                    && tp_reads > non_reads + fp_reads)
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
            let (circ_type, gene_id) =
                Self::annotate_circ(annotation, p_key[0], start, p_key[2].parse().unwrap_or(0));
            let mut ids: Vec<String> = circ_id_set3.into_iter().collect();
            ids.sort();
            for id in &ids {
                final_read_ids.insert(id.clone());
            }
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
            final_results
                .entry(p_key[0].to_string())
                .or_default()
                .push(CircSortItem {
                    start_site: start,
                    line,
                });
        }

        let out_file = File::create(result_path)?;
        let mut writer = BufWriter::with_capacity(1024 * 1024, out_file);
        writeln!(writer, "circRNA_ID\tchr\tcircRNA_start\tcircRNA_end\t#junction_reads\tSM_MS_SMS\t#non_junction_reads\tjunction_reads_ratio\tcircRNA_type\tgene_id\tstrand\tjunction_reads_ID\tScore")?;
        let mut chrs: Vec<String> = final_results.keys().cloned().collect();
        chrs.sort_by(|a, b| Self::compare_chr_names(a, b));
        for chr in chrs {
            let items = final_results.get_mut(&chr).expect("chromosome key exists");
            items.sort_by_key(|x| x.start_site);
            for item in items {
                writeln!(writer, "{}", item.line)?;
                self.circ_count += 1;
            }
        }
        writer.flush()?;
        self.final_bsj_reads = final_read_ids.len();
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::Summary;

    #[test]
    fn test_compare_chr_names_main_chromosomes_first() {
        let mut chrs = vec![
            "GL000225.1".to_string(),
            "chr10".to_string(),
            "chr2".to_string(),
            "KI270442.1".to_string(),
            "chrX".to_string(),
            "chr1".to_string(),
            "chrM".to_string(),
            "chrUn_KI270442v1".to_string(),
        ];
        chrs.sort_by(|a, b| Summary::compare_chr_names(a, b));
        assert_eq!(
            chrs,
            vec![
                "chr1",
                "chr2",
                "chr10",
                "chrX",
                "chrM",
                "chrUn_KI270442v1",
                "GL000225.1",
                "KI270442.1",
            ]
        );
    }
}
