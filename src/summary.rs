/// Summary module: Implements 1:1 pixel-level logic with global Read-ID uniqueness.
/// Ensures each read is assigned to exactly one circRNA based on a rigorous priority hierarchy.

use std::collections::{HashMap, HashSet, BTreeMap};
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};
use anyhow::Result;
use crate::annotation::Annotation;

pub struct Summary {
    pub stringency: i32,
}

struct CircSortItem {
    start_site: i32,
    data: CircOutputData,
}

#[derive(Clone)]
struct CircOutputData {
    id: String,
    chr: String,
    start: i32,
    end: i32,
    tp_reads_3: i32,
    cigar_count_3: [i32; 3],
    fsj: i32,
    ratio: f64,
    circ_type: String,
    gene_id: String,
    strand: String,
    junction_reads_id: String,
    tag_val: i32,
}

/// Internal structure to rank multiple evidences for the same Read ID
#[derive(Debug, Clone)]
struct ReadEvidence {
    line: String,
    _site_key: String,
    tag_type: i32,    // 1 (Scan1), 3 (Scan2)
    signal_type: i32, // 1 (AG/GT), 2 (CT/AC), 0 (Other)
    sum_q: i32,
    m_len: i32,
}

impl Summary {
    pub fn new(stringency: i32) -> Self {
        Self {
            stringency,
        }
    }

    pub fn run(&mut self, bsj1_file: &str, out_prefix: &str, fsj_map: &HashMap<String, i32>, _chr_tcga_map: &HashMap<String, String>, annotation: &Annotation) -> Result<()> {
        let mut read_evidences: HashMap<String, Vec<ReadEvidence>> = HashMap::new();

        // 1. Load all evidence and group by Read ID
        let file = File::open(bsj1_file)?;
        let reader = BufReader::new(file);
        for line_res in reader.lines() {
            let line = line_res?;
            let p: Vec<&str> = line.split('\t').collect();
            if p.len() < 10 { continue; }
            
            let read_id = p[0].to_string();
            let site_key = format!("{}\t{}\t{}", p[3], p[4], p[5]);
            let tag_type = p[2].parse::<i32>().unwrap_or(0);
            
            // For Scan1 (tag_type 1), sumQ is in p[10] if available, otherwise default to 0.
            // For Scan2 (tag_type 3), sumQ is not explicitly passed in BSJ1, default to 0.
            let sum_q = if p.len() > 10 { p[10].parse::<i32>().unwrap_or(0) } else { 0 };
            
            // signal_type is in p[9]
            let signal_type = p[9].parse::<i32>().unwrap_or(0);
            
            // Calculate total M length as a tie-breaker
            let m_len: i32 = p[1].chars().filter(|c| c.is_digit(10)).collect::<String>().parse().unwrap_or(0); 

            read_evidences.entry(read_id).or_insert_with(Vec::new).push(ReadEvidence {
                line, _site_key: site_key, tag_type, signal_type, sum_q, m_len
            });
        }

        // 2. Resolve Ambiguity: Each Read ID gets only ONE "Best" candidate site
        let mut unique_bsj_lines: Vec<String> = Vec::new();
        for (_, mut evs) in read_evidences {
            // Rank: 1. Scan1 > Scan2, 2. Signal AG/GT > CT/AC, 3. SumQ=1 > 0, 4. Longer M
            evs.sort_by(|a, b| {
                let a_type_score = if a.tag_type == 1 { 2 } else { 1 };
                let b_type_score = if b.tag_type == 1 { 2 } else { 1 };
                
                b_type_score.cmp(&a_type_score)
                    .then(b.signal_type.cmp(&a.signal_type))
                    .then(b.sum_q.cmp(&a.sum_q))
                    .then(b.m_len.cmp(&a.m_len))
            });
            unique_bsj_lines.push(evs[0].line.clone());
        }

        // 3. Populate statistics based on uniquely assigned reads
        let mut circ_map: HashMap<String, HashSet<String>> = HashMap::new();
        for line in unique_bsj_lines {
            let p: Vec<&str> = line.split('\t').collect();
            let site_key = format!("{}\t{}\t{}", p[3], p[4], p[5]);
            circ_map.entry(site_key).or_insert_with(HashSet::new).insert(line);
        }

        // 4. Final Filtering (Mirroring Summary.java logic)
        let mut final_results: BTreeMap<String, Vec<CircSortItem>> = BTreeMap::new();
        for (chr_start_end, lines) in &circ_map {
            let p_key: Vec<&str> = chr_start_end.split('\t').collect();
            
            let mut circ_id_set = HashSet::new();
            let mut circ_id_set_3 = HashSet::new();
            let mut cigar_set = HashSet::new();
            let mut cigar_set_3 = HashSet::new();
            let mut false_cigar_set = HashSet::new();
            let mut false_cigar_set_3 = HashSet::new();
            let mut fp_reads = 0;
            let mut non_reads = 0;
            let mut tag = 0;
            let mut strand = " ";

            for line in lines {
                let p: Vec<&str> = line.split('\t').collect();
                strand = p[6];
                let tag_type = p[2];
                let rescue_flag = p[9];
                
                if tag_type == "1" { // Scan1 TP
                    circ_id_set.insert(p[0].to_string());
                    circ_id_set_3.insert(p[0].to_string());
                    for c in p[1].split(';') { cigar_set.insert(c.to_string()); cigar_set_3.insert(c.to_string()); }
                    if rescue_flag == "1" { tag += 1; }
                } else if tag_type == "3" { // Scan2 Rescue TP
                    circ_id_set_3.insert(p[0].to_string());
                    cigar_set_3.insert(p[1].to_string());
                    if rescue_flag == "1" { tag += 1; }
                } else if tag_type == "-1" {
                    non_reads += 1;
                    for c in p[1].split(';') { false_cigar_set.insert(c.to_string()); }
                } else if tag_type == "-2" {
                    fp_reads += 1;
                    for c in p[1].split(';') { false_cigar_set.insert(c.to_string()); }
                } else if tag_type.starts_with("-") {
                    false_cigar_set_3.insert(p[1].to_string());
                }
            }

            let tp_reads = circ_id_set.len() as i32;
            let tp_reads_3 = circ_id_set_3.len() as i32;
            
            let passed = match self.stringency {
                2 => ((tp_reads > 19 * fp_reads || fp_reads <= 1) && tp_reads > (non_reads + fp_reads) && cigar_set.len() >= 3 && tp_reads >= 2) ||
                     (tag > 0 && false_cigar_set_3.is_empty() && cigar_set_3.len() >= 3 && tp_reads_3 >= 2),
                1 => ((tp_reads > 19 * fp_reads || false_cigar_set.len() <= 2) && tp_reads > (non_reads + fp_reads) && tp_reads >= 2) ||
                     (tag > 0 && false_cigar_set_3.is_empty() && tp_reads_3 >= 2),
                _ => ((tp_reads > 19 * fp_reads || false_cigar_set.len() <= 2) && tp_reads > (non_reads + fp_reads)) ||
                     (tag > 0 && false_cigar_set_3.is_empty() && tp_reads_3 >= 2),
            };

            if passed {
                let mut cigar_count_3 = [0, 0, 0];
                for c_str in &cigar_set_3 { self.update_cigar_counts(c_str, &mut cigar_count_3); }
                let fsj = *fsj_map.get(chr_start_end).unwrap_or(&0);
                let ratio = if (tp_reads_3 * 2 + fsj) > 0 { (tp_reads_3 * 2) as f64 / (tp_reads_3 * 2 + fsj) as f64 } else { 0.0 };
                let start_p = p_key[1].parse::<i32>().unwrap_or(0);
                let end_p = p_key[2].parse::<i32>().unwrap_or(0);
                
                let mut circ_type = "intergenic_region";
                let mut gene_id = "NA";
                if let (Some(v1), Some(v2)) = (annotation.chr_exon_start_map.get(&format!("{}\t{}", p_key[0], p_key[1])), 
                                               annotation.chr_exon_end_map.get(&format!("{}\t{}", p_key[0], p_key[2]))) {
                    let (g1, _) = v1.split_once('\t').unwrap();
                    let (g2, _) = v2.split_once('\t').unwrap();
                    if g1 == g2 { circ_type = "exon"; gene_id = g1; }
                }
                let mut ids: Vec<String> = circ_id_set_3.into_iter().collect();
                ids.sort();
                final_results.entry(p_key[0].to_string()).or_insert_with(Vec::new).push(CircSortItem {
                    start_site: start_p,
                    data: CircOutputData {
                        id: format!("{}:{}|{}", p_key[0], p_key[1], p_key[2]), chr: p_key[0].to_string(), start: start_p, end: end_p,
                        tp_reads_3, cigar_count_3, fsj, ratio, circ_type: circ_type.to_string(), gene_id: gene_id.to_string(),
                        strand: strand.to_string(), junction_reads_id: ids.join(","), tag_val: tag
                    }
                });
            }
        }

        // 5. Final Output
        let out_file = File::create(format!("{}.result", out_prefix))?;
        let mut writer = BufWriter::new(out_file);
        writeln!(writer, "circRNA_ID\tchr\tcircRNA_start\tcircRNA_end\t#junction_reads\tSM_MS_SMS\t#non_junction_reads\tjunction_reads_ratio\tcircRNA_type\tgene_id\tstrand\tjunction_reads_ID\tScore")?;
        for items in final_results.values_mut() {
            items.sort_by_key(|it| it.start_site);
            for it in items {
                let d = &it.data;
                writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}_{}_{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}", d.id, d.chr, d.start, d.end, d.tp_reads_3, d.cigar_count_3[0], d.cigar_count_3[1], d.cigar_count_3[2], d.fsj, d.ratio, d.circ_type, d.gene_id, d.strand, d.junction_reads_id, d.tag_val)?;
            }
        }
        Ok(())
    }

    fn update_cigar_counts(&self, cigar: &str, counts: &mut [i32; 3]) {
        let t: String = cigar.chars().filter(|x| x.is_alphabetic()).collect();
        if t == "M" { return; }
        else if t.find('M') > Some(0) && t.rfind('M') < Some(t.len() - 1) { counts[2] += 1; }
        else if t.find('M') == Some(0) { counts[1] += 1; } else { counts[0] += 1; }
    }
}
