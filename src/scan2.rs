/// Scan 2 module: Parallel sharded implementation for PEM rescue and FSJ counting.
/// Uses Mmap and zero-copy view structs to eliminate heap contention.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};
use anyhow::Result;
use rayon::prelude::*;
use memmap2::Mmap;
use memchr::memchr;
use crate::misd::misd;
use crate::is_bsj_hg2::IsBSJHg2;
use crate::utils::reverse_complement;

pub struct Scan2 {
    pub min_mapq_uni: i32,
    pub linear_range_size_min: i32,
    pub index1: HashMap<String, Vec<CandidateBreakpoint>>,
    pub index2: HashMap<String, Vec<CandidateBreakpoint>>,
    pub fsj_map: HashMap<String, i32>,
}

#[derive(Clone)]
pub struct CandidateBreakpoint {
    pub site: i32,
    pub data: Vec<String>,
}

#[derive(Clone)]
struct AlignmentView<'a> {
    flag: i32, chrom: &'a str, pos: i32, mapq: i32, cigar: &'a str, seq: &'a str,
}

#[inline]
fn fast_parse_i32(bytes: &[u8]) -> i32 {
    let mut res = 0;
    for &b in bytes { if b >= b'0' && b <= b'9' { res = res * 10 + (b - b'0') as i32; } }
    res
}

impl Scan2 {
    pub fn new(min_mapq_uni: i32, linear_range_size_min: i32, _seq_len: i32) -> Self {
        Self { min_mapq_uni, linear_range_size_min, index1: HashMap::new(), index2: HashMap::new(), fsj_map: HashMap::new() }
    }

    pub fn build_index(&mut self, bsj1_file: &str) -> Result<()> {
        let file = File::open(bsj1_file)?;
        let reader = BufReader::new(file);
        for line_res in reader.lines() {
            let line = line_res?;
            let p: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
            if p.len() < 10 { continue; }
            let (chr, site1, site2) = (p[3].clone(), p[4].parse().unwrap_or(0), p[5].parse().unwrap_or(0));
            self.fsj_map.entry(format!("{}\t{}\t{}", p[3], p[4], p[5])).or_insert(0);
            self.index1.entry(chr.clone()).or_insert_with(Vec::new).push(CandidateBreakpoint { site: site1, data: p.clone() });
            self.index2.entry(chr).or_insert_with(Vec::new).push(CandidateBreakpoint { site: site2, data: p });
        }
        Ok(())
    }

    pub fn run(&mut self, sam_file: &str, output_bsj2: &str, chr_tcga_map: &HashMap<String, String>) -> Result<()> {
        let file = File::open(sam_file)?;
        let mmap = unsafe { Mmap::map(&file)? };
        let num_threads = rayon::current_num_threads().max(1);
        let shard_size = mmap.len() / num_threads;

        let shard_results: Vec<(Vec<String>, HashMap<String, i32>)> = (0..num_threads).into_par_iter().map(|i| {
            let start = i * shard_size;
            let end = if i == num_threads - 1 { mmap.len() } else { (i + 1) * shard_size };
            self.process_shard_view(&mmap, start, end, chr_tcga_map).unwrap_or_default()
        }).collect();

        let out_file = std::fs::OpenOptions::new().append(true).open(output_bsj2)?;
        let mut writer = BufWriter::with_capacity(1024 * 1024, out_file);
        for (lines, partial_fsj) in shard_results {
            for line in lines { writeln!(writer, "{}", line)?; }
            for (key, count) in partial_fsj { *self.fsj_map.entry(key).or_insert(0) += count; }
        }
        writer.flush()?;
        Ok(())
    }

    fn process_shard_view(&self, mmap: &[u8], start: usize, end: usize, chr_tcga_map: &HashMap<String, String>) -> Result<(Vec<String>, HashMap<String, i32>)> {
        let mut results = Vec::new();
        let mut local_fsj = HashMap::new();
        let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
        let mut pos = if start == 0 { 0 } else { memchr(b'\n', &mmap[start..]).map(|p| start + p + 1).unwrap_or(end) };
        if pos >= end && start != 0 { return Ok((results, local_fsj)); }

        let mut current_id: &[u8] = &[];
        let mut alignments: Vec<AlignmentView> = Vec::with_capacity(16);
        let mut stand_map: HashMap<i32, (char, &str)> = HashMap::with_capacity(4);

        while pos < end || !current_id.is_empty() {
            let line_end = memchr(b'\n', &mmap[pos..]).map(|p| pos + p).unwrap_or(mmap.len());
            let line = &mmap[pos..line_end];
            if line.is_empty() { if pos >= mmap.len() { break; } pos += 1; continue; }
            if line[0] == b'@' { pos = line_end + 1; continue; }

            let mut cols = line.split(|&b| b == b'\t');
            let read_id = cols.next().unwrap();
            
            if read_id != current_id {
                if !current_id.is_empty() {
                    let id_str = unsafe { std::str::from_utf8_unchecked(current_id) };
                    self.process_group_view(id_str, &alignments, &stand_map, &mut results, &mut local_fsj, chr_tcga_map, &mut validator)?;
                }
                if pos >= end { break; }
                current_id = read_id; alignments.clear(); stand_map.clear();
            }
            
            let flag = fast_parse_i32(cols.next().unwrap_or(b"0"));
            let chrom = unsafe { std::str::from_utf8_unchecked(cols.next().unwrap_or(b"*")) };
            let start_pos = fast_parse_i32(cols.next().unwrap_or(b"0"));
            let mapq = fast_parse_i32(cols.next().unwrap_or(b"0"));
            let cigar = unsafe { std::str::from_utf8_unchecked(cols.next().unwrap_or(b"*")) };
            cols.next(); cols.next(); cols.next();
            let seq_bytes = cols.next().unwrap_or(b"*");
            let seq = unsafe { std::str::from_utf8_unchecked(seq_bytes) }.trim();
            
            if seq != "*" {
                let s_idx = if flag & 0x40 != 0 { 0 } else { 1 };
                let strand_char = if flag & 0x10 != 0 { '1' } else { '0' };
                let entry = stand_map.entry(s_idx).or_insert((strand_char, seq));
                if seq.len() > entry.1.len() { *entry = (strand_char, seq); }
            }
            alignments.push(AlignmentView { flag, chrom, pos: start_pos, mapq, cigar, seq });
            pos = line_end + 1;
        }
        Ok((results, local_fsj))
    }

    fn process_group_view(&self, id: &str, alignments: &[AlignmentView], stand_map: &HashMap<i32, (char, &str)>, results: &mut Vec<String>, local_fsj: &mut HashMap<String, i32>, chr_tcga_map: &HashMap<String, String>, is_bsj_hg2: &mut IsBSJHg2) -> Result<()> {
        let mut segments: HashMap<i32, Vec<&AlignmentView>> = HashMap::new();
        for aln in alignments {
            segments.entry(if aln.flag & 0x40 != 0 { 0 } else { 1 }).or_insert_with(Vec::new).push(aln);
        }
        let mut tem_fsj_keys = HashSet::new();
        for (&seg_idx, seg_alns) in &segments {
            let (read_strand, seq) = match stand_map.get(&seg_idx) { Some(&(st, s)) => (st, s), None => continue };
            let slen = seq.len() as i32;
            let mut p_str = String::new(); let mut s2_ok = "0";
            if let Some(&(p_strand, p_seq)) = stand_map.get(&(1 - seg_idx)) {
                if p_strand != read_strand { p_str = p_seq.to_string(); } else { p_str = reverse_complement(p_seq); }
                s2_ok = "1";
            }
            for aln in seg_alns {
                let chr = aln.chrom; if !self.index1.contains_key(chr) { continue; }
                let c = misd(aln.cigar, slen);
                if c[0] == -1 || c[0] == 10 {
                    if let Some(list) = self.index1.get(chr) {
                        for cand in list {
                            if (cand.site - aln.pos).abs() <= 6 {
                                let str_e = if aln.flag & 0x10 != 0 { if read_strand == '0' { reverse_complement(seq) } else { seq.to_string() } } else { if read_strand == '0' { seq.to_string() } else { reverse_complement(seq) } };
                                let e_idx = c[1] + (cand.site - aln.pos);
                                if e_idx > 0 && e_idx <= slen {
                                    let str_f = &str_e[0..e_idx as usize];
                                    let str3 = if c[0] == 10 { let si = slen - c[2]; if si >= 0 && si <= slen { &str_e[si as usize..] } else { "" } } else { "*" };
                                    let circ_c = vec![if aln.flag & 0x10 != 0 { "1" } else { "0" }.to_string(), chr.to_string(), "sm".to_string(), cand.data[4].clone(), cand.data[5].clone(), str_f.to_string(), p_str.clone(), str3.to_string(), s2_ok.to_string(), cand.data[6].clone(), cand.data[7].clone(), cand.data[8].clone(), aln.mapq.to_string()];
                                    let tag = is_bsj_hg2.is_bsj_hg2(&circ_c, chr_tcga_map.get(chr).unwrap());
                                    if tag == "0" { tem_fsj_keys.insert(format!("{}\t{}\t{}", chr, cand.data[4], cand.data[5])); }
                                    else if tag != "2" { results.push(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", id, aln.cigar, tag.chars().next_back().unwrap(), chr, cand.data[4], cand.data[5], cand.data[6], cand.data[7], cand.data[8], &tag[0..tag.len()-1])); return Ok(()); }
                                }
                            }
                        }
                    }
                }
                if c[0] == 1 || c[0] == 10 {
                    let new_site = aln.pos + c[3] - 1;
                    if let Some(list) = self.index2.get(chr) {
                        for cand in list {
                            if (cand.site - new_site).abs() <= 6 {
                                let str_e = if aln.flag & 0x10 != 0 { if read_strand == '0' { reverse_complement(seq) } else { seq.to_string() } } else { if read_strand == '0' { seq.to_string() } else { reverse_complement(seq) } };
                                let s_idx = if c[0] == 10 { c[1] + c[3] + (cand.site - new_site) } else { c[1] + (cand.site - new_site) };
                                if s_idx >= 0 && s_idx < slen {
                                    let str_f = &str_e[s_idx as usize..];
                                    let str3 = if c[0] == 10 { &str_e[0..c[1] as usize] } else { "*" };
                                    let circ_c = vec![if aln.flag & 0x10 != 0 { "1" } else { "0" }.to_string(), chr.to_string(), "ms".to_string(), cand.data[4].clone(), cand.data[5].clone(), str_f.to_string(), p_str.clone(), str3.to_string(), s2_ok.to_string(), cand.data[6].clone(), cand.data[7].clone(), cand.data[8].clone(), aln.mapq.to_string()];
                                    let tag = is_bsj_hg2.is_bsj_hg2(&circ_c, chr_tcga_map.get(chr).unwrap());
                                    if tag == "0" { tem_fsj_keys.insert(format!("{}\t{}\t{}", chr, cand.data[4], cand.data[5])); }
                                    else if tag != "2" { results.push(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", id, aln.cigar, tag.chars().next_back().unwrap(), chr, cand.data[4], cand.data[5], cand.data[6], cand.data[7], cand.data[8], &tag[0..tag.len()-1])); return Ok(()); }
                                }
                            }
                        }
                    }
                }
            }
        }
        for k in tem_fsj_keys { *local_fsj.entry(k).or_insert(0) += 1; }
        Ok(())
    }
}
