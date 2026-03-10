//! Scan 2 module: Parallel multi-format sequence rescue and FSJ counting.
//!
//! This module implements the second pass of identification, focusing on rescuing 
//! PEM (Paired-End Mapping) and SMS signals that were missed in Scan 1.
//! It also quantifies Forward-Spliced Junctions (FSJ) for ratio calculation.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};
use std::borrow::Cow;
use anyhow::Result;
use rayon::prelude::*;
use memmap2::Mmap;
use memchr::memchr;
use noodles::sam::{self, alignment::Record as _};
use indicatif::{ProgressBar, ProgressStyle};
use crate::misd::misd;
use crate::is_bsj_hg2::IsBSJHg2;
use crate::utils::{reverse_complement, AlignmentRecord};

/// Core logic for the second scan pass.
pub struct Scan2 {
    /// Minimum mapping quality for a read.
    pub min_mapq_uni: i32,
    /// Minimum linear range size for competition check.
    pub linear_range_size_min: i32,
    /// Index for BSJ candidates (indexed by Site1).
    pub index1: HashMap<String, Vec<CandidateBreakpoint>>,
    /// Index for BSJ candidates (indexed by Site2).
    pub index2: HashMap<String, Vec<CandidateBreakpoint>>,
    /// Map to store Forward Spliced Junction (FSJ) counts.
    pub fsj_map: HashMap<String, i32>,
}

/// Stores a potential breakpoint site and its associated Scan 1 metadata.
#[derive(Clone)]
pub struct CandidateBreakpoint {
    /// Genomic coordinate of the breakpoint.
    pub site: i32,
    /// Full tab-separated metadata from BSJ1 file.
    pub data: Vec<String>,
}

/// Fast integer parser from bytes.
#[inline]
fn fast_parse_i32(bytes: &[u8]) -> i32 {
    let mut res = 0;
    for &b in bytes { if b >= b'0' && b <= b'9' { res = res * 10 + (b - b'0') as i32; } }
    res
}

impl Scan2 {
    /// Creates a new `Scan2` instance with specified thresholds.
    pub fn new(min_mapq_uni: i32, linear_range_size_min: i32, _seq_len: i32) -> Self {
        Self { min_mapq_uni, linear_range_size_min, index1: HashMap::new(), index2: HashMap::new(), fsj_map: HashMap::new() }
    }

    /// Loads BSJ1 candidates and initializes the FSJ map.
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

    /// Entry point for Scan 2. Automatically detects file format.
    pub fn run(&mut self, sam_file: &str, output_bsj2: &str, chr_tcga_map: &HashMap<String, String>) -> Result<()> {
        use crate::sam_bam::{detect_format, InputFormat};
        let format = detect_format(sam_file)?;
        
        match format {
            InputFormat::Sam => self.run_sam(sam_file, output_bsj2, chr_tcga_map),
            InputFormat::Bam => self.run_bam(sam_file, output_bsj2, chr_tcga_map),
        }
    }

    /// Parallel runner for SAM files using Mmap.
    pub fn run_sam(&mut self, sam_file: &str, output_bsj2: &str, chr_tcga_map: &HashMap<String, String>) -> Result<()> {
        let file = File::open(sam_file)?;
        let mmap = unsafe { Mmap::map(&file)? };
        let file_size = mmap.len();
        let num_threads = rayon::current_num_threads().max(1);
        let shard_size = mmap.len() / num_threads;

        let pb = ProgressBar::new(file_size as u64);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta}) {msg}")?
            .progress_chars("#>-"));
        pb.set_message("Scan 2: Rescuing signals & counting FSJ");

        let shard_results: Vec<(Vec<String>, HashMap<String, i32>)> = (0..num_threads).into_par_iter().map(|i| {
            let start = i * shard_size;
            let end = if i == num_threads - 1 { mmap.len() } else { (i + 1) * shard_size };
            self.process_sam_shard(&mmap, start, end, chr_tcga_map, &pb).unwrap_or_default()
        }).collect();

        pb.finish_with_message("Scan 2: Completed");

        let out_file = std::fs::OpenOptions::new().create(true).append(true).open(output_bsj2)?;
        let mut writer = BufWriter::with_capacity(1024 * 1024, out_file);
        for (lines, partial_fsj) in shard_results {
            for line in lines { writeln!(writer, "{}", line)?; }
            for (key, count) in partial_fsj { *self.fsj_map.entry(key).or_insert(0) += count; }
        }
        writer.flush()?;
        Ok(())
    }

    /// Parallel runner for BAM files.
    pub fn run_bam(&mut self, bam_file: &str, output_bsj2: &str, chr_tcga_map: &HashMap<String, String>) -> Result<()> {
        use noodles::bam;
        let file = File::open(bam_file)?;
        let mmap = unsafe { Mmap::map(&file)? };
        let file_size = mmap.len();
        
        let mut header_reader = bam::io::Reader::new(&mmap[..]);
        let header = header_reader.read_header()?;
        
        let num_threads = rayon::current_num_threads().max(1);
        let shard_size = mmap.len() / num_threads;

        let pb = ProgressBar::new(file_size as u64);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta}) {msg}")?
            .progress_chars("#>-"));
        pb.set_message("Scan 2: Rescuing signals & counting FSJ (BAM)");

        let shard_results: Vec<(Vec<String>, HashMap<String, i32>)> = (0..num_threads).into_par_iter().map(|i| {
            let start = i * shard_size;
            let end = if i == num_threads - 1 { mmap.len() } else { (i + 1) * shard_size };
            self.process_bam_shard(&mmap, start, end, &header, chr_tcga_map, &pb).unwrap_or_default()
        }).collect();

        pb.finish_with_message("Scan 2: Completed");

        let out_file = std::fs::OpenOptions::new().create(true).append(true).open(output_bsj2)?;
        let mut writer = BufWriter::with_capacity(1024 * 1024, out_file);
        for (lines, partial_fsj) in shard_results {
            for line in lines { writeln!(writer, "{}", line)?; }
            for (key, count) in partial_fsj { *self.fsj_map.entry(key).or_insert(0) += count; }
        }
        writer.flush()?;
        Ok(())
    }

    /// Processes a BAM shard.
    fn process_bam_shard(&self, mmap: &[u8], start: usize, end: usize, header: &sam::Header, chr_tcga_map: &HashMap<String, String>, pb: &ProgressBar) -> Result<(Vec<String>, HashMap<String, i32>)> {
        use noodles::bam;
        let mut results = Vec::new();
        let mut local_fsj = HashMap::new();
        let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
        
        let pos = if start == 0 { 0 } else {
            let mut found = None;
            for i in start..end {
                if i + 3 < mmap.len() && &mmap[i..i+4] == b"\x1f\x8b\x08\x04" {
                    found = Some(i); break;
                }
            }
            found.unwrap_or(mmap.len())
        };
        if pos >= mmap.len() { return Ok((results, local_fsj)); }

        let mut reader = bam::io::Reader::new(&mmap[pos..]);
        let mut record = bam::Record::default();
        
        let mut current_id: Vec<u8> = Vec::new();
        let mut alignments: Vec<AlignmentRecord> = Vec::with_capacity(16);
        let mut stand_map: HashMap<i32, (char, Cow<str>)> = HashMap::with_capacity(4);
        let mut first_id_skipped = start == 0;
        let mut last_compressed_pos = 0;

        while reader.read_record(&mut record)? != 0 {
            let current_compressed_pos = reader.get_ref().virtual_position().compressed() as usize;
            pb.inc((current_compressed_pos - last_compressed_pos) as u64);
            last_compressed_pos = current_compressed_pos;

            let read_id = record.name().ok_or_else(|| anyhow::anyhow!("Missing read name"))?;
            
            if !first_id_skipped {
                if current_id.is_empty() { current_id = read_id.to_vec(); continue; }
                if read_id.to_vec() != current_id { first_id_skipped = true; current_id = read_id.to_vec(); } else { continue; }
            }

            if read_id.to_vec() != current_id {
                if !current_id.is_empty() {
                    let id_str = String::from_utf8_lossy(&current_id);
                    self.process_group_view(&id_str, &alignments, &stand_map, &mut results, &mut local_fsj, chr_tcga_map, &mut validator)?;
                    if pos + current_compressed_pos > end { current_id.clear(); break; }
                }
                current_id = read_id.to_vec(); alignments.clear(); stand_map.clear();
            }

            let flag = i32::from(u16::from(record.flags()));
            let chrom = match record.reference_sequence(header) { Some(Ok((name, _))) => String::from_utf8_lossy(name).to_string(), _ => "*".to_string() };
            let start_pos = record.alignment_start().transpose()?.map(|p| p.get() as i32).unwrap_or(0);
            let mapq = record.mapping_quality().map(u8::from).unwrap_or(0) as i32;
            
            let mut cigar = String::new();
            for res in record.cigar().iter() {
                let op = res?; use noodles::sam::alignment::record::cigar::op::Kind;
                let op_c = match op.kind() {
                    Kind::Match => 'M', Kind::Insertion => 'I', Kind::Deletion => 'D',
                    Kind::Skip => 'N', Kind::SoftClip => 'S', Kind::HardClip => 'H',
                    Kind::Pad => 'P', Kind::SequenceMatch => '=', Kind::SequenceMismatch => 'X',
                };
                cigar.push_str(&format!("{}{}", op.len(), op_c));
            }
            let mut seq = String::new(); for b in record.sequence().iter() { seq.push(char::from(b)); }

            if !seq.is_empty() && seq != "*" {
                let s_idx = if flag & 0x40 != 0 { 0 } else { 1 };
                let st_c = if flag & 0x10 != 0 { '1' } else { '0' };
                let entry = stand_map.entry(s_idx).or_insert((st_c, Cow::Owned(seq.clone())));
                if seq.len() > entry.1.len() { *entry = (st_c, Cow::Owned(seq.clone())); }
            }
            alignments.push(AlignmentRecord { flag, chrom: Cow::Owned(chrom), pos: start_pos, mapq, cigar: Cow::Owned(cigar), seq: Cow::Owned(seq) });
        }

        if !current_id.is_empty() {
            let id_str = String::from_utf8_lossy(&current_id);
            self.process_group_view(&id_str, &alignments, &stand_map, &mut results, &mut local_fsj, chr_tcga_map, &mut validator)?;
        }
        Ok((results, local_fsj))
    }

    /// Processes a SAM shard using Mmap.
    fn process_sam_shard<'a>(&self, mmap: &'a [u8], start: usize, end: usize, chr_tcga_map: &HashMap<String, String>, pb: &ProgressBar) -> Result<(Vec<String>, HashMap<String, i32>)> {
        let mut results = Vec::new();
        let mut local_fsj = HashMap::new();
        let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
        let mut pos = if start == 0 { 0 } else { memchr(b'\n', &mmap[start..]).map(|p| start + p + 1).unwrap_or(end) };
        if pos >= end && start != 0 { return Ok((results, local_fsj)); }

        let mut current_id: &[u8] = &[];
        let mut alignments: Vec<AlignmentRecord<'a>> = Vec::with_capacity(16);
        let mut stand_map: HashMap<i32, (char, Cow<'a, str>)> = HashMap::with_capacity(4);

        while pos < mmap.len() {
            let line_end = memchr(b'\n', &mmap[pos..]).map(|p| pos + p).unwrap_or(mmap.len());
            let line = &mmap[pos..line_end];
            if line.is_empty() { 
                pos = line_end + 1; 
                pb.inc(1);
                continue; 
            }
            if line[0] == b'@' { 
                pos = line_end + 1; 
                pb.inc((line_end - pos + 1) as u64);
                continue; 
            }

            let mut cols = line.split(|&b| b == b'\t');
            let read_id = cols.next().unwrap();
            
            if read_id != current_id {
                if !current_id.is_empty() { self.process_group_view(unsafe { std::str::from_utf8_unchecked(current_id) }, &alignments, &stand_map, &mut results, &mut local_fsj, chr_tcga_map, &mut validator)?; }
                if pos >= end { break; }
                current_id = read_id; alignments.clear(); stand_map.clear();
            }
            
            let flag = fast_parse_i32(cols.next().unwrap_or(b"0"));
            let chrom = unsafe { std::str::from_utf8_unchecked(cols.next().unwrap_or(b"*")) };
            let start_pos = fast_parse_i32(cols.next().unwrap_or(b"0"));
            let mapq = fast_parse_i32(cols.next().unwrap_or(b"0"));
            let cigar = unsafe { std::str::from_utf8_unchecked(cols.next().unwrap_or(b"*")) };
            cols.next(); cols.next(); cols.next();
            let seq = unsafe { std::str::from_utf8_unchecked(cols.next().unwrap_or(b"*")).trim() };
            
            if seq != "*" {
                let s_idx = if flag & 0x40 != 0 { 0 } else { 1 };
                let st_c = if flag & 0x10 != 0 { '1' } else { '0' };
                let entry = stand_map.entry(s_idx).or_insert((st_c, Cow::Borrowed(seq)));
                if seq.len() > entry.1.len() { *entry = (st_c, Cow::Borrowed(seq)); }
            }
            alignments.push(AlignmentRecord { flag, chrom: Cow::Borrowed(chrom), pos: start_pos, mapq, cigar: Cow::Borrowed(cigar), seq: Cow::Borrowed(seq) });
            pb.inc((line_end - pos + 1) as u64);
            pos = line_end + 1;
        }
        Ok((results, local_fsj))
    }

    /// Common identification logic for a group of alignments in Scan 2.
    pub(crate) fn process_group_view<'a>(&self, id: &str, alignments: &[AlignmentRecord<'a>], stand_map: &HashMap<i32, (char, Cow<'a, str>)>, results: &mut Vec<String>, local_fsj: &mut HashMap<String, i32>, chr_tcga_map: &HashMap<String, String>, is_bsj_hg2: &mut IsBSJHg2) -> Result<()> {
        let mut segments: HashMap<i32, Vec<&AlignmentRecord<'a>>> = HashMap::new();
        for aln in alignments {
            segments.entry(if aln.flag & 0x40 != 0 { 0 } else { 1 }).or_insert_with(Vec::new).push(aln);
        }
        let mut tem_fsj_keys = HashSet::new();
        for (&seg_idx, seg_alns) in &segments {
            let (read_strand, seq) = match stand_map.get(&seg_idx) { Some(&(st, ref s)) => (st, s.as_ref()), None => continue };
            let slen = seq.len() as i32;
            let mut p_str = String::new(); let mut s2_ok = "0";
            if let Some(&(p_strand, ref p_seq)) = stand_map.get(&(1 - seg_idx)) {
                if p_strand != read_strand { p_str = p_seq.to_string(); } else { p_str = reverse_complement(p_seq); }
                s2_ok = "1";
            }
            for aln in seg_alns {
                let chr = aln.chrom.as_ref(); if !self.index1.contains_key(chr) { continue; }
                let c = misd(&aln.cigar, slen);
                let (aln_flag, aln_pos) = (aln.flag, aln.pos);
                if c[0] == -1 || c[0] == 10 {
                    if let Some(list) = self.index1.get(chr) {
                        for cand in list {
                            if (cand.site - aln_pos).abs() <= 6 {
                                let str_e = if aln_flag & 0x10 != 0 { if read_strand == '0' { reverse_complement(seq) } else { seq.to_string() } } else { if read_strand == '0' { seq.to_string() } else { reverse_complement(seq) } };
                                let e_idx = c[1] + (cand.site - aln_pos);
                                if e_idx > 0 && e_idx <= slen {
                                    let str_f = &str_e[0..e_idx as usize];
                                    let str3 = if c[0] == 10 { let si = slen - c[2]; if si >= 0 && si <= slen { &str_e[si as usize..] } else { "" } } else { "*" };
                                    let circ_c = vec![(if aln_flag & 0x10 != 0 { "1" } else { "0" }).to_string(), chr.to_string(), "sm".to_string(), cand.data[4].clone(), cand.data[5].clone(), str_f.to_string(), p_str.clone(), str3.to_string(), s2_ok.to_string(), cand.data[6].clone(), cand.data[7].clone(), cand.data[8].clone(), aln.mapq.to_string()];
                                    let tag = is_bsj_hg2.is_bsj_hg2(&circ_c, chr_tcga_map.get(chr).unwrap());
                                    if tag == "0" { tem_fsj_keys.insert(format!("{}\t{}\t{}", chr, cand.data[4], cand.data[5])); }
                                    else if tag != "2" { results.push(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", id, aln.cigar, tag.chars().next_back().unwrap(), chr, cand.data[4], cand.data[5], cand.data[6], cand.data[7], cand.data[8], &tag[0..tag.len()-1])); return Ok(()); }
                                }
                            }
                        }
                    }
                }
                if c[0] == 1 || c[0] == 10 {
                    let new_site = aln_pos + c[3] - 1;
                    if let Some(list) = self.index2.get(chr) {
                        for cand in list {
                            if (cand.site - new_site).abs() <= 6 {
                                let str_e = if aln_flag & 0x10 != 0 { if read_strand == '0' { reverse_complement(seq) } else { seq.to_string() } } else { if read_strand == '0' { seq.to_string() } else { reverse_complement(seq) } };
                                let s_idx = if c[0] == 10 { c[1] + c[3] + (cand.site - new_site) } else { c[1] + (cand.site - new_site) };
                                if s_idx >= 0 && s_idx < slen {
                                    let str_f = &str_e[s_idx as usize..];
                                    let str3 = if c[0] == 10 { &str_e[0..c[1] as usize] } else { "*" };
                                    let circ_c = vec![(if aln_flag & 0x10 != 0 { "1" } else { "0" }).to_string(), chr.to_string(), "ms".to_string(), cand.data[4].clone(), cand.data[5].clone(), str_f.to_string(), p_str.clone(), str3.to_string(), s2_ok.to_string(), cand.data[6].clone(), cand.data[7].clone(), cand.data[8].clone(), aln.mapq.to_string()];
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
