/// Scan 1 module: Ultra-high performance Mmap implementation with TB-scale scalability.
/// Implements robust shard synchronization to handle Read IDs spanning shard boundaries.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{Write, BufWriter};
use std::borrow::Cow;
use anyhow::Result;
use rayon::prelude::*;
use memmap2::Mmap;
use memchr::memchr;
use noodles::sam;
use crate::misd::misd;
use crate::is_bsj_hg2::{IsBSJHg2, java_substring};
use crate::annotation::Annotation;
use crate::utils::AlignmentRecord;

pub struct Scan1 {
    pub min_mapq_uni: i32,
    pub max_circle: i32,
    pub min_circle: i32,
    pub linear_range_size_min: i32,
}

#[derive(Debug, Clone)]
struct BSJCandidate {
    result_str: String, tag: i32, sum_q: i32, total_mq: i32, total_m_len: i32, cigar_pair: String,
}

#[inline]
fn fast_parse_i32(bytes: &[u8]) -> i32 {
    let mut res = 0;
    for &b in bytes { if b >= b'0' && b <= b'9' { res = res * 10 + (b - b'0') as i32; } }
    res
}

impl Scan1 {
    pub fn new(min_mapq_uni: i32, min_circle: i32, max_circle: i32, linear_range_size_min: i32) -> Self {
        Self { min_mapq_uni, max_circle, min_circle, linear_range_size_min }
    }

    pub fn run(&mut self, sam_file: &str, out_prefix: &str, fasta_map: &HashMap<String, String>, annotation: &Annotation) -> Result<HashSet<String>> {
        use crate::sam_bam::{detect_format, InputFormat};
        let format = detect_format(sam_file)?;
        
        match format {
            InputFormat::Sam => self.run_sam(sam_file, out_prefix, fasta_map, annotation),
            InputFormat::Bam => self.run_bam(sam_file, out_prefix, fasta_map, annotation),
        }
    }

    pub fn run_sam(&mut self, sam_file: &str, out_prefix: &str, fasta_map: &HashMap<String, String>, annotation: &Annotation) -> Result<HashSet<String>> {
        let file = File::open(sam_file)?;
        let mmap = unsafe { Mmap::map(&file)? };
        let file_size = mmap.len();
        let num_threads = rayon::current_num_threads().max(1);
        let shard_size = file_size / num_threads;

        let shard_results: Vec<Vec<(String, String)>> = (0..num_threads).into_par_iter().map(|i| {
            let start = i * shard_size;
            let end = if i == num_threads - 1 { file_size } else { (i + 1) * shard_size };
            self.process_shard_robust(&mmap, start, end, fasta_map, annotation).unwrap_or_default()
        }).collect();

        let bsj1_path = format!("{}.BSJ1", out_prefix);
        let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(&bsj1_path)?);
        let mut scan1_id_map = HashSet::new();
        for shard in shard_results {
            for (id, line) in shard {
                writeln!(writer, "{}", line)?;
                scan1_id_map.insert(id);
            }
        }
        writer.flush()?;
        Ok(scan1_id_map)
    }

    pub fn run_bam(&mut self, bam_file: &str, out_prefix: &str, fasta_map: &HashMap<String, String>, annotation: &Annotation) -> Result<HashSet<String>> {
        use noodles::bam;
        let file = File::open(bam_file)?;
        let mmap = unsafe { Mmap::map(&file)? };
        
        let mut header_reader = bam::io::Reader::new(&mmap[..]);
        let header = header_reader.read_header()?;
        
        let file_size = mmap.len();
        let num_threads = rayon::current_num_threads().max(1);
        let shard_size = file_size / num_threads;

        let shard_results: Vec<Vec<(String, String)>> = (0..num_threads).into_par_iter().map(|i| {
            let start = i * shard_size;
            let end = if i == num_threads - 1 { file_size } else { (i + 1) * shard_size };
            self.process_bam_shard(&mmap, start, end, &header, fasta_map, annotation).unwrap_or_default()
        }).collect();

        let bsj1_path = format!("{}.BSJ1", out_prefix);
        let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(&bsj1_path)?);
        let mut scan1_id_map = HashSet::new();
        for shard in shard_results {
            for (id, line) in shard {
                writeln!(writer, "{}", line)?;
                scan1_id_map.insert(id);
            }
        }
        writer.flush()?;
        Ok(scan1_id_map)
    }

    fn process_bam_shard(&self, mmap: &[u8], start: usize, end: usize, header: &sam::Header, fasta_map: &HashMap<String, String>, annotation: &Annotation) -> Result<Vec<(String, String)>> {
        use noodles::bam;
        use noodles::sam::alignment::Record;
        let mut results = Vec::new();
        let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
        
        // 1. Find the first BGZF block start in this shard
        let pos = if start == 0 { 0 } else {
            let mut found = None;
            for i in start..end {
                if i + 3 < mmap.len() && &mmap[i..i+4] == b"\x1f\x8b\x08\x04" {
                    found = Some(i);
                    break;
                }
            }
            found.unwrap_or(mmap.len())
        };
        if pos >= mmap.len() { return Ok(results); }

        let mut reader = bam::io::Reader::new(&mmap[pos..]);
        let mut record = bam::Record::default();
        
        let mut current_id: Vec<u8> = Vec::new();
        let mut group: [Vec<AlignmentRecord>; 2] = [Vec::with_capacity(8), Vec::with_capacity(8)];
        let mut first_id_skipped = start == 0;

        while reader.read_record(&mut record)? != 0 {
            let read_id = record.name().ok_or_else(|| anyhow::anyhow!("Missing read name"))?;
            
            if !first_id_skipped {
                if current_id.is_empty() {
                    current_id = read_id.to_vec();
                    continue;
                }
                if read_id.to_vec() != current_id {
                    first_id_skipped = true;
                    current_id = read_id.to_vec();
                } else {
                    continue;
                }
            }

            if read_id.to_vec() != current_id {
                if !current_id.is_empty() {
                    let id_str = String::from_utf8_lossy(&current_id);
                    if let Some(res) = self.process_group_view(&id_str, &group, fasta_map, annotation, &mut validator) { results.push(res); }
                    
                    // Stop condition: We've passed the physical shard end AND finished the current Read ID group.
                    let current_compressed_pos = pos + reader.get_ref().virtual_position().compressed() as usize;
                    if current_compressed_pos > end {
                        current_id.clear(); 
                        break;
                    }
                }
                current_id = read_id.to_vec();
                group[0].clear(); group[1].clear();
            }

            // Convert BAM record to AlignmentRecord
            let chrom = match record.reference_sequence(header) {
                Some(Ok((name, _))) => String::from_utf8_lossy(name.as_ref()).to_string(),
                _ => "*".to_string(),
            };
            
            let flag = i32::from(u16::from(record.flags()));
            let start_pos = record.alignment_start().transpose()?.map(|p| p.get() as i32).unwrap_or(0);
            let mapq = record.mapping_quality().map(u8::from).unwrap_or(0) as i32;
            
            let mut cigar = String::new();
            for result in record.cigar().iter() {
                let op = result?;
                use noodles::sam::alignment::record::cigar::op::Kind;
                let op_char = match op.kind() {
                    Kind::Match => 'M',
                    Kind::Insertion => 'I',
                    Kind::Deletion => 'D',
                    Kind::Skip => 'N',
                    Kind::SoftClip => 'S',
                    Kind::HardClip => 'H',
                    Kind::Pad => 'P',
                    Kind::SequenceMatch => '=',
                    Kind::SequenceMismatch => 'X',
                };
                cigar.push_str(&format!("{}{}", op.len(), op_char));
            }
            
            let mut seq = String::new();
            for b in record.sequence().iter() {
                seq.push(char::from(b));
            }

            group[if flag & 0x40 != 0 { 0 } else { 1 }].push(AlignmentRecord { 
                flag, 
                chrom: Cow::Owned(chrom), 
                pos: start_pos, 
                mapq, 
                cigar: Cow::Owned(cigar), 
                seq: Cow::Owned(seq) 
            });
        }
        
        // Final flush
        if !current_id.is_empty() {
            let id_str = String::from_utf8_lossy(&current_id);
            if let Some(res) = self.process_group_view(&id_str, &group, fasta_map, annotation, &mut validator) { results.push(res); }
        }

        Ok(results)
    }

    fn process_shard_robust<'a>(&self, mmap: &'a [u8], start: usize, end: usize, fasta_map: &HashMap<String, String>, annotation: &Annotation) -> Result<Vec<(String, String)>> {
        let mut results = Vec::new();
        let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
        
        // 1. Shard Synchronization: Find the actual start line
        let mut pos = if start == 0 { 0 } else { 
            // Skip the first partial line
            let mut next_line = memchr(b'\n', &mmap[start..]).map(|p| start + p + 1).unwrap_or(mmap.len());
            if next_line >= mmap.len() { return Ok(results); }
            
            // Critical: If we start in the middle, we must skip ALL records of the first Read ID 
            // because they are being handled by the tail of the previous shard.
            let first_tab = memchr(b'\t', &mmap[next_line..]).map(|p| next_line + p).unwrap_or(mmap.len());
            let first_id = &mmap[next_line..first_tab];
            
            while next_line < mmap.len() {
                let line_end = memchr(b'\n', &mmap[next_line..]).map(|p| next_line + p).unwrap_or(mmap.len());
                let line_tab = memchr(b'\t', &mmap[next_line..line_end]).map(|p| next_line + p).unwrap_or(line_end);
                if &mmap[next_line..line_tab] != first_id {
                    break; // Found the next ID group
                }
                next_line = line_end + 1;
            }
            next_line
        };

        let mut current_id: &[u8] = &[];
        let mut group: [Vec<AlignmentRecord<'a>>; 2] = [Vec::with_capacity(8), Vec::with_capacity(8)];

        while pos < mmap.len() {
            let line_end = memchr(b'\n', &mmap[pos..]).map(|p| pos + p).unwrap_or(mmap.len());
            let line = &mmap[pos..line_end];
            if line.is_empty() { pos = line_end + 1; continue; }
            if line[0] == b'@' { pos = line_end + 1; continue; }

            let mut cols = line.split(|&b| b == b'\t');
            let read_id = cols.next().unwrap();
            
            if read_id != current_id {
                if !current_id.is_empty() {
                    let id_str = unsafe { std::str::from_utf8_unchecked(current_id) };
                    if let Some(res) = self.process_group_view(id_str, &group, fasta_map, annotation, &mut validator) { results.push(res); }
                    // Stop condition: We've passed the physical shard end AND finished the current Read ID group.
                    if pos > end { break; }
                }
                current_id = read_id;
                group[0].clear(); group[1].clear();
            }

            let flag_bytes = cols.next().unwrap_or(b"0");
            let flag = fast_parse_i32(flag_bytes);
            let chrom = unsafe { std::str::from_utf8_unchecked(cols.next().unwrap_or(b"*")) };
            let start_pos = fast_parse_i32(cols.next().unwrap_or(b"0"));
            let mapq = fast_parse_i32(cols.next().unwrap_or(b"0"));
            let cigar = unsafe { std::str::from_utf8_unchecked(cols.next().unwrap_or(b"*")) };
            cols.next(); cols.next(); cols.next();
            let seq = unsafe { std::str::from_utf8_unchecked(cols.next().unwrap_or(b"*")).trim() };

            group[if flag & 0x40 != 0 { 0 } else { 1 }].push(AlignmentRecord { 
                flag, 
                chrom: Cow::Borrowed(chrom), 
                pos: start_pos, 
                mapq, 
                cigar: Cow::Borrowed(cigar), 
                seq: Cow::Borrowed(seq) 
            });
            pos = line_end + 1;
        }
        
        // Final flush for the last ID in the shard
        if !current_id.is_empty() && pos >= mmap.len() {
            let id_str = unsafe { std::str::from_utf8_unchecked(current_id) };
            if let Some(res) = self.process_group_view(id_str, &group, fasta_map, annotation, &mut validator) { results.push(res); }
        }

        Ok(results)
    }

    fn process_group_view(&self, read_id: &str, group: &[Vec<AlignmentRecord>; 2], fasta_map: &HashMap<String, String>, annotation: &Annotation, validator: &mut IsBSJHg2) -> Option<(String, String)> {
        let mut candidates: Vec<BSJCandidate> = Vec::new();
        let [pair1, pair2] = group;
        for n in 0..=1 {
            let segments = if n == 0 { pair1 } else { pair2 };
            let mate_segments = if n == 0 { pair2 } else { pair1 };
            if segments.len() < 2 { continue; }
            let seq_len = segments[0].seq.len() as i32;
            for i in 0..segments.len() {
                for j in i + 1..segments.len() {
                    let (mut al1, mut al2) = (&segments[i], &segments[j]);
                    if al1.chrom != al2.chrom { continue; }
                    let mut c1 = misd(&al1.cigar, seq_len);
                    let mut c2 = misd(&al2.cigar, seq_len);
                    if c1[3] == -2 || c2[3] == -2 { continue; }
                    if c1[0] > c2[0] { std::mem::swap(&mut c1, &mut c2); std::mem::swap(&mut al1, &mut al2); }
                    let mut identified = false;
                    let (mut s1_n, mut s2_n, mut adj1, mut adj2): (i32, i32, i32, i32) = (0, 0, 0, 0);
                    let (mut str1, mut str2, mut str3, mut str4) = (String::new(), String::new(), String::new(), String::new());
                    let (mut q1, mut q2, mut sum_q) = (0, 0, 0);
                    if c1[0] * c2[0] == -1 {
                        let scale = c1[0] * al1.pos + c1[2] + c2[0] * al2.pos + c2[2];
                        if scale > 0 && (c1[1] - c2[1]).abs() <= 6 && scale <= self.max_circle && scale >= self.min_circle {
                            adj1 = (c1[1] * c1[0] + c2[1] * c2[0]) / 2; adj2 = (c1[1] * c1[0] + c2[1] * c2[0]) - adj1;
                            if adj1.abs() <= 4 {
                                identified = true; s1_n = al1.pos + adj1; s2_n = al2.pos + c2[3] - 1 - adj2;
                                if al1.flag & 0x40 != 0 { str2 = java_substring(&al1.seq, 0, c1[1] + adj1).to_string(); str1 = java_substring(&al1.seq, c1[1] + adj1, seq_len).to_string(); }
                                else { let cr = crate::utils::reverse_complement(&al1.seq); str2 = java_substring(&cr, 0, c1[1] + adj1).to_string(); str1 = java_substring(&cr, c1[1] + adj1, seq_len).to_string(); }
                                str3 = "*".to_string();
                                if al1.mapq >= self.min_mapq_uni && al2.mapq >= self.min_mapq_uni { q1 = 1; q2 = 1; sum_q = 1; }
                                else if al1.mapq >= self.min_mapq_uni { q1 = 1; } else if al2.mapq >= self.min_mapq_uni { q2 = 1; }
                            }
                        }
                    } else if (c1[0] * c2[0]).abs() == 10 {
                        if c1[0] == -1 {
                            let scale = al2.pos + c2[3] - 1 - al1.pos;
                            if scale > 0 && (seq_len - c2[2] - c1[1]).abs() <= 6 && scale <= self.max_circle && scale >= self.min_circle {
                                adj1 = (c2[1] + c2[3] - c1[1]) / 2; adj2 = (c2[1] + c2[3] - c1[1]) - adj1;
                                if adj1.abs() <= 4 {
                                    identified = true; s1_n = al1.pos + adj1; s2_n = al2.pos + c2[3] - 1 - adj2;
                                    if al1.flag & 0x40 != 0 { str1 = java_substring(&al1.seq, c1[1] + adj1, seq_len).to_string(); str2 = java_substring(&al1.seq, c2[1], c1[1] + adj1).to_string(); str3 = java_substring(&al1.seq, 0, c2[1]).to_string(); }
                                    else { let cr = crate::utils::reverse_complement(&al1.seq); str1 = java_substring(&cr, c1[1] + adj1, seq_len).to_string(); str2 = java_substring(&cr, c2[1], c1[1] + adj1).to_string(); str3 = java_substring(&cr, 0, c2[1]).to_string(); }
                                    if al1.mapq >= self.min_mapq_uni && al2.mapq >= self.min_mapq_uni { q1 = 1; q2 = 1; sum_q = 1; }
                                    else if al1.mapq >= self.min_mapq_uni { q1 = 1; } else if al2.mapq >= self.min_mapq_uni { q2 = 1; }
                                }
                            }
                        } else {
                            let scale = al1.pos + c1[3] - 1 - al2.pos;
                            if scale > 0 && (c1[1] - c2[1]).abs() <= 6 && scale <= self.max_circle && scale >= self.min_circle {
                                adj1 = (c1[1] - c2[1]) / 2; adj2 = (c1[1] - c2[1]) - adj1;
                                if adj1.abs() <= 4 {
                                    identified = true; s1_n = al2.pos + adj1; s2_n = al1.pos + c1[3] - 1 - adj2;
                                    if al1.flag & 0x40 != 0 { str2 = java_substring(&al1.seq, 0, c1[1] - adj2).to_string(); str1 = java_substring(&al1.seq, c1[1] - adj2, seq_len - c2[2]).to_string(); str3 = java_substring(&al1.seq, seq_len - c2[2], seq_len).to_string(); }
                                    else { let cr = crate::utils::reverse_complement(&al1.seq); str2 = java_substring(&cr, 0, c1[1] - adj2).to_string(); str1 = java_substring(&cr, c1[1] - adj2, seq_len - c2[2]).to_string(); str3 = java_substring(&cr, seq_len - c2[2], seq_len).to_string(); }
                                    if al1.mapq >= self.min_mapq_uni && al2.mapq >= self.min_mapq_uni { q1 = 1; q2 = 1; sum_q = 1; }
                                    else if al2.mapq >= self.min_mapq_uni { q1 = 1; } else if al1.mapq >= self.min_mapq_uni { q2 = 1; }
                                }
                            }
                        }
                    }
                    if identified {
                        let mut s4_ok = 0;
                        if !mate_segments.is_empty() {
                            str4 = if al1.flag & 0x10 != mate_segments[0].flag & 0x10 { mate_segments[0].seq.to_string() } else { crate::utils::reverse_complement(&mate_segments[0].seq) };
                            for m_aln in mate_segments {
                                let mc = misd(&m_aln.cigar, seq_len);
                                if m_aln.chrom == al1.chrom && m_aln.mapq >= self.min_mapq_uni {
                                    if m_aln.flag & 0x10 != al1.flag & 0x10 && m_aln.pos >= s1_n - 6 && m_aln.pos + mc[3] <= s2_n + 6 { s4_ok = 1; break; }
                                    else { s4_ok = -1; break; }
                                }
                            }
                        } else { s4_ok = 1; str4 = "".to_string(); }
                        if s1_n > s2_n { std::mem::swap(&mut s1_n, &mut s2_n); }
                        let mut line_arr = vec![if al1.flag & 0x10 != 0 { "1" } else { "0" }.to_string(), al1.chrom.to_string(), str1, str2, str3, str4, q1.to_string(), q2.to_string(), s4_ok.to_string(), s1_n.to_string(), s2_n.to_string(), adj1.to_string(), adj2.to_string()];
                        if let Some(chr_seq) = fasta_map.get(al1.chrom.as_ref()) {
                            if let Some(res) = validator.is_bsj_hg1(&mut line_arr, chr_seq, sum_q, "chrM", false, &annotation.chr_exon_start_map, &annotation.chr_exon_end_map) {
                                let t_v = res.split('\t').next().unwrap_or("0").parse().unwrap_or(0);
                                candidates.push(BSJCandidate { result_str: res, tag: t_v, sum_q, total_mq: al1.mapq + al2.mapq, total_m_len: c1[3] + c2[3], cigar_pair: format!("{};{}", al1.cigar, al2.cigar) });
                            }
                        }
                    }
                }
            }
        }
        if !candidates.is_empty() {
            candidates.sort_by(|a, b| b.tag.cmp(&a.tag).then(b.sum_q.cmp(&a.sum_q)).then(b.total_mq.cmp(&a.total_mq)).then(b.total_m_len.cmp(&a.total_m_len)));
            Some((read_id.to_string(), format!("{}\t{}\t{}", read_id, candidates[0].cigar_pair, candidates[0].result_str)))
        } else { None }
    }
}
