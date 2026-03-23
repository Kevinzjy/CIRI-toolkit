//! Scan 1 module: Ultra-high performance multi-format BSJ identification.
//! Optimized for TB-scale data using sharded I/O and strict Page Cache management.

use std::collections::{HashMap, HashSet};
use std::fs::{File};
use std::io::{Write, BufWriter, BufRead, BufReader};
use std::borrow::Cow;
use std::time::{Duration, Instant};
use anyhow::Result;
use rayon::prelude::*;
use memmap2::Mmap;
use memchr::memchr;
use noodles::sam::{self, alignment::Record as _};
use indicatif::{ProgressBar, ProgressStyle};
use crate::misd::misd;
use crate::is_bsj_hg2::{IsBSJHg2, java_substring};
use crate::annotation::Annotation;
use crate::utils::AlignmentRecord;

pub struct Scan1 {
    pub min_mapq_uni: i32,
    pub max_circle: i32,
    pub min_circle: i32,
    pub linear_range_size_min: i32,
    pub mem_limit: u64,
    pub read_len: i32,
}

#[inline]
fn fast_parse_i32(bytes: &[u8]) -> i32 {
    let mut res = 0;
    for &b in bytes { if b >= b'0' && b <= b'9' { res = res * 10 + (b - b'0') as i32; } }
    res
}

/// Active Page Cache eviction with strict alignment.
fn advise_dontneed_aligned(mmap: &Mmap, offset: usize, len: usize) {
    if len == 0 { return; }
    // Page alignment is mandatory for madvise (typically 4096 bytes)
    let page_size = 4096;
    let aligned_offset = (offset / page_size) * page_size;
    let aligned_len = ((offset + len + page_size - 1) / page_size) * page_size - aligned_offset;
    
    unsafe {
        let ptr = mmap.as_ptr().add(aligned_offset);
        libc::madvise(ptr as *mut libc::c_void, aligned_len, libc::MADV_DONTNEED);
    }
}

#[inline]
fn should_trace_read(read_id: &str) -> bool {
    // Optional targeted trace hook for parity debugging.
    // Enabled only when CIRI_TRACE_READS is explicitly set.
    if let Ok(raw) = std::env::var("CIRI_TRACE_READS") {
        for token in raw.split(',') {
            let t = token.trim();
            if !t.is_empty() && t == read_id {
                return true;
            }
        }
    }
    false
}

impl Scan1 {
    pub fn new(min_mapq_uni: i32, min_circle: i32, max_circle: i32, linear_range_size_min: i32) -> Self {
        Self {
            min_mapq_uni,
            max_circle,
            min_circle,
            linear_range_size_min,
            mem_limit: 2 * 1024 * 1024 * 1024,
            read_len: 0,
        }
    }

    pub fn set_mem_limit(&mut self, limit: u64) {
        self.mem_limit = limit;
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

        unsafe { libc::madvise(mmap.as_ptr() as *mut libc::c_void, file_size, libc::MADV_SEQUENTIAL); }

        let pb = ProgressBar::new(file_size as u64);
        pb.set_style(ProgressStyle::default_bar().template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta}) {msg}")?.progress_chars("#>-"));
        pb.set_message("Scan 1: Identifying BSJ candidates");

        (0..num_threads).into_par_iter().for_each(|i| {
            let start = i * shard_size;
            let end = if i == num_threads - 1 { file_size } else { (i + 1) * shard_size };
            let shard_out = format!("{}.BSJ1.shard_{}", out_prefix, i);
            let _ = self.process_shard_robust_to_file(&mmap, start, end, fasta_map, annotation, &pb, &shard_out);
        });

        pb.finish_with_message("Scan 1: Completed");
        self.merge_and_collect_ids(out_prefix, num_threads)
    }

    fn merge_and_collect_ids(&self, out_prefix: &str, num_threads: usize) -> Result<HashSet<String>> {
        let bsj1_path = format!("{}.BSJ1", out_prefix);
        let mut final_writer = BufWriter::with_capacity(1024 * 1024, File::create(&bsj1_path)?);
        let mut scan1_id_map = HashSet::new();
        for i in 0..num_threads {
            let shard_path = format!("{}.BSJ1.shard_{}", out_prefix, i);
            if let Ok(shard_file) = File::open(&shard_path) {
                let reader = BufReader::new(shard_file);
                for line_res in reader.lines() {
                    let line = line_res?;
                    let id = line.split('\t').next().unwrap_or("").to_string();
                    if !id.is_empty() {
                        writeln!(final_writer, "{}", line)?;
                        scan1_id_map.insert(id);
                    }
                }
            }
            let _ = std::fs::remove_file(shard_path);
        }
        final_writer.flush()?;
        Ok(scan1_id_map)
    }

    fn process_shard_robust_to_file<'a>(&self, mmap: &'a Mmap, start: usize, end: usize, fasta_map: &HashMap<String, String>, annotation: &Annotation, pb: &ProgressBar, out_path: &str) -> Result<()> {
        let mut writer = BufWriter::with_capacity(256 * 1024, File::create(out_path)?);
        let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
        let mut pos = if start == 0 { 0 } else { 
            let mut next_line = memchr(b'\n', &mmap[start..]).map(|p| start + p + 1).unwrap_or(mmap.len());
            if next_line >= mmap.len() { return Ok(()); }
            let first_tab = memchr(b'\t', &mmap[next_line..]).map(|p| next_line + p).unwrap_or(mmap.len());
            let first_id = &mmap[next_line..first_tab];
            while next_line < mmap.len() {
                let line_end = memchr(b'\n', &mmap[next_line..]).map(|p| next_line + p).unwrap_or(mmap.len());
                let line_tab = memchr(b'\t', &mmap[next_line..line_end]).map(|p| next_line + p).unwrap_or(line_end);
                if &mmap[next_line..line_tab] != first_id { break; }
                next_line = line_end + 1;
            }
            next_line
        };

        let mut current_id: &[u8] = &[];
        let mut group: [Vec<AlignmentRecord<'a>>; 2] = [Vec::with_capacity(8), Vec::with_capacity(8)];
        let mut stand_map: HashMap<i32, (char, Cow<'a, str>)> = HashMap::with_capacity(4);
        let mut align_num = 0usize;
        let mut one_read_key: i32 = -1;
        let mut last_evicted_pos = pos;
        
        // Aggressive eviction: use small fixed windows (e.g. 64MB) to keep RES low
        let eviction_threshold = 64 * 1024 * 1024; 

        while pos < mmap.len() {
            let line_end = memchr(b'\n', &mmap[pos..]).map(|p| pos + p).unwrap_or(mmap.len());
            let line = &mmap[pos..line_end];
            if line.is_empty() { pb.inc(1); pos = line_end + 1; continue; }
            
            if pos - last_evicted_pos > eviction_threshold {
                advise_dontneed_aligned(mmap, last_evicted_pos, pos - last_evicted_pos);
                last_evicted_pos = pos;
            }

            if line[0] == b'@' { pb.inc((line_end - pos + 1) as u64); pos = line_end + 1; continue; }
            let mut cols = line.split(|&b| b == b'\t');
            let read_id = cols.next().unwrap();
            if read_id != current_id {
                if !current_id.is_empty() {
                    let id_str = unsafe { std::str::from_utf8_unchecked(current_id) };
                    let non_empty_groups = group.iter().filter(|g| !g.is_empty()).count();
                    if align_num > 2 || non_empty_groups == 1 {
                        if let Some((_, res_line)) = self.process_group_view(id_str, &group, &stand_map, fasta_map, annotation, &mut validator) {
                            writeln!(writer, "{}", res_line)?;
                        }
                    }
                    if pos > end { break; }
                }
                current_id = read_id;
                group[0].clear();
                group[1].clear();
                stand_map.clear();
                align_num = 0;
                one_read_key = -1;
            }
            let flag = fast_parse_i32(cols.next().unwrap_or(b"0"));
            let chrom = unsafe { std::str::from_utf8_unchecked(cols.next().unwrap_or(b"*")) };
            let start_pos = fast_parse_i32(cols.next().unwrap_or(b"0"));
            let mapq = fast_parse_i32(cols.next().unwrap_or(b"0"));
            let cigar = unsafe { std::str::from_utf8_unchecked(cols.next().unwrap_or(b"*")) };
            cols.next(); cols.next(); cols.next();
            let seq = unsafe { std::str::from_utf8_unchecked(cols.next().unwrap_or(b"*")).trim() };
            let s_idx = if flag & 0x40 != 0 { 1 } else { 0 };
            if s_idx != one_read_key {
                group[s_idx as usize].clear();
                one_read_key = s_idx;
                let st_c = if flag & 0x10 != 0 { '1' } else { '0' };
                stand_map.insert(s_idx, (st_c, Cow::Borrowed(seq)));
            }
            group[s_idx as usize].push(AlignmentRecord { flag, chrom: Cow::Borrowed(chrom), pos: start_pos, mapq, cigar: Cow::Borrowed(cigar), seq: Cow::Borrowed(seq) });
            align_num += 1;
            pb.inc((line_end - pos + 1) as u64);
            pos = line_end + 1;
        }
        if !current_id.is_empty() && pos >= mmap.len() {
            let id_str = unsafe { std::str::from_utf8_unchecked(current_id) };
            if let Some((_, res_line)) = self.process_group_view(id_str, &group, &stand_map, fasta_map, annotation, &mut validator) {
                writeln!(writer, "{}", res_line)?;
            }
        }
        writer.flush()?;
        Ok(())
    }

    pub fn run_bam(&mut self, bam_file: &str, out_prefix: &str, fasta_map: &HashMap<String, String>, annotation: &Annotation) -> Result<HashSet<String>> {
        use noodles::bam;
        const WRITE_FLUSH_THRESHOLD: usize = 512 * 1024;
        const PROGRESS_UPDATE_INTERVAL: u64 = 4096;
        let file = File::open(bam_file)?;
        let file_size = std::fs::metadata(bam_file)?.len();
        let mut reader = bam::io::Reader::new(file);
        let header = reader.read_header()?;
        let mut record = bam::Record::default();
        let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);

        let bsj1_path = format!("{}.BSJ1", out_prefix);
        let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(&bsj1_path)?);
        let mut scan1_id_map = HashSet::new();
        let mut out_buf = String::with_capacity(WRITE_FLUSH_THRESHOLD);

        let mut current_id: Vec<u8> = Vec::new();
        let mut group: [Vec<AlignmentRecord>; 2] = [Vec::with_capacity(8), Vec::with_capacity(8)];
        let mut stand_map: HashMap<i32, (char, Cow<str>)> = HashMap::with_capacity(4);
        let mut align_num = 0usize;
        let mut one_read_key: i32 = -1;
        let mut total_records: u64 = 0;
        let mut total_read_groups: u64 = 0;
        let mut bsj_hits: u64 = 0;
        let mut cigar_buf = String::with_capacity(64);
        let mut seq_buf = String::with_capacity(256);
        let profile_scan1 = std::env::var("CIRI_PROFILE_SCAN1").ok().as_deref() == Some("1");
        let scan_start = Instant::now();
        let mut decode_time = Duration::ZERO;
        let mut group_eval_time = Duration::ZERO;
        let mut write_time = Duration::ZERO;

        let pb = ProgressBar::new(file_size);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {percent:>3}% ({eta}) {msg}")?
                .progress_chars("#>-"),
        );
        pb.set_message("Scan 1 (BAM): starting...");
        pb.enable_steady_tick(Duration::from_millis(120));

        while reader.read_record(&mut record)? != 0 {
            total_records += 1;
            let pos = reader.get_ref().virtual_position().compressed();
            if total_records % PROGRESS_UPDATE_INTERVAL == 0 {
                pb.set_position(pos.min(file_size));
            }
            let read_id = record.name().ok_or_else(|| anyhow::anyhow!("Missing read name"))?;
            let read_id_bytes = read_id.as_ref();
            if read_id_bytes != current_id.as_slice() {
                if !current_id.is_empty() {
                    total_read_groups += 1;
                    let id_str = String::from_utf8_lossy(&current_id);
                    let non_empty_groups = group.iter().filter(|g| !g.is_empty()).count();
                    if align_num > 2 || non_empty_groups == 1 {
                        let t0 = if profile_scan1 { Some(Instant::now()) } else { None };
                        if let Some((_, res_line)) = self.process_group_view(&id_str, &group, &stand_map, fasta_map, annotation, &mut validator) {
                            out_buf.push_str(&res_line);
                            out_buf.push('\n');
                            scan1_id_map.insert(id_str.to_string());
                            bsj_hits += 1;
                            if out_buf.len() >= WRITE_FLUSH_THRESHOLD {
                                let w0 = if profile_scan1 { Some(Instant::now()) } else { None };
                                writer.write_all(out_buf.as_bytes())?;
                                out_buf.clear();
                                if let Some(w0) = w0 {
                                    write_time += w0.elapsed();
                                }
                            }
                        }
                        if let Some(t0) = t0 {
                            group_eval_time += t0.elapsed();
                        }
                    }
                }
                current_id.clear();
                current_id.extend_from_slice(read_id_bytes);
                group[0].clear();
                group[1].clear();
                stand_map.clear();
                align_num = 0;
                one_read_key = -1;
            }

            let decode_start = if profile_scan1 { Some(Instant::now()) } else { None };
            let chrom = match record.reference_sequence(&header) {
                Some(Ok((name, _))) => String::from_utf8_lossy(name).to_string(),
                _ => "*".to_string(),
            };
            let flag = i32::from(u16::from(record.flags()));
            let start_pos = record.alignment_start().transpose()?.map(|p| p.get() as i32).unwrap_or(0);
            let mapq = record.mapping_quality().map(u8::from).unwrap_or(0) as i32;
            cigar_buf.clear();
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
                use std::fmt::Write as _;
                let _ = write!(&mut cigar_buf, "{}{}", op.len(), op_char);
            }
            seq_buf.clear();
            for b in record.sequence().iter() {
                seq_buf.push(char::from(b));
            }
            if (seq_buf.len() as i32) > self.read_len {
                self.read_len = seq_buf.len() as i32;
            }
            if let Some(decode_start) = decode_start {
                decode_time += decode_start.elapsed();
            }
            let seq = seq_buf.clone();
            let s_idx = if flag & 0x40 != 0 { 1 } else { 0 };
            if s_idx != one_read_key {
                group[s_idx as usize].clear();
                one_read_key = s_idx;
                let st_c = if flag & 0x10 != 0 { '1' } else { '0' };
                stand_map.insert(s_idx, (st_c, Cow::Owned(seq.clone())));
            }
            group[s_idx as usize].push(AlignmentRecord {
                flag,
                chrom: Cow::Owned(chrom),
                pos: start_pos,
                mapq,
                cigar: Cow::Owned(cigar_buf.clone()),
                seq: Cow::Owned(seq),
            });
            align_num += 1;
            if total_records % 100_000 == 0 {
                pb.set_message(format!(
                    "Scan 1 (BAM): records={} read_groups={} bsj_hits={}",
                    total_records, total_read_groups, bsj_hits
                ));
            }
        }

        if !current_id.is_empty() {
            total_read_groups += 1;
            let id_str = String::from_utf8_lossy(&current_id);
            let t0 = if profile_scan1 { Some(Instant::now()) } else { None };
            if let Some((_, res_line)) = self.process_group_view(&id_str, &group, &stand_map, fasta_map, annotation, &mut validator) {
                out_buf.push_str(&res_line);
                out_buf.push('\n');
                scan1_id_map.insert(id_str.to_string());
                bsj_hits += 1;
            }
            if let Some(t0) = t0 {
                group_eval_time += t0.elapsed();
            }
        }
        if !out_buf.is_empty() {
            let w0 = if profile_scan1 { Some(Instant::now()) } else { None };
            writer.write_all(out_buf.as_bytes())?;
            if let Some(w0) = w0 {
                write_time += w0.elapsed();
            }
        }
        writer.flush()?;
        pb.set_position(file_size);
        pb.finish_with_message(format!(
            "Scan 1 (BAM): done. records={} read_groups={} bsj_hits={}",
            total_records, total_read_groups, bsj_hits
        ));
        if profile_scan1 {
            let total = scan_start.elapsed();
            eprintln!(
                "[PROFILE_SCAN1_BAM] total={:.3}s decode={:.3}s group_eval={:.3}s write={:.3}s other={:.3}s records={} groups={} hits={}",
                total.as_secs_f64(),
                decode_time.as_secs_f64(),
                group_eval_time.as_secs_f64(),
                write_time.as_secs_f64(),
                (total
                    .saturating_sub(decode_time)
                    .saturating_sub(group_eval_time)
                    .saturating_sub(write_time))
                .as_secs_f64(),
                total_records,
                total_read_groups,
                bsj_hits
            );
        }
        Ok(scan1_id_map)
    }

    #[allow(dead_code)]
    fn process_bam_shard_to_file(&self, mmap: &Mmap, start: usize, end: usize, header: &sam::Header, fasta_map: &HashMap<String, String>, annotation: &Annotation, pb: &ProgressBar, out_path: &str) -> Result<()> {
        use noodles::bam;
        let mut writer = BufWriter::with_capacity(256 * 1024, File::create(out_path)?);
        let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
        let pos = if start == 0 { 0 } else {
            let mut found = None;
            for i in start..mmap.len() { if i + 3 < mmap.len() && &mmap[i..i+4] == b"\x1f\x8b\x08\x04" { found = Some(i); break; } }
            found.unwrap_or(mmap.len())
        };
        if pos >= mmap.len() { return Ok(()); }
        let mut reader = bam::io::Reader::new(&mmap[pos..]);
        let mut record = bam::Record::default();
        let mut current_id: Vec<u8> = Vec::new();
        let mut group: [Vec<AlignmentRecord>; 2] = [Vec::with_capacity(8), Vec::with_capacity(8)];
        let mut stand_map: HashMap<i32, (char, Cow<str>)> = HashMap::with_capacity(4);
        let mut align_num = 0usize;
        let mut first_id_skipped = start == 0;
        let mut last_compressed_pos = 0;
        let mut last_evicted_pos = pos;
        let eviction_threshold = 64 * 1024 * 1024;
        while reader.read_record(&mut record)? != 0 {
            let curr_c_pos = reader.get_ref().virtual_position().compressed() as usize;
            pb.inc((curr_c_pos - last_compressed_pos) as u64);
            last_compressed_pos = curr_c_pos;
            if curr_c_pos - (last_evicted_pos - pos) > eviction_threshold {
                advise_dontneed_aligned(mmap, last_evicted_pos, curr_c_pos - (last_evicted_pos - pos));
                last_evicted_pos = pos + curr_c_pos;
            }
            let read_id = record.name().ok_or_else(|| anyhow::anyhow!("Missing read name"))?;
            if !first_id_skipped {
                if current_id.is_empty() { current_id = read_id.to_vec(); continue; }
                if read_id.to_vec() != current_id { first_id_skipped = true; current_id = read_id.to_vec(); } else { continue; }
            }
            if read_id.to_vec() != current_id {
                if !current_id.is_empty() {
                    let id_str = String::from_utf8_lossy(&current_id);
                    let non_empty_groups = group.iter().filter(|g| !g.is_empty()).count();
                    if align_num > 2 || non_empty_groups == 1 {
                        if let Some((_, res_line)) = self.process_group_view(&id_str, &group, &stand_map, fasta_map, annotation, &mut validator) {
                            writeln!(writer, "{}", res_line)?;
                        }
                    }
                    if pos + curr_c_pos > end { current_id.clear(); break; }
                }
                current_id = read_id.to_vec();
                group[0].clear();
                group[1].clear();
                stand_map.clear();
                align_num = 0;
            }
            let chrom = match record.reference_sequence(header) { Some(Ok((name, _))) => String::from_utf8_lossy(name).to_string(), _ => "*".to_string() };
            let flag = i32::from(u16::from(record.flags()));
            let start_pos = record.alignment_start().transpose()?.map(|p| p.get() as i32).unwrap_or(0);
            let mapq = record.mapping_quality().map(u8::from).unwrap_or(0) as i32;
            let mut cigar = String::new();
            for result in record.cigar().iter() {
                let op = result?; use noodles::sam::alignment::record::cigar::op::Kind;
                let op_char = match op.kind() { Kind::Match => 'M', Kind::Insertion => 'I', Kind::Deletion => 'D', Kind::Skip => 'N', Kind::SoftClip => 'S', Kind::HardClip => 'H', Kind::Pad => 'P', Kind::SequenceMatch => '=', Kind::SequenceMismatch => 'X', };
                cigar.push_str(&format!("{}{}", op.len(), op_char));
            }
            let mut seq = String::new(); for b in record.sequence().iter() { seq.push(char::from(b)); }
            if !seq.is_empty() && seq != "*" {
                let s_idx = if flag & 0x40 != 0 { 1 } else { 0 };
                let st_c = if flag & 0x10 != 0 { '1' } else { '0' };
                stand_map.entry(s_idx).or_insert((st_c, Cow::Owned(seq.clone())));
            }
            group[if flag & 0x40 != 0 { 1 } else { 0 }].push(AlignmentRecord { flag, chrom: Cow::Owned(chrom), pos: start_pos, mapq, cigar: Cow::Owned(cigar), seq: Cow::Owned(seq) });
            align_num += 1;
        }
        if !current_id.is_empty() {
            let id_str = String::from_utf8_lossy(&current_id);
            if let Some((_, res_line)) = self.process_group_view(&id_str, &group, &stand_map, fasta_map, annotation, &mut validator) {
                writeln!(writer, "{}", res_line)?;
            }
        }
        writer.flush()?;
        Ok(())
    }

    fn process_group_view<'a>(&self, read_id: &str, group: &[Vec<AlignmentRecord<'a>>; 2], stand_map: &HashMap<i32, (char, Cow<'a, str>)>, fasta_map: &HashMap<String, String>, annotation: &Annotation, validator: &mut IsBSJHg2) -> Option<(String, String)> {
        // Debug trace is intentionally read-scoped to avoid overwhelming output.
        let trace_read = should_trace_read(read_id);
        let [pair1, pair2] = group;
        let group_count = (if !pair1.is_empty() { 1 } else { 0 }) + (if !pair2.is_empty() { 1 } else { 0 });
        if group_count == 0 {
            return None;
        }
        let is_paried = group_count - 1;
        for n in 0..=is_paried {
            let segments = if n == 0 { pair1 } else { pair2 };
            let mate_segments = if n == 0 { pair2 } else { pair1 };
            let (read_strand, read_seq) = match stand_map.get(&(n as i32)) {
                Some((st, seq)) => (*st, seq.as_ref()),
                None => continue,
            };
            if segments.len() < 2 { continue; }
            let seq_len = read_seq.len() as i32;
            for i in 0..segments.len() {
                for j in i + 1..segments.len() {
                    let (mut al1, mut al2) = (&segments[i], &segments[j]);
                    if al1.chrom != al2.chrom {
                        continue;
                    }
                    // Java parity: skip chrM unless explicitly enabled (currently always disabled).
                    if al1.chrom.as_ref() == "chrM" {
                        continue;
                    }
                    // Java parity: two candidate segments must share FLAG bit 0x10.
                    if (al1.flag & 0x10) != (al2.flag & 0x10) {
                        continue;
                    }
                    // Java parity: at least one segment MAPQ passes threshold.
                    if al1.mapq < self.min_mapq_uni && al2.mapq < self.min_mapq_uni {
                        continue;
                    }
                    if al1.cigar.as_ref() == "*" || al2.cigar.as_ref() == "*" {
                        continue;
                    }
                    let mut c1 = misd(&al1.cigar, seq_len);
                    let mut c2 = misd(&al2.cigar, seq_len);
                    if c1[0] > c2[0] {
                        std::mem::swap(&mut c1, &mut c2);
                        std::mem::swap(&mut al1, &mut al2);
                    }
                    let mut identified = false;
                    let (mut s1_n, mut s2_n, mut adj1, mut adj2): (i32, i32, i32, i32) = (0, 0, 0, 0);
                    let (mut str1, mut str2, mut str3, mut str4) = (String::new(), String::new(), String::new(), String::new());
                    let (mut q1, mut q2, mut sum_q) = (0, 0, 0);
                    let al1_strand = if al1.flag & 0x10 != 0 { '1' } else { '0' };
                    let seq_oriented = if al1_strand == read_strand { read_seq.to_string() } else { crate::utils::reverse_complement(read_seq) };
                    if c1[0] * c2[0] == -1 {
                        let scale = c1[0] * al1.pos + c1[2] + c2[0] * al2.pos + c2[2];
                        if scale > 0 && (c1[1] - c2[1]).abs() <= 6 && scale <= self.max_circle && scale >= self.min_circle {
                            adj1 = (c1[1] * c1[0] + c2[1] * c2[0]) / 2; adj2 = (c1[1] * c1[0] + c2[1] * c2[0]) - adj1;
                            if adj1.abs() <= 4 {
                                identified = true; s1_n = al1.pos + adj1; s2_n = al2.pos + c2[3] - 1 - adj2;
                                str2 = java_substring(&seq_oriented, 0, c1[1] + adj1).to_string();
                                str1 = java_substring(&seq_oriented, c1[1] + adj1, seq_len).to_string();
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
                                    str1 = java_substring(&seq_oriented, c1[1] + adj1, seq_len).to_string();
                                    str2 = java_substring(&seq_oriented, c2[1], c1[1] + adj1).to_string();
                                    str3 = java_substring(&seq_oriented, 0, c2[1]).to_string();
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
                                    str2 = java_substring(&seq_oriented, 0, c1[1] - adj2).to_string();
                                    str1 = java_substring(&seq_oriented, c1[1] - adj2, seq_len - c2[2]).to_string();
                                    str3 = java_substring(&seq_oriented, seq_len - c2[2], seq_len).to_string();
                                    if al1.mapq >= self.min_mapq_uni && al2.mapq >= self.min_mapq_uni { q1 = 1; q2 = 1; sum_q = 1; }
                                    else if al2.mapq >= self.min_mapq_uni { q1 = 1; } else if al1.mapq >= self.min_mapq_uni { q2 = 1; }
                                }
                            }
                        }
                    }
                    if identified {
                        let mut s4_ok = 0;
                        if is_paried == 1 {
                            if let Some((mate_stand, mate_seq)) = stand_map.get(&((1 - n) as i32)) {
                                str4 = if al1_strand != *mate_stand { mate_seq.to_string() } else { crate::utils::reverse_complement(mate_seq) };
                            }
                            for m_aln in mate_segments {
                                let mc = misd(&m_aln.cigar, seq_len);
                                if m_aln.chrom == al1.chrom && m_aln.mapq >= self.min_mapq_uni {
                                    if m_aln.flag & 0x10 != al1.flag & 0x10 && m_aln.pos >= s1_n - 6 && m_aln.pos + mc[3] - 1 <= s2_n + 6 { s4_ok = 1; break; }
                                    else { s4_ok = -1; break; }
                                }
                            }
                        } else { s4_ok = 1; str4 = "".to_string(); }
                        let mut line_arr = vec![if al1.flag & 0x10 != 0 { "1" } else { "0" }.to_string(), al1.chrom.to_string(), str1, str2, str3, str4, q1.to_string(), q2.to_string(), s4_ok.to_string(), s1_n.to_string(), s2_n.to_string(), adj1.to_string(), adj2.to_string()];
                        if trace_read {
                            eprintln!(
                                "[TRACE_SCAN1_CAND] id={} n={} al1=({}, {}, {}, {}) al2=({}, {}, {}, {}) c1={:?} c2={:?} s1_n={} s2_n={} adj1={} adj2={} q=({},{},{}) s4_ok={} line_arr6_8={}/{}/{}",
                                read_id,
                                n,
                                al1.chrom,
                                al1.pos,
                                al1.mapq,
                                al1.cigar,
                                al2.chrom,
                                al2.pos,
                                al2.mapq,
                                al2.cigar,
                                c1,
                                c2,
                                s1_n,
                                s2_n,
                                adj1,
                                adj2,
                                q1,
                                q2,
                                sum_q,
                                s4_ok,
                                line_arr[6],
                                line_arr[7],
                                line_arr[8]
                            );
                        }
                        if let Some(chr_seq) = fasta_map.get(al1.chrom.as_ref()) {
                            let res = validator.is_bsj_hg1(&mut line_arr, chr_seq, sum_q, "chrM", false, &annotation.chr_exon_start_map, &annotation.chr_exon_end_map);
                            if trace_read {
                                eprintln!(
                                    "[TRACE_SCAN1_HG1] id={} result={} post_line_arr6_8={}/{}/{} sites={}->{}",
                                    read_id,
                                    if res.is_some() { "Some" } else { "None" },
                                    line_arr[6],
                                    line_arr[7],
                                    line_arr[8],
                                    line_arr[9],
                                    line_arr[10]
                                );
                            }
                            if let Some(res) = res {
                                return Some((read_id.to_string(), format!("{}\t{}\t{}", read_id, format!("{};{}", al1.cigar, al2.cigar), res)));
                            }
                        }
                    }
                }
            }
        }
        None
    }
}
