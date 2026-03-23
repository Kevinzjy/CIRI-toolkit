//! Scan 1 module: Ultra-high performance multi-format BSJ identification.
//! Optimized for TB-scale data using sharded I/O and strict Page Cache management.

use std::collections::{HashMap, HashSet};
use std::fs::{File};
use std::io::{Write, BufWriter, BufRead, BufReader, Seek};
use std::borrow::Cow;
use std::fmt::Write as FmtWrite;
use std::sync::mpsc;
use std::thread;
use std::time::{Duration, Instant};
use std::sync::atomic::{AtomicU64, Ordering};
use anyhow::Result;
use rayon::prelude::*;
use memmap2::Mmap;
use noodles::sam::{self, alignment::{Record as _, record::Sequence as _}};
use indicatif::{ProgressBar, ProgressStyle};
use crate::misd::misd;
use crate::is_bsj_hg2::{IsBSJHg2, java_substring, report_scan1_hg_profile};
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

type OwnedAlignmentRecord = AlignmentRecord<'static>;
type OwnedStandMap = HashMap<i32, (char, Cow<'static, str>)>;

struct SamOwnedGroup {
    read_id: String,
    group: [Vec<OwnedAlignmentRecord>; 2],
    stand_map: OwnedStandMap,
    align_num: usize,
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

#[derive(Default)]
struct Scan1Profile {
    wall_total_ns: AtomicU64,
    shard_total_ns: AtomicU64,
    group_process_ns: AtomicU64,
    bsj_judge_ns: AtomicU64,
    write_ns: AtomicU64,
    merge_ns: AtomicU64,
    records: AtomicU64,
    groups: AtomicU64,
    hg1_calls: AtomicU64,
    hg1_hits: AtomicU64,
}

impl Scan1Profile {
    fn enabled_from_env() -> bool {
        matches!(std::env::var("CIRI_PROFILE_SCAN1"), Ok(v) if !v.is_empty() && v != "0")
    }

    fn report(&self) {
        let wall_total_ns = self.wall_total_ns.load(Ordering::Relaxed);
        let shard_total_ns = self.shard_total_ns.load(Ordering::Relaxed);
        let group_process_ns = self.group_process_ns.load(Ordering::Relaxed);
        let bsj_judge_ns = self.bsj_judge_ns.load(Ordering::Relaxed);
        let write_ns = self.write_ns.load(Ordering::Relaxed);
        let merge_ns = self.merge_ns.load(Ordering::Relaxed);
        let records = self.records.load(Ordering::Relaxed);
        let groups = self.groups.load(Ordering::Relaxed);
        let hg1_calls = self.hg1_calls.load(Ordering::Relaxed);
        let hg1_hits = self.hg1_hits.load(Ordering::Relaxed);
        let other_shard_ns = shard_total_ns.saturating_sub(group_process_ns.saturating_add(write_ns));
        let pct = |part: u64, whole: u64| -> f64 {
            if whole == 0 { 0.0 } else { part as f64 * 100.0 / whole as f64 }
        };
        eprintln!(
            "[PROFILE_SCAN1] wall_ms={:.3} shard_work_ms={:.3} merge_ms={:.3} records={} groups={} hg1_calls={} hg1_hits={}",
            wall_total_ns as f64 / 1_000_000.0,
            shard_total_ns as f64 / 1_000_000.0,
            merge_ns as f64 / 1_000_000.0,
            records,
            groups,
            hg1_calls,
            hg1_hits
        );
        eprintln!(
            "[PROFILE_SCAN1] shard_breakdown_ms group_process={:.3} ({:.1}%) bsj_judge={:.3} ({:.1}% of group) write={:.3} ({:.1}%) other={:.3} ({:.1}%)",
            group_process_ns as f64 / 1_000_000.0,
            pct(group_process_ns, shard_total_ns),
            bsj_judge_ns as f64 / 1_000_000.0,
            pct(bsj_judge_ns, group_process_ns),
            write_ns as f64 / 1_000_000.0,
            pct(write_ns, shard_total_ns),
            other_shard_ns as f64 / 1_000_000.0,
            pct(other_shard_ns, shard_total_ns),
        );
        report_scan1_hg_profile();
    }
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
        let file_size = std::fs::metadata(sam_file)?.len();
        let pb = ProgressBar::new(file_size);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {percent:>3}% ({eta}) {msg}")?
                .progress_chars("#>-"),
        );
        pb.set_message("Scan 1 (SAM): reading groups...");
        pb.enable_steady_tick(Duration::from_millis(120));

        let shard_out = format!("{}.BSJ1.shard_0", out_prefix);
        self.read_len = self.process_sam_file_to_file(sam_file, fasta_map, annotation, &pb, &shard_out)?;

        pb.finish_with_message("Scan 1 (SAM): completed");
        self.merge_and_collect_ids(out_prefix, 1)
    }

    fn process_sam_file_to_file(&self, sam_file: &str, fasta_map: &HashMap<String, String>, annotation: &Annotation, pb: &ProgressBar, out_path: &str) -> Result<i32> {
        let mut writer = BufWriter::with_capacity(256 * 1024, File::create(out_path)?);
        let batch_size = (rayon::current_num_threads().max(1) * 256).max(1024);
        let (tx, rx) = mpsc::sync_channel::<Vec<SamOwnedGroup>>(rayon::current_num_threads().max(2));
        let pb_clone = pb.clone();

        let max_read_len = thread::scope(|scope| -> Result<i32> {
            let producer = scope.spawn(|| self.stream_sam_group_batches(sam_file, &pb_clone, tx, batch_size));

            for mut batch in rx {
                self.process_sam_group_batch(&mut writer, &mut batch, fasta_map, annotation)?;
            }

            producer
                .join()
                .map_err(|_| anyhow::anyhow!("SAM group producer thread panicked"))?
        })?;

        writer.flush()?;
        Ok(max_read_len)
    }

    fn stream_sam_group_batches(
        &self,
        sam_file: &str,
        pb: &ProgressBar,
        tx: mpsc::SyncSender<Vec<SamOwnedGroup>>,
        batch_size: usize,
    ) -> Result<i32> {
        let sam_reader = BufReader::with_capacity(1024 * 1024, File::open(sam_file)?);
        let mut reader = sam::io::Reader::new(sam_reader);
        let header = reader.read_header()?;
        let mut record = sam::Record::default();
        let mut current_id: Vec<u8> = Vec::new();
        let mut group: [Vec<OwnedAlignmentRecord>; 2] = [Vec::with_capacity(8), Vec::with_capacity(8)];
        let mut stand_map: OwnedStandMap = HashMap::with_capacity(4);
        let mut align_num = 0usize;
        let mut max_read_len = 0i32;
        let mut cigar_buf = String::with_capacity(64);
        let mut seq_buf = String::with_capacity(256);
        let mut batch = Vec::with_capacity(batch_size);
        let mut records_since_progress = 0usize;
        let mut last_progress_pos = reader.get_mut().stream_position()?;

        while reader.read_record(&mut record)? != 0 {
            let read_id = record.name().ok_or_else(|| anyhow::anyhow!("Missing read name"))?;
            if read_id.to_vec() != current_id {
                if !current_id.is_empty() {
                    self.push_sam_owned_group(&mut batch, &current_id, &mut group, &mut stand_map, align_num);
                    if batch.len() >= batch_size {
                        tx.send(std::mem::take(&mut batch))
                            .map_err(|_| anyhow::anyhow!("SAM group consumer dropped"))?;
                        batch = Vec::with_capacity(batch_size);
                    }
                }
                current_id = read_id.to_vec();
                align_num = 0;
            }

            let chrom = match record.reference_sequence(&header) {
                Some(Ok((name, _))) => String::from_utf8_lossy(name).to_string(),
                _ => "*".to_string(),
            };
            let flag = i32::from(u16::from(record.flags()?));
            let start_pos = record.alignment_start().transpose()?.map(|p| p.get() as i32).unwrap_or(0);
            let mapq = record.mapping_quality().transpose()?.map(u8::from).unwrap_or(0) as i32;

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
                let _ = write!(&mut cigar_buf, "{}{}", op.len(), op_char);
            }

            seq_buf.clear();
            for b in record.sequence().iter() {
                seq_buf.push(char::from(b));
            }
            max_read_len = max_read_len.max(seq_buf.len() as i32);
            let seq = seq_buf.clone();
            let s_idx = if flag & 0x40 != 0 { 1 } else { 0 };
            if !seq.is_empty() && seq != "*" {
                let st_c = if flag & 0x10 != 0 { '1' } else { '0' };
                stand_map.entry(s_idx).or_insert_with(|| (st_c, Cow::Owned(seq.clone())));
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
            records_since_progress += 1;
            if records_since_progress >= 4096 {
                let pos = reader.get_mut().stream_position()?;
                if pos > last_progress_pos {
                    pb.inc(pos - last_progress_pos);
                    last_progress_pos = pos;
                }
                records_since_progress = 0;
            }
        }

        if !current_id.is_empty() {
            self.push_sam_owned_group(&mut batch, &current_id, &mut group, &mut stand_map, align_num);
        }
        if !batch.is_empty() {
            tx.send(batch)
                .map_err(|_| anyhow::anyhow!("SAM group consumer dropped"))?;
        }

        let final_pos = reader.get_mut().stream_position()?;
        if final_pos > last_progress_pos {
            pb.inc(final_pos - last_progress_pos);
        }

        Ok(max_read_len)
    }

    fn push_sam_owned_group(
        &self,
        batch: &mut Vec<SamOwnedGroup>,
        current_id: &[u8],
        group: &mut [Vec<OwnedAlignmentRecord>; 2],
        stand_map: &mut OwnedStandMap,
        align_num: usize,
    ) {
        batch.push(SamOwnedGroup {
            read_id: String::from_utf8_lossy(current_id).into_owned(),
            group: [std::mem::take(&mut group[0]), std::mem::take(&mut group[1])],
            stand_map: std::mem::take(stand_map),
            align_num,
        });
        group[0] = Vec::with_capacity(8);
        group[1] = Vec::with_capacity(8);
        *stand_map = HashMap::with_capacity(4);
    }

    fn process_sam_group_batch(
        &self,
        writer: &mut BufWriter<File>,
        batch: &mut Vec<SamOwnedGroup>,
        fasta_map: &HashMap<String, String>,
        annotation: &Annotation,
    ) -> Result<()> {
        let groups = std::mem::take(batch);
        let results: Vec<Option<String>> = groups
            .into_par_iter()
            .map(|owned| {
                let non_empty_groups = owned.group.iter().filter(|g| !g.is_empty()).count();
                if owned.align_num <= 2 && non_empty_groups != 1 {
                    return None;
                }
                let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
                self.process_group_view(&owned.read_id, &owned.group, &owned.stand_map, fasta_map, annotation, &mut validator, None)
                    .map(|(_, res_line)| res_line)
            })
            .collect();

        for res_line in results.into_iter().flatten() {
            writeln!(writer, "{}", res_line)?;
        }

        Ok(())
    }

    fn merge_and_collect_ids(&self, out_prefix: &str, num_threads: usize) -> Result<HashSet<String>> {
        let bsj1_path = format!("{}.BSJ1", out_prefix);
        let mut final_writer = BufWriter::with_capacity(1024 * 1024, File::create(&bsj1_path)?);
        let mut scan1_id_map = HashSet::new();
        for i in 0..num_threads {
            let shard_path = format!("{}.BSJ1.shard_{}", out_prefix, i);
            if let Ok(shard_file) = File::open(&shard_path) {
                let mut reader = BufReader::new(shard_file);
                let mut line = String::new();
                loop {
                    line.clear();
                    if reader.read_line(&mut line)? == 0 {
                        break;
                    }
                    let tab_idx = line.find('\t').unwrap_or(0);
                    if tab_idx > 0 {
                        scan1_id_map.insert(line[..tab_idx].to_string());
                        final_writer.write_all(line.as_bytes())?;
                    }
                }
            }
            let _ = std::fs::remove_file(shard_path);
        }
        final_writer.flush()?;
        Ok(scan1_id_map)
    }

    pub fn run_bam(&mut self, bam_file: &str, out_prefix: &str, fasta_map: &HashMap<String, String>, annotation: &Annotation) -> Result<HashSet<String>> {
        use noodles::bam;
        let run_started = Instant::now();
        let profile = if Scan1Profile::enabled_from_env() { Some(Scan1Profile::default()) } else { None };
        let profile_ref = profile.as_ref();
        let file = File::open(bam_file)?;
        let file_size = std::fs::metadata(bam_file)?.len();
        let mmap = unsafe { Mmap::map(&file)? };
        let num_threads = rayon::current_num_threads().max(1);
        let shard_size = mmap.len() / num_threads;
        let mut reader = bam::io::Reader::new(file);
        let header = reader.read_header()?;

        let pb = ProgressBar::new(file_size);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {percent:>3}% ({eta}) {msg}")?
                .progress_chars("#>-"),
        );
        pb.set_message("Scan 1 (BAM): starting...");
        pb.enable_steady_tick(Duration::from_millis(120));

        unsafe { libc::madvise(mmap.as_ptr() as *mut libc::c_void, mmap.len(), libc::MADV_SEQUENTIAL); }

        let header_ref = &header;
        let shard_max_read_lens: Vec<i32> = (0..num_threads)
            .into_par_iter()
            .map(|i| {
                let start = i * shard_size;
                let end = if i == num_threads - 1 { mmap.len() } else { (i + 1) * shard_size };
                let shard_out = format!("{}.BSJ1.shard_{}", out_prefix, i);
                self.process_bam_shard_to_file(&mmap, start, end, header_ref, fasta_map, annotation, &pb, &shard_out, profile_ref)
                    .unwrap_or(0)
            })
            .collect();

        self.read_len = shard_max_read_lens.into_iter().max().unwrap_or(0);
        pb.set_position(file_size);
        pb.finish_with_message("Scan 1 (BAM): completed");
        let merge_started = Instant::now();
        let result = self.merge_and_collect_ids(out_prefix, num_threads);
        if let Some(profile) = profile_ref {
            profile.merge_ns.fetch_add(merge_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile.wall_total_ns.fetch_add(run_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile.report();
        }
        result
    }

    fn process_bam_shard_to_file(&self, mmap: &Mmap, start: usize, end: usize, header: &sam::Header, fasta_map: &HashMap<String, String>, annotation: &Annotation, pb: &ProgressBar, out_path: &str, profile: Option<&Scan1Profile>) -> Result<i32> {
        use noodles::bam;
        let shard_started = profile.map(|_| Instant::now());
        let mut writer = BufWriter::with_capacity(256 * 1024, File::create(out_path)?);
        let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
        let pos = if start == 0 { 0 } else {
            let mut found = None;
            for i in start..mmap.len() { if i + 3 < mmap.len() && &mmap[i..i+4] == b"\x1f\x8b\x08\x04" { found = Some(i); break; } }
            found.unwrap_or(mmap.len())
        };
        if pos >= mmap.len() { return Ok(0); }
        let mut reader = bam::io::Reader::new(&mmap[pos..]);
        if start == 0 {
            let _ = reader.read_header()?;
        }
        let mut record = bam::Record::default();
        let mut current_id: Vec<u8> = Vec::new();
        let mut group: [Vec<AlignmentRecord>; 2] = [Vec::with_capacity(8), Vec::with_capacity(8)];
        let mut stand_map: HashMap<i32, (char, Cow<str>)> = HashMap::with_capacity(4);
        let mut align_num = 0usize;
        let mut first_id_skipped = start == 0;
        let mut last_compressed_pos = 0;
        let mut last_evicted_pos = pos;
        let eviction_threshold = 64 * 1024 * 1024;
        let mut shard_max_read_len = 0i32;
        let mut cigar_buf = String::with_capacity(64);
        let mut seq_buf = String::with_capacity(256);
        while reader.read_record(&mut record)? != 0 {
            if let Some(profile) = profile {
                profile.records.fetch_add(1, Ordering::Relaxed);
            }
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
                        if let Some((_, res_line)) = self.process_group_view(&id_str, &group, &stand_map, fasta_map, annotation, &mut validator, profile) {
                            let write_started = profile.map(|_| Instant::now());
                            writeln!(writer, "{}", res_line)?;
                            if let (Some(profile), Some(write_started)) = (profile, write_started) {
                                profile.write_ns.fetch_add(write_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
                            }
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
            cigar_buf.clear();
            for result in record.cigar().iter() {
                let op = result?; use noodles::sam::alignment::record::cigar::op::Kind;
                let op_char = match op.kind() { Kind::Match => 'M', Kind::Insertion => 'I', Kind::Deletion => 'D', Kind::Skip => 'N', Kind::SoftClip => 'S', Kind::HardClip => 'H', Kind::Pad => 'P', Kind::SequenceMatch => '=', Kind::SequenceMismatch => 'X', };
                let _ = write!(&mut cigar_buf, "{}{}", op.len(), op_char);
            }
            seq_buf.clear();
            for b in record.sequence().iter() { seq_buf.push(char::from(b)); }
            shard_max_read_len = shard_max_read_len.max(seq_buf.len() as i32);
            let seq = seq_buf.clone();
            if !seq.is_empty() && seq != "*" {
                let s_idx = if flag & 0x40 != 0 { 1 } else { 0 };
                let st_c = if flag & 0x10 != 0 { '1' } else { '0' };
                stand_map.entry(s_idx).or_insert_with(|| (st_c, Cow::Owned(seq.clone())));
            }
            group[if flag & 0x40 != 0 { 1 } else { 0 }].push(AlignmentRecord { flag, chrom: Cow::Owned(chrom), pos: start_pos, mapq, cigar: Cow::Owned(cigar_buf.clone()), seq: Cow::Owned(seq) });
            align_num += 1;
        }
        if !current_id.is_empty() {
            let id_str = String::from_utf8_lossy(&current_id);
            if let Some((_, res_line)) = self.process_group_view(&id_str, &group, &stand_map, fasta_map, annotation, &mut validator, profile) {
                let write_started = profile.map(|_| Instant::now());
                writeln!(writer, "{}", res_line)?;
                if let (Some(profile), Some(write_started)) = (profile, write_started) {
                    profile.write_ns.fetch_add(write_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
                }
            }
        }
        writer.flush()?;
        if let (Some(profile), Some(shard_started)) = (profile, shard_started) {
            profile.shard_total_ns.fetch_add(shard_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
        Ok(shard_max_read_len)
    }

    fn process_group_view<'a>(&self, read_id: &str, group: &[Vec<AlignmentRecord<'a>>; 2], stand_map: &HashMap<i32, (char, Cow<'a, str>)>, fasta_map: &HashMap<String, String>, annotation: &Annotation, validator: &mut IsBSJHg2, profile: Option<&Scan1Profile>) -> Option<(String, String)> {
        let started = profile.map(|_| Instant::now());
        let result = self.process_group_view_impl(read_id, group, stand_map, fasta_map, annotation, validator, profile);
        if let (Some(profile), Some(started)) = (profile, started) {
            profile.group_process_ns.fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile.groups.fetch_add(1, Ordering::Relaxed);
        }
        result
    }

    fn process_group_view_impl<'a>(&self, read_id: &str, group: &[Vec<AlignmentRecord<'a>>; 2], stand_map: &HashMap<i32, (char, Cow<'a, str>)>, fasta_map: &HashMap<String, String>, annotation: &Annotation, validator: &mut IsBSJHg2, profile: Option<&Scan1Profile>) -> Option<(String, String)> {
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
            let segment_misds: Vec<[i32; 4]> = segments.iter().map(|aln| misd(&aln.cigar, seq_len)).collect();
            let mate_misds: Vec<[i32; 4]> = mate_segments.iter().map(|aln| misd(&aln.cigar, seq_len)).collect();
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
                    let mut c1 = segment_misds[i];
                    let mut c2 = segment_misds[j];
                    if c1[0] > c2[0] {
                        std::mem::swap(&mut c1, &mut c2);
                        std::mem::swap(&mut al1, &mut al2);
                    }
                    let mut identified = false;
                    let (mut s1_n, mut s2_n, mut adj1, mut adj2): (i32, i32, i32, i32) = (0, 0, 0, 0);
                    let (mut str1, mut str2, mut str3, mut str4) = (String::new(), String::new(), String::new(), String::new());
                    let (mut q1, mut q2, mut sum_q) = (0, 0, 0);
                    let al1_strand = if al1.flag & 0x10 != 0 { '1' } else { '0' };
                    let seq_oriented: Cow<'_, str> = if al1_strand == read_strand { Cow::Borrowed(read_seq) } else { Cow::Owned(crate::utils::reverse_complement(read_seq)) };
                    if c1[0] * c2[0] == -1 {
                        let scale = c1[0] * al1.pos + c1[2] + c2[0] * al2.pos + c2[2];
                        if scale > 0 && (c1[1] - c2[1]).abs() <= 6 && scale <= self.max_circle && scale >= self.min_circle {
                            adj1 = (c1[1] * c1[0] + c2[1] * c2[0]) / 2; adj2 = (c1[1] * c1[0] + c2[1] * c2[0]) - adj1;
                            if adj1.abs() <= 4 {
                                identified = true; s1_n = al1.pos + adj1; s2_n = al2.pos + c2[3] - 1 - adj2;
                                str2 = java_substring(seq_oriented.as_ref(), 0, c1[1] + adj1).to_string();
                                str1 = java_substring(seq_oriented.as_ref(), c1[1] + adj1, seq_len).to_string();
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
                                    str1 = java_substring(seq_oriented.as_ref(), c1[1] + adj1, seq_len).to_string();
                                    str2 = java_substring(seq_oriented.as_ref(), c2[1], c1[1] + adj1).to_string();
                                    str3 = java_substring(seq_oriented.as_ref(), 0, c2[1]).to_string();
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
                                    str2 = java_substring(seq_oriented.as_ref(), 0, c1[1] - adj2).to_string();
                                    str1 = java_substring(seq_oriented.as_ref(), c1[1] - adj2, seq_len - c2[2]).to_string();
                                    str3 = java_substring(seq_oriented.as_ref(), seq_len - c2[2], seq_len).to_string();
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
                            for (m_aln, mc) in mate_segments.iter().zip(mate_misds.iter()) {
                                if m_aln.chrom == al1.chrom && m_aln.mapq >= self.min_mapq_uni {
                                    if m_aln.flag & 0x10 != al1.flag & 0x10 && m_aln.pos >= s1_n - 6 && m_aln.pos + mc[3] - 1 <= s2_n + 6 { s4_ok = 1; break; }
                                    else { s4_ok = -1; break; }
                                }
                            }
                        } else { s4_ok = 1; str4 = "".to_string(); }
                        let mut line_arr = [if al1.flag & 0x10 != 0 { "1" } else { "0" }.to_string(), al1.chrom.to_string(), str1, str2, str3, str4, q1.to_string(), q2.to_string(), s4_ok.to_string(), s1_n.to_string(), s2_n.to_string(), adj1.to_string(), adj2.to_string()];
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
                            let hg1_started = profile.map(|_| Instant::now());
                            let res = validator.is_bsj_hg1(&mut line_arr, chr_seq, sum_q, "chrM", false, &annotation.chr_exon_start_map, &annotation.chr_exon_end_map);
                            if let Some(profile) = profile {
                                profile.hg1_calls.fetch_add(1, Ordering::Relaxed);
                                if let Some(hg1_started) = hg1_started {
                                    profile.bsj_judge_ns.fetch_add(hg1_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
                                }
                                if res.is_some() {
                                    profile.hg1_hits.fetch_add(1, Ordering::Relaxed);
                                }
                            }
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
                                return Some((read_id.to_string(), format!("{}\t{};{}\t{}", read_id, al1.cigar, al2.cigar, res)));
                            }
                        }
                    }
                }
            }
        }
        None
    }
}
