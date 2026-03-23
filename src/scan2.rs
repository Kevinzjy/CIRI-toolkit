//! Scan 2 module: Parallel multi-format sequence rescue and FSJ counting.
//! Optimized for TB-scale data using sharded I/O and active Page Cache eviction.

use std::collections::{HashMap, HashSet};
use std::fs::{File};
use std::io::{BufRead, BufReader, Write, BufWriter};
use std::borrow::Cow;
use std::time::Duration;
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
    ///
    /// Java parity note:
    /// - Backed by de-duplicated `siteInfor` records (site1..sum_q), not raw BSJ1 lines.
    /// - Candidate payload layout is `[site1, site2, strand, signal1, signal2, sum_q]`.
    pub index1: HashMap<String, Vec<CandidateBreakpoint>>,
    /// Index for BSJ candidates (indexed by Site2).
    pub index2: HashMap<String, Vec<CandidateBreakpoint>>,
    /// Java parity gate: existing buckets for Site1 index.
    pub site_array1: HashMap<String, HashSet<i32>>,
    /// Java parity gate: existing buckets for Site2 index.
    pub site_array2: HashMap<String, HashSet<i32>>,
    /// Map to store Forward Spliced Junction (FSJ) counts.
    pub fsj_map: HashMap<String, i32>,
    /// Memory limit per thread for Page Cache eviction.
    pub mem_limit: u64,
    /// Read IDs already assigned in Scan1 (Java parity: skip in Scan2).
    pub scan1_ids: HashSet<String>,
    /// Java Scan2 bucket size (`seqLen`) used for directional candidate traversal.
    pub seq_len: i32,
}

#[derive(Clone)]
pub struct CandidateBreakpoint {
    /// Genomic coordinate used for bucket/range lookup.
    pub site: i32,
    /// Insertion order after Java-style de-duplication.
    /// Used as a stable tie-breaker when two candidates share the same `site`.
    pub order: usize,
    /// Candidate payload in Java `siteInfor` layout:
    /// `[site1, site2, strand, signal1, signal2, sum_q]`.
    pub data: Vec<String>,
}

#[inline]
fn fast_parse_i32(bytes: &[u8]) -> i32 {
    let mut res = 0;
    for &b in bytes { if b >= b'0' && b <= b'9' { res = res * 10 + (b - b'0') as i32; } }
    res
}

/// Advisory eviction: Tells the kernel we no longer need these memory pages.
fn advise_dontneed(mmap: &Mmap, offset: usize, len: usize) {
    if len == 0 { return; }
    unsafe {
        let ptr = mmap.as_ptr().add(offset);
        libc::madvise(ptr as *mut libc::c_void, len, libc::MADV_DONTNEED);
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

impl Scan2 {
    #[inline]
    fn java_string_hash(s: &str) -> i32 {
        // Java String.hashCode(): h = 31*h + ch (32-bit overflow).
        let mut h: i32 = 0;
        for ch in s.chars() {
            h = h.wrapping_mul(31).wrapping_add(ch as i32);
        }
        h
    }

    #[inline]
    fn java_hash_spread(h: i32) -> u32 {
        let x = h as u32;
        x ^ (x >> 16)
    }

    #[inline]
    fn java_hashset_capacity(n: usize) -> usize {
        // Java HashMap defaults: initial capacity 16, load factor 0.75.
        // Capacity is power-of-two and grows so that size <= cap * 0.75.
        let mut cap = 16usize;
        while n > ((cap as f64) * 0.75f64) as usize {
            cap = cap.saturating_mul(2);
        }
        cap.max(1)
    }

    /// Creates a new `Scan2` instance.
    pub fn new(min_mapq_uni: i32, linear_range_size_min: i32, seq_len: i32) -> Self {
        Self { 
            min_mapq_uni, 
            linear_range_size_min, 
            index1: HashMap::new(), 
            index2: HashMap::new(), 
            site_array1: HashMap::new(),
            site_array2: HashMap::new(),
            fsj_map: HashMap::new(),
            mem_limit: 2 * 1024 * 1024 * 1024,
            scan1_ids: HashSet::new(),
            seq_len,
        }
    }

    fn lower_bound(list: &[CandidateBreakpoint], key: i32) -> usize {
        let mut lo = 0usize;
        let mut hi = list.len();
        while lo < hi {
            let mid = (lo + hi) / 2;
            if list[mid].site < key {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    }

    fn bucket_range(list: &[CandidateBreakpoint], bucket: i32, bucket_size: i32) -> (usize, usize) {
        if bucket < 0 || bucket_size <= 0 {
            return (0, 0);
        }
        let start = bucket.saturating_mul(bucket_size);
        let end_exclusive = (bucket + 1).saturating_mul(bucket_size);
        let l = Self::lower_bound(list, start);
        let r = Self::lower_bound(list, end_exclusive);
        (l, r)
    }

    /// Sets the per-thread memory limit.
    pub fn set_mem_limit(&mut self, limit: u64) {
        self.mem_limit = limit;
    }

    /// Loads BSJ1 candidates and initializes the FSJ map.
    pub fn build_index(&mut self, bsj1_file: &str) -> Result<()> {
        // Java parity: Scan2 index is built from unique circ sites, not all BSJ1 rows.
        // Equivalent Java flow:
        //   BSJ1 -> scan1IdMap + chrCircSiteMap(HashSet) -> siteArray/siteMap.
        let file = File::open(bsj1_file)?;
        let reader = BufReader::new(file);
        let mut chr_circ_site_seen: HashMap<String, HashSet<String>> = HashMap::new();
        let mut chr_circ_site_insertion: HashMap<String, Vec<String>> = HashMap::new();
        for line_res in reader.lines() {
            let line = line_res?;
            let p: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
            if p.len() < 10 { continue; }
            self.scan1_ids.insert(p[0].clone());
            let chr = p[3].clone();
            let site_infor = p[4..].join("\t");
            let seen = chr_circ_site_seen.entry(chr.clone()).or_insert_with(HashSet::new);
            if seen.insert(site_infor.clone()) {
                chr_circ_site_insertion
                    .entry(chr)
                    .or_insert_with(Vec::new)
                    .push(site_infor);
            }
        }

        let mut order_counter: usize = 0;
        for (chr, insertion_list) in chr_circ_site_insertion {
            // Java parity: emulate HashSet iteration order:
            // bucket index ascending, and insertion order within each bucket.
            let cap = Self::java_hashset_capacity(insertion_list.len());
            let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); cap];
            for (idx, s) in insertion_list.iter().enumerate() {
                let spread = Self::java_hash_spread(Self::java_string_hash(s));
                let bucket_idx = (spread as usize) & (cap - 1);
                buckets[bucket_idx].push(idx);
            }

            for bucket in buckets {
                for idx in bucket {
                    let arr: Vec<String> = insertion_list[idx].split('\t').map(|s| s.to_string()).collect();
                    if arr.len() < 6 {
                        continue;
                    }
                let site1 = arr[0].parse().unwrap_or(0);
                let site2 = arr[1].parse().unwrap_or(0);
                self.fsj_map.entry(format!("{}\t{}\t{}", chr, arr[0], arr[1])).or_insert(0);
                self.index1
                    .entry(chr.clone())
                    .or_insert_with(Vec::new)
                    .push(CandidateBreakpoint { site: site1, order: order_counter, data: arr.clone() });
                self.index2
                    .entry(chr.clone())
                    .or_insert_with(Vec::new)
                    .push(CandidateBreakpoint { site: site2, order: order_counter, data: arr });
                let bucket_size = self.seq_len.max(1);
                self.site_array1
                    .entry(chr.clone())
                    .or_insert_with(HashSet::new)
                    .insert(site1 / bucket_size);
                self.site_array2
                    .entry(chr.clone())
                    .or_insert_with(HashSet::new)
                    .insert(site2 / bucket_size);
                order_counter += 1;
                }
            }
        }
        for list in self.index1.values_mut() {
            list.sort_by_key(|x| (x.site, x.order));
        }
        for list in self.index2.values_mut() {
            list.sort_by_key(|x| (x.site, x.order));
        }
        Ok(())
    }

    fn collect_fsj_keys_in_range(&self, chr: &str, start_tem: i32, end_tem: i32, style: i32, out: &mut HashSet<String>) {
        let lb = |list: &Vec<CandidateBreakpoint>, key: i32| -> usize {
            let mut lo = 0usize;
            let mut hi = list.len();
            while lo < hi {
                let mid = (lo + hi) / 2;
                if list[mid].site < key {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }
            lo
        };
        if let Some(list1) = self.index1.get(chr) {
            let mut i = lb(list1, start_tem);
            while i < list1.len() && list1[i].site <= end_tem {
                let cand = &list1[i];
                out.insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                i += 1;
            }
            if style == 10 || style == 1 {
                let mut j = i;
                while j < list1.len() && list1[j].site <= end_tem + 6 {
                    let cand = &list1[j];
                    out.insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                    j += 1;
                }
            }
        }
        if let Some(list2) = self.index2.get(chr) {
            let mut i = lb(list2, start_tem);
            while i < list2.len() && list2[i].site <= end_tem {
                let cand = &list2[i];
                out.insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                i += 1;
            }
            if style == 10 || style == -1 {
                let mut j = i;
                while j < list2.len() && list2[j].site <= end_tem + 6 {
                    let cand = &list2[j];
                    out.insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                    j += 1;
                }
            }
        }
    }

    /// Entry point for Scan 2.
    pub fn run(&mut self, sam_file: &str, output_bsj2: &str, chr_tcga_map: &HashMap<String, String>) -> Result<()> {
        use crate::sam_bam::{detect_format, InputFormat};
        let format = detect_format(sam_file)?;
        match format {
            InputFormat::Sam => self.run_sam(sam_file, output_bsj2, chr_tcga_map),
            InputFormat::Bam => self.run_bam(sam_file, output_bsj2, chr_tcga_map),
        }
    }

    pub fn run_sam(&mut self, sam_file: &str, output_bsj2: &str, chr_tcga_map: &HashMap<String, String>) -> Result<()> {
        let file = File::open(sam_file)?;
        let mmap = unsafe { Mmap::map(&file)? };
        let file_size = mmap.len();
        let num_threads = rayon::current_num_threads().max(1);
        let shard_size = file_size / num_threads;

        unsafe { libc::madvise(mmap.as_ptr() as *mut libc::c_void, file_size, libc::MADV_SEQUENTIAL); }

        let pb = ProgressBar::new(file_size as u64);
        pb.set_style(ProgressStyle::default_bar().template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta}) {msg}")?.progress_chars("#>-"));
        pb.set_message("Scan 2: Rescuing signals & counting FSJ");

        let shard_results: Vec<HashMap<String, i32>> = (0..num_threads).into_par_iter().map(|i| {
            let start = i * shard_size;
            let end = if i == num_threads - 1 { file_size } else { (i + 1) * shard_size };
            let shard_out = format!("{}.shard_{}", output_bsj2, i);
            let (_, partial_fsj) = self.process_sam_shard_to_file(&mmap, start, end, chr_tcga_map, &pb, &shard_out).unwrap_or_default();
            partial_fsj
        }).collect();

        pb.finish_with_message("Scan 2: Completed");
        self.merge_shards_and_fsj(output_bsj2, shard_results, num_threads)
    }

    fn merge_shards_and_fsj(&mut self, output_bsj2: &str, shard_fsjs: Vec<HashMap<String, i32>>, num_threads: usize) -> Result<()> {
        let mut writer = BufWriter::with_capacity(1024 * 1024, std::fs::OpenOptions::new().create(true).append(true).open(output_bsj2)?);
        for i in 0..num_threads {
            let shard_path = format!("{}.shard_{}", output_bsj2, i);
            if let Ok(shard_file) = File::open(&shard_path) {
                let mut shard_reader = BufReader::new(shard_file);
                let mut line = String::new();
                while shard_reader.read_line(&mut line)? != 0 {
                    writer.write_all(line.as_bytes())?;
                    line.clear();
                }
            }
            let _ = std::fs::remove_file(shard_path);
        }
        for partial_fsj in shard_fsjs {
            for (key, count) in partial_fsj { *self.fsj_map.entry(key).or_insert(0) += count; }
        }
        writer.flush()?;
        Ok(())
    }

    pub fn run_bam(&mut self, bam_file: &str, output_bsj2: &str, chr_tcga_map: &HashMap<String, String>) -> Result<()> {
        use noodles::bam;
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
        pb.set_message("Scan 2 (BAM): starting...");
        pb.enable_steady_tick(Duration::from_millis(120));

        unsafe { libc::madvise(mmap.as_ptr() as *mut libc::c_void, mmap.len(), libc::MADV_SEQUENTIAL); }

        let header_ref = &header;
        let shard_results: Vec<HashMap<String, i32>> = (0..num_threads)
            .into_par_iter()
            .map(|i| {
                let start = i * shard_size;
                let end = if i == num_threads - 1 { mmap.len() } else { (i + 1) * shard_size };
                let shard_out = format!("{}.shard_{}", output_bsj2, i);
                let (_, partial_fsj) = self
                    .process_bam_shard_to_file(&mmap, start, end, header_ref, chr_tcga_map, &pb, &shard_out)
                    .unwrap_or_default();
                partial_fsj
            })
            .collect();

        pb.set_position(file_size);
        pb.finish_with_message("Scan 2 (BAM): completed");
        self.merge_shards_and_fsj(output_bsj2, shard_results, num_threads)
    }

    fn process_bam_shard_to_file(&self, mmap: &Mmap, start: usize, end: usize, header: &sam::Header, chr_tcga_map: &HashMap<String, String>, pb: &ProgressBar, out_path: &str) -> Result<(Vec<String>, HashMap<String, i32>)> {
        use noodles::bam;
        let mut writer = BufWriter::with_capacity(256 * 1024, File::create(out_path)?);
        let mut local_fsj = HashMap::new();
        let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
        let pos = if start == 0 { 0 } else {
            let mut found = None;
            for i in start..mmap.len() { if i + 3 < mmap.len() && &mmap[i..i+4] == b"\x1f\x8b\x08\x04" { found = Some(i); break; } }
            found.unwrap_or(mmap.len())
        };
        if pos >= mmap.len() { return Ok((Vec::new(), local_fsj)); }

        let mut reader = bam::io::Reader::new(&mmap[pos..]);
        if start == 0 {
            let _ = reader.read_header()?;
        }
        let mut record = bam::Record::default();
        let mut current_id: Vec<u8> = Vec::new();
        let mut alignments: Vec<AlignmentRecord> = Vec::with_capacity(16);
        let mut stand_map: HashMap<i32, (char, Cow<str>)> = HashMap::with_capacity(4);
        let mut first_id_skipped = start == 0;
        let mut last_compressed_pos = 0;
        let mut last_evicted_pos = pos;
        let eviction_threshold = (self.mem_limit as f64 * 0.8) as usize;

        while reader.read_record(&mut record)? != 0 {
            let curr_c_pos = reader.get_ref().virtual_position().compressed() as usize;
            pb.inc((curr_c_pos - last_compressed_pos) as u64);
            last_compressed_pos = curr_c_pos;

            if curr_c_pos - (last_evicted_pos - pos) > eviction_threshold {
                advise_dontneed(mmap, last_evicted_pos, curr_c_pos - (last_evicted_pos - pos));
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
                    let mut res_batch = Vec::new();
                    self.process_group_view(&id_str, &alignments, &stand_map, &mut res_batch, &mut local_fsj, chr_tcga_map, &mut validator)?;
                    for line in res_batch { writeln!(writer, "{}", line)?; }
                    if pos + curr_c_pos > end { current_id.clear(); break; }
                }
                current_id = read_id.to_vec(); alignments.clear(); stand_map.clear();
            }

            let flag = i32::from(u16::from(record.flags()));
            let chrom = match record.reference_sequence(header) { Some(Ok((name, _))) => String::from_utf8_lossy(name).to_string(), _ => "*".to_string() };
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
            alignments.push(AlignmentRecord { flag, chrom: Cow::Owned(chrom), pos: start_pos, mapq, cigar: Cow::Owned(cigar), seq: Cow::Owned(seq) });
        }
        if !current_id.is_empty() {
            let id_str = String::from_utf8_lossy(&current_id);
            let mut res_batch = Vec::new();
            self.process_group_view(&id_str, &alignments, &stand_map, &mut res_batch, &mut local_fsj, chr_tcga_map, &mut validator)?;
            for line in res_batch { writeln!(writer, "{}", line)?; }
        }
        writer.flush()?;
        Ok((Vec::new(), local_fsj))
    }

    fn process_sam_shard_to_file<'a>(&self, mmap: &'a Mmap, start: usize, end: usize, chr_tcga_map: &HashMap<String, String>, pb: &ProgressBar, out_path: &str) -> Result<(Vec<String>, HashMap<String, i32>)> {
        let mut writer = BufWriter::with_capacity(256 * 1024, File::create(out_path)?);
        let mut local_fsj = HashMap::new();
        let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
        let mut pos = if start == 0 { 0 } else { memchr(b'\n', &mmap[start..]).map(|p| start + p + 1).unwrap_or(mmap.len()) };
        if pos >= mmap.len() { return Ok((Vec::new(), local_fsj)); }

        let mut current_id: &[u8] = &[];
        let mut alignments: Vec<AlignmentRecord<'a>> = Vec::with_capacity(16);
        let mut stand_map: HashMap<i32, (char, Cow<'a, str>)> = HashMap::with_capacity(4);
        let mut one_read_key: i32 = -1;
        let mut last_evicted_pos = pos;
        let eviction_threshold = (self.mem_limit as f64 * 0.8) as usize;

        while pos < mmap.len() {
            let line_end = memchr(b'\n', &mmap[pos..]).map(|p| pos + p).unwrap_or(mmap.len());
            let line = &mmap[pos..line_end];
            if line.is_empty() { pb.inc(1); pos = line_end + 1; continue; }
            if pos - last_evicted_pos > eviction_threshold {
                advise_dontneed(mmap, last_evicted_pos, pos - last_evicted_pos);
                last_evicted_pos = pos;
            }
            if line[0] == b'@' { pb.inc((line_end - pos + 1) as u64); pos = line_end + 1; continue; }

            let mut cols = line.split(|&b| b == b'\t');
            let read_id = cols.next().unwrap();
            if read_id != current_id {
                if !current_id.is_empty() {
                    let mut res_batch = Vec::new();
                    self.process_group_view(unsafe { std::str::from_utf8_unchecked(current_id) }, &alignments, &stand_map, &mut res_batch, &mut local_fsj, chr_tcga_map, &mut validator)?;
                    for l in res_batch { writeln!(writer, "{}", l)?; }
                }
                if pos >= end { break; }
                current_id = read_id; alignments.clear(); stand_map.clear(); one_read_key = -1;
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
                one_read_key = s_idx;
                alignments.retain(|a| {
                    let idx = if a.flag & 0x40 != 0 { 1 } else { 0 };
                    idx != s_idx
                });
                let st_c = if flag & 0x10 != 0 { '1' } else { '0' };
                stand_map.insert(s_idx, (st_c, Cow::Borrowed(seq)));
            }
            alignments.push(AlignmentRecord { flag, chrom: Cow::Borrowed(chrom), pos: start_pos, mapq, cigar: Cow::Borrowed(cigar), seq: Cow::Borrowed(seq) });
            pb.inc((line_end - pos + 1) as u64);
            pos = line_end + 1;
        }
        writer.flush()?;
        Ok((Vec::new(), local_fsj))
    }

    pub(crate) fn process_group_view<'a>(&self, id: &str, alignments: &[AlignmentRecord<'a>], stand_map: &HashMap<i32, (char, Cow<'a, str>)>, results: &mut Vec<String>, local_fsj: &mut HashMap<String, i32>, chr_tcga_map: &HashMap<String, String>, is_bsj_hg2: &mut IsBSJHg2) -> Result<()> {
        let trace_read = should_trace_read(id);
        let trace_all_candidates = trace_read && std::env::var("CIRI_TRACE_ALL_CANDS").ok().as_deref() == Some("1");
        if self.scan1_ids.contains(id) {
            return Ok(());
        }
        let mut segments: HashMap<i32, Vec<&AlignmentRecord<'a>>> = HashMap::new();
        for aln in alignments {
            segments.entry(if aln.flag & 0x40 != 0 { 1 } else { 0 }).or_insert_with(Vec::new).push(aln);
        }
        let mut tem_fsj_keys = HashSet::new();
        for seg_idx in [0_i32, 1_i32] {
            let Some(seg_alns) = segments.get(&seg_idx) else { continue };
            let (read_strand, seq) = match stand_map.get(&seg_idx) { Some(&(st, ref s)) => (st, s.as_ref()), None => continue };
            let slen = seq.len() as i32;
            for aln in seg_alns {
                let chr = aln.chrom.as_ref(); if !self.index1.contains_key(chr) { continue; }
                let cigar_ref: Cow<'_, str> = if aln.cigar.contains('H') {
                    Cow::Owned(aln.cigar.replace('H', "S"))
                } else {
                    Cow::Borrowed(aln.cigar.as_ref())
                };
                if cigar_ref == "*" { continue; }
                let c = misd(cigar_ref.as_ref(), slen);
                if cigar_ref == format!("{}M", slen) {
                    let start_tem = aln.pos + 6;
                    let end_tem = aln.pos + slen - 7;
                    self.collect_fsj_keys_in_range(chr, start_tem, end_tem, 0, &mut tem_fsj_keys);
                    continue;
                }
                if c[0] == -1 || c[0] == 10 {
                    if let Some(list) = self.index1.get(chr) {
                        let new_num_site = aln.pos;
                        // Java parity: fixed seqLen from constructor.
                        let bucket_size = self.seq_len.max(1);
                        let num1 = (new_num_site - 6) / bucket_size;
                        let num2 = (new_num_site + 6) / bucket_size;

                        let mut eval_candidate = |cand: &CandidateBreakpoint| -> Result<bool> {
                            let str_e = if aln.flag & 0x10 != 0 {
                                if read_strand == '0' { reverse_complement(seq) } else { seq.to_string() }
                            } else if read_strand == '0' {
                                seq.to_string()
                            } else {
                                reverse_complement(seq)
                            };
                            let e_idx = c[1] + (cand.site - new_num_site);
                            if e_idx <= 0 || e_idx > slen {
                                return Ok(false);
                            }
                            let str_f = &str_e[0..e_idx as usize];
                            let str3 = if c[0] == 10 {
                                let si = slen - c[2];
                                if si >= 0 && si <= slen { &str_e[si as usize..] } else { "" }
                            } else {
                                "*"
                            };
                            let site1 = cand.data[0].parse::<i32>().unwrap_or(0);
                            let site2 = cand.data[1].parse::<i32>().unwrap_or(0);
                            let curr_strand = if aln.flag & 0x10 != 0 { '1' } else { '0' };
                            let p_str = if let Some(&(p_strand, ref p_seq)) = stand_map.get(&(1 - seg_idx)) {
                                if p_strand != curr_strand { p_seq.to_string() } else { reverse_complement(p_seq) }
                            } else {
                                String::new()
                            };
                            let mut s2_ok = if segments.get(&(1 - seg_idx)).is_some() { 0 } else { 1 };
                            if let Some(mate_seg_alns) = segments.get(&(1 - seg_idx)) {
                                for mate_aln in mate_seg_alns {
                                    if mate_aln.chrom.as_ref() == chr && mate_aln.mapq >= self.min_mapq_uni {
                                        let c_ano = misd(&mate_aln.cigar, slen);
                                        let mate_strand = if mate_aln.flag & 0x10 != 0 { '1' } else { '0' };
                                        if mate_strand != curr_strand
                                            && mate_aln.pos >= site1 - 6
                                            && mate_aln.pos + c_ano[3] - 1 <= site2 + 6
                                        {
                                            s2_ok = 1;
                                        } else {
                                            s2_ok = -1;
                                        }
                                        break;
                                    }
                                }
                            }
                            let circ_c = vec![
                                (if aln.flag & 0x10 != 0 { "1" } else { "0" }).to_string(),
                                chr.to_string(),
                                "sm".to_string(),
                                cand.data[0].clone(),
                                cand.data[1].clone(),
                                str_f.to_string(),
                                p_str,
                                str3.to_string(),
                                s2_ok.to_string(),
                                cand.data[2].clone(),
                                cand.data[3].clone(),
                                cand.data[4].clone(),
                                aln.mapq.to_string(),
                            ];
                            let tag = is_bsj_hg2.is_bsj_hg2(&circ_c, chr_tcga_map.get(chr).unwrap());
                            if trace_read {
                                eprintln!(
                                    "[TRACE_SCAN2_CAND] id={} type=sm seg={} aln_pos={} chr={} site1={} site2={} cand_site={} cigar={} mapq={} s2_ok={} str_len={} pair_len={} tag={}",
                                    id,
                                    seg_idx,
                                    aln.pos,
                                    chr,
                                    cand.data[0],
                                    cand.data[1],
                                    cand.site,
                                    cigar_ref.as_ref(),
                                    aln.mapq,
                                    circ_c[8],
                                    circ_c[5].len(),
                                    circ_c[6].len(),
                                    tag
                                );
                            }
                            if tag == "0" {
                                tem_fsj_keys.insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                                Ok(false)
                            } else if tag != "2" {
                                let tag_body = &tag[0..tag.len() - 1];
                                results.push(format!(
                                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                    id,
                                    cigar_ref.as_ref(),
                                    tag_body,
                                    chr,
                                    cand.data[0],
                                    cand.data[1],
                                    cand.data[2],
                                    cand.data[3],
                                    cand.data[4],
                                    tag.chars().next_back().unwrap()
                                ));
                                Ok(true)
                            } else {
                                Ok(false)
                            }
                        };

                        // Java parity: num1 bucket reverse order.
                        if num1 > 0 && self.site_array1.get(chr).is_some_and(|s| s.contains(&num1)) {
                            let (l1, r1) = Self::bucket_range(list, num1, bucket_size);
                            for idx in (l1..r1).rev() {
                                let cand = &list[idx];
                                let bias = cand.site - new_num_site;
                                if bias >= -6 {
                                    if bias <= 6 && eval_candidate(cand)? && !trace_all_candidates {
                                        return Ok(());
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                        // Java parity: num2 bucket forward order.
                        if num2 != num1 && self.site_array1.get(chr).is_some_and(|s| s.contains(&num2)) {
                            let (l2, r2) = Self::bucket_range(list, num2, bucket_size);
                            for idx in l2..r2 {
                                let cand = &list[idx];
                                let bias = cand.site - new_num_site;
                                if bias <= 6 {
                                    if eval_candidate(cand)? && !trace_all_candidates {
                                        return Ok(());
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                }
                if c[0] == 1 || c[0] == 10 {
                    let new_site = aln.pos + c[3] - 1;
                    if let Some(list) = self.index2.get(chr) {
                        // Java parity: fixed seqLen from constructor.
                        let bucket_size = self.seq_len.max(1);
                        let num1 = (new_site - 6) / bucket_size;
                        let num2 = (new_site + 6) / bucket_size;

                        let mut eval_candidate = |cand: &CandidateBreakpoint| -> Result<bool> {
                            let str_e = if aln.flag & 0x10 != 0 {
                                if read_strand == '0' { reverse_complement(seq) } else { seq.to_string() }
                            } else if read_strand == '0' {
                                seq.to_string()
                            } else {
                                reverse_complement(seq)
                            };
                            let s_idx = if c[0] == 10 {
                                c[1] + c[3] + (cand.site - new_site)
                            } else {
                                c[1] + (cand.site - new_site)
                            };
                            if s_idx < 0 || s_idx >= slen {
                                return Ok(false);
                            }
                            let str_f = &str_e[s_idx as usize..];
                            let str3 = if c[0] == 10 { &str_e[0..c[1] as usize] } else { "*" };
                            let site1 = cand.data[0].parse::<i32>().unwrap_or(0);
                            let site2 = cand.data[1].parse::<i32>().unwrap_or(0);
                            let curr_strand = if aln.flag & 0x10 != 0 { '1' } else { '0' };
                            let p_str = if let Some(&(p_strand, ref p_seq)) = stand_map.get(&(1 - seg_idx)) {
                                if p_strand != curr_strand { p_seq.to_string() } else { reverse_complement(p_seq) }
                            } else {
                                String::new()
                            };
                            let mut s2_ok = if segments.get(&(1 - seg_idx)).is_some() { 0 } else { 1 };
                            if let Some(mate_seg_alns) = segments.get(&(1 - seg_idx)) {
                                for mate_aln in mate_seg_alns {
                                    if mate_aln.chrom.as_ref() == chr && mate_aln.mapq >= self.min_mapq_uni {
                                        let c_ano = misd(&mate_aln.cigar, slen);
                                        let mate_strand = if mate_aln.flag & 0x10 != 0 { '1' } else { '0' };
                                        if mate_strand != curr_strand
                                            && mate_aln.pos >= site1 - 6
                                            && mate_aln.pos + c_ano[3] - 1 <= site2 + 6
                                        {
                                            s2_ok = 1;
                                        } else {
                                            s2_ok = -1;
                                        }
                                        break;
                                    }
                                }
                            }
                            let circ_c = vec![
                                (if aln.flag & 0x10 != 0 { "1" } else { "0" }).to_string(),
                                chr.to_string(),
                                "ms".to_string(),
                                cand.data[0].clone(),
                                cand.data[1].clone(),
                                str_f.to_string(),
                                p_str,
                                str3.to_string(),
                                s2_ok.to_string(),
                                cand.data[2].clone(),
                                cand.data[3].clone(),
                                cand.data[4].clone(),
                                aln.mapq.to_string(),
                            ];
                            let tag = is_bsj_hg2.is_bsj_hg2(&circ_c, chr_tcga_map.get(chr).unwrap());
                            if trace_read {
                                eprintln!(
                                    "[TRACE_SCAN2_CAND] id={} type=ms seg={} aln_pos={} chr={} site1={} site2={} cand_site={} cigar={} mapq={} s2_ok={} str_len={} pair_len={} tag={}",
                                    id,
                                    seg_idx,
                                    aln.pos,
                                    chr,
                                    cand.data[0],
                                    cand.data[1],
                                    cand.site,
                                    cigar_ref.as_ref(),
                                    aln.mapq,
                                    circ_c[8],
                                    circ_c[5].len(),
                                    circ_c[6].len(),
                                    tag
                                );
                            }
                            if tag == "0" {
                                tem_fsj_keys.insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                                Ok(false)
                            } else if tag != "2" {
                                let tag_body = &tag[0..tag.len() - 1];
                                results.push(format!(
                                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                    id,
                                    cigar_ref.as_ref(),
                                    tag_body,
                                    chr,
                                    cand.data[0],
                                    cand.data[1],
                                    cand.data[2],
                                    cand.data[3],
                                    cand.data[4],
                                    tag.chars().next_back().unwrap()
                                ));
                                Ok(true)
                            } else {
                                Ok(false)
                            }
                        };

                        // Java parity: num1 bucket reverse order.
                        if self.site_array2.get(chr).is_some_and(|s| s.contains(&num1)) {
                            let (l1, r1) = Self::bucket_range(list, num1, bucket_size);
                            for idx in (l1..r1).rev() {
                                let cand = &list[idx];
                                let bias = cand.site - new_site;
                                if bias >= -6 {
                                    if bias <= 6 && eval_candidate(cand)? && !trace_all_candidates {
                                        return Ok(());
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                        // Java parity: num2 bucket forward order.
                        if num2 != num1 && self.site_array2.get(chr).is_some_and(|s| s.contains(&num2)) {
                            let (l2, r2) = Self::bucket_range(list, num2, bucket_size);
                            for idx in l2..r2 {
                                let cand = &list[idx];
                                let bias = cand.site - new_site;
                                if bias <= 6 {
                                    if eval_candidate(cand)? && !trace_all_candidates {
                                        return Ok(());
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                }
                let start_tem = aln.pos + 6;
                let end_tem = aln.pos + c[3] - 7;
                self.collect_fsj_keys_in_range(chr, start_tem, end_tem, c[0], &mut tem_fsj_keys);
            }
        }
        for k in tem_fsj_keys { *local_fsj.entry(k).or_insert(0) += 1; }
        Ok(())
    }
}
