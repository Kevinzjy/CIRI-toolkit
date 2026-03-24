//! Scan2: second-pass rescue, candidate validation, and FSJ counting.
//!
//! This stage consumes the Java-compatible BSJ output from Scan1, rebuilds
//! the de-duplicated candidate indexes expected by Java CIRI3, and then revisits
//! the input alignments to rescue additional support while counting FSJ evidence.

use crate::is_bsj_hg2::{report_scan2_hg_profile, IsBSJHg2};
use crate::misd::misd;
use crate::utils::{bam_shard_count, part_path, reverse_complement, AlignmentRecord};
use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use memchr::memchr;
use memmap2::Mmap;
use noodles::sam::{self, alignment::Record as _};
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Duration;
use std::time::Instant;

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
    /// Unique read IDs rescued by Scan2 and appended to the final BSJ file.
    pub rescued_reads: usize,
    /// Total unique BSJ-supporting read IDs in the final merged BSJ file.
    pub final_bsj_reads: usize,
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

/// Fast parser for positive SAM integer fields.
///
/// This avoids `str::parse` in the SAM text path, where flag/position/MAPQ are
/// decoded for every record and profiling showed the generic parser had no value.
#[inline]
fn fast_parse_i32(bytes: &[u8]) -> i32 {
    let mut res = 0;
    for &b in bytes {
        if b >= b'0' && b <= b'9' {
            res = res * 10 + (b - b'0') as i32;
        }
    }
    res
}

/// Advises the kernel that a processed `mmap` slice can be evicted from cache.
///
/// The offset is aligned manually because `madvise` requires page alignment on
/// Linux. Without that alignment the call can fail silently from the pipeline's
/// point of view and the full BAM/SAM page cache footprint keeps growing.
fn advise_dontneed(mmap: &Mmap, offset: usize, len: usize) {
    if len == 0 {
        return;
    }
    let page_size = 4096usize;
    let aligned_offset = (offset / page_size) * page_size;
    let aligned_len = ((offset + len + page_size - 1) / page_size) * page_size - aligned_offset;
    unsafe {
        let ptr = mmap.as_ptr().add(aligned_offset);
        libc::madvise(ptr as *mut libc::c_void, aligned_len, libc::MADV_DONTNEED);
    }
}

/// Returns whether a read is selected for targeted Scan2 tracing.
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

/// Optional Scan2 profiler used only when `CIRI_PROFILE_SCAN2` is enabled.
///
/// The counters stay intentionally coarse so the profiler can be left inside the
/// candidate hot path without materially perturbing release timings.
#[derive(Default)]
pub(crate) struct Scan2Profile {
    wall_total_ns: AtomicU64,
    shard_total_ns: AtomicU64,
    group_process_ns: AtomicU64,
    validator_call_ns: AtomicU64,
    write_ns: AtomicU64,
    merge_ns: AtomicU64,
    records: AtomicU64,
    groups: AtomicU64,
    candidate_checks: AtomicU64,
    candidate_hits: AtomicU64,
}

impl Scan2Profile {
    /// Checks whether release profiling is enabled for the current Scan2 run.
    fn enabled_from_env() -> bool {
        matches!(std::env::var("CIRI_PROFILE_SCAN2"), Ok(v) if !v.is_empty() && v != "0")
    }

    /// Emits the aggregated Scan2 timing summary.
    fn report(&self) {
        let wall_total_ns = self.wall_total_ns.load(Ordering::Relaxed);
        let shard_total_ns = self.shard_total_ns.load(Ordering::Relaxed);
        let group_process_ns = self.group_process_ns.load(Ordering::Relaxed);
        let validator_call_ns = self.validator_call_ns.load(Ordering::Relaxed);
        let write_ns = self.write_ns.load(Ordering::Relaxed);
        let merge_ns = self.merge_ns.load(Ordering::Relaxed);
        let records = self.records.load(Ordering::Relaxed);
        let groups = self.groups.load(Ordering::Relaxed);
        let candidate_checks = self.candidate_checks.load(Ordering::Relaxed);
        let candidate_hits = self.candidate_hits.load(Ordering::Relaxed);
        let other_shard_ns =
            shard_total_ns.saturating_sub(group_process_ns.saturating_add(write_ns));
        let pct = |part: u64, whole: u64| -> f64 {
            if whole == 0 {
                0.0
            } else {
                part as f64 * 100.0 / whole as f64
            }
        };
        eprintln!(
            "[PROFILE_SCAN2] wall_ms={:.3} shard_work_ms={:.3} merge_ms={:.3} records={} groups={} candidate_checks={} candidate_hits={}",
            wall_total_ns as f64 / 1_000_000.0,
            shard_total_ns as f64 / 1_000_000.0,
            merge_ns as f64 / 1_000_000.0,
            records,
            groups,
            candidate_checks,
            candidate_hits,
        );
        eprintln!(
            "[PROFILE_SCAN2] shard_breakdown_ms group_process={:.3} ({:.1}%) validator={:.3} ({:.1}% of group) write={:.3} ({:.1}%) other={:.3} ({:.1}%)",
            group_process_ns as f64 / 1_000_000.0,
            pct(group_process_ns, shard_total_ns),
            validator_call_ns as f64 / 1_000_000.0,
            pct(validator_call_ns, group_process_ns),
            write_ns as f64 / 1_000_000.0,
            pct(write_ns, shard_total_ns),
            other_shard_ns as f64 / 1_000_000.0,
            pct(other_shard_ns, shard_total_ns),
        );
    }
}

impl Scan2 {
    /// Mirrors Java `String.hashCode` so Scan2 can rebuild Java-like `HashSet`
    /// iteration order when constructing the candidate index.
    #[inline]
    fn java_string_hash(s: &str) -> i32 {
        // Java String.hashCode(): h = 31*h + ch (32-bit overflow).
        let mut h: i32 = 0;
        for ch in s.chars() {
            h = h.wrapping_mul(31).wrapping_add(ch as i32);
        }
        h
    }

    /// Applies the same spread step Java `HashMap` uses before bucket selection.
    #[inline]
    fn java_hash_spread(h: i32) -> u32 {
        let x = h as u32;
        x ^ (x >> 16)
    }

    /// Computes the minimum Java-style `HashSet` backing capacity for `n` items.
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

    /// Creates a Scan2 runner with CIRI3-compatible thresholds and bucket size.
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
            rescued_reads: 0,
            final_bsj_reads: 0,
        }
    }

    /// Binary-search lower bound used by the sorted candidate buckets.
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

    /// Returns the candidate slice covering one Java bucket.
    ///
    /// The index vectors are sorted once during `build_index`, so hot-path bucket
    /// scans can stay branch-light and avoid re-filtering the whole chromosome list.
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

    /// Overrides the per-thread page-cache eviction threshold.
    pub fn set_mem_limit(&mut self, limit: u64) {
        self.mem_limit = limit;
    }

    /// Drops Scan2-only indexing state once rescue is complete.
    ///
    /// Summary only needs the merged `fsj_map`, so the large candidate indexes and
    /// Scan1 read-id filter can be released explicitly before the final stage to
    /// keep peak RSS lower on whole-genome inputs.
    pub fn release_working_set(&mut self) {
        self.index1.clear();
        self.index1.shrink_to_fit();
        self.index2.clear();
        self.index2.shrink_to_fit();
        self.site_array1.clear();
        self.site_array1.shrink_to_fit();
        self.site_array2.clear();
        self.site_array2.shrink_to_fit();
        self.scan1_ids.clear();
        self.scan1_ids.shrink_to_fit();
    }

    /// Returns the temporary path used for one shard-local FSJ spill file.
    #[inline]
    fn shard_fsj_path(output_fsj: &str, shard_idx: usize) -> String {
        part_path(output_fsj, shard_idx)
    }

    /// Writes one shard-local FSJ map to disk so the main thread can merge it
    /// later without retaining every shard map in memory at once.
    fn write_fsj_shard(path: &str, local_fsj: &HashMap<String, i32>) -> Result<()> {
        let mut writer = BufWriter::with_capacity(256 * 1024, File::create(path)?);
        for (key, count) in local_fsj {
            writeln!(writer, "{}\t{}", key, count)?;
        }
        writer.flush()?;
        Ok(())
    }

    /// Loads Scan1 BSJ1 output and builds the de-duplicated Scan2 candidate index.
    ///
    /// This function intentionally mirrors Java's unique-site indexing semantics:
    /// the traversal order of candidates must match Java, because Scan2 returns on
    /// the first valid non-`2` tag it encounters.
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
            if p.len() < 10 {
                continue;
            }
            self.scan1_ids.insert(p[0].clone());
            let chr = p[3].clone();
            let site_infor = p[4..].join("\t");
            let seen = chr_circ_site_seen
                .entry(chr.clone())
                .or_insert_with(HashSet::new);
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
                    let arr: Vec<String> = insertion_list[idx]
                        .split('\t')
                        .map(|s| s.to_string())
                        .collect();
                    if arr.len() < 6 {
                        continue;
                    }
                    let site1 = arr[0].parse().unwrap_or(0);
                    let site2 = arr[1].parse().unwrap_or(0);
                    self.fsj_map
                        .entry(format!("{}\t{}\t{}", chr, arr[0], arr[1]))
                        .or_insert(0);
                    self.index1
                        .entry(chr.clone())
                        .or_insert_with(Vec::new)
                        .push(CandidateBreakpoint {
                            site: site1,
                            order: order_counter,
                            data: arr.clone(),
                        });
                    self.index2
                        .entry(chr.clone())
                        .or_insert_with(Vec::new)
                        .push(CandidateBreakpoint {
                            site: site2,
                            order: order_counter,
                            data: arr,
                        });
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

    /// Collects FSJ keys overlapped by a linear alignment span.
    ///
    /// This mirrors Java `GetFSJClass.getFSJ(...)`: only the two buckets touched
    /// by `[start_tem, end_tem]` are scanned, with reverse traversal on `num1`
    /// and forward traversal on `num2`. The broader lower-bound scan is cheaper
    /// to write but not Java-compatible, and it over-counts FSJs on chr1 while
    /// leaving BSJ rescue unchanged.
    fn collect_fsj_keys_in_range(
        &self,
        chr: &str,
        start_tem: i32,
        end_tem: i32,
        style: i32,
        out: &mut HashSet<String>,
    ) {
        let bucket_size = self.seq_len.max(1);
        let num1 = start_tem / bucket_size;
        let num2 = end_tem / bucket_size;
        if let Some(list1) = self.index1.get(chr) {
            if self.site_array1.get(chr).is_some_and(|s| s.contains(&num1)) {
                let (l1, r1) = Self::bucket_range(list1, num1, bucket_size);
                for idx in (l1..r1).rev() {
                    let cand = &list1[idx];
                    if cand.site >= start_tem {
                        if cand.site <= end_tem || style == 10 || style == 1 {
                            out.insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                        }
                    } else {
                        break;
                    }
                }
            }
            if num2 != num1 && self.site_array1.get(chr).is_some_and(|s| s.contains(&num2)) {
                let (l2, r2) = Self::bucket_range(list1, num2, bucket_size);
                for idx in l2..r2 {
                    let cand = &list1[idx];
                    if cand.site <= end_tem {
                        if cand.site >= start_tem || style == 10 || style == -1 {
                            out.insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                        }
                    } else {
                        break;
                    }
                }
            }
        }
        if let Some(list2) = self.index2.get(chr) {
            if self.site_array2.get(chr).is_some_and(|s| s.contains(&num1)) {
                let (l1, r1) = Self::bucket_range(list2, num1, bucket_size);
                for idx in (l1..r1).rev() {
                    let cand = &list2[idx];
                    if cand.site >= start_tem {
                        if cand.site <= end_tem || style == 10 || style == 1 {
                            out.insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                        }
                    } else {
                        break;
                    }
                }
            }
            if num2 != num1 && self.site_array2.get(chr).is_some_and(|s| s.contains(&num2)) {
                let (l2, r2) = Self::bucket_range(list2, num2, bucket_size);
                for idx in l2..r2 {
                    let cand = &list2[idx];
                    if cand.site <= end_tem {
                        if cand.site >= start_tem || style == 10 || style == -1 {
                            out.insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                        }
                    } else {
                        break;
                    }
                }
            }
        }
    }

    /// Dispatches to the SAM or BAM Scan2 implementation.
    pub fn run(
        &mut self,
        sam_file: &str,
        input_bsj1: &str,
        output_bsj: &str,
        output_bsj2: &str,
        output_fsj: &str,
        chr_tcga_map: &HashMap<String, String>,
    ) -> Result<()> {
        use crate::sam_bam::{detect_format, InputFormat};
        let format = detect_format(sam_file)?;
        match format {
            InputFormat::Sam => self.run_sam(
                sam_file,
                input_bsj1,
                output_bsj,
                output_bsj2,
                output_fsj,
                chr_tcga_map,
            ),
            InputFormat::Bam => self.run_bam(
                sam_file,
                input_bsj1,
                output_bsj,
                output_bsj2,
                output_fsj,
                chr_tcga_map,
            ),
        }
    }

    /// Runs Scan2 on SAM input.
    ///
    /// The SAM path is sharded by byte range, but each shard still goes through
    /// the same `process_group_view` logic as BAM so rescue behavior stays unified.
    pub fn run_sam(
        &mut self,
        sam_file: &str,
        input_bsj1: &str,
        output_bsj: &str,
        output_bsj2: &str,
        output_fsj: &str,
        chr_tcga_map: &HashMap<String, String>,
    ) -> Result<()> {
        let run_started = Instant::now();
        let profile = if Scan2Profile::enabled_from_env() {
            Some(Scan2Profile::default())
        } else {
            None
        };
        let profile_ref = profile.as_ref();
        let file = File::open(sam_file)?;
        let mmap = unsafe { Mmap::map(&file)? };
        let file_size = mmap.len();
        let num_threads = rayon::current_num_threads().max(1);
        let shard_size = file_size / num_threads;

        unsafe {
            libc::madvise(
                mmap.as_ptr() as *mut libc::c_void,
                file_size,
                libc::MADV_SEQUENTIAL,
            );
        }

        let pb = ProgressBar::new(file_size as u64);
        pb.set_style(ProgressStyle::default_bar().template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {percent:>3}% ({eta}) {msg}")?.progress_chars("#>-"));
        pb.set_message("");

        (0..num_threads).into_par_iter().try_for_each(|i| {
            let start = i * shard_size;
            let end = if i == num_threads - 1 {
                file_size
            } else {
                (i + 1) * shard_size
            };
            let shard_out = part_path(output_bsj2, i);
            let fsj_out = Self::shard_fsj_path(output_fsj, i);
            self.process_sam_shard_to_file(
                &mmap,
                start,
                end,
                chr_tcga_map,
                &pb,
                &shard_out,
                &fsj_out,
                profile_ref,
            )
        })?;

        pb.set_style(ProgressStyle::default_bar().template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {msg}")?.progress_chars("#>-"));
        pb.finish_with_message("Completed");
        let merge_started = Instant::now();
        let result =
            self.merge_shards_and_fsj(input_bsj1, output_bsj, output_bsj2, output_fsj, num_threads);
        if let Some(profile) = profile_ref {
            profile
                .merge_ns
                .fetch_add(merge_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile
                .wall_total_ns
                .fetch_add(run_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile.report();
            report_scan2_hg_profile();
        }
        result
    }

    /// Concatenates Scan1 and Scan2 BSJ outputs and merges per-shard FSJ spill files.
    ///
    /// The final `<prefix>.bsj` keeps one extra trailing column describing where
    /// each line came from (`scan1` or `scan2`). Summary ignores the extra field
    /// because all behaviorally relevant columns stay in their original positions.
    fn merge_shards_and_fsj(
        &mut self,
        input_bsj1: &str,
        output_bsj: &str,
        output_bsj2: &str,
        output_fsj: &str,
        num_threads: usize,
    ) -> Result<()> {
        let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(output_bsj)?);
        if let Ok(bsj1_file) = File::open(input_bsj1) {
            let mut bsj1_reader = BufReader::new(bsj1_file);
            let mut line = String::new();
            while bsj1_reader.read_line(&mut line)? != 0 {
                let trimmed = line.trim_end();
                if !trimmed.is_empty() {
                    writeln!(writer, "{}\tscan1", trimmed)?;
                }
                line.clear();
            }
        }
        let mut rescued_ids = HashSet::new();
        for i in 0..num_threads {
            let shard_path = part_path(output_bsj2, i);
            if let Ok(shard_file) = File::open(&shard_path) {
                let mut shard_reader = BufReader::new(shard_file);
                let mut line = String::new();
                while shard_reader.read_line(&mut line)? != 0 {
                    let trimmed = line.trim_end();
                    if !trimmed.is_empty() {
                        if let Some(tab_idx) = trimmed.find('\t') {
                            rescued_ids.insert(trimmed[..tab_idx].to_string());
                        }
                        writeln!(writer, "{}\tscan2", trimmed)?;
                    }
                    line.clear();
                }
            }
            let _ = std::fs::remove_file(shard_path);
        }
        for i in 0..num_threads {
            let fsj_path = Self::shard_fsj_path(output_fsj, i);
            if let Ok(fsj_file) = File::open(&fsj_path) {
                let mut fsj_reader = BufReader::new(fsj_file);
                let mut line = String::new();
                while fsj_reader.read_line(&mut line)? != 0 {
                    let trimmed = line.trim_end();
                    if !trimmed.is_empty() {
                        let mut parts = trimmed.rsplitn(2, '\t');
                        let count = parts
                            .next()
                            .and_then(|s| s.parse::<i32>().ok())
                            .unwrap_or(0);
                        if let Some(key) = parts.next() {
                            *self.fsj_map.entry(key.to_string()).or_insert(0) += count;
                        }
                    }
                    line.clear();
                }
            }
            let _ = std::fs::remove_file(fsj_path);
        }
        writer.flush()?;
        self.rescued_reads = rescued_ids.len();
        self.final_bsj_reads = self.scan1_ids.len() + self.rescued_reads;
        let _ = std::fs::remove_file(input_bsj1);
        Ok(())
    }

    /// Runs Scan2 on BAM input using BGZF-aware sharding.
    pub fn run_bam(
        &mut self,
        bam_file: &str,
        input_bsj1: &str,
        output_bsj: &str,
        output_bsj2: &str,
        output_fsj: &str,
        chr_tcga_map: &HashMap<String, String>,
    ) -> Result<()> {
        use noodles::bam;
        let run_started = Instant::now();
        let profile = if Scan2Profile::enabled_from_env() {
            Some(Scan2Profile::default())
        } else {
            None
        };
        let profile_ref = profile.as_ref();
        let file = File::open(bam_file)?;
        let file_size = std::fs::metadata(bam_file)?.len();
        let mmap = unsafe { Mmap::map(&file)? };
        let num_threads = bam_shard_count(mmap.len(), rayon::current_num_threads());
        let shard_size = mmap.len() / num_threads;
        let mut reader = bam::io::Reader::new(file);
        let header = reader.read_header()?;

        let pb = ProgressBar::new(file_size);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {percent:>3}% ({eta}) {msg}")?
                .progress_chars("#>-"),
        );
        pb.set_message("");
        pb.enable_steady_tick(Duration::from_millis(120));

        unsafe {
            libc::madvise(
                mmap.as_ptr() as *mut libc::c_void,
                mmap.len(),
                libc::MADV_SEQUENTIAL,
            );
        }

        let header_ref = &header;
        (0..num_threads).into_par_iter().try_for_each(|i| {
            let start = i * shard_size;
            let end = if i == num_threads - 1 {
                mmap.len()
            } else {
                (i + 1) * shard_size
            };
            let shard_out = part_path(output_bsj2, i);
            let fsj_out = Self::shard_fsj_path(output_fsj, i);
            self.process_bam_shard_to_file(
                &mmap,
                start,
                end,
                header_ref,
                chr_tcga_map,
                &pb,
                &shard_out,
                &fsj_out,
                profile_ref,
            )
        })?;

        pb.set_position(file_size);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {msg}")?
                .progress_chars("#>-"),
        );
        pb.finish_with_message("Completed");
        let merge_started = Instant::now();
        let result =
            self.merge_shards_and_fsj(input_bsj1, output_bsj, output_bsj2, output_fsj, num_threads);
        if let Some(profile) = profile_ref {
            profile
                .merge_ns
                .fetch_add(merge_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile
                .wall_total_ns
                .fetch_add(run_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile.report();
            report_scan2_hg_profile();
        }
        result
    }

    /// Processes one BAM shard and emits rescued BSJ2 lines plus one local FSJ spill file.
    ///
    /// The small buffer reuses inside this function are intentional low-risk
    /// optimizations: they reduce per-record allocation churn without changing the
    /// order or content of any candidate checks.
    fn process_bam_shard_to_file(
        &self,
        mmap: &Mmap,
        start: usize,
        end: usize,
        header: &sam::Header,
        chr_tcga_map: &HashMap<String, String>,
        pb: &ProgressBar,
        out_path: &str,
        fsj_path: &str,
        profile: Option<&Scan2Profile>,
    ) -> Result<()> {
        use noodles::bam;
        let shard_started = profile.map(|_| Instant::now());
        let mut writer = BufWriter::with_capacity(256 * 1024, File::create(out_path)?);
        let mut local_fsj = HashMap::new();
        let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
        let block_start = if start == 0 {
            0
        } else {
            let mut found = None;
            for i in start..mmap.len().saturating_sub(3) {
                if &mmap[i..i + 4] == b"\x1f\x8b\x08\x04" {
                    found = Some(i);
                    break;
                }
            }
            found.unwrap_or(mmap.len())
        };
        let pos = if start == 0 {
            0
        } else {
            // Decode from the previous BGZF block so we can distinguish a true
            // new read group from one that started in the earlier shard.
            let mut found = None;
            let mut i = block_start
                .saturating_sub(1)
                .min(mmap.len().saturating_sub(4));
            loop {
                if i + 3 < mmap.len() && &mmap[i..i + 4] == b"\x1f\x8b\x08\x04" {
                    found = Some(i);
                    break;
                }
                if i == 0 {
                    break;
                }
                i -= 1;
            }
            found.unwrap_or(0)
        };
        if pos >= mmap.len() {
            return Self::write_fsj_shard(fsj_path, &local_fsj);
        }

        let mut reader = bam::io::Reader::new(&mmap[pos..]);
        if start == 0 {
            let _ = reader.read_header()?;
        }
        let mut record = bam::Record::default();
        let mut current_id: Vec<u8> = Vec::new();
        let mut alignments: Vec<AlignmentRecord> = Vec::with_capacity(16);
        let mut stand_map: HashMap<i32, (char, Cow<str>)> = HashMap::with_capacity(4);
        let mut one_read_key: i32 = -1;
        let mut crossed_start = start == 0;
        let mut leading_partial_id: Option<Vec<u8>> = None;
        let mut last_compressed_pos = block_start.saturating_sub(pos);
        let mut last_evicted_pos = pos;
        // Keep page-cache eviction on a small sliding window. Large thresholds let
        // the whole mapped BAM drift into RSS before any MADV_DONTNEED call fires.
        let eviction_threshold =
            (self.mem_limit as usize / 8).clamp(8 * 1024 * 1024, 64 * 1024 * 1024);
        let mut cigar_buf = String::with_capacity(64);
        let mut seq_buf = String::with_capacity(256);
        let mut res_batch = Vec::new();

        while reader.read_record(&mut record)? != 0 {
            if let Some(profile) = profile {
                profile.records.fetch_add(1, Ordering::Relaxed);
            }
            let curr_c_pos = reader.get_ref().virtual_position().compressed() as usize;
            pb.inc((curr_c_pos - last_compressed_pos) as u64);
            last_compressed_pos = curr_c_pos;

            if curr_c_pos - (last_evicted_pos - pos) > eviction_threshold {
                advise_dontneed(
                    mmap,
                    last_evicted_pos,
                    curr_c_pos - (last_evicted_pos - pos),
                );
                last_evicted_pos = pos + curr_c_pos;
            }

            let read_id = record
                .name()
                .ok_or_else(|| anyhow::anyhow!("Missing read name"))?;
            let abs_c_pos = pos + curr_c_pos;
            if !crossed_start {
                if abs_c_pos < block_start {
                    leading_partial_id = Some(read_id.to_vec());
                    continue;
                }
                crossed_start = true;
            }
            if let Some(partial_id) = &leading_partial_id {
                if read_id.to_vec() == *partial_id {
                    // Let the earlier shard own the cross-boundary read group
                    // completely; otherwise Scan2 can silently drop the first
                    // full read group visible in this shard.
                    continue;
                }
                leading_partial_id = None;
            }

            if read_id.to_vec() != current_id {
                if !current_id.is_empty() {
                    let id_str = String::from_utf8_lossy(&current_id);
                    res_batch.clear();
                    let write_started = profile.map(|_| Instant::now());
                    self.process_group_view(
                        &id_str,
                        &alignments,
                        &stand_map,
                        &mut res_batch,
                        &mut local_fsj,
                        chr_tcga_map,
                        &mut validator,
                        profile,
                    )?;
                    for line in &res_batch {
                        writeln!(writer, "{}", line)?;
                    }
                    if let (Some(profile), Some(write_started)) = (profile, write_started) {
                        profile.write_ns.fetch_add(
                            write_started.elapsed().as_nanos() as u64,
                            Ordering::Relaxed,
                        );
                    }
                    if abs_c_pos > end {
                        current_id.clear();
                        break;
                    }
                }
                current_id = read_id.to_vec();
                alignments.clear();
                stand_map.clear();
                one_read_key = -1;
            }

            let flag = i32::from(u16::from(record.flags()));
            let chrom = match record.reference_sequence(header) {
                Some(Ok((name, _))) => String::from_utf8_lossy(name).to_string(),
                _ => "*".to_string(),
            };
            let start_pos = record
                .alignment_start()
                .transpose()?
                .map(|p| p.get() as i32)
                .unwrap_or(0);
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
                let _ = write!(&mut cigar_buf, "{}{}", op.len(), op_char);
            }
            seq_buf.clear();
            for b in record.sequence().iter() {
                seq_buf.push(char::from(b));
            }
            let seq = seq_buf.clone();
            let s_idx = if flag & 0x40 != 0 { 1 } else { 0 };
            let st_c = if flag & 0x10 != 0 { '1' } else { '0' };
            if s_idx != one_read_key {
                one_read_key = s_idx;
                // Java parity: BAM Scan2 rebuilds the per-mate alignment list each
                // time the iterator switches between R1 and R2 within the same
                // read group. Keeping every earlier alignment for that mate looks
                // harmless, but it changes which candidate payloads contribute to
                // `temFSJId` on complex supplementary-heavy reads.
                alignments.retain(|a| {
                    let idx = if a.flag & 0x40 != 0 { 1 } else { 0 };
                    idx != s_idx
                });
                // Java parity: `standMap` is overwritten on mate switches with the
                // current record's sequence; it is not a "longest-sequence wins"
                // cache in Scan2. This only affects Scan2 candidate validation,
                // not Scan1 representative-sequence handling.
                if !seq.is_empty() && seq != "*" {
                    stand_map.insert(s_idx, (st_c, Cow::Owned(seq.clone())));
                }
            }
            alignments.push(AlignmentRecord {
                flag,
                chrom: Cow::Owned(chrom),
                pos: start_pos,
                mapq,
                cigar: Cow::Owned(cigar_buf.clone()),
                seq: Cow::Owned(seq),
            });
        }
        if !current_id.is_empty() {
            let id_str = String::from_utf8_lossy(&current_id);
            res_batch.clear();
            let write_started = profile.map(|_| Instant::now());
            self.process_group_view(
                &id_str,
                &alignments,
                &stand_map,
                &mut res_batch,
                &mut local_fsj,
                chr_tcga_map,
                &mut validator,
                profile,
            )?;
            for line in &res_batch {
                writeln!(writer, "{}", line)?;
            }
            if let (Some(profile), Some(write_started)) = (profile, write_started) {
                profile
                    .write_ns
                    .fetch_add(write_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            }
        }
        if last_compressed_pos > (last_evicted_pos - pos) {
            advise_dontneed(
                mmap,
                last_evicted_pos,
                last_compressed_pos - (last_evicted_pos - pos),
            );
        }
        writer.flush()?;
        if let (Some(profile), Some(shard_started)) = (profile, shard_started) {
            profile
                .shard_total_ns
                .fetch_add(shard_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
        Self::write_fsj_shard(fsj_path, &local_fsj)
    }

    /// Processes one SAM shard into BSJ2 lines plus one local FSJ spill file.
    ///
    /// This path stays text-based for SAM, but reuses the exact same read-group
    /// rescue logic as BAM after records have been assembled.
    fn process_sam_shard_to_file<'a>(
        &self,
        mmap: &'a Mmap,
        start: usize,
        end: usize,
        chr_tcga_map: &HashMap<String, String>,
        pb: &ProgressBar,
        out_path: &str,
        fsj_path: &str,
        profile: Option<&Scan2Profile>,
    ) -> Result<()> {
        let shard_started = profile.map(|_| Instant::now());
        let mut writer = BufWriter::with_capacity(256 * 1024, File::create(out_path)?);
        let mut local_fsj = HashMap::new();
        let mut validator = IsBSJHg2::new(self.linear_range_size_min, self.min_mapq_uni);
        let mut pos = if start == 0 {
            0
        } else {
            memchr(b'\n', &mmap[start..])
                .map(|p| start + p + 1)
                .unwrap_or(mmap.len())
        };
        if pos >= mmap.len() {
            return Self::write_fsj_shard(fsj_path, &local_fsj);
        }

        let mut current_id: &[u8] = &[];
        let mut alignments: Vec<AlignmentRecord<'a>> = Vec::with_capacity(16);
        let mut stand_map: HashMap<i32, (char, Cow<'a, str>)> = HashMap::with_capacity(4);
        let mut one_read_key: i32 = -1;
        let mut last_evicted_pos = pos;
        let eviction_threshold =
            (self.mem_limit as usize / 8).clamp(8 * 1024 * 1024, 64 * 1024 * 1024);
        let mut res_batch = Vec::new();

        while pos < mmap.len() {
            let line_end = memchr(b'\n', &mmap[pos..])
                .map(|p| pos + p)
                .unwrap_or(mmap.len());
            let line = &mmap[pos..line_end];
            if line.is_empty() {
                pb.inc(1);
                pos = line_end + 1;
                continue;
            }
            if pos - last_evicted_pos > eviction_threshold {
                advise_dontneed(mmap, last_evicted_pos, pos - last_evicted_pos);
                last_evicted_pos = pos;
            }
            if line[0] == b'@' {
                pb.inc((line_end - pos + 1) as u64);
                pos = line_end + 1;
                continue;
            }

            let mut cols = line.split(|&b| b == b'\t');
            let read_id = cols.next().unwrap();
            if let Some(profile) = profile {
                profile.records.fetch_add(1, Ordering::Relaxed);
            }
            if read_id != current_id {
                if !current_id.is_empty() {
                    res_batch.clear();
                    let write_started = profile.map(|_| Instant::now());
                    self.process_group_view(
                        unsafe { std::str::from_utf8_unchecked(current_id) },
                        &alignments,
                        &stand_map,
                        &mut res_batch,
                        &mut local_fsj,
                        chr_tcga_map,
                        &mut validator,
                        profile,
                    )?;
                    for l in &res_batch {
                        writeln!(writer, "{}", l)?;
                    }
                    if let (Some(profile), Some(write_started)) = (profile, write_started) {
                        profile.write_ns.fetch_add(
                            write_started.elapsed().as_nanos() as u64,
                            Ordering::Relaxed,
                        );
                    }
                }
                if pos >= end {
                    break;
                }
                current_id = read_id;
                alignments.clear();
                stand_map.clear();
                one_read_key = -1;
            }
            let flag = fast_parse_i32(cols.next().unwrap_or(b"0"));
            let chrom = unsafe { std::str::from_utf8_unchecked(cols.next().unwrap_or(b"*")) };
            let start_pos = fast_parse_i32(cols.next().unwrap_or(b"0"));
            let mapq = fast_parse_i32(cols.next().unwrap_or(b"0"));
            let cigar = unsafe { std::str::from_utf8_unchecked(cols.next().unwrap_or(b"*")) };
            cols.next();
            cols.next();
            cols.next();
            let seq = unsafe { std::str::from_utf8_unchecked(cols.next().unwrap_or(b"*")).trim() };
            let s_idx = if flag & 0x40 != 0 { 1 } else { 0 };
            let st_c = if flag & 0x10 != 0 { '1' } else { '0' };
            if s_idx != one_read_key {
                one_read_key = s_idx;
                alignments.retain(|a| {
                    let idx = if a.flag & 0x40 != 0 { 1 } else { 0 };
                    idx != s_idx
                });
                if !seq.is_empty() && seq != "*" {
                    stand_map.insert(s_idx, (st_c, Cow::Borrowed(seq)));
                }
            }
            alignments.push(AlignmentRecord {
                flag,
                chrom: Cow::Borrowed(chrom),
                pos: start_pos,
                mapq,
                cigar: Cow::Borrowed(cigar),
                seq: Cow::Borrowed(seq),
            });
            pb.inc((line_end - pos + 1) as u64);
            pos = line_end + 1;
        }
        if pos > last_evicted_pos {
            advise_dontneed(mmap, last_evicted_pos, pos - last_evicted_pos);
        }
        writer.flush()?;
        if let (Some(profile), Some(shard_started)) = (profile, shard_started) {
            profile
                .shard_total_ns
                .fetch_add(shard_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
        Self::write_fsj_shard(fsj_path, &local_fsj)
    }

    /// Rescues BSJ candidates and counts overlapping FSJs for one read group.
    ///
    /// This is the Scan2 hot path. The preserved optimizations are intentionally
    /// modest:
    /// - mate-segment lookup is hoisted once per segment side
    /// - oriented read/mate strings are computed once per alignment, not once per
    ///   candidate
    /// - bucket traversal order remains Java-compatible
    pub(crate) fn process_group_view<'a>(
        &self,
        id: &str,
        alignments: &[AlignmentRecord<'a>],
        stand_map: &HashMap<i32, (char, Cow<'a, str>)>,
        results: &mut Vec<String>,
        local_fsj: &mut HashMap<String, i32>,
        chr_tcga_map: &HashMap<String, String>,
        is_bsj_hg2: &mut IsBSJHg2,
        profile: Option<&Scan2Profile>,
    ) -> Result<()> {
        let started = profile.map(|_| Instant::now());
        let trace_read = should_trace_read(id);
        let trace_all_candidates =
            trace_read && std::env::var("CIRI_TRACE_ALL_CANDS").ok().as_deref() == Some("1");
        if self.scan1_ids.contains(id) {
            return Ok(());
        }
        let mut segments: HashMap<i32, Vec<&AlignmentRecord<'a>>> = HashMap::new();
        for aln in alignments {
            segments
                .entry(if aln.flag & 0x40 != 0 { 1 } else { 0 })
                .or_insert_with(Vec::new)
                .push(aln);
        }
        let mut tem_fsj_keys = HashSet::new();
        for seg_idx in [0_i32, 1_i32] {
            let Some(seg_alns) = segments.get(&seg_idx) else {
                continue;
            };
            let (read_strand, seq) = match stand_map.get(&seg_idx) {
                Some(&(st, ref s)) => (st, s.as_ref()),
                None => continue,
            };
            let slen = seq.len() as i32;
            // Mate alignment lookup is invariant for all candidates generated from
            // this read side, so keep the reference once and reuse it below.
            let mate_seg_alns = segments.get(&(1 - seg_idx));
            for aln in seg_alns {
                let chr = aln.chrom.as_ref();
                if !self.index1.contains_key(chr) {
                    continue;
                }
                let cigar_ref: Cow<'_, str> = if aln.cigar.contains('H') {
                    Cow::Owned(aln.cigar.replace('H', "S"))
                } else {
                    Cow::Borrowed(aln.cigar.as_ref())
                };
                if cigar_ref == "*" {
                    continue;
                }
                let c = misd(cigar_ref.as_ref(), slen);
                if cigar_ref == format!("{}M", slen) {
                    let start_tem = aln.pos + 6;
                    let end_tem = aln.pos + slen - 7;
                    self.collect_fsj_keys_in_range(chr, start_tem, end_tem, 0, &mut tem_fsj_keys);
                    continue;
                }
                if c[0] == -1 || c[0] == 10 {
                    if let Some(list) = self.index1.get(chr) {
                        let curr_strand = if aln.flag & 0x10 != 0 { '1' } else { '0' };
                        // Orientation depends on the alignment strand of the
                        // current record, so the normalized strings are hoisted to
                        // alignment scope but not further.
                        let str_e_owned = if aln.flag & 0x10 != 0 {
                            if read_strand == '0' {
                                reverse_complement(seq)
                            } else {
                                seq.to_string()
                            }
                        } else if read_strand == '0' {
                            seq.to_string()
                        } else {
                            reverse_complement(seq)
                        };
                        let p_str_owned =
                            if let Some(&(p_strand, ref p_seq)) = stand_map.get(&(1 - seg_idx)) {
                                if p_strand != curr_strand {
                                    p_seq.to_string()
                                } else {
                                    reverse_complement(p_seq)
                                }
                            } else {
                                String::new()
                            };
                        let new_num_site = aln.pos;
                        // Java parity: fixed seqLen from constructor.
                        let bucket_size = self.seq_len.max(1);
                        let num1 = (new_num_site - 6) / bucket_size;
                        let num2 = (new_num_site + 6) / bucket_size;

                        let mut eval_candidate = |cand: &CandidateBreakpoint| -> Result<bool> {
                            let e_idx = c[1] + (cand.site - new_num_site);
                            if e_idx <= 0 || e_idx > slen {
                                return Ok(false);
                            }
                            let str_f = &str_e_owned[0..e_idx as usize];
                            let str3 = if c[0] == 10 {
                                let si = slen - c[2];
                                if si >= 0 && si <= slen {
                                    &str_e_owned[si as usize..]
                                } else {
                                    ""
                                }
                            } else {
                                "*"
                            };
                            let site1 = cand.data[0].parse::<i32>().unwrap_or(0);
                            let site2 = cand.data[1].parse::<i32>().unwrap_or(0);
                            let mut s2_ok = if mate_seg_alns.is_some() { 0 } else { 1 };
                            if let Some(mate_seg_alns) = mate_seg_alns {
                                for mate_aln in mate_seg_alns {
                                    if mate_aln.chrom.as_ref() == chr
                                        && mate_aln.mapq >= self.min_mapq_uni
                                    {
                                        let c_ano = misd(&mate_aln.cigar, slen);
                                        let mate_strand =
                                            if mate_aln.flag & 0x10 != 0 { '1' } else { '0' };
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
                                p_str_owned.clone(),
                                str3.to_string(),
                                s2_ok.to_string(),
                                cand.data[2].clone(),
                                cand.data[3].clone(),
                                cand.data[4].clone(),
                                aln.mapq.to_string(),
                            ];
                            let validator_started = profile.map(|_| Instant::now());
                            let tag =
                                is_bsj_hg2.is_bsj_hg2(&circ_c, chr_tcga_map.get(chr).unwrap());
                            if let Some(profile) = profile {
                                profile.candidate_checks.fetch_add(1, Ordering::Relaxed);
                                if let Some(validator_started) = validator_started {
                                    profile.validator_call_ns.fetch_add(
                                        validator_started.elapsed().as_nanos() as u64,
                                        Ordering::Relaxed,
                                    );
                                }
                            }
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
                                tem_fsj_keys
                                    .insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                                Ok(false)
                            } else if tag != "2" {
                                if let Some(profile) = profile {
                                    profile.candidate_hits.fetch_add(1, Ordering::Relaxed);
                                }
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
                        if num1 > 0 && self.site_array1.get(chr).is_some_and(|s| s.contains(&num1))
                        {
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
                        if num2 != num1
                            && self.site_array1.get(chr).is_some_and(|s| s.contains(&num2))
                        {
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
                        let curr_strand = if aln.flag & 0x10 != 0 { '1' } else { '0' };
                        let str_e_owned = if aln.flag & 0x10 != 0 {
                            if read_strand == '0' {
                                reverse_complement(seq)
                            } else {
                                seq.to_string()
                            }
                        } else if read_strand == '0' {
                            seq.to_string()
                        } else {
                            reverse_complement(seq)
                        };
                        let p_str_owned =
                            if let Some(&(p_strand, ref p_seq)) = stand_map.get(&(1 - seg_idx)) {
                                if p_strand != curr_strand {
                                    p_seq.to_string()
                                } else {
                                    reverse_complement(p_seq)
                                }
                            } else {
                                String::new()
                            };
                        // Java parity: fixed seqLen from constructor.
                        let bucket_size = self.seq_len.max(1);
                        let num1 = (new_site - 6) / bucket_size;
                        let num2 = (new_site + 6) / bucket_size;

                        let mut eval_candidate = |cand: &CandidateBreakpoint| -> Result<bool> {
                            let s_idx = if c[0] == 10 {
                                // Java parity: SMS right-anchor slicing uses
                                // `seqLen - tail_soft_clip + bias`, not
                                // `start_soft_clip + mapped_len + bias`.
                                // The latter drifts whenever the middle CIGAR
                                // contains insertions or deletions.
                                slen - c[2] + (cand.site - new_site)
                            } else {
                                c[1] + (cand.site - new_site)
                            };
                            if s_idx < 0 || s_idx >= slen {
                                return Ok(false);
                            }
                            let str_f = &str_e_owned[s_idx as usize..];
                            let str3 = if c[0] == 10 {
                                &str_e_owned[0..c[1] as usize]
                            } else {
                                "*"
                            };
                            let site1 = cand.data[0].parse::<i32>().unwrap_or(0);
                            let site2 = cand.data[1].parse::<i32>().unwrap_or(0);
                            let mut s2_ok = if mate_seg_alns.is_some() { 0 } else { 1 };
                            if let Some(mate_seg_alns) = mate_seg_alns {
                                for mate_aln in mate_seg_alns {
                                    if mate_aln.chrom.as_ref() == chr
                                        && mate_aln.mapq >= self.min_mapq_uni
                                    {
                                        let c_ano = misd(&mate_aln.cigar, slen);
                                        let mate_strand =
                                            if mate_aln.flag & 0x10 != 0 { '1' } else { '0' };
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
                                p_str_owned.clone(),
                                str3.to_string(),
                                s2_ok.to_string(),
                                cand.data[2].clone(),
                                cand.data[3].clone(),
                                cand.data[4].clone(),
                                aln.mapq.to_string(),
                            ];
                            let validator_started = profile.map(|_| Instant::now());
                            let tag =
                                is_bsj_hg2.is_bsj_hg2(&circ_c, chr_tcga_map.get(chr).unwrap());
                            if let Some(profile) = profile {
                                profile.candidate_checks.fetch_add(1, Ordering::Relaxed);
                                if let Some(validator_started) = validator_started {
                                    profile.validator_call_ns.fetch_add(
                                        validator_started.elapsed().as_nanos() as u64,
                                        Ordering::Relaxed,
                                    );
                                }
                            }
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
                                tem_fsj_keys
                                    .insert(format!("{}\t{}\t{}", chr, cand.data[0], cand.data[1]));
                                Ok(false)
                            } else if tag != "2" {
                                if let Some(profile) = profile {
                                    profile.candidate_hits.fetch_add(1, Ordering::Relaxed);
                                }
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
                        if num2 != num1
                            && self.site_array2.get(chr).is_some_and(|s| s.contains(&num2))
                        {
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
        for k in tem_fsj_keys {
            *local_fsj.entry(k).or_insert(0) += 1;
        }
        if let (Some(profile), Some(started)) = (profile, started) {
            profile
                .group_process_ns
                .fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile.groups.fetch_add(1, Ordering::Relaxed);
        }
        Ok(())
    }
}
