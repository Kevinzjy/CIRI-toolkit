//! Scan1: first-pass BSJ candidate discovery.
//!
//! This module owns the initial read-group scan over SAM/BAM input and emits the
//! Java-compatible BSJ intermediate file consumed by Scan2 and Summary.
//! Performance work here is deliberately limited to ingestion, sharding, and
//! temporary allocation control; the candidate semantics still follow Java CIRI3.

use crate::annotation::Annotation;
use crate::is_bsj_hg2::{java_substring, report_scan1_hg_profile, IsBSJHg2};
use crate::misd::misd;
use crate::utils::{bam_shard_count, part_path, AlignmentRecord};
use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use memmap2::Mmap;
use noodles::sam::{
    self,
    alignment::{record::Sequence as _, Record as _},
};
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Seek, Write};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc;
use std::thread;
use std::time::{Duration, Instant};

/// First-pass BSJ detector.
///
/// Design constraints:
/// - Behavior is anchored to Java CIRI3, so this type keeps the same two-pass
///   boundaries and candidate layout.
/// - Performance work is intentionally concentrated in record ingestion and
///   `process_group_view`, because profiling showed merge and final write are not
///   meaningful bottlenecks.
pub struct Scan1 {
    pub min_mapq_uni: i32,
    pub max_circle: i32,
    pub min_circle: i32,
    pub linear_range_size_min: i32,
    pub mem_limit: u64,
    pub read_len: i32,
    /// Total read groups observed during Scan1 input traversal.
    pub mapped_reads: u64,
    /// Unique read IDs emitted into the merged Scan1 BSJ file.
    pub bsj1_reads: usize,
}

type OwnedAlignmentRecord = AlignmentRecord<'static>;
type OwnedStandMap = HashMap<i32, (char, Cow<'static, str>)>;

/// Owned SAM read-group payload passed from the sequential parser thread to the
/// parallel candidate-evaluation workers.
///
/// The SAM path uses this owned form so parsing stays single-source-of-truth
/// through `noodles`, while the expensive BSJ judgment can still run in parallel.
struct SamOwnedGroup {
    read_id: String,
    group: [Vec<OwnedAlignmentRecord>; 2],
    stand_map: OwnedStandMap,
    align_num: usize,
}

/// Lightweight Scan1 traversal summary used for user-facing stage logs.
struct Scan1TraversalStats {
    max_read_len: i32,
    mapped_reads: u64,
}

/// One shard-local Scan1 summary.
struct Scan1ShardStats {
    max_read_len: i32,
    mapped_reads: u64,
}

/// Advises the kernel that an already-processed BAM byte range can be evicted.
///
/// The offset is aligned manually because `madvise` requires page alignment.
/// This is a memory-footprint optimization, not a correctness requirement.
fn advise_dontneed_aligned(mmap: &Mmap, offset: usize, len: usize) {
    if len == 0 {
        return;
    }
    // Page alignment is mandatory for madvise (typically 4096 bytes)
    let page_size = 4096;
    let aligned_offset = (offset / page_size) * page_size;
    let aligned_len = ((offset + len + page_size - 1) / page_size) * page_size - aligned_offset;

    unsafe {
        let ptr = mmap.as_ptr().add(aligned_offset);
        libc::madvise(ptr as *mut libc::c_void, aligned_len, libc::MADV_DONTNEED);
    }
}

/// Returns whether a read is selected for targeted parity tracing.
///
/// The hook is intentionally cheap so it can stay in hot paths without affecting
/// normal runs when `CIRI_TRACE_READS` is unset.
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

/// Emits shard-boundary trace lines for targeted BAM parity debugging.
///
/// This stays behind `CIRI_TRACE_READS` so it can remain in the codebase without
/// affecting normal runs. The output is intentionally narrow: only shard
/// ownership decisions around read-group boundaries are logged.
fn trace_bam_shard_event(
    shard_idx: usize,
    read_id: &str,
    stage: &str,
    start: usize,
    block_start: usize,
    end: usize,
    abs_c_pos: usize,
) {
    if should_trace_read(read_id) {
        eprintln!(
            "[TRACE_SCAN1_SHARD] shard={} stage={} id={} start={} block_start={} end={} abs_c_pos={}",
            shard_idx, stage, read_id, start, block_start, end, abs_c_pos
        );
    }
}

/// Keeps the best representative sequence for one mate inside a read group.
///
/// BAM/SAM records that carry hard clips expose a shorter `SEQ` than the
/// corresponding soft-clipped alignment of the same template. BSJ judgment later
/// slices the representative read sequence using clip-derived offsets, so taking
/// the first record verbatim can accidentally anchor Scan1 on a truncated
/// hard-clipped sequence. To stay compatible with CIRI3's expectation of using
/// the full read sequence, we retain the longest observed `SEQ` for each mate.
#[inline]
fn update_best_stand_seq<'a>(
    stand_map: &mut HashMap<i32, (char, Cow<'a, str>)>,
    s_idx: i32,
    st_c: char,
    seq: Cow<'a, str>,
) {
    match stand_map.get(&s_idx) {
        Some((_, existing_seq)) if existing_seq.len() >= seq.len() => {}
        _ => {
            stand_map.insert(s_idx, (st_c, seq));
        }
    }
}

/// Optional Scan1 profiler used only when `CIRI_PROFILE_SCAN1` is enabled.
///
/// The counters are coarse on purpose: they are cheap enough to leave in the hot
/// path, but still isolate whether time is spent in group assembly, BSJ judgment,
/// or write/merge scaffolding.
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
    /// Checks whether release profiling is enabled for the current run.
    fn enabled_from_env() -> bool {
        matches!(std::env::var("CIRI_PROFILE_SCAN1"), Ok(v) if !v.is_empty() && v != "0")
    }

    /// Emits the aggregated Scan1 timing summary.
    ///
    /// This keeps the profiler output stable across optimization rounds so A/B
    /// runs remain comparable.
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
    /// Normalizes the split-CIGAR text that gets written into the persisted BSJ file.
    ///
    /// Java's downstream `Summary` and `Misd` code effectively reasons about
    /// hard clips as soft clips (`H -> S`). Keeping raw `H` in the persisted
    /// BSJ1 text makes Summary-side CIGAR parsing diverge on real BAMs that
    /// contain hard-clipped supplementary alignments, even though the hot-path
    /// validator already treated those clips as soft clips during Scan1. The
    /// on-disk BSJ1 protocol therefore needs the same normalization.
    #[inline]
    fn normalize_bsj_cigar(cigar: &str) -> std::borrow::Cow<'_, str> {
        if cigar.contains('H') {
            std::borrow::Cow::Owned(cigar.replace('H', "S"))
        } else {
            std::borrow::Cow::Borrowed(cigar)
        }
    }

    /// Creates a Scan1 runner with CIRI3-compatible thresholds.
    pub fn new(
        min_mapq_uni: i32,
        min_circle: i32,
        max_circle: i32,
        linear_range_size_min: i32,
    ) -> Self {
        Self {
            min_mapq_uni,
            max_circle,
            min_circle,
            linear_range_size_min,
            mem_limit: 2 * 1024 * 1024 * 1024,
            read_len: 0,
            mapped_reads: 0,
            bsj1_reads: 0,
        }
    }

    /// Overrides the per-thread page-cache eviction threshold.
    pub fn set_mem_limit(&mut self, limit: u64) {
        self.mem_limit = limit;
    }

    /// Dispatches to the SAM or BAM implementation while preserving one Scan1 API.
    pub fn run(
        &mut self,
        sam_file: &str,
        bsj_path: &str,
        fasta_map: &HashMap<String, String>,
        annotation: &Annotation,
    ) -> Result<HashSet<String>> {
        use crate::sam_bam::{detect_format, InputFormat};
        let format = detect_format(sam_file)?;
        match format {
            InputFormat::Sam => self.run_sam(sam_file, bsj_path, fasta_map, annotation),
            InputFormat::Bam => self.run_bam(sam_file, bsj_path, fasta_map, annotation),
        }
    }

    /// Runs Scan1 on SAM input.
    ///
    /// The current SAM implementation keeps parsing sequential and parallelizes at
    /// the read-group level. This was chosen deliberately because an earlier
    /// hand-written text fast path was faster but produced SAM/BAM mismatches.
    pub fn run_sam(
        &mut self,
        sam_file: &str,
        bsj_path: &str,
        fasta_map: &HashMap<String, String>,
        annotation: &Annotation,
    ) -> Result<HashSet<String>> {
        let file_size = std::fs::metadata(sam_file)?.len();
        let pb = ProgressBar::new(file_size);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {percent:>3}% ({eta}) {msg}")?
                .progress_chars("#>-"),
        );
        pb.set_message("");
        pb.enable_steady_tick(Duration::from_millis(120));

        let shard_out = part_path(bsj_path, 0);
        let stats =
            self.process_sam_file_to_file(sam_file, fasta_map, annotation, &pb, &shard_out)?;
        self.read_len = stats.max_read_len;
        self.mapped_reads = stats.mapped_reads;

        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {msg}")?
                .progress_chars("#>-"),
        );
        pb.finish_with_message("");
        let scan1_ids = self.merge_and_collect_ids(bsj_path, 1)?;
        self.bsj1_reads = scan1_ids.len();
        Ok(scan1_ids)
    }

    /// Streams a SAM file into one shard output while batching read groups for
    /// parallel candidate evaluation.
    ///
    /// There is only one shard on the SAM path because correctness depends on
    /// preserving full read-group boundaries from the sequential parser.
    fn process_sam_file_to_file(
        &self,
        sam_file: &str,
        fasta_map: &HashMap<String, String>,
        annotation: &Annotation,
        pb: &ProgressBar,
        out_path: &str,
    ) -> Result<Scan1TraversalStats> {
        let mut writer = BufWriter::with_capacity(256 * 1024, File::create(out_path)?);
        let batch_size = (rayon::current_num_threads().max(1) * 256).max(1024);
        let (tx, rx) =
            mpsc::sync_channel::<Vec<SamOwnedGroup>>(rayon::current_num_threads().max(2));
        let pb_clone = pb.clone();

        let stats = thread::scope(|scope| -> Result<Scan1TraversalStats> {
            let producer =
                scope.spawn(|| self.stream_sam_group_batches(sam_file, &pb_clone, tx, batch_size));

            for mut batch in rx {
                self.process_sam_group_batch(&mut writer, &mut batch, fasta_map, annotation)?;
            }

            producer
                .join()
                .map_err(|_| anyhow::anyhow!("SAM group producer thread panicked"))?
        })?;

        writer.flush()?;
        Ok(stats)
    }

    /// Sequentially parses SAM records and sends owned read groups to workers.
    ///
    /// This split keeps all SAM semantics in one place while still allowing the
    /// heavier BSJ judgment to scale with Rayon.
    fn stream_sam_group_batches(
        &self,
        sam_file: &str,
        pb: &ProgressBar,
        tx: mpsc::SyncSender<Vec<SamOwnedGroup>>,
        batch_size: usize,
    ) -> Result<Scan1TraversalStats> {
        let sam_reader = BufReader::with_capacity(1024 * 1024, File::open(sam_file)?);
        let mut reader = sam::io::Reader::new(sam_reader);
        let header = reader.read_header()?;
        let mut record = sam::Record::default();
        let mut current_id: Vec<u8> = Vec::new();
        let mut group: [Vec<OwnedAlignmentRecord>; 2] =
            [Vec::with_capacity(8), Vec::with_capacity(8)];
        let mut stand_map: OwnedStandMap = HashMap::with_capacity(4);
        let mut align_num = 0usize;
        let mut max_read_len = 0i32;
        let mut cigar_buf = String::with_capacity(64);
        let mut seq_buf = String::with_capacity(256);
        let mut batch = Vec::with_capacity(batch_size);
        let mut records_since_progress = 0usize;
        let mut last_progress_pos = reader.get_mut().stream_position()?;
        let mut mapped_reads = 0u64;

        while reader.read_record(&mut record)? != 0 {
            let read_id = record
                .name()
                .ok_or_else(|| anyhow::anyhow!("Missing read name"))?;
            if read_id.to_vec() != current_id {
                if !current_id.is_empty() {
                    self.push_sam_owned_group(
                        &mut batch,
                        &current_id,
                        &mut group,
                        &mut stand_map,
                        align_num,
                    );
                    mapped_reads += 1;
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
            let start_pos = record
                .alignment_start()
                .transpose()?
                .map(|p| p.get() as i32)
                .unwrap_or(0);
            let mapq = record
                .mapping_quality()
                .transpose()?
                .map(u8::from)
                .unwrap_or(0) as i32;

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
                update_best_stand_seq(&mut stand_map, s_idx, st_c, Cow::Owned(seq.clone()));
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
            self.push_sam_owned_group(
                &mut batch,
                &current_id,
                &mut group,
                &mut stand_map,
                align_num,
            );
            mapped_reads += 1;
        }
        if !batch.is_empty() {
            tx.send(batch)
                .map_err(|_| anyhow::anyhow!("SAM group consumer dropped"))?;
        }

        let final_pos = reader.get_mut().stream_position()?;
        if final_pos > last_progress_pos {
            pb.inc(final_pos - last_progress_pos);
        }

        Ok(Scan1TraversalStats {
            max_read_len,
            mapped_reads,
        })
    }

    /// Moves the currently accumulated SAM read group into the worker batch.
    ///
    /// `mem::take` is used to avoid cloning alignment payloads on the hot path.
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

    /// Evaluates one batch of owned SAM read groups in parallel and writes BSJ1 lines.
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
                self.process_group_view(
                    &owned.read_id,
                    &owned.group,
                    &owned.stand_map,
                    fasta_map,
                    annotation,
                    &mut validator,
                    None,
                )
                .map(|(_, res_line)| res_line)
            })
            .collect();

        for res_line in results.into_iter().flatten() {
            writeln!(writer, "{}", res_line)?;
        }

        Ok(())
    }

    /// Concatenates per-shard BSJ1 outputs and returns the read IDs already claimed
    /// by Scan1.
    ///
    /// The final `HashSet` is later used by Scan2 to implement Java's "skip
    /// already-assigned reads" behavior.
    fn merge_and_collect_ids(&self, bsj_path: &str, num_threads: usize) -> Result<HashSet<String>> {
        let mut final_writer = BufWriter::with_capacity(1024 * 1024, File::create(bsj_path)?);
        let mut scan1_id_map = HashSet::new();
        for i in 0..num_threads {
            let shard_path = part_path(bsj_path, i);
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

    /// Runs Scan1 on BAM input using independent compressed-byte shards.
    ///
    /// Unlike SAM, BAM can be safely split near BGZF boundaries, so the parallel
    /// strategy here is shard-based instead of producer/consumer batching.
    pub fn run_bam(
        &mut self,
        bam_file: &str,
        bsj_path: &str,
        fasta_map: &HashMap<String, String>,
        annotation: &Annotation,
    ) -> Result<HashSet<String>> {
        use noodles::bam;
        let run_started = Instant::now();
        let profile = if Scan1Profile::enabled_from_env() {
            Some(Scan1Profile::default())
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
        let shard_stats: Vec<Scan1ShardStats> = (0..num_threads)
            .into_par_iter()
            .map(|i| {
                let start = i * shard_size;
                let end = if i == num_threads - 1 {
                    mmap.len()
                } else {
                    (i + 1) * shard_size
                };
                let shard_out = part_path(bsj_path, i);
                self.process_bam_shard_to_file(
                    i,
                    &mmap,
                    start,
                    end,
                    header_ref,
                    fasta_map,
                    annotation,
                    &pb,
                    &shard_out,
                    profile_ref,
                )
            })
            .collect::<Result<Vec<_>>>()?;

        self.read_len = shard_stats
            .iter()
            .map(|s| s.max_read_len)
            .max()
            .unwrap_or(0);
        self.mapped_reads = shard_stats.iter().map(|s| s.mapped_reads).sum();
        pb.set_position(file_size);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {msg}")?
                .progress_chars("#>-"),
        );
        pb.finish_with_message("Completed");
        let merge_started = Instant::now();
        let result = self.merge_and_collect_ids(bsj_path, num_threads);
        if let Some(profile) = profile_ref {
            profile
                .merge_ns
                .fetch_add(merge_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile
                .wall_total_ns
                .fetch_add(run_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile.report();
        }
        let scan1_ids = result?;
        self.bsj1_reads = scan1_ids.len();
        Ok(scan1_ids)
    }

    /// Processes one BAM shard into a temporary BSJ1 shard file.
    ///
    /// The shard starts at the first BGZF header found at or after `start`; the
    /// leading partial read group is skipped for non-zero shards so every read is
    /// judged exactly once.
    fn process_bam_shard_to_file(
        &self,
        shard_idx: usize,
        mmap: &Mmap,
        start: usize,
        end: usize,
        header: &sam::Header,
        fasta_map: &HashMap<String, String>,
        annotation: &Annotation,
        pb: &ProgressBar,
        out_path: &str,
        profile: Option<&Scan1Profile>,
    ) -> Result<Scan1ShardStats> {
        use noodles::bam;
        let shard_started = profile.map(|_| Instant::now());
        let mut writer = BufWriter::with_capacity(256 * 1024, File::create(out_path)?);
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
            // Decode from the previous BGZF block so we can decide whether the
            // first read group in `block_start` is a continuation from the prior
            // shard or a complete new group that should be kept.
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
            return Ok(Scan1ShardStats {
                max_read_len: 0,
                mapped_reads: 0,
            });
        }
        let mut reader = bam::io::Reader::new(&mmap[pos..]);
        if start == 0 {
            let _ = reader.read_header()?;
        }
        let mut record = bam::Record::default();
        let mut current_id: Vec<u8> = Vec::new();
        let mut group: [Vec<AlignmentRecord>; 2] = [Vec::with_capacity(8), Vec::with_capacity(8)];
        let mut stand_map: HashMap<i32, (char, Cow<str>)> = HashMap::with_capacity(4);
        let mut align_num = 0usize;
        let mut crossed_start = start == 0;
        let mut leading_partial_id: Option<Vec<u8>> = None;
        let mut last_progress_pos = block_start.saturating_sub(pos);
        let mut last_evicted_pos = pos;
        let eviction_threshold = 64 * 1024 * 1024;
        let mut shard_max_read_len = 0i32;
        let mut mapped_reads = 0u64;
        let mut cigar_buf = String::with_capacity(64);
        let mut seq_buf = String::with_capacity(256);
        while reader.read_record(&mut record)? != 0 {
            if let Some(profile) = profile {
                profile.records.fetch_add(1, Ordering::Relaxed);
            }
            let curr_c_pos = reader.get_ref().virtual_position().compressed() as usize;
            if curr_c_pos > last_progress_pos {
                pb.inc((curr_c_pos - last_progress_pos) as u64);
                last_progress_pos = curr_c_pos;
            }
            if curr_c_pos - (last_evicted_pos - pos) > eviction_threshold {
                advise_dontneed_aligned(
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
            let trace_read_id = String::from_utf8_lossy(read_id);
            if !crossed_start {
                if abs_c_pos < block_start {
                    trace_bam_shard_event(
                        shard_idx,
                        &trace_read_id,
                        "before_block_start",
                        start,
                        block_start,
                        end,
                        abs_c_pos,
                    );
                    leading_partial_id = Some(read_id.to_vec());
                    continue;
                }
                crossed_start = true;
            }
            if let Some(partial_id) = &leading_partial_id {
                if read_id.to_vec() == *partial_id {
                    // Non-zero shards start decoding from the previous BGZF block,
                    // so the first logical read group may already have records in
                    // the prior shard. Skip that entire group here and let the
                    // earlier shard own it completely.
                    trace_bam_shard_event(
                        shard_idx,
                        &trace_read_id,
                        "skip_leading_partial_group",
                        start,
                        block_start,
                        end,
                        abs_c_pos,
                    );
                    continue;
                }
                trace_bam_shard_event(
                    shard_idx,
                    &trace_read_id,
                    "leading_partial_cleared",
                    start,
                    block_start,
                    end,
                    abs_c_pos,
                );
                leading_partial_id = None;
            }
            if read_id.to_vec() != current_id {
                trace_bam_shard_event(
                    shard_idx,
                    &trace_read_id,
                    "group_boundary",
                    start,
                    block_start,
                    end,
                    abs_c_pos,
                );
                if !current_id.is_empty() {
                    mapped_reads += 1;
                    let id_str = String::from_utf8_lossy(&current_id);
                    let non_empty_groups = group.iter().filter(|g| !g.is_empty()).count();
                    if align_num > 2 || non_empty_groups == 1 {
                        if let Some((_, res_line)) = self.process_group_view(
                            &id_str,
                            &group,
                            &stand_map,
                            fasta_map,
                            annotation,
                            &mut validator,
                            profile,
                        ) {
                            trace_bam_shard_event(
                                shard_idx,
                                &id_str,
                                "write_bsj1",
                                start,
                                block_start,
                                end,
                                abs_c_pos,
                            );
                            let write_started = profile.map(|_| Instant::now());
                            writeln!(writer, "{}", res_line)?;
                            if let (Some(profile), Some(write_started)) = (profile, write_started) {
                                profile.write_ns.fetch_add(
                                    write_started.elapsed().as_nanos() as u64,
                                    Ordering::Relaxed,
                                );
                            }
                        }
                    }
                    if pos + curr_c_pos > end {
                        current_id.clear();
                        break;
                    }
                }
                current_id = read_id.to_vec();
                group[0].clear();
                group[1].clear();
                stand_map.clear();
                align_num = 0;
            }
            let chrom = match record.reference_sequence(header) {
                Some(Ok((name, _))) => String::from_utf8_lossy(name).to_string(),
                _ => "*".to_string(),
            };
            let flag = i32::from(u16::from(record.flags()));
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
            shard_max_read_len = shard_max_read_len.max(seq_buf.len() as i32);
            let seq = seq_buf.clone();
            if !seq.is_empty() && seq != "*" {
                let s_idx = if flag & 0x40 != 0 { 1 } else { 0 };
                let st_c = if flag & 0x10 != 0 { '1' } else { '0' };
                update_best_stand_seq(&mut stand_map, s_idx, st_c, Cow::Owned(seq.clone()));
            }
            group[if flag & 0x40 != 0 { 1 } else { 0 }].push(AlignmentRecord {
                flag,
                chrom: Cow::Owned(chrom),
                pos: start_pos,
                mapq,
                cigar: Cow::Owned(cigar_buf.clone()),
                seq: Cow::Owned(seq),
            });
            align_num += 1;
        }
        if !current_id.is_empty() {
            mapped_reads += 1;
            let id_str = String::from_utf8_lossy(&current_id);
            if let Some((_, res_line)) = self.process_group_view(
                &id_str,
                &group,
                &stand_map,
                fasta_map,
                annotation,
                &mut validator,
                profile,
            ) {
                trace_bam_shard_event(
                    shard_idx,
                    &id_str,
                    "write_bsj1_tail",
                    start,
                    block_start,
                    end,
                    pos + last_progress_pos,
                );
                let write_started = profile.map(|_| Instant::now());
                writeln!(writer, "{}", res_line)?;
                if let (Some(profile), Some(write_started)) = (profile, write_started) {
                    profile
                        .write_ns
                        .fetch_add(write_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
                }
            }
        }
        writer.flush()?;
        if let (Some(profile), Some(shard_started)) = (profile, shard_started) {
            profile
                .shard_total_ns
                .fetch_add(shard_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
        Ok(Scan1ShardStats {
            max_read_len: shard_max_read_len,
            mapped_reads,
        })
    }

    /// Thin profiling wrapper around the actual read-group judge.
    fn process_group_view<'a>(
        &self,
        read_id: &str,
        group: &[Vec<AlignmentRecord<'a>>; 2],
        stand_map: &HashMap<i32, (char, Cow<'a, str>)>,
        fasta_map: &HashMap<String, String>,
        annotation: &Annotation,
        validator: &mut IsBSJHg2,
        profile: Option<&Scan1Profile>,
    ) -> Option<(String, String)> {
        let started = profile.map(|_| Instant::now());
        let result = self.process_group_view_impl(
            read_id, group, stand_map, fasta_map, annotation, validator, profile,
        );
        if let (Some(profile), Some(started)) = (profile, started) {
            profile
                .group_process_ns
                .fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            profile.groups.fetch_add(1, Ordering::Relaxed);
        }
        result
    }

    /// Reproduces Java Scan1 BSJ candidate assembly for one read group.
    ///
    /// Most Scan1 CPU time ends up here and inside `IsBSJHg2::is_bsj_hg1`, so the
    /// low-risk optimizations in this function are limited to reusing derived
    /// values such as `misd` outputs and oriented read sequence strings. The
    /// branch structure itself intentionally stays close to Java for parity.
    fn process_group_view_impl<'a>(
        &self,
        read_id: &str,
        group: &[Vec<AlignmentRecord<'a>>; 2],
        stand_map: &HashMap<i32, (char, Cow<'a, str>)>,
        fasta_map: &HashMap<String, String>,
        annotation: &Annotation,
        validator: &mut IsBSJHg2,
        profile: Option<&Scan1Profile>,
    ) -> Option<(String, String)> {
        // Debug trace is intentionally read-scoped to avoid overwhelming output.
        let trace_read = should_trace_read(read_id);
        let [pair1, pair2] = group;
        let group_count =
            (if !pair1.is_empty() { 1 } else { 0 }) + (if !pair2.is_empty() { 1 } else { 0 });
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
            if segments.len() < 2 {
                continue;
            }
            let seq_len = read_seq.len() as i32;
            // `misd` is deterministic for a given alignment and sequence length, so
            // computing it once per segment is a safe hot-path optimization.
            let segment_misds: Vec<[i32; 4]> = segments
                .iter()
                .map(|aln| misd(&aln.cigar, seq_len))
                .collect();
            let mate_misds: Vec<[i32; 4]> = mate_segments
                .iter()
                .map(|aln| misd(&aln.cigar, seq_len))
                .collect();
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
                    let (mut s1_n, mut s2_n, mut adj1, mut adj2): (i32, i32, i32, i32) =
                        (0, 0, 0, 0);
                    let (mut str1, mut str2, mut str3, mut str4) =
                        (String::new(), String::new(), String::new(), String::new());
                    let (mut q1, mut q2, mut sum_q) = (0, 0, 0);
                    let al1_strand = if al1.flag & 0x10 != 0 { '1' } else { '0' };
                    // Sequence orientation depends on the current alignment strand,
                    // not just on the original read slot. Hoisting it here avoids
                    // rebuilding the same reverse-complement per candidate branch.
                    let seq_oriented: Cow<'_, str> = if al1_strand == read_strand {
                        Cow::Borrowed(read_seq)
                    } else {
                        Cow::Owned(crate::utils::reverse_complement(read_seq))
                    };
                    if c1[0] * c2[0] == -1 {
                        let scale = c1[0] * al1.pos + c1[2] + c2[0] * al2.pos + c2[2];
                        if scale > 0
                            && (c1[1] - c2[1]).abs() <= 6
                            && scale <= self.max_circle
                            && scale >= self.min_circle
                        {
                            adj1 = (c1[1] * c1[0] + c2[1] * c2[0]) / 2;
                            adj2 = (c1[1] * c1[0] + c2[1] * c2[0]) - adj1;
                            if adj1.abs() <= 4 {
                                identified = true;
                                s1_n = al1.pos + adj1;
                                s2_n = al2.pos + c2[3] - 1 - adj2;
                                str2 = java_substring(seq_oriented.as_ref(), 0, c1[1] + adj1)
                                    .to_string();
                                str1 = java_substring(seq_oriented.as_ref(), c1[1] + adj1, seq_len)
                                    .to_string();
                                str3 = "*".to_string();
                                if al1.mapq >= self.min_mapq_uni && al2.mapq >= self.min_mapq_uni {
                                    q1 = 1;
                                    q2 = 1;
                                    sum_q = 1;
                                } else if al1.mapq >= self.min_mapq_uni {
                                    q1 = 1;
                                } else if al2.mapq >= self.min_mapq_uni {
                                    q2 = 1;
                                }
                            }
                        }
                    } else if (c1[0] * c2[0]).abs() == 10 {
                        if c1[0] == -1 {
                            let scale = al2.pos + c2[3] - 1 - al1.pos;
                            if scale > 0
                                && (seq_len - c2[2] - c1[1]).abs() <= 6
                                && scale <= self.max_circle
                                && scale >= self.min_circle
                            {
                                adj1 = (c2[1] + c2[3] - c1[1]) / 2;
                                adj2 = (c2[1] + c2[3] - c1[1]) - adj1;
                                if adj1.abs() <= 4 {
                                    identified = true;
                                    s1_n = al1.pos + adj1;
                                    s2_n = al2.pos + c2[3] - 1 - adj2;
                                    str1 = java_substring(
                                        seq_oriented.as_ref(),
                                        c1[1] + adj1,
                                        seq_len,
                                    )
                                    .to_string();
                                    str2 =
                                        java_substring(seq_oriented.as_ref(), c2[1], c1[1] + adj1)
                                            .to_string();
                                    str3 =
                                        java_substring(seq_oriented.as_ref(), 0, c2[1]).to_string();
                                    if al1.mapq >= self.min_mapq_uni
                                        && al2.mapq >= self.min_mapq_uni
                                    {
                                        q1 = 1;
                                        q2 = 1;
                                        sum_q = 1;
                                    } else if al1.mapq >= self.min_mapq_uni {
                                        q1 = 1;
                                    } else if al2.mapq >= self.min_mapq_uni {
                                        q2 = 1;
                                    }
                                }
                            }
                        } else {
                            let scale = al1.pos + c1[3] - 1 - al2.pos;
                            if scale > 0
                                && (c1[1] - c2[1]).abs() <= 6
                                && scale <= self.max_circle
                                && scale >= self.min_circle
                            {
                                adj1 = (c1[1] - c2[1]) / 2;
                                adj2 = (c1[1] - c2[1]) - adj1;
                                if adj1.abs() <= 4 {
                                    identified = true;
                                    s1_n = al2.pos + adj1;
                                    s2_n = al1.pos + c1[3] - 1 - adj2;
                                    str2 = java_substring(seq_oriented.as_ref(), 0, c1[1] - adj2)
                                        .to_string();
                                    str1 = java_substring(
                                        seq_oriented.as_ref(),
                                        c1[1] - adj2,
                                        seq_len - c2[2],
                                    )
                                    .to_string();
                                    str3 = java_substring(
                                        seq_oriented.as_ref(),
                                        seq_len - c2[2],
                                        seq_len,
                                    )
                                    .to_string();
                                    if al1.mapq >= self.min_mapq_uni
                                        && al2.mapq >= self.min_mapq_uni
                                    {
                                        q1 = 1;
                                        q2 = 1;
                                        sum_q = 1;
                                    } else if al2.mapq >= self.min_mapq_uni {
                                        q1 = 1;
                                    } else if al1.mapq >= self.min_mapq_uni {
                                        q2 = 1;
                                    }
                                }
                            }
                        }
                    }
                    if identified {
                        let mut s4_ok = 0;
                        if is_paried == 1 {
                            if let Some((mate_stand, mate_seq)) = stand_map.get(&((1 - n) as i32)) {
                                str4 = if al1_strand != *mate_stand {
                                    mate_seq.to_string()
                                } else {
                                    crate::utils::reverse_complement(mate_seq)
                                };
                            }
                            for (m_aln, mc) in mate_segments.iter().zip(mate_misds.iter()) {
                                if m_aln.chrom == al1.chrom && m_aln.mapq >= self.min_mapq_uni {
                                    if m_aln.flag & 0x10 != al1.flag & 0x10
                                        && m_aln.pos >= s1_n - 6
                                        && m_aln.pos + mc[3] - 1 <= s2_n + 6
                                    {
                                        s4_ok = 1;
                                        break;
                                    } else {
                                        s4_ok = -1;
                                        break;
                                    }
                                }
                            }
                        } else {
                            s4_ok = 1;
                            str4 = "".to_string();
                        }
                        let mut line_arr = [
                            if al1.flag & 0x10 != 0 { "1" } else { "0" }.to_string(),
                            al1.chrom.to_string(),
                            str1,
                            str2,
                            str3,
                            str4,
                            q1.to_string(),
                            q2.to_string(),
                            s4_ok.to_string(),
                            s1_n.to_string(),
                            s2_n.to_string(),
                            adj1.to_string(),
                            adj2.to_string(),
                        ];
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
                            let res = validator.is_bsj_hg1(
                                &mut line_arr,
                                chr_seq,
                                sum_q,
                                "chrM",
                                false,
                                &annotation.chr_exon_start_map,
                                &annotation.chr_exon_end_map,
                            );
                            if let Some(profile) = profile {
                                profile.hg1_calls.fetch_add(1, Ordering::Relaxed);
                                if let Some(hg1_started) = hg1_started {
                                    profile.bsj_judge_ns.fetch_add(
                                        hg1_started.elapsed().as_nanos() as u64,
                                        Ordering::Relaxed,
                                    );
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
                                let cigar1 = Self::normalize_bsj_cigar(al1.cigar.as_ref());
                                let cigar2 = Self::normalize_bsj_cigar(al2.cigar.as_ref());
                                return Some((
                                    read_id.to_string(),
                                    format!("{}\t{};{}\t{}", read_id, cigar1, cigar2, res),
                                ));
                            }
                        }
                    }
                }
            }
        }
        None
    }
}
