//! Validator module: Reproduces the core validation logic of Java CIRI3.
//!
//! This module includes comprehensive sequence validation, canonical splice signal
//! identification, and linear competition checks to distinguish BSJs from linear splicing noise.

use crate::index_compare::IndexCompare;
use crate::runtime::{
    emit_debug_line, emit_perf_line, scan1_profile_enabled, scan2_profile_enabled,
    trace_hg2_enabled,
};
use std::collections::HashMap;
use std::fmt::Write as FmtWrite;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

/// Aggregated profiling counters for the Scan1 validator hot path.
///
/// This stays private because the profiler is a maintenance aid, not part of the
/// validation API surface.
#[derive(Default)]
struct HgProfile {
    hg1_total_ns: AtomicU64,
    hg1_calls: AtomicU64,
    index_compare_ns: AtomicU64,
    index_compare_calls: AtomicU64,
    exon_lookup_ns: AtomicU64,
    exon_lookup_calls: AtomicU64,
    shift_prepare_ns: AtomicU64,
    shift_prepare_calls: AtomicU64,
    shift_entries: AtomicU64,
    linear_check_ns: AtomicU64,
    linear_check_calls: AtomicU64,
    circ2_ns: AtomicU64,
    circ2_calls: AtomicU64,
    circ3_ns: AtomicU64,
    circ3_calls: AtomicU64,
    hg1_hits: AtomicU64,
    sw_ns: AtomicU64,
    sw_calls: AtomicU64,
}

static HG_PROFILE: HgProfile = HgProfile {
    hg1_total_ns: AtomicU64::new(0),
    hg1_calls: AtomicU64::new(0),
    index_compare_ns: AtomicU64::new(0),
    index_compare_calls: AtomicU64::new(0),
    exon_lookup_ns: AtomicU64::new(0),
    exon_lookup_calls: AtomicU64::new(0),
    shift_prepare_ns: AtomicU64::new(0),
    shift_prepare_calls: AtomicU64::new(0),
    shift_entries: AtomicU64::new(0),
    linear_check_ns: AtomicU64::new(0),
    linear_check_calls: AtomicU64::new(0),
    circ2_ns: AtomicU64::new(0),
    circ2_calls: AtomicU64::new(0),
    circ3_ns: AtomicU64::new(0),
    circ3_calls: AtomicU64::new(0),
    hg1_hits: AtomicU64::new(0),
    sw_ns: AtomicU64::new(0),
    sw_calls: AtomicU64::new(0),
};

/// Aggregated profiling counters for the Scan2 validator hot path.
#[derive(Default)]
struct Hg2Profile {
    hg2_total_ns: AtomicU64,
    hg2_calls: AtomicU64,
    sw_ns: AtomicU64,
    sw_calls: AtomicU64,
    linear11_ns: AtomicU64,
    linear11_calls: AtomicU64,
    linear12_ns: AtomicU64,
    linear12_calls: AtomicU64,
    circ2_ns: AtomicU64,
    circ2_calls: AtomicU64,
    circ3_ns: AtomicU64,
    circ3_calls: AtomicU64,
}

static HG2_PROFILE: Hg2Profile = Hg2Profile {
    hg2_total_ns: AtomicU64::new(0),
    hg2_calls: AtomicU64::new(0),
    sw_ns: AtomicU64::new(0),
    sw_calls: AtomicU64::new(0),
    linear11_ns: AtomicU64::new(0),
    linear11_calls: AtomicU64::new(0),
    linear12_ns: AtomicU64::new(0),
    linear12_calls: AtomicU64::new(0),
    circ2_ns: AtomicU64::new(0),
    circ2_calls: AtomicU64::new(0),
    circ3_ns: AtomicU64::new(0),
    circ3_calls: AtomicU64::new(0),
};

/// Checks whether verbose `is_bsj_hg2` stage tracing is enabled.
///
/// This trace is intentionally separate from `CIRI_TRACE_ALL_CANDS`: the latter
/// changes candidate traversal so developers can inspect the full search space,
/// while `--debug` / `CIRI_TRACE_HG2` keep the original traversal and only
/// annotate which validator branch accepted or rejected a traced candidate.
fn trace_scan2_hg2_enabled() -> bool {
    trace_hg2_enabled()
}

/// Emits one targeted `is_bsj_hg2` branch trace line.
///
/// The payload mirrors the Java-shaped `circ_line_arr` contract so a developer
/// can compare Rust and Java branch decisions without first reconstructing the
/// candidate fields by hand. This stays off by default because the validator hot
/// path calls it frequently once tracing is enabled.
fn trace_scan2_hg2(stage: &str, circ_line_arr: &[String], extra: &str) {
    if !trace_scan2_hg2_enabled() {
        return;
    }
    emit_debug_line(&format!(
        "[TRACE_SCAN2_HG2] stage={} type={} strand={} chr={} site1={} site2={} mapq={} s2_ok={} str_len={} pair_len={} str3_len={} {}",
        stage,
        circ_line_arr[2],
        circ_line_arr[0],
        circ_line_arr[1],
        circ_line_arr[3],
        circ_line_arr[4],
        circ_line_arr[12],
        circ_line_arr[8],
        circ_line_arr[5].len(),
        circ_line_arr[6].len(),
        circ_line_arr[7].len(),
        extra,
    ));
}

/// Adds elapsed nanoseconds to one profiling counter.
#[inline]
fn add_ns(counter: &AtomicU64, started: Instant) {
    counter.fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
}

/// Finalizes one profiled Scan2 validator call and returns its result.
#[inline]
fn finish_scan2_hg2(result: String, started: Option<Instant>) -> String {
    if let Some(started) = started {
        add_ns(&HG2_PROFILE.hg2_total_ns, started);
    }
    result
}

/// Prints the collected Scan1 validator profile.
///
/// The report is intentionally stable across optimization rounds so release A/B
/// measurements stay easy to compare.
pub fn report_scan1_hg_profile() {
    if !scan1_profile_enabled() {
        return;
    }
    let hg1_total_ns = HG_PROFILE.hg1_total_ns.load(Ordering::Relaxed);
    let hg1_calls = HG_PROFILE.hg1_calls.load(Ordering::Relaxed);
    let index_compare_ns = HG_PROFILE.index_compare_ns.load(Ordering::Relaxed);
    let index_compare_calls = HG_PROFILE.index_compare_calls.load(Ordering::Relaxed);
    let exon_lookup_ns = HG_PROFILE.exon_lookup_ns.load(Ordering::Relaxed);
    let exon_lookup_calls = HG_PROFILE.exon_lookup_calls.load(Ordering::Relaxed);
    let shift_prepare_ns = HG_PROFILE.shift_prepare_ns.load(Ordering::Relaxed);
    let shift_prepare_calls = HG_PROFILE.shift_prepare_calls.load(Ordering::Relaxed);
    let shift_entries = HG_PROFILE.shift_entries.load(Ordering::Relaxed);
    let linear_check_ns = HG_PROFILE.linear_check_ns.load(Ordering::Relaxed);
    let linear_check_calls = HG_PROFILE.linear_check_calls.load(Ordering::Relaxed);
    let circ2_ns = HG_PROFILE.circ2_ns.load(Ordering::Relaxed);
    let circ2_calls = HG_PROFILE.circ2_calls.load(Ordering::Relaxed);
    let circ3_ns = HG_PROFILE.circ3_ns.load(Ordering::Relaxed);
    let circ3_calls = HG_PROFILE.circ3_calls.load(Ordering::Relaxed);
    let hg1_hits = HG_PROFILE.hg1_hits.load(Ordering::Relaxed);
    let sw_ns = HG_PROFILE.sw_ns.load(Ordering::Relaxed);
    let sw_calls = HG_PROFILE.sw_calls.load(Ordering::Relaxed);
    let known_ns = index_compare_ns
        .saturating_add(exon_lookup_ns)
        .saturating_add(shift_prepare_ns)
        .saturating_add(linear_check_ns)
        .saturating_add(circ2_ns)
        .saturating_add(circ3_ns);
    let other_ns = hg1_total_ns.saturating_sub(known_ns);
    let pct = |part: u64, whole: u64| -> f64 {
        if whole == 0 {
            0.0
        } else {
            part as f64 * 100.0 / whole as f64
        }
    };
    emit_perf_line(&format!(
        "[PROFILE_SCAN1_HG1] calls={} hits={} shift_entries={} sw_calls={}",
        hg1_calls, hg1_hits, shift_entries, sw_calls
    ));
    emit_perf_line(&format!(
        "[PROFILE_SCAN1_HG1] total_ms={:.3} index_compare_ms={:.3} ({:.1}%, calls={}) exon_lookup_ms={:.3} ({:.1}%, calls={}) shift_prepare_ms={:.3} ({:.1}%, calls={}) linear_check_ms={:.3} ({:.1}%, calls={}) circ2_ms={:.3} ({:.1}%, calls={}) circ3_ms={:.3} ({:.1}%, calls={}) other_ms={:.3} ({:.1}%) sw_ms={:.3}",
        hg1_total_ns as f64 / 1_000_000.0,
        index_compare_ns as f64 / 1_000_000.0,
        pct(index_compare_ns, hg1_total_ns),
        index_compare_calls,
        exon_lookup_ns as f64 / 1_000_000.0,
        pct(exon_lookup_ns, hg1_total_ns),
        exon_lookup_calls,
        shift_prepare_ns as f64 / 1_000_000.0,
        pct(shift_prepare_ns, hg1_total_ns),
        shift_prepare_calls,
        linear_check_ns as f64 / 1_000_000.0,
        pct(linear_check_ns, hg1_total_ns),
        linear_check_calls,
        circ2_ns as f64 / 1_000_000.0,
        pct(circ2_ns, hg1_total_ns),
        circ2_calls,
        circ3_ns as f64 / 1_000_000.0,
        pct(circ3_ns, hg1_total_ns),
        circ3_calls,
        other_ns as f64 / 1_000_000.0,
        pct(other_ns, hg1_total_ns),
        sw_ns as f64 / 1_000_000.0,
    ));
}

/// Prints the collected Scan2 validator profile.
pub fn report_scan2_hg_profile() {
    if !scan2_profile_enabled() {
        return;
    }
    let hg2_total_ns = HG2_PROFILE.hg2_total_ns.load(Ordering::Relaxed);
    let hg2_calls = HG2_PROFILE.hg2_calls.load(Ordering::Relaxed);
    let sw_ns = HG2_PROFILE.sw_ns.load(Ordering::Relaxed);
    let sw_calls = HG2_PROFILE.sw_calls.load(Ordering::Relaxed);
    let linear11_ns = HG2_PROFILE.linear11_ns.load(Ordering::Relaxed);
    let linear11_calls = HG2_PROFILE.linear11_calls.load(Ordering::Relaxed);
    let linear12_ns = HG2_PROFILE.linear12_ns.load(Ordering::Relaxed);
    let linear12_calls = HG2_PROFILE.linear12_calls.load(Ordering::Relaxed);
    let circ2_ns = HG2_PROFILE.circ2_ns.load(Ordering::Relaxed);
    let circ2_calls = HG2_PROFILE.circ2_calls.load(Ordering::Relaxed);
    let circ3_ns = HG2_PROFILE.circ3_ns.load(Ordering::Relaxed);
    let circ3_calls = HG2_PROFILE.circ3_calls.load(Ordering::Relaxed);
    let known_ns = sw_ns
        .saturating_add(linear11_ns)
        .saturating_add(linear12_ns)
        .saturating_add(circ2_ns)
        .saturating_add(circ3_ns);
    let other_ns = hg2_total_ns.saturating_sub(known_ns);
    let pct = |part: u64, whole: u64| -> f64 {
        if whole == 0 {
            0.0
        } else {
            part as f64 * 100.0 / whole as f64
        }
    };
    emit_perf_line(&format!(
        "[PROFILE_SCAN2_HG2] calls={} sw_calls={} linear11_calls={} linear12_calls={} circ2_calls={} circ3_calls={}",
        hg2_calls, sw_calls, linear11_calls, linear12_calls, circ2_calls, circ3_calls
    ));
    emit_perf_line(&format!(
        "[PROFILE_SCAN2_HG2] total_ms={:.3} sw_ms={:.3} ({:.1}%) linear11_ms={:.3} ({:.1}%) linear12_ms={:.3} ({:.1}%) circ2_ms={:.3} ({:.1}%) circ3_ms={:.3} ({:.1}%) other_ms={:.3} ({:.1}%)",
        hg2_total_ns as f64 / 1_000_000.0,
        sw_ns as f64 / 1_000_000.0,
        pct(sw_ns, hg2_total_ns),
        linear11_ns as f64 / 1_000_000.0,
        pct(linear11_ns, hg2_total_ns),
        linear12_ns as f64 / 1_000_000.0,
        pct(linear12_ns, hg2_total_ns),
        circ2_ns as f64 / 1_000_000.0,
        pct(circ2_ns, hg2_total_ns),
        circ3_ns as f64 / 1_000_000.0,
        pct(circ3_ns, hg2_total_ns),
        other_ns as f64 / 1_000_000.0,
        pct(other_ns, hg2_total_ns),
    ));
}

/// Smith-Waterman local alignment implementation with traceback.
///
/// This is kept close to Java's scoring semantics even though more specialized
/// implementations could be faster, because Scan2 parity depends on matching the
/// original acceptance thresholds.
pub struct SmithWaterman {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_penalty: i32,
    pub score: i32,
    /// Number of steps in the traceback, equivalent to `alignment[1].length()` in Java.
    pub aligned_len: i32,
    seq1: String,
    seq2: String,
}

impl SmithWaterman {
    /// Initializes an aligner with Java-compatible scoring parameters.
    pub fn new(m: i32, mis: i32, gap: i32) -> Self {
        Self {
            match_score: m,
            mismatch_penalty: mis,
            gap_penalty: gap,
            score: 0,
            aligned_len: 0,
            seq1: String::new(),
            seq2: String::new(),
        }
    }

    /// Sets the sequences to be aligned.
    ///
    /// Uppercasing here mirrors Java behavior and keeps later scoring branches
    /// case-insensitive without adding checks inside the DP loop.
    pub fn set_seq(&mut self, s1: &str, s2: &str) {
        self.seq1 = s1.to_uppercase();
        self.seq2 = s2.to_uppercase();
    }

    /// Performs local alignment and traceback to calculate the optimal score and
    /// aligned length.
    ///
    /// The implementation favors parity and traceability over algorithmic novelty:
    /// it uses an explicit DP table and recomputes the final score from traceback,
    /// matching Java's `getAlignmentScore` behavior.
    pub fn align(&mut self) {
        let sw_started = scan1_profile_enabled().then(Instant::now);
        // Java parity (`smith` package): table rows = seq2, cols = seq1.
        let cols = self.seq1.len();
        let rows = self.seq2.len();
        if cols == 0 || rows == 0 {
            self.score = 0;
            self.aligned_len = 0;
            if let Some(sw_started) = sw_started {
                HG_PROFILE.sw_calls.fetch_add(1, Ordering::Relaxed);
                add_ns(&HG_PROFILE.sw_ns, sw_started);
            }
            return;
        }

        let s1 = self.seq1.as_bytes();
        let s2 = self.seq2.as_bytes();
        let mut score_table = vec![vec![0i32; cols + 1]; rows + 1];
        let mut prev: Vec<Vec<Option<(usize, usize)>>> = vec![vec![None; cols + 1]; rows + 1];
        let mut high_row = 0usize;
        let mut high_col = 0usize;

        for row in 1..=rows {
            for col in 1..=cols {
                let row_space_score = score_table[row - 1][col] + self.gap_penalty;
                let col_space_score = score_table[row][col - 1] + self.gap_penalty;
                let mut match_or_mismatch_score = score_table[row - 1][col - 1];
                if s2[row - 1] == s1[col - 1] {
                    match_or_mismatch_score += self.match_score;
                } else {
                    match_or_mismatch_score += self.mismatch_penalty;
                }

                let mut cell_score = 0i32;
                let mut cell_prev: Option<(usize, usize)> = None;
                if row_space_score >= col_space_score {
                    if match_or_mismatch_score >= row_space_score {
                        if match_or_mismatch_score > 0 {
                            cell_score = match_or_mismatch_score;
                            cell_prev = Some((row - 1, col - 1));
                        }
                    } else if row_space_score > 0 {
                        cell_score = row_space_score;
                        cell_prev = Some((row - 1, col));
                    }
                } else if match_or_mismatch_score >= col_space_score {
                    if match_or_mismatch_score > 0 {
                        cell_score = match_or_mismatch_score;
                        cell_prev = Some((row - 1, col - 1));
                    }
                } else if col_space_score > 0 {
                    cell_score = col_space_score;
                    cell_prev = Some((row, col - 1));
                }

                score_table[row][col] = cell_score;
                prev[row][col] = cell_prev;
                if cell_score > score_table[high_row][high_col] {
                    high_row = row;
                    high_col = col;
                }
            }
        }

        // Java traceback stops at score == 0.
        let mut row = high_row;
        let mut col = high_col;
        let mut align1: Vec<u8> = Vec::new();
        let mut align2: Vec<u8> = Vec::new();
        while score_table[row][col] != 0 {
            let (pr, pc) = match prev[row][col] {
                Some(p) => p,
                None => break,
            };
            if row - pr == 1 {
                align2.push(s2[row - 1]);
            } else {
                align2.push(b'-');
            }
            if col - pc == 1 {
                align1.push(s1[col - 1]);
            } else {
                align1.push(b'-');
            }
            row = pr;
            col = pc;
        }
        self.aligned_len = align2.len() as i32;

        // Java getAlignmentScore recomputes score from traceback alignments.
        let mut total = 0i32;
        for i in 0..align1.len() {
            let c1 = align1[i];
            let c2 = align2[i];
            if c1 == b'-' || c2 == b'-' {
                total += self.gap_penalty;
            } else if c1 == c2 {
                total += self.match_score;
            } else {
                total += self.mismatch_penalty;
            }
        }
        self.score = total;
        if let Some(sw_started) = sw_started {
            HG_PROFILE.sw_calls.fetch_add(1, Ordering::Relaxed);
            add_ns(&HG_PROFILE.sw_ns, sw_started);
        }
    }
}

/// Shared validator for Scan1 and Scan2 BSJ candidate verification.
///
/// The scratch buffers are a retained low-risk optimization from Scan1 tuning:
/// `is_in_circ_rna_3` is called often enough that reusing window-search buffers
/// measurably reduces allocation pressure without changing decision logic.
pub struct IsBSJHg2 {
    pub linear_range_size_min: i32,
    pub min_mapq_uni: i32,
    pub initial_size1: i32,
    pub aligner: SmithWaterman,
    window_unit: [i32; 5],
    query_code_buf: Vec<u64>,
    circ_pos_buf: Vec<i32>,
    pem_pos_buf: Vec<i32>,
}

/// Mirrors Java `String.substring` with explicit bound checks.
///
/// Keeping this helper centralizes the end-exclusive indexing contract used
/// throughout the port and makes parity bugs easier to audit.
pub fn java_substring(s: &str, start: i32, end: i32) -> &str {
    let len = s.len() as i32;
    assert!(
        start >= 0,
        "java_substring start<0: start={}, end={}, len={}",
        start,
        end,
        len
    );
    assert!(
        end >= 0,
        "java_substring end<0: start={}, end={}, len={}",
        start,
        end,
        len
    );
    assert!(
        start <= end,
        "java_substring start>end: start={}, end={}, len={}",
        start,
        end,
        len
    );
    assert!(
        end <= len,
        "java_substring end>len: start={}, end={}, len={}",
        start,
        end,
        len
    );
    &s[start as usize..end as usize]
}

impl IsBSJHg2 {
    /// Creates a new validator with Java-compatible defaults.
    pub fn new(linear_range_size_min: i32, min_mapq_uni: i32) -> Self {
        Self {
            linear_range_size_min,
            min_mapq_uni,
            initial_size1: 7,
            aligner: SmithWaterman::new(1, -1, -1),
            window_unit: [9, 7, 5, 4, 3],
            query_code_buf: Vec::with_capacity(64),
            circ_pos_buf: Vec::with_capacity(64),
            pem_pos_buf: Vec::with_capacity(64),
        }
    }

    /// Shared distance heuristic used by the linear-competition window scans.
    fn distance_loci_stats(
        &self,
        locus_count: usize,
        locus_sum: i32,
        locus2_count: usize,
        locus2_sum: i32,
        window_step: i32,
    ) -> i32 {
        if locus_count < 2 && locus2_count < 2 {
            return 0;
        }
        if locus_sum <= window_step * locus_count as i32 && locus_sum * 20 < locus2_sum {
            1
        } else {
            0
        }
    }

    /// Returns the left-side linear competitor window used by Scan1 checks.
    ///
    /// This helper exists mostly to document and reuse Java's boundary arithmetic;
    /// the tiny speedup from avoiding duplicate branch trees is only a secondary
    /// benefit.
    #[inline]
    fn linear_range_left<'a>(&self, chr_taga: &'a str, site1_new: i32, site2_new: i32) -> &'a str {
        if site2_new - site1_new + 5 >= self.linear_range_size_min {
            if 2 * site1_new >= site2_new + 6 {
                java_substring(chr_taga, 2 * site1_new - site2_new - 6, site1_new - 1)
            } else {
                java_substring(chr_taga, 0, site1_new - 1)
            }
        } else if site1_new >= self.linear_range_size_min + 1 {
            java_substring(
                chr_taga,
                site1_new - self.linear_range_size_min - 1,
                site1_new - 1,
            )
        } else {
            java_substring(chr_taga, 0, site1_new - 1)
        }
    }

    /// Returns the right-side linear competitor window used by Scan1 checks.
    #[inline]
    fn linear_range_right<'a>(
        &self,
        chr_taga: &'a str,
        site1_new: i32,
        site2_new: i32,
        chr_taga_len: i32,
    ) -> &'a str {
        if site2_new - site1_new + 5 >= self.linear_range_size_min {
            if 2 * site2_new - site1_new + 5 > chr_taga_len {
                java_substring(chr_taga, site2_new, chr_taga_len)
            } else {
                java_substring(chr_taga, site2_new, 2 * site2_new - site1_new + 5)
            }
        } else if site2_new + self.linear_range_size_min > chr_taga_len {
            java_substring(chr_taga, site2_new, chr_taga_len)
        } else {
            java_substring(chr_taga, site2_new, site2_new + self.linear_range_size_min)
        }
    }

    /// Linear competition check for the left-anchored junction fragment
    /// (`Phase 1.1` in the original implementation).
    pub fn is_in_circ_rna_1_1(
        &mut self,
        len_str: i32,
        str_val: &str,
        circ_range_seq: &str,
        linear_range: &str,
    ) -> i32 {
        for &step in &self.window_unit {
            if len_str < step * 2 {
                continue;
            }
            let window_size = step * 2;
            let mut trial = (len_str - window_size) / step;
            let mut window_count = trial.max(0) as usize + 1;
            if len_str % step != 0 {
                window_count += 1;
            }
            self.query_code_buf.clear();
            self.query_code_buf
                .reserve(window_count.saturating_sub(self.query_code_buf.capacity()));
            for j in 0..=trial {
                let s_idx = len_str - j * step - window_size;
                let e_idx = len_str - j * step;
                self.query_code_buf.push(Self::encode_window(
                    &str_val.as_bytes()[s_idx as usize..e_idx as usize],
                ));
            }
            if len_str % step != 0 {
                self.query_code_buf.push(Self::encode_window(
                    &str_val.as_bytes()[0..window_size as usize],
                ));
            }
            let mut hash_keys = [u64::MAX; 128];
            let mut hash_masks = [0u64; 128];
            let use_hash_fast_path =
                Self::build_query_hash(&self.query_code_buf, &mut hash_keys, &mut hash_masks);
            self.circ_pos_buf.clear();
            self.circ_pos_buf.resize(window_count, -1);
            Self::scan_window_positions(
                &self.query_code_buf,
                circ_range_seq.as_bytes(),
                &mut self.circ_pos_buf,
                window_size as usize,
                true,
                use_hash_fast_path,
                &hash_keys,
                &hash_masks,
            );
            self.pem_pos_buf.clear();
            self.pem_pos_buf.resize(window_count, -1);
            Self::scan_window_positions(
                &self.query_code_buf,
                linear_range.as_bytes(),
                &mut self.pem_pos_buf,
                window_size as usize,
                true,
                use_hash_fast_path,
                &hash_keys,
                &hash_masks,
            );
            let mut locus_count = 0usize;
            let mut locus2_count = 0usize;
            let mut locus_sum = 0i32;
            let mut locus2_sum = 0i32;
            let mut prev_locus: Option<i32> = None;
            let mut prev_locus2: Option<i32> = None;
            let mut miss_count = [0, 0, 0];
            let mut miss_count2_total = 0;
            if len_str % step != 0 {
                trial += 1;
            }
            for (&pos, &pos2) in self.circ_pos_buf.iter().zip(self.pem_pos_buf.iter()) {
                if pos >= 0 {
                    if let Some(prev) = prev_locus {
                        locus_sum += (pos - prev).abs();
                    }
                    prev_locus = Some(pos);
                    locus_count += 1;
                    miss_count[0] = 0;
                } else {
                    miss_count[1] += 1;
                    miss_count[0] += 1;
                    if miss_count[0] > miss_count[2] {
                        miss_count[2] = miss_count[0];
                    }
                }
                if pos2 >= 0 {
                    if let Some(prev) = prev_locus2 {
                        locus2_sum += (pos2 - prev).abs();
                    }
                    prev_locus2 = Some(pos2);
                    locus2_count += 1;
                } else {
                    miss_count2_total += 1;
                }
            }
            if miss_count2_total == 0 && miss_count[1] == 0 {
                if self.distance_loci_stats(locus_count, locus_sum, locus2_count, locus2_sum, step)
                    == 1
                {
                    return 1;
                } else {
                    return 0;
                }
            } else if miss_count2_total <= miss_count[1] {
                if locus2_count != 0 {
                    return 0;
                }
            } else if miss_count[2] > 5 || miss_count[1] * 2 > trial {
                continue;
            } else {
                return 1;
            }
        }
        0
    }

    /// Linear competition check for the right-anchored junction fragment
    /// (`Phase 1.2` in the original implementation).
    pub fn is_in_circ_rna_1_2(
        &mut self,
        len_str: i32,
        str_val: &str,
        circ_range_seq: &str,
        linear_range: &str,
    ) -> i32 {
        for &step in &self.window_unit {
            if len_str < step * 2 {
                continue;
            }
            let window_size = step * 2;
            let mut trial = (len_str - window_size) / step;
            let mut window_count = trial.max(0) as usize + 1;
            if len_str % step != 0 {
                window_count += 1;
            }
            self.query_code_buf.clear();
            self.query_code_buf
                .reserve(window_count.saturating_sub(self.query_code_buf.capacity()));
            for j in 0..=trial {
                let s_idx = j * step;
                let e_idx = j * step + window_size;
                self.query_code_buf.push(Self::encode_window(
                    &str_val.as_bytes()[s_idx as usize..e_idx as usize],
                ));
            }
            if len_str % step != 0 {
                self.query_code_buf.push(Self::encode_window(
                    &str_val.as_bytes()[(len_str - window_size) as usize..len_str as usize],
                ));
            }
            let mut hash_keys = [u64::MAX; 128];
            let mut hash_masks = [0u64; 128];
            let use_hash_fast_path =
                Self::build_query_hash(&self.query_code_buf, &mut hash_keys, &mut hash_masks);
            self.circ_pos_buf.clear();
            self.circ_pos_buf.resize(window_count, -1);
            Self::scan_window_positions(
                &self.query_code_buf,
                circ_range_seq.as_bytes(),
                &mut self.circ_pos_buf,
                window_size as usize,
                false,
                use_hash_fast_path,
                &hash_keys,
                &hash_masks,
            );
            self.pem_pos_buf.clear();
            self.pem_pos_buf.resize(window_count, -1);
            Self::scan_window_positions(
                &self.query_code_buf,
                linear_range.as_bytes(),
                &mut self.pem_pos_buf,
                window_size as usize,
                false,
                use_hash_fast_path,
                &hash_keys,
                &hash_masks,
            );
            let mut locus_count = 0usize;
            let mut locus2_count = 0usize;
            let mut locus_sum = 0i32;
            let mut locus2_sum = 0i32;
            let mut prev_locus: Option<i32> = None;
            let mut prev_locus2: Option<i32> = None;
            let mut miss_count = [0, 0, 0];
            let mut miss_count2_total = 0;
            if len_str % step != 0 {
                trial += 1;
            }
            for (&pos, &pos2) in self.circ_pos_buf.iter().zip(self.pem_pos_buf.iter()) {
                if pos >= 0 {
                    if let Some(prev) = prev_locus {
                        locus_sum += (pos - prev).abs();
                    }
                    prev_locus = Some(pos);
                    locus_count += 1;
                    miss_count[0] = 0;
                } else {
                    miss_count[1] += 1;
                    miss_count[0] += 1;
                    if miss_count[0] > miss_count[2] {
                        miss_count[2] = miss_count[0];
                    }
                }
                if pos2 >= 0 {
                    if let Some(prev) = prev_locus2 {
                        locus2_sum += (pos2 - prev).abs();
                    }
                    prev_locus2 = Some(pos2);
                    locus2_count += 1;
                } else {
                    miss_count2_total += 1;
                }
            }
            if miss_count2_total == 0 && miss_count[1] == 0 {
                if self.distance_loci_stats(locus_count, locus_sum, locus2_count, locus2_sum, step)
                    == 1
                {
                    return 1;
                } else {
                    return 0;
                }
            } else if miss_count2_total <= miss_count[1] {
                if locus2_count != 0 {
                    return 0;
                }
            } else if miss_count[2] > 5 || miss_count[1] * 2 > trial {
                continue;
            } else {
                return 1;
            }
        }
        0
    }

    /// Validates whether the rescued soft/unmapped sequence is still consistent
    /// with the circ region (`Phase 2`).
    pub fn is_in_circ_rna_2(&self, unmap_seq: &str, circ_range_seq: &str) -> i32 {
        let seq_len = unmap_seq.len() as i32;
        let window_step = 5;
        let window_size = 10;
        let mut miss_count3 = [0, 0, 0];
        let trial = (seq_len - window_size) / window_step;
        if trial < 0 {
            return 1;
        }
        for j in 0..=trial {
            let seq = if seq_len <= 10 {
                unmap_seq
            } else {
                &unmap_seq[(j * window_step) as usize..(j * window_step + window_size) as usize]
            };
            if circ_range_seq.contains(seq) {
                miss_count3[0] = 0;
            } else {
                miss_count3[1] += 1;
                miss_count3[0] += 1;
                if miss_count3[0] > miss_count3[2] {
                    miss_count3[2] = miss_count3[0];
                }
            }
        }
        if miss_count3[2] > 5 || (miss_count3[1] - 1) * 2 > trial {
            return 0;
        }
        1
    }

    /// Encodes one base into the compact 3-bit alphabet used by the rolling
    /// window matcher in `is_in_circ_rna_3`.
    #[inline]
    fn encode_window_base(base: u8) -> u64 {
        match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            b'N' | b'n' => 4,
            _ => 7,
        }
    }

    /// Encodes a fixed-size sequence window for fast equality checks.
    #[inline]
    fn encode_window(seq: &[u8]) -> u64 {
        let mut code = 0u64;
        for &base in seq {
            code = (code << 3) | Self::encode_window_base(base);
        }
        code
    }

    /// Builds the tiny query-code hash table used by the rolling matcher.
    ///
    /// The table only helps for the common small-query case; larger query sets
    /// still fall back to direct comparison to keep the helper simple.
    fn build_query_hash(
        query_codes: &[u64],
        hash_keys: &mut [u64; 128],
        hash_masks: &mut [u64; 128],
    ) -> bool {
        const EMPTY_CODE: u64 = u64::MAX;

        if query_codes.len() > 64 {
            return false;
        }

        hash_keys.fill(EMPTY_CODE);
        hash_masks.fill(0);
        for (idx, &code) in query_codes.iter().enumerate() {
            let mut slot = (code as usize).wrapping_mul(0x9E37_79B1usize) & (128 - 1);
            loop {
                if hash_keys[slot] == EMPTY_CODE {
                    hash_keys[slot] = code;
                    hash_masks[slot] = 1u64 << idx;
                    break;
                }
                if hash_keys[slot] == code {
                    hash_masks[slot] |= 1u64 << idx;
                    break;
                }
                slot = (slot + 1) & (128 - 1);
            }
        }
        true
    }

    /// Scans one haystack with a prebuilt matcher and records first/last hits.
    fn scan_window_positions(
        query_codes: &[u64],
        haystack: &[u8],
        out: &mut [i32],
        window_size: usize,
        keep_last: bool,
        use_hash_fast_path: bool,
        hash_keys: &[u64; 128],
        hash_masks: &[u64; 128],
    ) {
        const EMPTY_CODE: u64 = u64::MAX;

        if haystack.len() < window_size || query_codes.is_empty() {
            return;
        }

        let mut remaining = out.len();
        let keep_last_bases_mask = (1u64 << ((window_size as u32 - 1) * 3)) - 1;
        let mut rolling_code = Self::encode_window(&haystack[..window_size]);
        let last_start = haystack.len() - window_size;

        for start in 0..=last_start {
            if start != 0 {
                rolling_code = ((rolling_code & keep_last_bases_mask) << 3)
                    | Self::encode_window_base(haystack[start + window_size - 1]);
            }

            if use_hash_fast_path {
                let mut slot = (rolling_code as usize).wrapping_mul(0x9E37_79B1usize) & (128 - 1);
                let mut mask = 0u64;
                loop {
                    let key = hash_keys[slot];
                    if key == EMPTY_CODE {
                        break;
                    }
                    if key == rolling_code {
                        mask = hash_masks[slot];
                        break;
                    }
                    slot = (slot + 1) & (128 - 1);
                }
                while mask != 0 {
                    let idx = mask.trailing_zeros() as usize;
                    if keep_last {
                        out[idx] = start as i32;
                    } else if out[idx] == -1 {
                        out[idx] = start as i32;
                        remaining -= 1;
                    }
                    mask &= mask - 1;
                }
            } else {
                for (idx, &query_code) in query_codes.iter().enumerate() {
                    if query_code == rolling_code {
                        if keep_last {
                            out[idx] = start as i32;
                        } else if out[idx] == -1 {
                            out[idx] = start as i32;
                            remaining -= 1;
                        }
                    }
                }
            }

            if !keep_last && remaining == 0 {
                break;
            }
        }
    }

    /// Finds the first or last occurrence of each encoded query window in `haystack`.
    ///
    /// This is the core retained optimization from Scan1 tuning: the old version
    /// called `find`/`rfind` repeatedly for each window, while the current
    /// version scans the haystack once with a rolling code. For the common
    /// small-query case, a tiny fixed-capacity hash table avoids a nested linear
    /// scan without changing Java-visible decisions.
    fn find_window_positions(
        query_codes: &[u64],
        haystack: &[u8],
        out: &mut [i32],
        window_size: usize,
        keep_last: bool,
    ) {
        let mut hash_keys = [u64::MAX; 128];
        let mut hash_masks = [0u64; 128];
        let use_hash_fast_path =
            Self::build_query_hash(query_codes, &mut hash_keys, &mut hash_masks);
        Self::scan_window_positions(
            query_codes,
            haystack,
            out,
            window_size,
            keep_last,
            use_hash_fast_path,
            &hash_keys,
            &hash_masks,
        );
    }

    /// Validates the mate-support sequence against circ and linear competitor
    /// regions (`Phase 3`).
    ///
    /// The branch logic is still Java-shaped; the only substantive optimization is
    /// the buffer-reusing rolling-window matcher above.
    pub fn is_in_circ_rna_3(
        &mut self,
        ano_read: &str,
        pre_judge: &str,
        circ_range_seq: &str,
        pem_null_range_seq: &str,
    ) -> i32 {
        const WINDOW_STEP: usize = 5;
        const WINDOW_SIZE: usize = 10;
        let ano_read_len = ano_read.len() as i32;
        let trial = (ano_read_len - WINDOW_SIZE as i32) / WINDOW_STEP as i32;
        if trial < 0 {
            if !pem_null_range_seq.is_empty() {
                return -2;
            }
            return 1;
        }
        let mut miss_count = [0, 0, 0];
        let mut miss_count2_total = 0;
        let mut locus_count = 0usize;
        let mut locus2_count = 0usize;
        let mut locus_sum = 0i32;
        let mut locus2_sum = 0i32;
        let mut prev_locus: Option<i32> = None;
        let mut prev_locus2: Option<i32> = None;
        let ano_read_bytes = ano_read.as_bytes();
        let window_count = trial.max(0) as usize + 1;
        self.query_code_buf.clear();
        self.query_code_buf
            .reserve(window_count.saturating_sub(self.query_code_buf.capacity()));
        for j in 0..window_count {
            let start = j * WINDOW_STEP;
            let end = start + WINDOW_SIZE;
            self.query_code_buf
                .push(Self::encode_window(&ano_read_bytes[start..end]));
        }

        self.circ_pos_buf.clear();
        self.circ_pos_buf.resize(window_count, -1);
        Self::find_window_positions(
            &self.query_code_buf,
            circ_range_seq.as_bytes(),
            &mut self.circ_pos_buf,
            WINDOW_SIZE,
            false,
        );

        let has_pem_positions = !pem_null_range_seq.is_empty();
        if has_pem_positions {
            self.pem_pos_buf.clear();
            self.pem_pos_buf.resize(window_count, -1);
            Self::find_window_positions(
                &self.query_code_buf,
                pem_null_range_seq.as_bytes(),
                &mut self.pem_pos_buf,
                WINDOW_SIZE,
                false,
            );
        }

        for j in 0..window_count {
            let pos = self.circ_pos_buf[j];
            if pos >= 0 {
                if let Some(prev) = prev_locus {
                    locus_sum += (pos - prev).abs();
                }
                prev_locus = Some(pos);
                locus_count += 1;
                miss_count[0] = 0;
            } else {
                miss_count[1] += 1;
                miss_count[0] += 1;
                if miss_count[0] > miss_count[2] {
                    miss_count[2] = miss_count[0];
                }
            }

            if has_pem_positions {
                let pos2 = self.pem_pos_buf[j];
                if pos2 >= 0 {
                    if let Some(prev) = prev_locus2 {
                        locus2_sum += (pos2 - prev).abs();
                    }
                    prev_locus2 = Some(pos2);
                    locus2_count += 1;
                } else {
                    miss_count2_total += 1;
                }
            }
        }
        if !pem_null_range_seq.is_empty() {
            if miss_count2_total == 0 && miss_count[1] == 0 {
                if self.distance_loci_stats(
                    locus_count,
                    locus_sum,
                    locus2_count,
                    locus2_sum,
                    WINDOW_STEP as i32,
                ) == 1
                {
                    return 1;
                } else {
                    return -2;
                }
            } else if miss_count2_total <= miss_count[1] {
                if locus2_count != 0 {
                    return -2;
                } else {
                    return -1;
                }
            } else if miss_count[1] * 4 > trial * 3 && pre_judge == "0" {
                return -2;
            } else if miss_count[2] > 5 || miss_count[1] * 2 > trial {
                return -1;
            } else {
                return 1;
            }
        } else {
            if miss_count[1] * 4 > trial * 3 && pre_judge == "0" {
                return -2;
            } else if miss_count[2] > 5 || miss_count[1] * 2 > trial {
                return -1;
            }
        }
        1
    }

    /// High-level BSJ identification from Scan1 metadata.
    ///
    /// This is the main Scan1 validator entry. The function remains intentionally
    /// close to Java branch order, because any "cleaner" rewrite here tends to
    /// introduce subtle parity drift. Retained optimizations are limited to helper
    /// extraction and scratch reuse around the hottest inner checks.
    pub fn is_bsj_hg1(
        &mut self,
        circ_line_arr: &mut [String],
        chr_taga: &str,
        sum_q: i32,
        mitochondrion: &str,
        sp_label: bool,
        chr_exon_start_map: &HashMap<String, String>,
        chr_exon_end_map: &HashMap<String, String>,
    ) -> Option<String> {
        let hg1_started = scan1_profile_enabled().then(Instant::now);
        HG_PROFILE.hg1_calls.fetch_add(1, Ordering::Relaxed);
        let site1 = circ_line_arr[9].parse::<i32>().unwrap_or(0);
        let site2 = circ_line_arr[10].parse::<i32>().unwrap_or(0);
        let end_adjt1 = circ_line_arr[11].parse::<i32>().unwrap_or(0);
        let end_adjt2 = circ_line_arr[12].parse::<i32>().unwrap_or(0);
        let total_adjustment = end_adjt1 + end_adjt2;
        let chr_taga_len = chr_taga.len() as i32;

        let (end_string1, end_string2, tmp_site1, tmp_site2, adjt_bp);
        if end_adjt2 >= 0 {
            tmp_site1 = site1 - end_adjt1 - 1;
            tmp_site2 = site2 - end_adjt1 - 1;
            adjt_bp = 2 + total_adjustment;
            if site1 - end_adjt1 - 4 >= 0 {
                end_string1 =
                    java_substring(chr_taga, site1 - end_adjt1 - 4, end_adjt2 + site1).to_string();
            } else {
                let n_pad = (0 - (site1 - end_adjt1 - 4)) as usize;
                end_string1 = format!(
                    "{}{}",
                    "N".repeat(n_pad),
                    java_substring(chr_taga, 0, end_adjt2 + site1)
                );
            }
            if 3 + end_adjt2 + site2 > chr_taga_len {
                end_string2 =
                    java_substring(chr_taga, site2 - end_adjt1 - 1, chr_taga_len).to_string();
            } else {
                end_string2 =
                    java_substring(chr_taga, site2 - end_adjt1 - 1, 3 + end_adjt2 + site2)
                        .to_string();
            }
        } else {
            tmp_site1 = site1 + end_adjt2 - 1;
            tmp_site2 = site2 + end_adjt2 - 1;
            adjt_bp = 2 - total_adjustment;
            if site1 + end_adjt2 - 4 >= 0 {
                end_string1 =
                    java_substring(chr_taga, site1 + end_adjt2 - 4, site1 - end_adjt1).to_string();
            } else {
                let n_pad = (0 - (site1 + end_adjt2 - 4)) as usize;
                end_string1 = format!(
                    "{}{}",
                    "N".repeat(n_pad),
                    java_substring(chr_taga, 0, site1 - end_adjt1)
                );
            }
            if 3 + site2 - end_adjt1 > chr_taga_len {
                end_string2 =
                    java_substring(chr_taga, site2 + end_adjt2 - 1, chr_taga_len).to_string();
            } else {
                end_string2 =
                    java_substring(chr_taga, site2 + end_adjt2 - 1, 3 + site2 - end_adjt1)
                        .to_string();
            }
        }

        let index_compare_started = scan1_profile_enabled().then(Instant::now);
        let mut index_strand_map = if circ_line_arr[1] == mitochondrion || sp_label {
            IndexCompare::index_compare_chrm(&end_string1, &end_string2)
        } else {
            IndexCompare::index_compare(&end_string1, &end_string2)
        };
        if let Some(index_compare_started) = index_compare_started {
            HG_PROFILE
                .index_compare_calls
                .fetch_add(1, Ordering::Relaxed);
            add_ns(&HG_PROFILE.index_compare_ns, index_compare_started);
        }
        let chr = circ_line_arr[1].as_str();
        let mut start_key = String::with_capacity(chr.len() + 16);
        let mut end_key = String::with_capacity(chr.len() + 16);
        let exon_lookup_started = scan1_profile_enabled().then(Instant::now);
        for i in 0..=adjt_bp {
            start_key.clear();
            start_key.push_str(chr);
            start_key.push('\t');
            let _ = write!(&mut start_key, "{}", tmp_site1 + i);
            end_key.clear();
            end_key.push_str(chr);
            end_key.push('\t');
            let _ = write!(&mut end_key, "{}", tmp_site2 + i);
            if !index_strand_map.contains_key(&i) {
                if let (Some(gene_start), Some(gene_end)) = (
                    chr_exon_start_map.get(&start_key),
                    chr_exon_end_map.get(&end_key),
                ) {
                    if gene_start == gene_end {
                        let mut parts = gene_start.split('\t');
                        let _ = parts.next();
                        let strand = parts.next().unwrap_or("");
                        index_strand_map.insert(
                            i,
                            format!(
                                "{}\t{}\t{}\t{}",
                                i,
                                strand,
                                java_substring(chr_taga, tmp_site1 + i - 3, tmp_site1 + i - 1),
                                java_substring(chr_taga, tmp_site2 + i, tmp_site2 + i + 2)
                            ),
                        );
                    }
                }
            }
        }
        if let Some(exon_lookup_started) = exon_lookup_started {
            HG_PROFILE.exon_lookup_calls.fetch_add(1, Ordering::Relaxed);
            add_ns(&HG_PROFILE.exon_lookup_ns, exon_lookup_started);
        }

        if !index_strand_map.is_empty() {
            for (&shift, sig) in &index_strand_map {
                HG_PROFILE.shift_entries.fetch_add(1, Ordering::Relaxed);
                let mut shift_parts = sig.split('\t');
                let _shift_idx = shift_parts.next().unwrap_or("");
                let shift_strand = shift_parts.next().unwrap_or("");
                let shift_left = shift_parts.next().unwrap_or("");
                let shift_right = shift_parts.next().unwrap_or("");
                let diff_adjt = if end_adjt2 >= 0 {
                    shift - 1 - end_adjt1
                } else {
                    shift - 1 + total_adjustment - end_adjt1
                };
                let site1_new = site1 + diff_adjt;
                let site2_new = site2 + diff_adjt;
                let shift_prepare_started = scan1_profile_enabled().then(Instant::now);
                let mut str_new = ["".to_string(), "".to_string()];
                if diff_adjt >= 0 {
                    let str_adj = java_substring(&circ_line_arr[2], 0, diff_adjt);
                    str_new[1] = String::with_capacity(circ_line_arr[3].len() + str_adj.len());
                    str_new[1].push_str(&circ_line_arr[3]);
                    str_new[1].push_str(str_adj);
                    str_new[0] =
                        java_substring(&circ_line_arr[2], diff_adjt, circ_line_arr[2].len() as i32)
                            .to_string();
                } else {
                    let str_adj = java_substring(
                        &circ_line_arr[3],
                        circ_line_arr[3].len() as i32 + diff_adjt,
                        circ_line_arr[3].len() as i32,
                    );
                    str_new[0] = String::with_capacity(str_adj.len() + circ_line_arr[2].len());
                    str_new[0].push_str(str_adj);
                    str_new[0].push_str(&circ_line_arr[2]);
                    str_new[1] = java_substring(
                        &circ_line_arr[3],
                        0,
                        circ_line_arr[3].len() as i32 + diff_adjt,
                    )
                    .to_string();
                }
                {
                    let mut merged = String::with_capacity(shift_left.len() + str_new[0].len());
                    merged.push_str(shift_left);
                    merged.push_str(&str_new[0]);
                    str_new[0] = merged;
                }
                str_new[1].push_str(shift_right);
                let initial_seq1 = java_substring(&str_new[0], 0, self.initial_size1);
                let initial_seq2 = java_substring(
                    &str_new[1],
                    str_new[1].len() as i32 - self.initial_size1,
                    str_new[1].len() as i32,
                );
                let circ_range_seq = if site1_new - 3 < 0 && site2_new + 2 > chr_taga_len {
                    &chr_taga[..]
                } else if site1_new - 3 < 0 {
                    java_substring(chr_taga, 0, site2_new + 2)
                } else if site2_new + 2 > chr_taga_len {
                    java_substring(chr_taga, site1_new - 3, chr_taga_len)
                } else {
                    java_substring(chr_taga, site1_new - 3, site2_new + 2)
                };
                if let Some(shift_prepare_started) = shift_prepare_started {
                    HG_PROFILE
                        .shift_prepare_calls
                        .fetch_add(1, Ordering::Relaxed);
                    add_ns(&HG_PROFILE.shift_prepare_ns, shift_prepare_started);
                }
                if circ_range_seq.starts_with(initial_seq1)
                    && circ_range_seq.ends_with(initial_seq2)
                {
                    let right_linear_range =
                        self.linear_range_right(chr_taga, site1_new, site2_new, chr_taga_len);
                    let left_linear_range = self.linear_range_left(chr_taga, site1_new, site2_new);
                    if circ_line_arr[6] != "1" {
                        let linear_started = scan1_profile_enabled().then(Instant::now);
                        circ_line_arr[6] = self
                            .is_in_circ_rna_1_2(
                                str_new[0].len() as i32,
                                &str_new[0],
                                circ_range_seq,
                                right_linear_range,
                            )
                            .to_string();
                        if let Some(linear_started) = linear_started {
                            HG_PROFILE
                                .linear_check_calls
                                .fetch_add(1, Ordering::Relaxed);
                            add_ns(&HG_PROFILE.linear_check_ns, linear_started);
                        }
                    }
                    if circ_line_arr[7] != "1" {
                        let linear_started = scan1_profile_enabled().then(Instant::now);
                        circ_line_arr[7] = self
                            .is_in_circ_rna_1_2(
                                str_new[1].len() as i32,
                                &str_new[1],
                                circ_range_seq,
                                left_linear_range,
                            )
                            .to_string();
                        if let Some(linear_started) = linear_started {
                            HG_PROFILE
                                .linear_check_calls
                                .fetch_add(1, Ordering::Relaxed);
                            add_ns(&HG_PROFILE.linear_check_ns, linear_started);
                        }
                    }
                    if circ_line_arr[6] == "1" && circ_line_arr[7] == "1" {
                        if circ_line_arr[4] != "*" {
                            let circ2_started = scan1_profile_enabled().then(Instant::now);
                            let circ2_ok = self.is_in_circ_rna_2(&circ_line_arr[4], circ_range_seq);
                            if let Some(circ2_started) = circ2_started {
                                HG_PROFILE.circ2_calls.fetch_add(1, Ordering::Relaxed);
                                add_ns(&HG_PROFILE.circ2_ns, circ2_started);
                            }
                            if circ2_ok == 0 {
                                if let Some(hg1_started) = hg1_started {
                                    add_ns(&HG_PROFILE.hg1_total_ns, hg1_started);
                                }
                                return None;
                            }
                        }
                        let mut tag = 1;
                        if circ_line_arr[5].len() > 5 {
                            // Java parity: preserve original branch order/conditions.
                            let pem_null = if circ_line_arr[0] == "1"
                                && site1_new - site1_new + 5 >= self.linear_range_size_min
                            {
                                if 2 * site1_new >= site2_new + 6 {
                                    java_substring(
                                        chr_taga,
                                        2 * site1_new - site2_new - 6,
                                        site1_new - 1,
                                    )
                                } else {
                                    java_substring(chr_taga, 0, site1_new - 1)
                                }
                            } else if circ_line_arr[0] == "1" {
                                if site1_new >= self.linear_range_size_min + 1 {
                                    java_substring(
                                        chr_taga,
                                        site1_new - self.linear_range_size_min - 1,
                                        site1_new - 1,
                                    )
                                } else {
                                    java_substring(chr_taga, 0, site1_new - 1)
                                }
                            } else if circ_line_arr[0] == "0"
                                && site2_new - site1_new + 5 > self.linear_range_size_min
                            {
                                if 2 * site2_new - site1_new + 5 > chr_taga_len {
                                    java_substring(chr_taga, site2_new, chr_taga_len)
                                } else {
                                    java_substring(
                                        chr_taga,
                                        site2_new,
                                        2 * site2_new - site1_new + 5,
                                    )
                                }
                            } else {
                                if site2_new + self.linear_range_size_min > chr_taga_len {
                                    java_substring(chr_taga, site2_new, chr_taga_len)
                                } else {
                                    java_substring(
                                        chr_taga,
                                        site2_new,
                                        site2_new + self.linear_range_size_min,
                                    )
                                }
                            };
                            let circ3_started = scan1_profile_enabled().then(Instant::now);
                            tag = self.is_in_circ_rna_3(
                                &circ_line_arr[5],
                                &circ_line_arr[8],
                                circ_range_seq,
                                pem_null,
                            );
                            if let Some(circ3_started) = circ3_started {
                                HG_PROFILE.circ3_calls.fetch_add(1, Ordering::Relaxed);
                                add_ns(&HG_PROFILE.circ3_ns, circ3_started);
                            }
                        }
                        HG_PROFILE.hg1_hits.fetch_add(1, Ordering::Relaxed);
                        if let Some(hg1_started) = hg1_started {
                            add_ns(&HG_PROFILE.hg1_total_ns, hg1_started);
                        }
                        return Some(format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            tag,
                            circ_line_arr[1],
                            site1_new,
                            site2_new,
                            shift_strand,
                            shift_left,
                            shift_right,
                            sum_q
                        ));
                    }
                }
            }
        }
        if let Some(hg1_started) = hg1_started {
            add_ns(&HG_PROFILE.hg1_total_ns, hg1_started);
        }
        None
    }

    /// Low-level BSJ identification from Scan2 rescue metadata.
    ///
    /// Scan2 still passes a string-array payload here to match the historical
    /// Java data layout. Several attempts to replace that protocol with more typed
    /// structures did not deliver reliable speedups, so the current design favors
    /// auditability over abstraction.
    pub fn is_bsj_hg2(&mut self, circ_line_arr: &[String], chr_taga: &str) -> String {
        let hg2_started = scan2_profile_enabled().then(Instant::now);
        if hg2_started.is_some() {
            HG2_PROFILE.hg2_calls.fetch_add(1, Ordering::Relaxed);
        }
        trace_scan2_hg2("enter", circ_line_arr, "");
        let mut judge_tag = "3".to_string();
        let start_site = circ_line_arr[3].parse::<i32>().unwrap_or(0);
        let end_site = circ_line_arr[4].parse::<i32>().unwrap_or(0);
        let quant = circ_line_arr[12].parse::<i32>().unwrap_or(0);
        let chr_taga_len = chr_taga.len() as i32;
        let circ_range_seq = if start_site - 3 < 0 && end_site + 2 > chr_taga_len {
            &chr_taga[..]
        } else if start_site - 3 < 0 {
            java_substring(chr_taga, 0, end_site + 2)
        } else if end_site + 2 > chr_taga_len {
            java_substring(chr_taga, start_site - 3, chr_taga_len)
        } else {
            java_substring(chr_taga, start_site - 3, end_site + 2)
        };
        let circ_range_len = circ_range_seq.len() as i32;
        let str_full;
        let mut pem_null_range_seq = "";
        if circ_line_arr[2] == "sm" {
            str_full = format!("{}{}", circ_line_arr[5], circ_line_arr[11]);
            let len_str = str_full.len() as i32;
            if len_str < 7 {
                if java_substring(circ_range_seq, circ_range_len - len_str, circ_range_len)
                    != str_full
                {
                    trace_scan2_hg2("sm_short_mismatch", circ_line_arr, "");
                    return finish_scan2_hg2("0".to_string(), hg2_started);
                } else {
                    trace_scan2_hg2("sm_short_exact", circ_line_arr, "");
                    return finish_scan2_hg2("2".to_string(), hg2_started);
                }
            } else {
                let s1_part1 =
                    java_substring(circ_range_seq, circ_range_len - 4, circ_range_len - 2);
                let s1_part2 = java_substring(&str_full, len_str - 4, len_str - 2);
                let s2_part1 =
                    java_substring(circ_range_seq, circ_range_len - 7, circ_range_len - 4);
                let s2_part2 = java_substring(&str_full, len_str - 7, len_str - 4);
                if s1_part1 == s1_part2 || s2_part1 == s2_part2 {
                    let mut label = true;
                    let linear_range;
                    if end_site - start_site + 5 >= self.linear_range_size_min {
                        if 2 * start_site >= end_site + 6 {
                            linear_range = java_substring(
                                chr_taga,
                                2 * start_site - end_site - 6,
                                start_site - 1,
                            );
                        } else {
                            linear_range = java_substring(chr_taga, 0, start_site - 1);
                        }
                    } else if start_site - 1 >= self.linear_range_size_min {
                        linear_range = java_substring(
                            chr_taga,
                            start_site - self.linear_range_size_min - 1,
                            start_site - 1,
                        );
                    } else {
                        linear_range = java_substring(chr_taga, 0, start_site - 1);
                    }
                    if circ_range_len < len_str {
                        return finish_scan2_hg2("2".to_string(), hg2_started);
                    }
                    let sw_started = scan2_profile_enabled().then(Instant::now);
                    self.aligner.set_seq(
                        java_substring(circ_range_seq, circ_range_len - len_str, circ_range_len),
                        &str_full,
                    );
                    self.aligner.align();
                    if let Some(sw_started) = sw_started {
                        HG2_PROFILE.sw_calls.fetch_add(1, Ordering::Relaxed);
                        add_ns(&HG2_PROFILE.sw_ns, sw_started);
                    }
                    if self.aligner.score >= (len_str - (len_str - 2) / 10 * 2) {
                        if quant >= self.min_mapq_uni {
                            let initial_seq =
                                java_substring(&str_full, len_str - self.initial_size1, len_str);
                            if java_substring(
                                circ_range_seq,
                                circ_range_len - self.initial_size1,
                                circ_range_len,
                            ) == initial_seq
                            {
                                let linear11_started = scan2_profile_enabled().then(Instant::now);
                                let linear11_ok = self.is_in_circ_rna_1_1(
                                    len_str,
                                    &str_full,
                                    circ_range_seq,
                                    linear_range,
                                ) == 1;
                                if let Some(linear11_started) = linear11_started {
                                    HG2_PROFILE.linear11_calls.fetch_add(1, Ordering::Relaxed);
                                    add_ns(&HG2_PROFILE.linear11_ns, linear11_started);
                                }
                                if linear11_ok {
                                    judge_tag = "2".to_string();
                                } else {
                                    judge_tag = "3".to_string();
                                }
                            }
                        }
                        if circ_line_arr[0] == "1" {
                            pem_null_range_seq = linear_range;
                        }
                        label = false;
                    }
                    if label {
                        let initial_seq =
                            java_substring(&str_full, len_str - self.initial_size1, len_str);
                        if java_substring(
                            circ_range_seq,
                            circ_range_len - self.initial_size1,
                            circ_range_len,
                        ) == initial_seq
                        {
                            let linear11_started = scan2_profile_enabled().then(Instant::now);
                            let linear11_ok = self.is_in_circ_rna_1_1(
                                len_str,
                                &str_full,
                                circ_range_seq,
                                linear_range,
                            ) == 1;
                            if let Some(linear11_started) = linear11_started {
                                HG2_PROFILE.linear11_calls.fetch_add(1, Ordering::Relaxed);
                                add_ns(&HG2_PROFILE.linear11_ns, linear11_started);
                            }
                            if linear11_ok {
                                if quant >= self.min_mapq_uni {
                                    judge_tag = "2".to_string();
                                    if circ_line_arr[0] == "1" {
                                        pem_null_range_seq = linear_range;
                                    }
                                } else {
                                    trace_scan2_hg2("sm_linear11_low_mapq", circ_line_arr, "");
                                    return finish_scan2_hg2("2".to_string(), hg2_started);
                                }
                            } else {
                                trace_scan2_hg2("sm_linear11_fail", circ_line_arr, "");
                                return finish_scan2_hg2("0".to_string(), hg2_started);
                            }
                        } else {
                            trace_scan2_hg2("sm_initial_fail", circ_line_arr, "");
                            return finish_scan2_hg2("0".to_string(), hg2_started);
                        }
                    }
                } else {
                    trace_scan2_hg2("sm_seed_fail", circ_line_arr, "");
                    return finish_scan2_hg2("0".to_string(), hg2_started);
                }
            }
        } else {
            str_full = format!("{}{}", circ_line_arr[10], circ_line_arr[5]);
            let len_str = str_full.len() as i32;
            if len_str < 7 {
                if java_substring(circ_range_seq, 0, len_str) != str_full {
                    trace_scan2_hg2("ms_short_mismatch", circ_line_arr, "");
                    return finish_scan2_hg2("0".to_string(), hg2_started);
                } else {
                    trace_scan2_hg2("ms_short_exact", circ_line_arr, "");
                    return finish_scan2_hg2("2".to_string(), hg2_started);
                }
            } else {
                let s1_part1 = java_substring(circ_range_seq, 2, 4);
                let s1_part2 = java_substring(&str_full, 2, 4);
                let s2_part1 = java_substring(circ_range_seq, 4, 7);
                let s2_part2 = java_substring(&str_full, 4, 7);
                if s1_part1 == s1_part2 || s2_part1 == s2_part2 {
                    let mut label = true;
                    let linear_range;
                    if end_site - start_site + 5 >= self.linear_range_size_min {
                        if 2 * end_site - start_site + 5 > chr_taga_len {
                            linear_range = java_substring(chr_taga, end_site, chr_taga_len);
                        } else {
                            linear_range =
                                java_substring(chr_taga, end_site, 2 * end_site - start_site + 5);
                        }
                    } else if end_site + self.linear_range_size_min > chr_taga_len {
                        linear_range = java_substring(chr_taga, end_site, chr_taga_len);
                    } else {
                        linear_range = java_substring(
                            chr_taga,
                            end_site,
                            end_site + self.linear_range_size_min,
                        );
                    }
                    if circ_range_len < len_str {
                        return finish_scan2_hg2("2".to_string(), hg2_started);
                    }
                    let sw_started = scan2_profile_enabled().then(Instant::now);
                    self.aligner
                        .set_seq(java_substring(circ_range_seq, 0, len_str), &str_full);
                    self.aligner.align();
                    if let Some(sw_started) = sw_started {
                        HG2_PROFILE.sw_calls.fetch_add(1, Ordering::Relaxed);
                        add_ns(&HG2_PROFILE.sw_ns, sw_started);
                    }
                    if self.aligner.score >= (len_str - (len_str - 2) / 10 * 2) {
                        let initial_seq = java_substring(&str_full, 0, self.initial_size1);
                        if java_substring(circ_range_seq, 0, self.initial_size1) == initial_seq {
                            let linear12_started = scan2_profile_enabled().then(Instant::now);
                            let linear12_ok = self.is_in_circ_rna_1_2(
                                len_str,
                                &str_full,
                                circ_range_seq,
                                linear_range,
                            ) == 1;
                            if let Some(linear12_started) = linear12_started {
                                HG2_PROFILE.linear12_calls.fetch_add(1, Ordering::Relaxed);
                                add_ns(&HG2_PROFILE.linear12_ns, linear12_started);
                            }
                            if linear12_ok {
                                if quant >= self.min_mapq_uni {
                                    judge_tag = "2".to_string();
                                }
                            } else {
                                judge_tag = "3".to_string();
                            }
                        }
                        if circ_line_arr[0] == "0" {
                            pem_null_range_seq = linear_range;
                        }
                        label = false;
                    }
                    if label {
                        let initial_seq = java_substring(&str_full, 0, self.initial_size1);
                        if java_substring(circ_range_seq, 0, self.initial_size1) == initial_seq {
                            let linear12_started = scan2_profile_enabled().then(Instant::now);
                            let linear12_ok = self.is_in_circ_rna_1_2(
                                len_str,
                                &str_full,
                                circ_range_seq,
                                linear_range,
                            ) == 1;
                            if let Some(linear12_started) = linear12_started {
                                HG2_PROFILE.linear12_calls.fetch_add(1, Ordering::Relaxed);
                                add_ns(&HG2_PROFILE.linear12_ns, linear12_started);
                            }
                            if linear12_ok {
                                if quant >= self.min_mapq_uni {
                                    judge_tag = "2".to_string();
                                    if circ_line_arr[0] == "0" {
                                        pem_null_range_seq = linear_range;
                                    }
                                } else {
                                    trace_scan2_hg2("ms_linear12_low_mapq", circ_line_arr, "");
                                    return finish_scan2_hg2("2".to_string(), hg2_started);
                                }
                            } else {
                                trace_scan2_hg2("ms_linear12_fail", circ_line_arr, "");
                                return finish_scan2_hg2("0".to_string(), hg2_started);
                            }
                        } else {
                            trace_scan2_hg2("ms_initial_fail", circ_line_arr, "");
                            return finish_scan2_hg2("0".to_string(), hg2_started);
                        }
                    }
                } else {
                    trace_scan2_hg2("ms_seed_fail", circ_line_arr, "");
                    return finish_scan2_hg2("0".to_string(), hg2_started);
                }
            }
        }
        if circ_line_arr[7] != "*" {
            let circ2_started = scan2_profile_enabled().then(Instant::now);
            let circ2_ok = self.is_in_circ_rna_2(&circ_line_arr[7], circ_range_seq) != 0;
            if let Some(circ2_started) = circ2_started {
                HG2_PROFILE.circ2_calls.fetch_add(1, Ordering::Relaxed);
                add_ns(&HG2_PROFILE.circ2_ns, circ2_started);
            }
            if !circ2_ok {
                trace_scan2_hg2("circ2_fail", circ_line_arr, "");
                return finish_scan2_hg2("0".to_string(), hg2_started);
            }
        }
        if circ_line_arr[6].len() > 5 {
            let circ3_started = scan2_profile_enabled().then(Instant::now);
            let res = self.is_in_circ_rna_3(
                &circ_line_arr[6],
                &circ_line_arr[8],
                circ_range_seq,
                pem_null_range_seq,
            );
            if let Some(circ3_started) = circ3_started {
                HG2_PROFILE.circ3_calls.fetch_add(1, Ordering::Relaxed);
                add_ns(&HG2_PROFILE.circ3_ns, circ3_started);
            }
            trace_scan2_hg2(
                "final",
                circ_line_arr,
                &format!("judge_tag={} circ3_res={}", judge_tag, res),
            );
            return finish_scan2_hg2(format!("{}{}", res, judge_tag), hg2_started);
        }
        trace_scan2_hg2(
            "final_short_pair",
            circ_line_arr,
            &format!("judge_tag={}", judge_tag),
        );
        finish_scan2_hg2(format!("1{}", judge_tag), hg2_started)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_java_substring() {
        assert_eq!(java_substring("ABCDE", 1, 3), "BC");
        assert_eq!(java_substring("ABCDE", 0, 5), "ABCDE");
    }

    #[test]
    fn test_smith_waterman_exact() {
        let mut sw = SmithWaterman::new(1, -1, -1);
        sw.set_seq("ATGC", "ATGC");
        sw.align();
        assert_eq!(sw.score, 4);
        assert_eq!(sw.aligned_len, 4);
    }

    #[test]
    fn test_smith_waterman_mismatch() {
        let mut sw = SmithWaterman::new(1, -1, -1);
        sw.set_seq("ATGC", "ATGG");
        sw.align();
        assert_eq!(sw.score, 3);
    }

    #[test]
    fn test_smith_waterman_short() {
        let mut sw = SmithWaterman::new(1, -1, -1);
        sw.set_seq("GCAT", "GC");
        sw.align();
        assert_eq!(sw.score, 2);
    }
}
