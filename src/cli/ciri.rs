//! User-facing CLI entry point for the `ciri` binary.
//!
//! This module owns argument parsing and orchestration for the user pipeline.
//! Keeping it separate from `src/bin/ciri.rs` lets the binary remain a thin
//! wrapper while the behaviorally sensitive logic stays testable in the library.

use anyhow::{bail, Result};
use chrono::Local;
use clap::Parser;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use crate::annotation::Annotation;
use crate::ciri_as::{rebuild_major_isoforms_from_segments, run_ciri_as, AsConfig};
use crate::fasta::FastaReader;
use crate::runtime::init_runtime;
use crate::sam_bam::{check_bam_sorting, detect_format, InputFormat};
use crate::scan1::Scan1;
use crate::scan2::Scan2;
use crate::summary::Summary;
use crate::utils::{
    bsj1_path_for_output, bsj2_path_for_output, bsj_path_for_output, fsj_path_for_output,
    log_path_for_output, parse_mem_str, perf_path_for_output, result_path_for_output,
    segments1_path_for_output, segments2_path_for_output, segments_non_bsj_path_for_output,
    trace_path_for_output,
};

/// Parsed command-line arguments for the end-to-end pipeline.
///
/// Defaults are kept aligned with CIRI3 unless there is explicit evidence that a
/// different value is required for parity.
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Run the CIRI Rust analysis pipeline",
    long_about = None
)]
struct Args {
    /// Path to the input SAM/BAM file.
    #[arg(short = 'i', long = "in")]
    in_sam: String,

    /// Output prefix; final outputs include `.out`, `.bsj`, `.segments`,
    /// `.isoforms.gtf`, and `.isoforms.fa` sidecars.
    #[arg(short = 'o', long = "out")]
    out_prefix: String,

    /// Path to the reference genome FASTA file
    #[arg(short = 'r', long = "ref")]
    ref_fasta: String,

    /// (Optional) Path to the GTF annotation file
    #[arg(short = 'a', long = "anno")]
    gtf: Option<String>,

    /// Minimum Mapping Quality (MAPQ) for candidate BSJ reads
    #[arg(short = 'm', long = "mapq", default_value_t = 10)]
    min_mapq: i32,

    /// Stringency level (0, 1, 2), Java default is 2
    #[arg(short = 's', long = "stringency", default_value_t = 2)]
    stringency: i32,

    /// Max spanning distance of circRNAs (Java -Max, default 200000)
    #[arg(long = "max-span", default_value_t = 200000)]
    max_span: i32,

    /// Min spanning distance of circRNAs (Java -Min, default 140)
    #[arg(long = "min-span", default_value_t = 140)]
    min_span: i32,

    /// Linear competition search range size (Java internal default 50000)
    #[arg(long = "linear-range-size-min", default_value_t = 50000)]
    linear_range_size_min: i32,

    /// Number of threads to use (default: auto)
    #[arg(short = 't', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Maximum memory per thread (e.g., 512M, 2G)
    #[arg(short = 'M', long = "mem-per-thread", default_value = "512M")]
    mem_per_thread: String,

    /// Comma-separated read IDs to trace through Scan1/Scan2.
    ///
    /// When set, detailed trace lines are written to `<prefix>.trace.log`.
    #[arg(long = "trace", value_name = "READS")]
    trace_reads: Option<String>,

    /// Keep internal pipeline temporary files for debugging.
    #[arg(long = "debug", default_value_t = false)]
    debug: bool,

    /// Enable profiling and write the report to `<prefix>.perf.log`.
    #[arg(long = "perf", default_value_t = false)]
    perf: bool,

    /// Resume from the latest completed merged checkpoint.
    ///
    /// The resume logic intentionally ignores shard-local `.part_*.tmp` files.
    /// If `<prefix>.segments` exists, only isoforms are rebuilt. Otherwise, if
    /// `<prefix>.out` and `<prefix>.bsj` exist, segments and isoforms are rebuilt
    /// from the input BAM/SAM. All normal required inputs still must be provided
    /// so `--continue` stays an execution-mode switch, not a separate CLI shape.
    #[arg(long = "continue", default_value_t = false)]
    continue_run: bool,
}

/// Emits one aligned, timestamped progress line.
///
/// Logs are mirrored to stdout and to `<prefix>.log` so the user can inspect a
/// durable run record after the process exits.
fn log_info(log_writer: &mut BufWriter<File>, label: &str, msg: &str) -> Result<()> {
    let now = Local::now().format("%Y-%m-%d %H:%M:%S");
    let line = format!("{} [INFO] {:<18}: {}", now, label, msg);
    println!("{}", line);
    writeln!(log_writer, "{}", line)?;
    log_writer.flush()?;
    Ok(())
}

/// Sorts the mate-level BSJ rows and writes the final user-facing `.bsj`.
///
/// This runs strictly after Summary has consumed only `priority=1` rows, so the
/// user-facing sort order and `priority=0` evidence cannot perturb `.out`.
fn write_display_bsj(final_bsj: &str, bsj1_path: &str, bsj2_path: &str) -> Result<()> {
    let mut rows: Vec<(String, usize, usize, usize, String)> = Vec::new();
    for (path, stage_rank, stage_name) in
        [(bsj1_path, 0usize, "scan1"), (bsj2_path, 1usize, "scan2")]
    {
        if let Ok(file) = File::open(path) {
            let reader = BufReader::new(file);
            for line_res in reader.lines() {
                let line = line_res?;
                if line.is_empty() {
                    continue;
                }
                let mut parts = line.splitn(4, '\t');
                if let (Some(read_id), Some(mate_label), Some(priority), Some(_)) =
                    (parts.next(), parts.next(), parts.next(), parts.next())
                {
                    let mate_rank = if mate_label == "R1" { 0 } else { 1 };
                    let priority_rank = if priority == "1" { 0 } else { 1 };
                    rows.push((
                        read_id.to_string(),
                        mate_rank,
                        priority_rank,
                        stage_rank,
                        format!("{}\t{}", line, stage_name),
                    ));
                }
            }
        }
    }
    rows.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| a.1.cmp(&b.1))
            .then_with(|| a.2.cmp(&b.2))
            .then_with(|| a.3.cmp(&b.3))
    });
    let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(final_bsj)?);
    for (_, _, _, _, line) in rows {
        writeln!(writer, "{}", line)?;
    }
    writer.flush()?;
    Ok(())
}

/// Loads inputs, runs Scan1 -> Scan2 -> Summary, and writes outputs.
///
/// The current `ciri` entry keeps the historical direct CIRI3-style arguments
/// instead of forcing a `detect` subcommand. After Summary finishes, the CLI now
/// always performs the read-level circRNA segments pass that feeds future
/// full-length reconstruction.
pub fn main() -> Result<()> {
    let run_started = Instant::now();
    let args = Args::parse();
    let mem_limit = parse_mem_str(&args.mem_per_thread);
    let result_output = result_path_for_output(&args.out_prefix);
    let log_output = log_path_for_output(&args.out_prefix);
    let trace_output = trace_path_for_output(&args.out_prefix);
    let perf_output = perf_path_for_output(&args.out_prefix);
    let bsj1_output = bsj1_path_for_output(&args.out_prefix);
    let bsj_output = bsj_path_for_output(&args.out_prefix);
    let bsj2_output = bsj2_path_for_output(&args.out_prefix);
    let segments1_output = segments1_path_for_output(&args.out_prefix);
    let segments2_output = segments2_path_for_output(&args.out_prefix);
    let segments_non_bsj_output = segments_non_bsj_path_for_output(&args.out_prefix);
    let fsj_output = fsj_path_for_output(&args.out_prefix);
    let mut log_writer = BufWriter::new(File::create(&log_output)?);

    init_runtime(
        args.trace_reads.as_deref(),
        args.trace_reads.as_ref().map(|_| trace_output.as_str()),
        args.perf.then_some(perf_output.as_str()),
    )?;

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()?;
    }

    log_info(&mut log_writer, "Reference FASTA", &args.ref_fasta)?;
    let mut fasta = FastaReader::new();
    fasta.read_fasta(&args.ref_fasta)?;

    log_info(
        &mut log_writer,
        "Annotation GTF",
        args.gtf.as_deref().unwrap_or("None provided"),
    )?;
    let mut annotation = Annotation::new();
    if let Some(gtf_path) = &args.gtf {
        annotation.read_gtf(gtf_path)?;
    }

    if args.trace_reads.is_some() {
        log_info(&mut log_writer, "Read trace", &trace_output)?;
    }
    if args.debug {
        log_info(
            &mut log_writer,
            "Debug temp",
            "Keeping internal temporary files after successful completion",
        )?;
    }
    if args.perf {
        log_info(&mut log_writer, "Perf report", &perf_output)?;
    }

    let input_path = args.in_sam.as_str();
    let format = detect_format(input_path)?;
    let format_str = match format {
        InputFormat::Bam => {
            check_bam_sorting(input_path)?;
            "BAM (queryname-sorted)"
        }
        InputFormat::Sam => "SAM (text-based)",
    };
    log_info(&mut log_writer, "Input format", format_str)?;

    let segments_output = format!("{}.segments", args.out_prefix);
    if args.continue_run {
        if Path::new(&segments_output).is_file() {
            if !Path::new(&result_output).is_file() {
                bail!(
                    "--continue found {} but missing required circRNA table {}",
                    segments_output,
                    result_output
                );
            }
            log_info(
                &mut log_writer,
                "Continue",
                "Resuming from completed .segments; rebuilding isoforms only...",
            )?;
            let major_isoforms = rebuild_major_isoforms_from_segments(
                &result_output,
                &segments_output,
                &args.out_prefix,
                &fasta.chr_tcga_map,
                args.gtf.as_ref().map(|_| &annotation),
            )?;
            log_info(
                &mut log_writer,
                "Isoforms output",
                &format!(
                    "{}.isoforms.gtf, {}.isoforms.fa",
                    args.out_prefix, args.out_prefix
                ),
            )?;
            log_info(
                &mut log_writer,
                "Isoforms summary",
                &format!("{} major isoforms", major_isoforms),
            )?;
            log_info(
                &mut log_writer,
                "Total runtime",
                &format!("{:.2} seconds", run_started.elapsed().as_secs_f64()),
            )?;
            return Ok(());
        }

        if Path::new(&result_output).is_file() && Path::new(&bsj_output).is_file() {
            log_info(
                &mut log_writer,
                "Continue",
                "Resuming from completed .out/.bsj; rebuilding segments and isoforms...",
            )?;
            let mut segment_progress_log =
                |label: &str, message: &str| log_info(&mut log_writer, label, message);
            let segment_summary = run_ciri_as(AsConfig {
                input_path,
                circ_path: &result_output,
                bsj_path: Some(&bsj_output),
                segment_evidence_paths: Vec::new(),
                non_bsj_segment_evidence_paths: Vec::new(),
                out_prefix: &args.out_prefix,
                keep_temp_files: args.debug,
                reference: &fasta.chr_tcga_map,
                annotation: args.gtf.as_ref().map(|_| &annotation),
                min_mapq: args.min_mapq,
                progress_log: Some(&mut segment_progress_log),
            })?;
            log_info(&mut log_writer, "Segments output", &segments_output)?;
            log_info(
                &mut log_writer,
                "Segments summary",
                &format!(
                    "{} rows ({} BSJ, {} backward, {} outward)",
                    segment_summary.total_segments,
                    segment_summary.bsj_segments,
                    segment_summary.backward_segments,
                    segment_summary.outward_segments
                ),
            )?;
            log_info(
                &mut log_writer,
                "Isoforms output",
                &format!(
                    "{}.isoforms.gtf, {}.isoforms.fa",
                    args.out_prefix, args.out_prefix
                ),
            )?;
            log_info(
                &mut log_writer,
                "Isoforms summary",
                &format!("{} major isoforms", segment_summary.major_isoforms),
            )?;
            log_info(
                &mut log_writer,
                "Total runtime",
                &format!("{:.2} seconds", run_started.elapsed().as_secs_f64()),
            )?;
            return Ok(());
        }

        bail!(
            "--continue did not find a resumable checkpoint: expected {} or both {} and {}",
            segments_output,
            result_output,
            bsj_output
        );
    }

    // Stage boundaries are logged explicitly because most benchmarking and parity
    // work is reasoned about in terms of Scan1 / Scan2 / Summary timings.
    // 3. Scan 1
    log_info(
        &mut log_writer,
        "Running scan 1",
        "Identifying back-spliced junctions...",
    )?;
    let mut scan1 = Scan1::new(
        args.min_mapq,
        args.min_span,
        args.max_span,
        args.linear_range_size_min,
    );
    scan1.set_mem_limit(mem_limit);
    scan1.run_with_priority_and_segments(
        input_path,
        &bsj1_output,
        Some(&segments1_output),
        &fasta.chr_tcga_map,
        &annotation,
    )?;
    log_info(
        &mut log_writer,
        "Scan 1 summary",
        &format!(
            "{} mapped reads, {} BSJ1 reads",
            scan1.mapped_reads, scan1.bsj1_reads
        ),
    )?;

    // 4. Indexing
    log_info(
        &mut log_writer,
        "Loading BSJ sites",
        "Generating candidate BSJ index...",
    )?;
    let scan2_seq_len = (scan1.read_len - 12).max(1);
    // Scan1 is no longer needed once its output file and derived read length have
    // been materialized, so drop it before Scan2 to keep whole-genome RSS lower.
    drop(scan1);
    let mut scan2 = Scan2::new(args.min_mapq, args.linear_range_size_min, scan2_seq_len);
    scan2.set_mem_limit(mem_limit);
    scan2.build_index(&bsj1_output)?;

    // 5. Scan 2
    log_info(
        &mut log_writer,
        "Running scan 2",
        "Curating splicing signals & counting FSJs...",
    )?;
    let scan2_segment_artifacts = scan2.run_with_display_and_segments(
        input_path,
        &bsj2_output,
        &fsj_output,
        Some(&bsj1_output),
        None,
        Some(&segments2_output),
        Some(&segments_non_bsj_output),
        &fasta.chr_tcga_map,
    )?;
    log_info(
        &mut log_writer,
        "Scan 2 summary",
        &format!("{} BSJ2 reads rescued", scan2.rescued_reads),
    )?;
    scan2.release_working_set();

    // 6. Finalization
    log_info(
        &mut log_writer,
        "Post-processing",
        "Clustering sites and filtering results...",
    )?;
    let mut summary = Summary::new(args.stringency);
    summary.run_from_bsj_files(
        &[&bsj1_output, &bsj2_output],
        &result_output,
        &scan2.fsj_map,
        &fasta.chr_tcga_map,
        &annotation,
    )?;
    log_info(
        &mut log_writer,
        "Final summary",
        &format!(
            // This count comes from Summary's retained circ/read assignments,
            // not from `Scan1 + Scan2` raw BSJ accumulation. Keeping the final
            // user-facing metric here prevents stage-local bookkeeping from
            // being mistaken for the clustered output size.
            "{} circRNAs, {} BSJ reads detected",
            summary.circ_count, summary.final_bsj_reads
        ),
    )?;
    scan2.release_working_set();

    log_info(
        &mut log_writer,
        "Formatting BSJ",
        "Sorting mate-level BSJ display...",
    )?;
    write_display_bsj(&bsj_output, &bsj1_output, &bsj2_output)?;

    log_info(&mut log_writer, "Output file", &result_output)?;
    log_info(
        &mut log_writer,
        "Running segments",
        "Reconstructing circRNA read-level segments...",
    )?;
    let mut segment_progress_log =
        |label: &str, message: &str| log_info(&mut log_writer, label, message);
    let non_bsj_segment_evidence_paths: Vec<&str> = scan2_segment_artifacts
        .non_bsj_segment_evidence_paths
        .iter()
        .map(String::as_str)
        .collect();
    let segment_summary_result = run_ciri_as(AsConfig {
        input_path,
        circ_path: &result_output,
        bsj_path: Some(&bsj_output),
        segment_evidence_paths: vec![&segments1_output, &segments2_output],
        non_bsj_segment_evidence_paths,
        out_prefix: &args.out_prefix,
        keep_temp_files: args.debug,
        reference: &fasta.chr_tcga_map,
        annotation: args.gtf.as_ref().map(|_| &annotation),
        min_mapq: args.min_mapq,
        progress_log: Some(&mut segment_progress_log),
    });
    let segment_summary = segment_summary_result?;
    log_info(
        &mut log_writer,
        "Segments output",
        &format!("{}.segments", args.out_prefix),
    )?;
    log_info(
        &mut log_writer,
        "Segments summary",
        &format!(
            "{} rows ({} BSJ, {} backward, {} outward)",
            segment_summary.total_segments,
            segment_summary.bsj_segments,
            segment_summary.backward_segments,
            segment_summary.outward_segments
        ),
    )?;
    log_info(
        &mut log_writer,
        "Isoforms output",
        &format!(
            "{}.isoforms.gtf, {}.isoforms.fa",
            args.out_prefix, args.out_prefix
        ),
    )?;
    log_info(
        &mut log_writer,
        "Isoforms summary",
        &format!("{} major isoforms", segment_summary.major_isoforms),
    )?;

    log_info(
        &mut log_writer,
        "Total runtime",
        &format!("{:.2} seconds", run_started.elapsed().as_secs_f64()),
    )?;

    if !args.debug {
        cleanup_pipeline_temp_files(
            &[
                &bsj1_output,
                &bsj2_output,
                &segments1_output,
                &segments2_output,
                &segments_non_bsj_output,
                &fsj_output,
            ],
            &scan2_segment_artifacts.non_bsj_segment_evidence_paths,
        );
    }

    Ok(())
}

/// Removes internal stage sidecars after the final user-facing outputs exist.
///
/// Official CLI outputs should remain focused on the final result files. The
/// intermediate `.bsj1/.bsj2/.segments1/.segments2` stems are implementation
/// details kept only when `--debug` is set, while shard paths are removed here
/// because Scan1/Scan2/segments may leave either merged stems or `.part_*.tmp`
/// files depending on which fast path was active.
fn cleanup_pipeline_temp_files(stems: &[&str], extra_paths: &[String]) {
    for path in extra_paths {
        let _ = std::fs::remove_file(path);
    }
    for stem in stems {
        let _ = std::fs::remove_file(stem);
        for idx in 0..rayon::current_num_threads().max(1) {
            let _ = std::fs::remove_file(crate::utils::part_path(stem, idx));
        }
    }
}
