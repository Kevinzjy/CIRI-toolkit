//! CLI entry point for the `ciri-toolkit` binary.
//!
//! The binary is intentionally thin: it wires together the reference loaders and
//! the three pipeline stages, while the behaviorally sensitive logic stays in the
//! library modules for easier testing and parity verification.

use anyhow::Result;
use chrono::Local;
use clap::Parser;
use mimalloc::MiMalloc;
use std::time::Instant;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use ciri_toolkit::annotation::Annotation;
use ciri_toolkit::fasta::FastaReader;
use ciri_toolkit::sam_bam::{check_bam_sorting, detect_format, InputFormat};
use ciri_toolkit::scan1::Scan1;
use ciri_toolkit::scan2::Scan2;
use ciri_toolkit::summary::Summary;
use ciri_toolkit::utils::{
    bsj1_path_for_output, bsj2_path_for_output, bsj_path_for_output, fsj_path_for_output,
    parse_mem_str, result_path_for_output,
};

/// Parsed command-line arguments for the end-to-end pipeline.
///
/// Defaults are kept aligned with CIRI3 unless there is explicit evidence that a
/// different value is required for parity.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the input SAM/BAM file
    #[arg(short = 'i', long = "in")]
    in_sam: String,

    /// Output prefix; the pipeline writes `<prefix>.out/.bsj1/.bsj`
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
}

/// Emits one aligned, timestamped progress line.
fn log_info(label: &str, msg: &str) {
    let now = Local::now().format("%Y-%m-%d %H:%M:%S");
    println!("{} [INFO] {:<18}: {}", now, label, msg);
}

/// Loads inputs, runs Scan1 -> Scan2 -> Summary, and writes the final report.
fn main() -> Result<()> {
    let run_started = Instant::now();
    let args = Args::parse();
    let mem_limit = parse_mem_str(&args.mem_per_thread);
    let result_output = result_path_for_output(&args.out_prefix);
    let bsj1_output = bsj1_path_for_output(&args.out_prefix);
    let bsj_output = bsj_path_for_output(&args.out_prefix);
    let bsj2_output = bsj2_path_for_output(&args.out_prefix);
    let fsj_output = fsj_path_for_output(&args.out_prefix);

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()?;
    }

    log_info("Reference FASTA", &args.ref_fasta);
    let mut fasta = FastaReader::new();
    fasta.read_fasta(&args.ref_fasta)?;

    log_info(
        "Annotation GTF",
        args.gtf.as_deref().unwrap_or("None provided"),
    );
    let mut annotation = Annotation::new();
    if let Some(gtf_path) = &args.gtf {
        annotation.read_gtf(gtf_path)?;
    }

    let format = detect_format(&args.in_sam)?;
    let format_str = match format {
        InputFormat::Bam => {
            check_bam_sorting(&args.in_sam)?;
            "BAM (queryname-sorted)"
        }
        InputFormat::Sam => "SAM (text-based)",
    };
    log_info("Input format", format_str);

    // Stage boundaries are logged explicitly because most benchmarking and parity
    // work is reasoned about in terms of Scan1 / Scan2 / Summary timings.
    // 3. Scan 1
    log_info("Running scan 1", "Identifying back-spliced junctions...");
    let mut scan1 = Scan1::new(
        args.min_mapq,
        args.min_span,
        args.max_span,
        args.linear_range_size_min,
    );
    scan1.set_mem_limit(mem_limit);
    scan1.run(&args.in_sam, &bsj1_output, &fasta.chr_tcga_map, &annotation)?;
    log_info(
        "Scan 1 summary",
        &format!(
            "{} mapped reads, {} BSJ1 reads",
            scan1.mapped_reads, scan1.bsj1_reads
        ),
    );

    // 4. Indexing
    log_info("Loading BSJ sites", "generating candidate BSJ index...");
    let scan2_seq_len = (scan1.read_len - 12).max(1);
    // Scan1 is no longer needed once its output file and derived read length have
    // been materialized, so drop it before Scan2 to keep whole-genome RSS lower.
    drop(scan1);
    let mut scan2 = Scan2::new(args.min_mapq, args.linear_range_size_min, scan2_seq_len);
    scan2.set_mem_limit(mem_limit);
    scan2.build_index(&bsj1_output)?;

    // 5. Scan 2
    log_info("Running scan 2", "Curating splicing signals & counting FSJs...");
    scan2.run(
        &args.in_sam,
        &bsj1_output,
        &bsj_output,
        &bsj2_output,
        &fsj_output,
        &fasta.chr_tcga_map,
    )?;
    log_info(
        "Scan 2 summary",
        &format!("{} BSJ2 reads rescued", scan2.rescued_reads),
    );
    scan2.release_working_set();

    // 6. Finalization
    log_info("Post-processing", "Clustering sites and filtering results...");
    let mut summary = Summary::new(args.stringency);
    summary.run(
        &bsj_output,
        &result_output,
        &scan2.fsj_map,
        &fasta.chr_tcga_map,
        &annotation,
    )?;
    log_info(
        "Final summary",
        &format!(
            "{} circRNAs, {} BSJ reads detected",
            summary.circ_count, summary.final_bsj_reads
        ),
    );

    log_info("Output file", &result_output);
    log_info(
        "Total runtime",
        &format!("{:.2} seconds", run_started.elapsed().as_secs_f64()),
    );

    Ok(())
}
