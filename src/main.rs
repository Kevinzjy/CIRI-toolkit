use anyhow::Result;
use clap::Parser;
use mimalloc::MiMalloc;
use chrono::Local;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use ciri_toolkit::fasta::FastaReader;
use ciri_toolkit::annotation::Annotation;
use ciri_toolkit::scan1::Scan1;
use ciri_toolkit::scan2::Scan2;
use ciri_toolkit::summary::Summary;
use ciri_toolkit::sam_bam::{detect_format, check_bam_sorting, InputFormat};
use ciri_toolkit::utils::parse_mem_str;

/// CIRI-toolkit: High-performance circular RNA identification.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the input SAM/BAM file
    #[arg(short = 'i', long = "in")]
    in_sam: String,
    
    /// Prefix for the output files
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

    /// Maximum memory per thread (e.g., 2G, 512M)
    #[arg(short = 'M', long = "mem-per-thread", default_value = "2G")]
    mem_per_thread: String,
}

/// Formatted logging with aligned labels and timestamps.
fn log_info(label: &str, msg: &str) {
    let now = Local::now().format("%Y-%m-%d %H:%M:%S");
    println!("{} [INFO] {:<18}: {}", now, label, msg);
}

fn main() -> Result<()> {
    let args = Args::parse();
    let mem_limit = parse_mem_str(&args.mem_per_thread);

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global()?;
    }

    log_info("Reference FASTA", &args.ref_fasta);
    let mut fasta = FastaReader::new();
    fasta.read_fasta(&args.ref_fasta)?;

    log_info("Annotation GTF", args.gtf.as_deref().unwrap_or("None provided"));
    let mut annotation = Annotation::new();
    if let Some(gtf_path) = &args.gtf {
        annotation.read_gtf(gtf_path)?;
    }

    let format = detect_format(&args.in_sam)?;
    let format_str = match format {
        InputFormat::Bam => {
            check_bam_sorting(&args.in_sam)?;
            "BAM (queryname-sorted)"
        },
        InputFormat::Sam => "SAM (text-based)",
    };
    log_info("Input format", format_str);

    // 3. Scan 1
    log_info("Processing Scan 1", "Identifying Back-Spliced Junctions...");
    let bsj1_output = format!("{}.BSJ1", args.out_prefix);
    let mut scan1 = Scan1::new(args.min_mapq, args.min_span, args.max_span, args.linear_range_size_min);
    scan1.set_mem_limit(mem_limit);
    scan1.run(&args.in_sam, &args.out_prefix, &fasta.chr_tcga_map, &annotation)?;
    println!();

    // 4. Indexing
    log_info("Index Mapping", "Constructing candidate site lookup...");
    let scan2_seq_len = (scan1.read_len - 12).max(1);
    let mut scan2 = Scan2::new(args.min_mapq, args.linear_range_size_min, scan2_seq_len);
    scan2.set_mem_limit(mem_limit);
    scan2.build_index(&bsj1_output)?;

    // 5. Scan 2
    log_info("Processing Scan 2", "Rescuing signals & quantifying FSJ...");
    scan2.run(&args.in_sam, &bsj1_output, &fasta.chr_tcga_map)?;
    println!();

    // 6. Finalization
    log_info("Summarizing", "Clustering sites and filtering results...");
    let mut summary = Summary::new(args.stringency);
    summary.run(&bsj1_output, &args.out_prefix, &scan2.fsj_map, &fasta.chr_tcga_map, &annotation)?;

    log_info("Final Report", &format!("{}.result", args.out_prefix));
    
    Ok(())
}
