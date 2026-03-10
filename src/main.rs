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

/// CIRI-toolkit: High-performance circular RNA identification.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the input SAM file
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

    /// Stringency level (0, 1, 2)
    #[arg(short = 's', long = "stringency", default_value_t = 1)]
    stringency: i32,

    /// Number of threads to use (default: auto)
    #[arg(short = 't', long = "threads", default_value_t = 0)]
    threads: usize,
}

/// Formatted logging with aligned labels and timestamps.
fn log_info(label: &str, msg: &str) {
    let now = Local::now().format("%Y-%m-%d %H:%M:%S");
    println!("{} [INFO] {:<18}: {}", now, label, msg);
}

fn main() -> Result<()> {
    let args = Args::parse();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global()?;
    }

    // 1. Reference Loading
    log_info("Reference FASTA", &args.ref_fasta);
    let mut fasta = FastaReader::new();
    fasta.read_fasta(&args.ref_fasta)?;

    log_info("Annotation GTF", args.gtf.as_deref().unwrap_or("None provided"));
    let mut annotation = Annotation::new();
    if let Some(gtf_path) = &args.gtf {
        annotation.read_gtf(gtf_path)?;
    }

    // 2. Format Discovery
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
    let mut scan1 = Scan1::new(args.min_mapq, 140, 200000, 5);
    scan1.run(&args.in_sam, &args.out_prefix, &fasta.chr_tcga_map, &annotation)?;
    println!(); // Clear line after progress bar

    // 4. Indexing
    log_info("Index Mapping", "Constructing candidate site lookup...");
    let mut scan2 = Scan2::new(args.min_mapq, 5, 100);
    scan2.build_index(&bsj1_output)?;

    // 5. Scan 2
    log_info("Processing Scan 2", "Rescuing signals & quantifying FSJ...");
    scan2.run(&args.in_sam, &bsj1_output, &fasta.chr_tcga_map)?;
    println!(); // Clear line after progress bar

    // 6. Finalization
    log_info("Summarizing", "Clustering sites and filtering results...");
    let mut summary = Summary::new(args.stringency);
    summary.run(&bsj1_output, &args.out_prefix, &scan2.fsj_map, &fasta.chr_tcga_map, &annotation)?;

    log_info("Final Report", &format!("{}.result", args.out_prefix));
    
    Ok(())
}
