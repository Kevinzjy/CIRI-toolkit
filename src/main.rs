use anyhow::Result;
use clap::Parser;
use mimalloc::MiMalloc;

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

fn main() -> Result<()> {
    let args = Args::parse();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global()?;
    }

    println!("Reading FASTA: {}", args.ref_fasta);
    let mut fasta = FastaReader::new();
    fasta.read_fasta(&args.ref_fasta)?;

    println!("Reading GTF: {:?}", args.gtf);
    let mut annotation = Annotation::new();
    if let Some(gtf_path) = args.gtf {
        annotation.read_gtf(&gtf_path)?;
    }

    let bsj1_output = format!("{}.BSJ1", args.out_prefix);
    
    // 4. Input Format Detection
    let format = detect_format(&args.in_sam)?;
    if matches!(format, InputFormat::Bam) {
        println!("[Input] BAM format detected. Checking sort order...");
        check_bam_sorting(&args.in_sam)?;
    }

    let mut scan1 = Scan1::new(args.min_mapq, 140, 200000, 5);
    println!("Running Scan 1: {} (Parallel)", args.in_sam);
    scan1.run(&args.in_sam, &args.out_prefix, &fasta.chr_tcga_map, &annotation)?;

    println!("Building Index for Scan 2...");
    let mut scan2 = Scan2::new(args.min_mapq, 5, 100);
    scan2.build_index(&bsj1_output)?;

    println!("Running Scan 2: {} (Parallel)", args.in_sam);
    scan2.run(&args.in_sam, &bsj1_output, &fasta.chr_tcga_map)?;

    println!("Summarizing Results...");
    let mut summary = Summary::new(args.stringency);
    summary.run(&bsj1_output, &args.out_prefix, &scan2.fsj_map, &fasta.chr_tcga_map, &annotation)?;

    println!("CIRI-toolkit finished. Final report generated at: {}.result", args.out_prefix);
    Ok(())
}
