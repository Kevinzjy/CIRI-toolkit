//! CLI implementation for `ciri-assemble` multi-sample isoform assembly.
//!
//! The command consumes completed second-pass sample prefixes rather than BAM
//! files. This keeps integrated assembly on the stable `.out + .segments` boundary
//! and lets users rerun usage/switching analysis without repeating expensive
//! read scanning.

use anyhow::{bail, Context, Result};
use clap::Parser;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::annotation::Annotation;
use crate::ciri_as::{run_ciri_assemble, CohortAssembleConfig, CohortSampleInput};
use crate::fasta::FastaReader;

/// Parsed command-line arguments for `ciri-assemble`.
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Integrative assemble of circRNA isoforms from multi-sample segments",
    long_about = None
)]
struct Args {
    /// Sample list TSV.
    #[arg(short = 'i', long = "in")]
    sample_list: String,

    /// Output prefix.
    #[arg(short = 'o', long = "out")]
    out_prefix: String,

    /// Reference genome FASTA.
    #[arg(short = 'r', long = "ref")]
    ref_fasta: String,

    /// GTF annotation file.
    #[arg(short = 'a', long = "anno")]
    gtf: String,
}

/// Parses arguments, loads reference resources, and runs multi-sample assembly.
pub fn main() -> Result<()> {
    let args = Args::parse();
    let samples = read_sample_list(&args.sample_list)?;
    let mut fasta = FastaReader::new();
    fasta.read_fasta(&args.ref_fasta)?;
    let mut annotation = Annotation::new();
    annotation.read_gtf(&args.gtf)?;
    let summary = run_ciri_assemble(CohortAssembleConfig {
        samples,
        out_prefix: &args.out_prefix,
        reference: &fasta.chr_tcga_map,
        annotation: Some(&annotation),
    })?;
    eprintln!(
        "Assembled {} isoforms from {} circRNAs across {} samples ({} isoform switching events)",
        summary.isoforms, summary.circ_rnas, summary.samples, summary.switching_circ_rnas
    );
    Ok(())
}

/// Reads a headerless two-column sample manifest.
///
/// The second column is a sample prefix. `ciri-assemble` derives
/// `<prefix>.out` and `<prefix>.segments` from it so the manifest stays compact
/// and matches the second-pass output contract.
fn read_sample_list(path: &str) -> Result<Vec<CohortSampleInput>> {
    let file = File::open(path).with_context(|| format!("open sample list {path}"))?;
    let reader = BufReader::new(file);
    let mut samples = Vec::new();
    let mut seen = HashSet::new();
    for (line_idx, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() != 2 {
            bail!(
                "{}:{} expected exactly two columns: sample_id and prefix",
                path,
                line_idx + 1
            );
        }
        if !seen.insert(fields[0].to_string()) {
            bail!(
                "{}:{} duplicate sample_id {}",
                path,
                line_idx + 1,
                fields[0]
            );
        }
        let prefix = fields[1].to_string();
        let out_path = format!("{}.out", prefix);
        let segments_path = format!("{}.segments", prefix);
        if !Path::new(&out_path).is_file() {
            bail!(
                "{}:{} missing derived .out file {}",
                path,
                line_idx + 1,
                out_path
            );
        }
        if !Path::new(&segments_path).is_file() {
            bail!(
                "{}:{} missing derived .segments file {}",
                path,
                line_idx + 1,
                segments_path
            );
        }
        samples.push(CohortSampleInput {
            sample_id: fields[0].to_string(),
            prefix,
        });
    }
    if samples.is_empty() {
        bail!("sample list {} did not contain any samples", path);
    }
    Ok(samples)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn sample_list_requires_headerless_sample_and_prefix() {
        let dir = tempdir().unwrap();
        let prefix = dir.path().join("s1");
        std::fs::write(prefix.with_extension("out"), "circRNA_ID\n").unwrap();
        std::fs::write(prefix.with_extension("segments"), "read_id\n").unwrap();
        let list = dir.path().join("samples.tsv");
        std::fs::write(&list, format!("sampleA\t{}\n", prefix.to_string_lossy())).unwrap();

        let samples = read_sample_list(list.to_str().unwrap()).unwrap();

        assert_eq!(samples.len(), 1);
        assert_eq!(samples[0].sample_id, "sampleA");
        assert_eq!(samples[0].prefix, prefix.to_string_lossy().to_string());
    }
}
