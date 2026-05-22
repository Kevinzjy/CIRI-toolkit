//! CLI implementation for `ciri-merge` cohort BSJ catalog construction.
//!
//! The command is intentionally narrower than the main `ciri` pipeline: it
//! consumes only completed first-pass `.out` files and writes a deterministic
//! BED6 cohort circRNA catalog. It does not read `.segments` or original
//! BAM/SAM files, so the output can be used as the external `--circ` catalog
//! for a later sample-specific second pass.

use anyhow::{anyhow, bail, Context, Result};
use clap::Parser;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Parsed command-line arguments for `ciri-merge`.
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Merge first-pass CIRI .out files into a circRNA bed file",
    long_about = None
)]
struct Args {
    /// First-pass CIRI `.out` files to merge.
    #[arg(short = 'i', long = "in", required = true, num_args = 1..)]
    inputs: Vec<String>,

    /// Output circRNA bed file path.
    #[arg(long = "out", short = 'o')]
    output: String,
}

/// One sample entry derived from a first-pass `.out` path.
#[derive(Debug, Clone)]
struct MergeSample {
    sample_id: String,
    out_path: String,
}

/// Cohort-level grouping key for one circRNA/BSJ boundary.
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
struct CircKey {
    chrom: String,
    start: i64,
    end: i64,
    strand: String,
}

/// Stable sort key that keeps common chromosome names in genomic order.
#[derive(Debug, Clone, Eq, PartialEq)]
struct ChromSortKey {
    group: u8,
    rank: u32,
    suffix: String,
}

impl Ord for ChromSortKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.group
            .cmp(&other.group)
            .then_with(|| self.rank.cmp(&other.rank))
            .then_with(|| self.suffix.cmp(&other.suffix))
    }
}

impl PartialOrd for ChromSortKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Per-sample evidence contributing to one cohort circRNA boundary.
#[derive(Debug, Clone)]
struct SampleCircSupport {
    source_circ_id: String,
    junction_reads: i64,
}

/// Aggregated cohort support for one [`CircKey`].
#[derive(Debug, Default)]
struct CohortCircSupport {
    by_sample: BTreeMap<String, SampleCircSupport>,
}

/// User-facing summary returned by the merge implementation.
#[derive(Debug, Default)]
struct MergeSummary {
    samples: usize,
    cohort_circs: usize,
    input_rows: usize,
    total_bsj_reads: i64,
}

/// Parses arguments, runs the merge, and prints a compact summary.
pub fn main() -> Result<()> {
    let args = Args::parse();
    let summary = run_merge(&args.inputs, &args.output)?;
    eprintln!(
        "Merged {} samples: {} circRNAs from {} rows, {} BSJ reads",
        summary.samples, summary.cohort_circs, summary.input_rows, summary.total_bsj_reads
    );
    Ok(())
}

/// Runs cohort circRNA catalog construction from first-pass `.out` files.
///
/// This function keeps only one compact support record per `(sample, circRNA
/// boundary)` in memory. It is therefore bounded by the number of detected
/// circRNA sites rather than by read count or segment rows.
fn run_merge(inputs: &[String], output_path: &str) -> Result<MergeSummary> {
    let samples = samples_from_inputs(inputs)?;
    if samples.is_empty() {
        bail!("no first-pass .out inputs were provided");
    }

    let mut cohort: HashMap<CircKey, CohortCircSupport> = HashMap::new();
    let mut input_rows = 0usize;
    let mut total_bsj_reads = 0i64;
    for sample in &samples {
        let rows = read_out_file(sample, &mut cohort)?;
        input_rows += rows.rows;
        total_bsj_reads += rows.total_bsj_reads;
    }

    let mut entries: Vec<_> = cohort.into_iter().collect();
    entries.sort_by(|(left, _), (right, _)| compare_circ_keys(left, right));

    write_circ_bed_catalog(output_path, &entries)?;

    Ok(MergeSummary {
        samples: samples.len(),
        cohort_circs: entries.len(),
        input_rows,
        total_bsj_reads,
    })
}

/// Builds merge sample entries from CLI input paths.
///
/// Sample IDs are derived from the basename after stripping one trailing `.out`
/// suffix. This keeps the common `sample1.out` and `sample1.ciri.out` cases
/// stable while avoiding a manifest file for the first implementation.
fn samples_from_inputs(inputs: &[String]) -> Result<Vec<MergeSample>> {
    let mut samples = Vec::with_capacity(inputs.len());
    let mut seen: HashMap<String, usize> = HashMap::new();
    for input in inputs {
        if !Path::new(input).is_file() {
            bail!("missing first-pass .out file: {input}");
        }
        let mut sample_id = sample_id_from_out_path(input)?;
        let count = seen.entry(sample_id.clone()).or_insert(0);
        if *count > 0 {
            sample_id = format!("{}_{}", sample_id, *count + 1);
        }
        *count += 1;
        samples.push(MergeSample {
            sample_id,
            out_path: input.clone(),
        });
    }
    Ok(samples)
}

/// Derives a stable sample ID from one `.out` path.
fn sample_id_from_out_path(path: &str) -> Result<String> {
    let name = Path::new(path)
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| anyhow!("cannot derive sample_id from path: {path}"))?;
    let sample = name.strip_suffix(".out").unwrap_or(name);
    if sample.is_empty() {
        bail!("empty sample_id derived from path: {path}");
    }
    Ok(sample.to_string())
}

#[derive(Debug, Default)]
struct OutReadSummary {
    rows: usize,
    total_bsj_reads: i64,
}

/// Streams one CIRI `.out` file into the cohort support map.
fn read_out_file(
    sample: &MergeSample,
    cohort: &mut HashMap<CircKey, CohortCircSupport>,
) -> Result<OutReadSummary> {
    let file = File::open(&sample.out_path)
        .with_context(|| format!("open first-pass .out {}", sample.out_path))?;
    let reader = BufReader::new(file);
    let mut summary = OutReadSummary::default();
    for (line_idx, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        if line.trim().is_empty() || line.starts_with('#') || line.starts_with("circRNA_ID\t") {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            bail!(
                "{}:{} has {} columns, expected at least 11",
                sample.out_path,
                line_idx + 1,
                fields.len()
            );
        }
        let circ_id = fields[0].to_string();
        let key = CircKey {
            chrom: fields[1].to_string(),
            start: fields[2].parse().with_context(|| {
                format!(
                    "parse circRNA_start at {}:{}",
                    sample.out_path,
                    line_idx + 1
                )
            })?,
            end: fields[3].parse().with_context(|| {
                format!("parse circRNA_end at {}:{}", sample.out_path, line_idx + 1)
            })?,
            strand: normalize_strand(fields[10]),
        };
        let junction_reads: i64 = fields[4].parse().with_context(|| {
            format!(
                "parse #junction_reads at {}:{}",
                sample.out_path,
                line_idx + 1
            )
        })?;
        if junction_reads < 0 {
            bail!(
                "{}:{} has negative #junction_reads: {}",
                sample.out_path,
                line_idx + 1,
                junction_reads
            );
        }
        let support = cohort.entry(key).or_default();
        support
            .by_sample
            .entry(sample.sample_id.clone())
            .and_modify(|existing| {
                existing.junction_reads += junction_reads;
                if existing.source_circ_id != circ_id {
                    existing.source_circ_id.push(',');
                    existing.source_circ_id.push_str(&circ_id);
                }
            })
            .or_insert(SampleCircSupport {
                source_circ_id: circ_id,
                junction_reads,
            });
        summary.rows += 1;
        summary.total_bsj_reads += junction_reads;
    }
    Ok(summary)
}

/// Normalizes missing or unknown CIRI strand labels for stable cohort keys.
fn normalize_strand(strand: &str) -> String {
    match strand {
        "+" | "-" => strand.to_string(),
        _ => "NA".to_string(),
    }
}

/// Writes the BED6 cohort circRNA catalog used by the future `--circ` second pass.
///
/// The first six columns follow BED6 conventions with 0-based half-open
/// coordinates. The score uses capped total BSJ support so the file remains
/// directly usable by BED-aware tools.
fn write_circ_bed_catalog(path: &str, entries: &[(CircKey, CohortCircSupport)]) -> Result<()> {
    let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(path)?);
    for (key, support) in entries {
        let total_bsj_reads = total_bsj_reads(support);
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            key.chrom,
            key.start - 1,
            key.end,
            cohort_circ_id(key),
            bed_score(total_bsj_reads),
            key.strand
        )?;
    }
    writer.flush()?;
    Ok(())
}

/// Converts total BSJ read support into the standard BED score range.
fn bed_score(total_bsj_reads: i64) -> i64 {
    total_bsj_reads.clamp(0, 1000)
}

/// Returns the stable cohort circRNA identifier for one merged BSJ boundary.
fn cohort_circ_id(key: &CircKey) -> String {
    format!("{}:{}|{}:{}", key.chrom, key.start, key.end, key.strand)
}

/// Sums BSJ read support across samples for one cohort circRNA.
fn total_bsj_reads(support: &CohortCircSupport) -> i64 {
    support
        .by_sample
        .values()
        .map(|support| support.junction_reads)
        .sum()
}

/// Compares cohort circRNA keys in coordinate order.
fn compare_circ_keys(left: &CircKey, right: &CircKey) -> Ordering {
    chrom_sort_key(&left.chrom)
        .cmp(&chrom_sort_key(&right.chrom))
        .then_with(|| left.start.cmp(&right.start))
        .then_with(|| left.end.cmp(&right.end))
        .then_with(|| left.strand.cmp(&right.strand))
}

/// Builds a human-genome-friendly chromosome key while keeping generic contigs deterministic.
fn chrom_sort_key(chrom: &str) -> ChromSortKey {
    let core = chrom.strip_prefix("chr").unwrap_or(chrom);
    if let Ok(rank) = core.parse::<u32>() {
        return ChromSortKey {
            group: 0,
            rank,
            suffix: String::new(),
        };
    }
    match core {
        "X" => ChromSortKey {
            group: 0,
            rank: 23,
            suffix: String::new(),
        },
        "Y" => ChromSortKey {
            group: 0,
            rank: 24,
            suffix: String::new(),
        },
        "M" | "MT" => ChromSortKey {
            group: 0,
            rank: 25,
            suffix: String::new(),
        },
        _ => ChromSortKey {
            group: 1,
            rank: 0,
            suffix: chrom.to_string(),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn merge_out_files_builds_bed_catalog() {
        let dir = tempdir().unwrap();
        let out1 = dir.path().join("s1.out");
        let out2 = dir.path().join("s2.out");
        std::fs::write(
            &out1,
            concat!(
                "circRNA_ID\tchr\tcircRNA_start\tcircRNA_end\t#junction_reads\tSM_MS_SMS\t#non_junction_reads\tjunction_reads_ratio\tcircRNA_type\tgene_id\tstrand\tjunction_reads_ID\tScore\n",
                "chr1:10|20\tchr1\t10\t20\t3\t0_0_0\t0\t1\texon\tG1\t+\tr1,r2,r3\t3\n",
                "chr2:30|40\tchr2\t30\t40\t5\t0_0_0\t0\t1\texon\tG2\t-\tr4\t5\n",
            ),
        )
        .unwrap();
        std::fs::write(
            &out2,
            concat!(
                "circRNA_ID\tchr\tcircRNA_start\tcircRNA_end\t#junction_reads\tSM_MS_SMS\t#non_junction_reads\tjunction_reads_ratio\tcircRNA_type\tgene_id\tstrand\tjunction_reads_ID\tScore\n",
                "chr1:10|20\tchr1\t10\t20\t7\t0_0_0\t0\t1\texon\tG1\t+\tr5\t7\n",
            ),
        )
        .unwrap();
        let output = dir.path().join("cohort.out");
        let summary = run_merge(
            &[
                out1.to_string_lossy().to_string(),
                out2.to_string_lossy().to_string(),
            ],
            output.to_str().unwrap(),
        )
        .unwrap();
        assert_eq!(summary.samples, 2);
        assert_eq!(summary.input_rows, 3);
        assert_eq!(summary.cohort_circs, 2);
        assert_eq!(summary.total_bsj_reads, 15);

        let circ = std::fs::read_to_string(&output).unwrap();
        assert!(!circ.contains("chrom\tchromStart\tchromEnd\tname\tscore\tstrand"));
        assert!(circ.contains("chr1\t9\t20\tchr1:10|20:+\t10\t+"));
        assert!(circ.contains("chr2\t29\t40\tchr2:30|40:-\t5\t-"));
        assert!(!output.with_extension("out.bsj_matrix.tsv").exists());
    }
}
