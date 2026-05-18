//! Standalone simulator for CIRI-toolkit full-length reconstruction fixtures.
//!
//! The simulator is intentionally separate from the user-facing `ciri`
//! pipeline. It reads FASTA/GTF directly, samples circular isoforms and linear
//! background reads, writes paired FASTQ, emits a linear-only annotation GTF, and
//! records two truth tables: one row per circRNA and one row per read pair. The
//! output contract follows `docs/06-simulation-truth-design.md` and is designed
//! for CIRI-AS-style internal-structure/full-length isoform validation rather
//! than for reproducing the historical Perl simulator files.

use crate::fasta::FastaReader;
use anyhow::{bail, Context, Result};
use clap::Parser;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::env;
use std::f64::consts::TAU;
use std::fs::{create_dir_all, remove_file, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// Probability used by the initial simulator model to add a second isoform.
const SECOND_ISOFORM_PROBABILITY: f64 = 0.5;

/// CircRNAs at or below this length use an insert distribution centered near
/// the circle length so small circles can generate full-length-supporting read
/// pairs at useful rates.
const SHORT_CIRC_FULL_LENGTH_MAX_LEN: usize = 300;

/// Fraction of short-circ fragments drawn from the full-length-enriched insert
/// distribution instead of the ordinary PE mixture.
const SHORT_CIRC_FULL_LENGTH_BIAS_FRACTION: f64 = 0.70;

/// Minimum coordinate displacement required before simulator truth calls a
/// non-BSJ circular pair `type=outward`.
///
/// The simulator does not have mapper MAPQ, alternative-hit, or soft-clip
/// evidence, so truth-side outward is intentionally a stricter pair-geometry
/// proxy: both outward-facing block boundaries must shift by at least the same
/// 19 bp anchor used by BSJ-style terminal evidence.
const OUTWARD_TRUTH_MIN_PAIR_OFFSET: usize = 19;

/// Command-line options for the standalone simulator.
///
/// Defaults approximate a compact PE150 RNA-seq fixture. The circular and
/// linear abundance knobs are coverage means; sampled values are converted into
/// paired-end read-pair counts by the fixed formulas in the design document.
#[derive(Parser, Debug, Clone)]
#[command(
    author,
    version,
    about = "Generate CIRI simulator FASTQ and truth tables",
    long_about = None
)]
pub struct SimulateArgs {
    /// Reference genome FASTA.
    #[arg(short = 'r', long = "ref")]
    pub ref_fasta: String,

    /// Source GTF annotation used to choose circular and linear transcript exon chains.
    #[arg(short = 'a', long = "anno")]
    pub gtf: String,

    /// Output prefix; files are written as `<prefix>.*`.
    #[arg(short = 'o', long = "out")]
    pub out_prefix: String,

    /// Optional chromosome filter, for example `chr1`.
    #[arg(long = "chrom")]
    pub chrom: Option<String>,

    /// Number of circular RNA loci to simulate.
    #[arg(long = "circ-count", default_value_t = 100)]
    pub circ_count: usize,

    /// Mean circRNA coverage; per-circ coverage is sampled with `--scale`.
    #[arg(long = "circ-coverage", default_value_t = 10.0)]
    pub circ_coverage: f64,

    /// Mean linear transcript coverage used to size the background read set.
    #[arg(long = "linear-coverage", default_value_t = 0.1)]
    pub linear_coverage: f64,

    /// Relative standard deviation for Gaussian coverage sampling.
    #[arg(long = "scale", default_value_t = 0.5)]
    pub scale: f64,

    /// Read length for ordinary paired-end reads.
    #[arg(long = "read-len", default_value_t = 150)]
    pub read_len: usize,

    /// Mean fragment span in transcript/circular coordinates.
    #[arg(long = "insert-len", default_value_t = 260)]
    pub insert_len: usize,

    /// Standard deviation of the primary fragment span distribution.
    #[arg(long = "insert-sd", default_value_t = 40.0)]
    pub insert_sd: f64,

    /// Mean fragment span for the minor long-fragment distribution.
    #[arg(long = "insert-len-minor", default_value_t = 420)]
    pub insert_len_minor: usize,

    /// Standard deviation of the minor long-fragment distribution.
    #[arg(long = "insert-sd-minor", default_value_t = 60.0)]
    pub insert_sd_minor: f64,

    /// Fraction of fragments sampled from the minor long-fragment distribution.
    #[arg(long = "minor-insert-fraction", default_value_t = 0.10)]
    pub minor_insert_fraction: f64,

    /// Per-base substitution probability for ordinary reads.
    #[arg(long = "error-rate", default_value_t = 0.002)]
    pub error_rate: f64,

    /// Per-unique-exon probability of marking an exon coordinate circRNA-exclusive.
    #[arg(long = "exon-exclusive-rate", default_value_t = 0.25)]
    pub exon_exclusive_rate: f64,

    /// Deterministic seed controlling annotation sampling, isoforms and reads.
    #[arg(long = "seed", default_value_t = 5)]
    pub seed: u64,
}

/// User-facing summary returned by one simulator run.
///
/// The display text is intentionally stable because development scripts use it
/// for quick fixture sanity checks. Count wording is explicit: total/circ/linear
/// values are paired-end read pairs, while the first BSJ value counts individual
/// mates crossing the circular boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SimulationSummary {
    /// Number of simulated circRNA loci.
    pub circ_count: usize,
    /// Number of circular isoforms across all simulated circRNA loci.
    pub isoform_count: usize,
    /// Number of paired-end read pairs written to FASTQ.
    pub total_read_pairs: usize,
    /// Number of read pairs sampled from circular isoforms.
    pub circ_read_pairs: usize,
    /// Number of read pairs sampled from the filtered linear annotation.
    pub linear_read_pairs: usize,
    /// Number of individual mates crossing a circular boundary.
    pub bsj_reads: usize,
    /// Number of read pairs where at least one mate crosses a circular boundary.
    pub bsj_read_pairs: usize,
}

impl std::fmt::Display for SimulationSummary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Simulated {} circRNAs, {} circular isoforms\nTotal {} read pairs, {} circRNA read pairs, {} linear read pairs\nBSJ feature: {} reads / {} read pairs",
            self.circ_count,
            self.isoform_count,
            self.total_read_pairs,
            self.circ_read_pairs,
            self.linear_read_pairs,
            self.bsj_reads,
            self.bsj_read_pairs
        )
    }
}

/// One deduplicated exon coordinate used as a sampling unit.
///
/// `gene_id` and `transcript_id` are deliberately excluded so shared genomic
/// exons across transcripts or genes are marked exclusive together.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct ExonKey {
    chrom: String,
    start: usize,
    end: usize,
    strand: char,
}

/// One exon interval from the source GTF.
#[derive(Clone, Debug)]
struct Exon {
    start: usize,
    end: usize,
    exclusive: bool,
}

impl Exon {
    /// Returns the 1-based closed interval length.
    fn len(&self) -> usize {
        self.end.saturating_sub(self.start).saturating_add(1)
    }
}

/// Transcript-level exon chain used by circular or linear read simulation.
#[derive(Clone, Debug)]
struct Transcript {
    chrom: String,
    gene_id: String,
    transcript_id: String,
    strand: char,
    exons: Vec<Exon>,
}

/// Parsed GTF record needed to write the filtered linear annotation.
///
/// Storing the original line lets the simulator preserve raw GTF order and
/// attributes for surviving `gene`, `transcript` and `exon` records.
#[derive(Clone, Debug)]
struct GtfRecord {
    line: String,
    feature: String,
    gene_id: Option<String>,
    transcript_id: Option<String>,
    exon_key: Option<ExonKey>,
}

/// Parsed annotation state after circRNA-exclusive exon sampling.
#[derive(Clone, Debug)]
struct AnnotationModel {
    records: Vec<GtfRecord>,
    transcripts: Vec<Transcript>,
    linear_transcripts: Vec<Transcript>,
    exclusive_exons: BTreeSet<ExonKey>,
    surviving_genes: HashSet<String>,
    surviving_transcripts: HashSet<String>,
}

/// One simulated circRNA locus.
#[derive(Clone, Debug)]
struct SimCirc {
    circ_id: String,
    chrom: String,
    start: usize,
    end: usize,
    strand: char,
    gene_id: String,
    transcript_id: String,
    isoforms: Vec<SimIsoform>,
    sampled_circ_coverage: f64,
    read_pair_count: usize,
    bsj_read_pair_count: usize,
    bsj_read_count: usize,
}

/// One transcript-derived circular isoform and its read counters.
#[derive(Clone, Debug)]
struct SimIsoform {
    isoform_id: String,
    exon_chain: Vec<Exon>,
    seq: String,
    read_pair_count: usize,
    bsj_read_pair_count: usize,
}

/// Per-base source coordinate in transcript-oriented sequence space.
#[derive(Clone, Debug)]
struct SourceBase {
    coord: usize,
    exon_idx: usize,
}

/// Segment protocol fields for one mate.
#[derive(Clone, Debug)]
struct MateSegments {
    text: String,
    is_bsj: bool,
}

/// One contiguous genomic segment produced while traversing simulated source bases.
#[derive(Clone, Debug)]
struct SegmentToken {
    start: usize,
    end: usize,
    exon_idx: usize,
    partition: usize,
}

/// Small deterministic RNG so the simulator does not add another crate-level dependency.
#[derive(Clone, Debug)]
struct Lcg64 {
    state: u64,
}

impl Lcg64 {
    /// Creates a deterministic generator from the CLI seed.
    fn new(seed: u64) -> Self {
        Self {
            state: seed ^ 0x9E37_79B9_7F4A_7C15,
        }
    }

    /// Advances the generator and returns the next pseudo-random value.
    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }

    /// Returns a value in `0..upper`.
    fn gen_range(&mut self, upper: usize) -> usize {
        if upper <= 1 {
            0
        } else {
            (self.next_u64() as usize) % upper
        }
    }

    /// Returns a floating point value in `[0, 1)`.
    fn gen_f64(&mut self) -> f64 {
        const SCALE: f64 = (1u64 << 53) as f64;
        ((self.next_u64() >> 11) as f64) / SCALE
    }

    /// Samples one standard normal value with Box-Muller.
    fn gen_standard_normal(&mut self) -> f64 {
        let u1 = self.gen_f64().clamp(f64::MIN_POSITIVE, 1.0 - f64::EPSILON);
        let u2 = self.gen_f64();
        (-2.0 * u1.ln()).sqrt() * (TAU * u2).cos()
    }
}

/// Creates an output path by appending a stable suffix to the user prefix.
fn prefixed_path(prefix: &str, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}{}", prefix, suffix))
}

/// Extracts a quoted GTF attribute such as `gene_id` or `transcript_id`.
fn attr_value(attrs: &str, key: &str) -> Option<String> {
    for item in attrs.split(';') {
        let trimmed = item.trim();
        if !trimmed.starts_with(key) {
            continue;
        }
        return trimmed.split('"').nth(1).map(str::to_string);
    }
    None
}

/// Parses the GTF, samples circRNA-exclusive unique exon coordinates, and builds models.
///
/// The exclusive sampling unit is `(chrom, start, end, strand)`. All records
/// sharing that coordinate are removed from the filtered linear annotation and
/// from linear read simulation, while circular simulation continues to use the
/// complete source transcript set.
fn read_annotation_model(
    gtf: &str,
    chrom_filter: Option<&str>,
    exon_exclusive_rate: f64,
    rng: &mut Lcg64,
) -> Result<AnnotationModel> {
    let file = File::open(gtf).with_context(|| format!("cannot open GTF: {gtf}"))?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut transcript_exons: BTreeMap<(String, String, String), (char, Vec<(ExonKey, Exon)>)> =
        BTreeMap::new();
    let mut unique_exons = BTreeSet::new();

    for line_res in reader.lines() {
        let line = line_res?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        if let Some(chrom) = chrom_filter {
            if fields[0] != chrom {
                continue;
            }
        }
        let feature = fields[2].to_string();
        if !matches!(feature.as_str(), "gene" | "transcript" | "exon") {
            continue;
        }
        let chrom = fields[0].to_string();
        let start = fields[3].parse::<usize>()?;
        let end = fields[4].parse::<usize>()?;
        let strand = fields[6].chars().next().unwrap_or('+');
        let gene_id = attr_value(fields[8], "gene_id");
        let transcript_id = attr_value(fields[8], "transcript_id");
        let exon_key = (feature == "exon").then(|| ExonKey {
            chrom: chrom.clone(),
            start,
            end,
            strand,
        });

        if let Some(key) = &exon_key {
            unique_exons.insert(key.clone());
            let gene = gene_id.clone().unwrap_or_else(|| "NA".to_string());
            let tx = transcript_id
                .clone()
                .unwrap_or_else(|| format!("{gene}:tx"));
            transcript_exons
                .entry((chrom.clone(), gene, tx))
                .or_insert_with(|| (strand, Vec::new()))
                .1
                .push((
                    key.clone(),
                    Exon {
                        start,
                        end,
                        exclusive: false,
                    },
                ));
        }

        records.push(GtfRecord {
            line,
            feature,
            gene_id,
            transcript_id,
            exon_key,
        });
    }

    let mut exclusive_exons = BTreeSet::new();
    for exon in &unique_exons {
        if rng.gen_f64() <= exon_exclusive_rate {
            exclusive_exons.insert(exon.clone());
        }
    }

    let mut transcripts = Vec::new();
    let mut linear_transcripts = Vec::new();
    let mut surviving_genes = HashSet::new();
    let mut surviving_transcripts = HashSet::new();

    for ((chrom, gene_id, transcript_id), (strand, exon_items)) in transcript_exons {
        let mut full_exons = Vec::new();
        let mut linear_exons = Vec::new();
        for (key, mut exon) in exon_items {
            exon.exclusive = exclusive_exons.contains(&key);
            full_exons.push(exon.clone());
            if !exon.exclusive {
                linear_exons.push(exon);
            }
        }
        full_exons.sort_by_key(|exon| exon.start);
        linear_exons.sort_by_key(|exon| exon.start);
        if !full_exons.is_empty() {
            transcripts.push(Transcript {
                chrom: chrom.clone(),
                gene_id: gene_id.clone(),
                transcript_id: transcript_id.clone(),
                strand,
                exons: full_exons,
            });
        }
        if !linear_exons.is_empty() {
            surviving_genes.insert(gene_id.clone());
            surviving_transcripts.insert(transcript_id.clone());
            linear_transcripts.push(Transcript {
                chrom,
                gene_id,
                transcript_id,
                strand,
                exons: linear_exons,
            });
        }
    }

    Ok(AnnotationModel {
        records,
        transcripts,
        linear_transcripts,
        exclusive_exons,
        surviving_genes,
        surviving_transcripts,
    })
}

/// Writes the linear annotation GTF after removing circRNA-exclusive exon coordinates.
///
/// Only `gene`, `transcript` and `exon` records are emitted, and the surviving
/// records are written in the same order as the raw input GTF so annotation
/// diffs are deterministic and easy to inspect.
fn write_linear_annotation(model: &AnnotationModel, prefix: &str) -> Result<()> {
    let mut writer = BufWriter::new(File::create(prefixed_path(prefix, ".annotation.gtf"))?);
    for record in &model.records {
        match record.feature.as_str() {
            "gene" => {
                let Some(gene_id) = &record.gene_id else {
                    continue;
                };
                if model.surviving_genes.contains(gene_id) {
                    writeln!(writer, "{}", record.line)?;
                }
            }
            "transcript" => {
                let Some(transcript_id) = &record.transcript_id else {
                    continue;
                };
                if model.surviving_transcripts.contains(transcript_id) {
                    writeln!(writer, "{}", record.line)?;
                }
            }
            "exon" => {
                let Some(exon_key) = &record.exon_key else {
                    continue;
                };
                let Some(transcript_id) = &record.transcript_id else {
                    continue;
                };
                if !model.exclusive_exons.contains(exon_key)
                    && model.surviving_transcripts.contains(transcript_id)
                {
                    writeln!(writer, "{}", record.line)?;
                }
            }
            _ => {}
        }
    }
    Ok(())
}

/// Returns reverse-complemented DNA, mapping non-ACGT bases to `N`.
fn revcomp(seq: &str) -> String {
    seq.as_bytes()
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => 'T',
            b'T' => 'A',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'N',
        })
        .collect()
}

/// Slices one 1-based closed genomic exon interval from the reference.
fn genomic_slice<'a>(
    reference: &'a HashMap<String, String>,
    chrom: &str,
    exon: &Exon,
) -> Result<&'a str> {
    let chr_seq = reference
        .get(chrom)
        .with_context(|| format!("reference FASTA does not contain chromosome {chrom}"))?;
    if exon.start == 0 || exon.end > chr_seq.len() || exon.start > exon.end {
        bail!("invalid exon interval {chrom}:{}-{}", exon.start, exon.end);
    }
    Ok(&chr_seq[exon.start - 1..exon.end])
}

/// Builds transcript-oriented sequence for a selected exon chain.
///
/// Negative-strand transcripts are emitted in reverse genomic exon order with
/// each exon reverse-complemented, so sampled offsets follow transcript/read
/// orientation while truth coordinates remain 1-based closed genomic intervals.
fn isoform_sequence(
    reference: &HashMap<String, String>,
    chrom: &str,
    strand: char,
    exons: &[Exon],
) -> Result<String> {
    let mut seq = String::new();
    if strand == '-' {
        for exon in exons.iter().rev() {
            seq.push_str(&revcomp(genomic_slice(reference, chrom, exon)?));
        }
    } else {
        for exon in exons {
            seq.push_str(genomic_slice(reference, chrom, exon)?);
        }
    }
    Ok(seq)
}

/// Builds a per-base coordinate map in transcript-oriented sequence order.
///
/// The map lets read truth convert simulated offsets back into compact genomic
/// segment tokens without guessing across exon boundaries.
fn source_map_for_exons(exons: &[Exon], strand: char) -> Vec<SourceBase> {
    let mut map = Vec::new();
    if strand == '-' {
        for (exon_idx, exon) in exons.iter().enumerate().rev() {
            for coord in (exon.start..=exon.end).rev() {
                map.push(SourceBase { coord, exon_idx });
            }
        }
    } else {
        for (exon_idx, exon) in exons.iter().enumerate() {
            for coord in exon.start..=exon.end {
                map.push(SourceBase { coord, exon_idx });
            }
        }
    }
    map
}

/// Formats exon coordinates in the compact `start-end:strand` protocol.
fn exon_chain_string(exons: &[Exon], strand: char) -> String {
    exons
        .iter()
        .map(|exon| {
            let suffix = if exon.exclusive { "*" } else { "" };
            format!("{}-{}:{}{}", exon.start, exon.end, strand, suffix)
        })
        .collect::<Vec<_>>()
        .join(",")
}

/// Takes a sequence slice in circular coordinate space.
fn circular_slice(seq: &str, start: usize, len: usize) -> String {
    let bytes = seq.as_bytes();
    let mut out = String::with_capacity(len);
    for offset in 0..len {
        out.push(bytes[(start + offset) % bytes.len()] as char);
    }
    out
}

/// Takes a linear slice if it fits entirely within the sequence.
fn linear_slice(seq: &str, start: usize, len: usize) -> Option<&str> {
    seq.get(start..start + len)
}

/// Returns whether an interval crosses the circular boundary in template space.
fn crosses_boundary(start: usize, len: usize, seq_len: usize) -> bool {
    start + len > seq_len
}

/// Formats mate segments for one simulated read.
///
/// The simulator truth emits segments in read-chain order.
///
/// This preserves circRNA topology. A plus-strand BSJ read is represented as
/// the circ-end-side source blocks followed by `<bsj>` and then the
/// circ-start-side blocks; negative-strand isoforms use the transcript-oriented
/// source map, so walking forward through `source_map` gives the comparable
/// read-chain order after reverse-complement normalization.
fn format_segments(
    source_map: &[SourceBase],
    strand: char,
    start: usize,
    len: usize,
    circular: bool,
) -> MateSegments {
    let seq_len = source_map.len();
    let positions: Vec<usize> = (0..len)
        .map(|offset| {
            if circular {
                (start + offset) % seq_len
            } else {
                start + offset
            }
        })
        .collect();
    let mut segments = Vec::new();
    let mut current: Option<(usize, usize, usize, usize)> = None;
    let mut is_bsj = false;
    let mut partition = 0usize;
    let flush_current =
        |segments: &mut Vec<SegmentToken>, current: &mut Option<(usize, usize, usize, usize)>| {
            if let Some((min_coord, max_coord, exon_idx, partition)) = current.take() {
                segments.push(SegmentToken {
                    start: min_coord,
                    end: max_coord,
                    exon_idx,
                    partition,
                });
            }
        };

    for (idx, pos) in positions.iter().copied().enumerate() {
        if idx > 0 && circular {
            let prev = positions[idx - 1];
            if (prev == seq_len - 1 && pos == 0) || (prev == 0 && pos == seq_len - 1) {
                flush_current(&mut segments, &mut current);
                partition += 1;
                is_bsj = true;
            }
        }

        let base = &source_map[pos];
        match current {
            Some((ref mut min_coord, ref mut max_coord, exon_idx, existing_partition))
                if exon_idx == base.exon_idx
                    && existing_partition == partition
                    && (base.coord.abs_diff(*min_coord) == 1
                        || base.coord.abs_diff(*max_coord) == 1) =>
            {
                *min_coord = (*min_coord).min(base.coord);
                *max_coord = (*max_coord).max(base.coord);
            }
            Some(_) => {
                flush_current(&mut segments, &mut current);
                current = Some((base.coord, base.coord, base.exon_idx, partition));
            }
            None => current = Some((base.coord, base.coord, base.exon_idx, partition)),
        }
    }
    flush_current(&mut segments, &mut current);

    let mut out = Vec::with_capacity(segments.len() + usize::from(is_bsj));
    let mut previous_partition = None;
    for segment in segments {
        if previous_partition.is_some_and(|partition| partition != segment.partition) {
            out.push("<bsj>".to_string());
        }
        out.push(format!("{}-{}:{strand}", segment.start, segment.end));
        previous_partition = Some(segment.partition);
    }
    let text = out.join("|");

    MateSegments { text, is_bsj }
}

/// Returns the genomic span of the first read-chain block for one mate.
///
/// This mirrors the first mapper-visible block used by the simulator's
/// pair-orientation truth rule. The first block is enough here because
/// `type=outward` is a pair-level orientation class; downstream exon-junction
/// support is still carried by the full segment chain.
fn first_read_chain_segment_span(
    source_map: &[SourceBase],
    start: usize,
    len: usize,
    circular: bool,
) -> Option<(usize, usize)> {
    let seq_len = source_map.len();
    let mut iter = (0..len).map(|offset| {
        if circular {
            (start + offset) % seq_len
        } else {
            start + offset
        }
    });
    let first_pos = iter.next()?;
    let first = &source_map[first_pos];
    let mut min_coord = first.coord;
    let mut max_coord = first.coord;
    let exon_idx = first.exon_idx;
    let mut prev_pos = first_pos;

    for pos in iter {
        if circular
            && ((prev_pos == seq_len - 1 && pos == 0) || (prev_pos == 0 && pos == seq_len - 1))
        {
            break;
        }
        let base = &source_map[pos];
        if base.exon_idx != exon_idx
            || (base.coord.abs_diff(min_coord) != 1 && base.coord.abs_diff(max_coord) != 1)
        {
            break;
        }
        min_coord = min_coord.min(base.coord);
        max_coord = max_coord.max(base.coord);
        prev_pos = pos;
    }

    Some((min_coord, max_coord))
}

/// Tests the simulator-side outward-facing pair truth geometry.
///
/// R1 maps in the transcript strand direction, while R2 maps in the opposite
/// direction because its emitted sequence is reverse-complemented. A simulator
/// `type=outward` pair requires the reverse-oriented block to sit left of the
/// forward-oriented block with a clear 19 bp displacement at both block
/// boundaries. Overlap is not required because the truth label describes
/// observable pair orientation, not a RO-overlap subtype; however tiny shifts
/// are left as `forward` so evaluation does not reward ambiguous outward calls.
fn is_outward_facing_truth(
    source_map: &[SourceBase],
    strand: char,
    r1_start: usize,
    r2_start: usize,
    read_len: usize,
) -> bool {
    let Some(r1_span) = first_read_chain_segment_span(source_map, r1_start, read_len, true) else {
        return false;
    };
    let Some(r2_span) = first_read_chain_segment_span(source_map, r2_start, read_len, true) else {
        return false;
    };
    let (reverse_span, forward_span) = if strand == '-' {
        (r1_span, r2_span)
    } else {
        (r2_span, r1_span)
    };
    let start_offset = forward_span.0.saturating_sub(reverse_span.0);
    let end_offset = forward_span.1.saturating_sub(reverse_span.1);
    start_offset >= OUTWARD_TRUTH_MIN_PAIR_OFFSET && end_offset >= OUTWARD_TRUTH_MIN_PAIR_OFFSET
}

/// Samples a non-negative coverage value from the configured Gaussian model.
fn sample_coverage(mean: f64, scale: f64, rng: &mut Lcg64) -> f64 {
    if mean <= 0.0 {
        return 0.0;
    }
    if scale <= 0.0 {
        return mean;
    }
    (mean + rng.gen_standard_normal() * mean * scale).max(0.0)
}

/// Resolves the number of circular read pairs for one circRNA.
///
/// The `max(1, ...)` rule makes `--circ-count` mean "circRNAs with at least one
/// circular read pair" and is part of the simulator truth contract.
fn circular_pair_count(circ_len: usize, coverage: f64, read_len: usize) -> usize {
    ((circ_len as f64 * coverage) / (read_len as f64 * 2.0))
        .round()
        .max(1.0) as usize
}

/// Resolves the linear background size from sampled transcript coverage.
fn linear_pair_target(
    transcripts: &[Transcript],
    args: &SimulateArgs,
    sampled_linear_coverage: f64,
) -> usize {
    let total_bases: usize = transcripts
        .iter()
        .flat_map(|tx| tx.exons.iter())
        .map(Exon::len)
        .sum();
    ((total_bases as f64 * sampled_linear_coverage) / (args.read_len as f64 * 2.0))
        .round()
        .max(0.0) as usize
}

/// Samples a PE-style fragment span from a primary/minor normal mixture.
///
/// The result is clamped to at least the read length because paired-end geometry
/// in this simulator represents an outer fragment span.
fn sample_insert_len(args: &SimulateArgs, rng: &mut Lcg64) -> usize {
    let use_minor = rng.gen_f64() < args.minor_insert_fraction;
    let (mean, sd) = if use_minor {
        (args.insert_len_minor as f64, args.insert_sd_minor)
    } else {
        (args.insert_len as f64, args.insert_sd)
    };
    let sampled = mean + rng.gen_standard_normal() * sd;
    sampled.round().max(args.read_len as f64) as usize
}

/// Samples a fragment span for one circular molecule.
fn sample_circular_insert_len(seq_len: usize, args: &SimulateArgs, rng: &mut Lcg64) -> usize {
    if seq_len > SHORT_CIRC_FULL_LENGTH_MAX_LEN
        || rng.gen_f64() >= SHORT_CIRC_FULL_LENGTH_BIAS_FRACTION
    {
        return sample_insert_len(args, rng);
    }

    let mean = seq_len as f64;
    let sd = ((seq_len as f64) * 0.06).max(8.0);
    let lower = args.read_len as f64;
    let upper = (seq_len + args.read_len / 3).max(args.read_len) as f64;
    let sampled = mean + rng.gen_standard_normal() * sd;
    sampled.round().clamp(lower, upper) as usize
}

/// Applies simple substitution errors to sampled reads.
fn mutate_read(seq: &str, error_rate: f64, rng: &mut Lcg64) -> String {
    if error_rate <= 0.0 {
        return seq.to_string();
    }
    let mut out = String::with_capacity(seq.len());
    for base in seq.bytes() {
        let upper = base.to_ascii_uppercase();
        if !matches!(upper, b'A' | b'C' | b'G' | b'T') || rng.gen_f64() >= error_rate {
            out.push(upper as char);
            continue;
        }
        let choices = match upper {
            b'A' => [b'C', b'G', b'T'],
            b'C' => [b'A', b'G', b'T'],
            b'G' => [b'A', b'C', b'T'],
            b'T' => [b'A', b'C', b'G'],
            _ => unreachable!(),
        };
        out.push(choices[rng.gen_range(choices.len())] as char);
    }
    out
}

/// Writes one FASTQ record with stable metadata in the header.
fn write_fastq(
    writer: &mut BufWriter<File>,
    read_id: &str,
    seq: &str,
    qual: &str,
    metadata: &str,
) -> Result<()> {
    writeln!(writer, "@{read_id} {metadata}")?;
    writeln!(writer, "{seq}")?;
    writeln!(writer, "+")?;
    writeln!(writer, "{qual}")?;
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FastqCompressor {
    Pigz,
    Gzip,
}

impl FastqCompressor {
    /// Returns the external program name used for compression.
    fn program(self) -> &'static str {
        match self {
            Self::Pigz => "pigz",
            Self::Gzip => "gzip",
        }
    }
}

/// Selects the external FASTQ compressor before simulation work starts.
fn select_fastq_compressor() -> Result<FastqCompressor> {
    let compressor = select_fastq_compressor_from_availability(
        executable_exists("pigz"),
        executable_exists("gzip"),
    )?;
    if compressor == FastqCompressor::Gzip {
        eprintln!("pigz not found; falling back to gzip for simulator FASTQ compression");
    }
    Ok(compressor)
}

/// Returns whether an executable can be found in `PATH`.
fn executable_exists(program: &str) -> bool {
    if program.is_empty() {
        return false;
    }
    let Some(paths) = env::var_os("PATH") else {
        return false;
    };
    env::split_paths(&paths).any(|dir| is_executable_file(&dir.join(program)))
}

#[cfg(unix)]
fn is_executable_file(path: &Path) -> bool {
    use std::os::unix::fs::PermissionsExt;

    path.is_file()
        && path
            .metadata()
            .map(|metadata| metadata.permissions().mode() & 0o111 != 0)
            .unwrap_or(false)
}

#[cfg(not(unix))]
fn is_executable_file(path: &Path) -> bool {
    path.is_file()
}

fn select_fastq_compressor_from_availability(
    has_pigz: bool,
    has_gzip: bool,
) -> Result<FastqCompressor> {
    if has_pigz {
        Ok(FastqCompressor::Pigz)
    } else if has_gzip {
        Ok(FastqCompressor::Gzip)
    } else {
        bail!("failed to compress simulator FASTQ: neither pigz nor gzip was found in PATH")
    }
}

/// Compresses one completed simulator FASTQ with the selected compressor.
fn compress_fastq_to_gz(
    input_path: &Path,
    output_path: &Path,
    compressor: FastqCompressor,
) -> Result<()> {
    run_fastq_compressor(compressor.program(), input_path, output_path)?;
    remove_file(input_path)
        .with_context(|| format!("remove temporary FASTQ {}", input_path.display()))?;
    Ok(())
}

/// Runs one external FASTQ compressor and fails if it cannot produce output.
fn run_fastq_compressor(program: &str, input_path: &Path, output_path: &Path) -> Result<()> {
    let output = File::create(output_path)
        .with_context(|| format!("create compressed FASTQ {}", output_path.display()))?;
    let status = match Command::new(program)
        .arg("-n")
        .arg("-c")
        .arg(input_path)
        .stdout(Stdio::from(output))
        .stderr(Stdio::null())
        .status()
    {
        Ok(status) => status,
        Err(err) if err.kind() == std::io::ErrorKind::NotFound => {
            let _ = remove_file(output_path);
            bail!("{program} disappeared before FASTQ compression could run");
        }
        Err(err) => {
            let _ = remove_file(output_path);
            return Err(err).with_context(|| format!("run {program}"));
        }
    };

    if status.success() {
        Ok(())
    } else {
        let _ = remove_file(output_path);
        bail!(
            "{program} exited with {status} while compressing {}",
            input_path.display()
        )
    }
}

/// Selects circular exon windows from full source transcripts.
fn build_circs(
    transcripts: &[Transcript],
    reference: &HashMap<String, String>,
    args: &SimulateArgs,
    rng: &mut Lcg64,
) -> Result<Vec<SimCirc>> {
    let min_len = args.read_len.max(80);
    let mut candidates: Vec<usize> = transcripts
        .iter()
        .enumerate()
        .filter_map(|(idx, tx)| (tx.exons.len() >= 2).then_some(idx))
        .collect();
    for i in (1..candidates.len()).rev() {
        let j = rng.gen_range(i + 1);
        candidates.swap(i, j);
    }

    let mut circs = Vec::new();
    for tx_idx in candidates {
        if circs.len() >= args.circ_count {
            break;
        }
        let tx = &transcripts[tx_idx];
        let max_window = tx.exons.len().min(8);
        let min_window = 2;
        if max_window < min_window {
            continue;
        }
        let window_len = min_window + rng.gen_range(max_window - min_window + 1);
        let window_start = rng.gen_range(tx.exons.len() - window_len + 1);
        let mut exon_window = tx.exons[window_start..window_start + window_len].to_vec();
        exon_window.sort_by_key(|exon| exon.start);

        let seq = isoform_sequence(reference, &tx.chrom, tx.strand, &exon_window)?;
        if seq.len() < min_len {
            continue;
        }
        let sampled_circ_coverage = sample_coverage(args.circ_coverage, args.scale, rng);
        let read_pair_count = circular_pair_count(seq.len(), sampled_circ_coverage, args.read_len);
        let start = exon_window.iter().map(|exon| exon.start).min().unwrap();
        let end = exon_window.iter().map(|exon| exon.end).max().unwrap();
        let circ_id = format!("{}:{}|{}", tx.chrom, start, end);
        let mut isoforms = vec![SimIsoform {
            isoform_id: "isoform1".to_string(),
            exon_chain: exon_window.clone(),
            seq,
            read_pair_count: 0,
            bsj_read_pair_count: 0,
        }];

        if exon_window.len() >= 3
            && read_pair_count >= 2
            && rng.gen_f64() < SECOND_ISOFORM_PROBABILITY
        {
            let skip_idx = 1 + rng.gen_range(exon_window.len() - 2);
            let skipped: Vec<Exon> = exon_window
                .iter()
                .enumerate()
                .filter_map(|(idx, exon)| (idx != skip_idx).then_some(exon.clone()))
                .collect();
            let skipped_seq = isoform_sequence(reference, &tx.chrom, tx.strand, &skipped)?;
            if skipped_seq.len() >= min_len {
                isoforms.push(SimIsoform {
                    isoform_id: "isoform2".to_string(),
                    exon_chain: skipped,
                    seq: skipped_seq,
                    read_pair_count: 0,
                    bsj_read_pair_count: 0,
                });
            }
        }

        circs.push(SimCirc {
            circ_id,
            chrom: tx.chrom.clone(),
            start,
            end,
            strand: tx.strand,
            gene_id: tx.gene_id.clone(),
            transcript_id: tx.transcript_id.clone(),
            isoforms,
            sampled_circ_coverage,
            read_pair_count,
            bsj_read_pair_count: 0,
            bsj_read_count: 0,
        });
    }
    Ok(circs)
}

/// Chooses an isoform while ensuring each output isoform receives at least one read.
fn choose_isoform_for_pair(circ: &SimCirc, pair_idx: usize, rng: &mut Lcg64) -> usize {
    if pair_idx < circ.isoforms.len() {
        pair_idx
    } else {
        rng.gen_range(circ.isoforms.len())
    }
}

/// Emits ordinary paired-end reads from circular isoforms and records read truth.
fn write_circular_reads(
    circs: &mut [SimCirc],
    args: &SimulateArgs,
    rng: &mut Lcg64,
    r1_writer: &mut BufWriter<File>,
    r2_writer: &mut BufWriter<File>,
    read_truth: &mut BufWriter<File>,
    next_read_index: &mut usize,
) -> Result<()> {
    let qual = "I".repeat(args.read_len);
    for circ in circs {
        for pair_idx in 0..circ.read_pair_count {
            let iso_idx = choose_isoform_for_pair(circ, pair_idx, rng);
            let iso = &mut circ.isoforms[iso_idx];
            let seq_len = iso.seq.len();
            let source_map = source_map_for_exons(&iso.exon_chain, circ.strand);
            let insert_len = sample_circular_insert_len(seq_len, args, rng);
            let r1_start = rng.gen_range(seq_len);
            let r2_start = (r1_start + insert_len.saturating_sub(args.read_len)) % seq_len;
            let r1_template = circular_slice(&iso.seq, r1_start, args.read_len);
            let r2_template = circular_slice(&iso.seq, r2_start, args.read_len);
            let r1_seq = mutate_read(&r1_template, args.error_rate, rng);
            let r2_seq = mutate_read(&revcomp(&r2_template), args.error_rate, rng);
            let r1_segments =
                format_segments(&source_map, circ.strand, r1_start, args.read_len, true);
            let r2_segments =
                format_segments(&source_map, circ.strand, r2_start, args.read_len, true);
            let is_bsj = r1_segments.is_bsj || r2_segments.is_bsj;
            let pair_crosses_boundary = crosses_boundary(r1_start, insert_len, seq_len);
            // Simulator truth is restricted to read-level observable evidence.
            // A read crossing the BSJ is already `bsj`; there is no separate
            // truth-side `backward` bucket. Non-BSJ circular-origin fragments
            // become `outward` only when their pair orientation is observable.
            let is_outward = !is_bsj
                && pair_crosses_boundary
                && is_outward_facing_truth(
                    &source_map,
                    circ.strand,
                    r1_start,
                    r2_start,
                    args.read_len,
                );
            let is_circular = is_bsj || is_outward;
            let truth_type = if is_bsj {
                "bsj"
            } else if is_outward {
                "outward"
            } else {
                "forward"
            };
            let read_id = format!("sim:{}", *next_read_index);
            let metadata = format!(
                "read_kind=circ circ_id={} isoform_id={} insert_len={}",
                circ.circ_id, iso.isoform_id, insert_len
            );
            write_fastq(
                r1_writer,
                &format!("{read_id}/1"),
                &r1_seq,
                &qual,
                &metadata,
            )?;
            write_fastq(
                r2_writer,
                &format!("{read_id}/2"),
                &r2_seq,
                &qual,
                &metadata,
            )?;
            writeln!(
                read_truth,
                "{read_id}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{truth_type}",
                circ.circ_id,
                circ.chrom,
                circ.start,
                circ.end,
                circ.strand,
                iso.isoform_id,
                usize::from(is_circular),
                usize::from(is_bsj),
                r1_segments.text,
                usize::from(r1_segments.is_bsj),
                r2_segments.text,
                usize::from(r2_segments.is_bsj)
            )?;
            iso.read_pair_count += 1;
            if is_bsj {
                iso.bsj_read_pair_count += 1;
                circ.bsj_read_pair_count += 1;
                circ.bsj_read_count +=
                    usize::from(r1_segments.is_bsj) + usize::from(r2_segments.is_bsj);
            }
            *next_read_index += 1;
        }
    }
    Ok(())
}

/// Emits linear background reads from the filtered annotation.
fn write_linear_reads(
    transcripts: &[Transcript],
    reference: &HashMap<String, String>,
    args: &SimulateArgs,
    sampled_linear_coverage: f64,
    rng: &mut Lcg64,
    r1_writer: &mut BufWriter<File>,
    r2_writer: &mut BufWriter<File>,
    read_truth: &mut BufWriter<File>,
    next_read_index: &mut usize,
) -> Result<()> {
    let target_pairs = linear_pair_target(transcripts, args, sampled_linear_coverage);
    if target_pairs == 0 {
        return Ok(());
    }
    let usable: Vec<&Transcript> = transcripts
        .iter()
        .filter(|tx| tx.exons.len() >= 2)
        .collect();
    if usable.is_empty() {
        bail!("no usable transcripts for linear read simulation");
    }
    let qual = "I".repeat(args.read_len);
    let mut emitted = 0usize;
    let mut attempts = 0usize;
    while emitted < target_pairs && attempts < target_pairs * 30 {
        attempts += 1;
        let tx = usable[rng.gen_range(usable.len())];
        let tx_seq = isoform_sequence(reference, &tx.chrom, tx.strand, &tx.exons)?;
        let source_map = source_map_for_exons(&tx.exons, tx.strand);
        let insert_len = sample_insert_len(args, rng);
        if tx_seq.len() <= insert_len || insert_len < args.read_len {
            continue;
        }
        let start = rng.gen_range(tx_seq.len() - insert_len + 1);
        let r2_start = start + insert_len - args.read_len;
        let Some(r1_template) = linear_slice(&tx_seq, start, args.read_len) else {
            continue;
        };
        let Some(r2_template) = linear_slice(&tx_seq, r2_start, args.read_len) else {
            continue;
        };
        let r1_seq = mutate_read(r1_template, args.error_rate, rng);
        let r2_seq = mutate_read(&revcomp(r2_template), args.error_rate, rng);
        let r1_segments = format_segments(&source_map, tx.strand, start, args.read_len, false);
        let r2_segments = format_segments(&source_map, tx.strand, r2_start, args.read_len, false);
        let read_id = format!("sim:{}", *next_read_index);
        let metadata = format!(
            "read_kind=linear gene_id={} transcript_id={} insert_len={}",
            tx.gene_id, tx.transcript_id, insert_len
        );
        write_fastq(
            r1_writer,
            &format!("{read_id}/1"),
            &r1_seq,
            &qual,
            &metadata,
        )?;
        write_fastq(
            r2_writer,
            &format!("{read_id}/2"),
            &r2_seq,
            &qual,
            &metadata,
        )?;
        writeln!(
            read_truth,
            "{read_id}\tNA\t{}\tNA\tNA\t{}\tNA\t0\t0\t{}\t0\t{}\t0\tforward",
            tx.chrom, tx.strand, r1_segments.text, r2_segments.text
        )?;
        *next_read_index += 1;
        emitted += 1;
    }
    if emitted < target_pairs {
        bail!(
            "only emitted {emitted} linear read pairs after {attempts} attempts; try reducing --insert-len or --linear-coverage"
        );
    }
    Ok(())
}

/// Writes the circRNA/isoform truth table after read counts are known.
fn write_isoform_truth(circs: &[SimCirc], prefix: &str) -> Result<()> {
    let mut writer = BufWriter::new(File::create(prefixed_path(prefix, ".isoforms.tsv"))?);
    writeln!(
        writer,
        "circ_id\tchrom\tstart\tend\tstrand\tgene_id\ttranscript_id\tcoverage\tread_cnt\tbsj_read_cnt\tisoform_cnt\tisoform_exons\tisoform_len\tisoform_read_cnt\tisoform_bsj_read_cnt"
    )?;
    for circ in circs {
        let isoform_exons = circ
            .isoforms
            .iter()
            .map(|iso| {
                format!(
                    "{}={}",
                    iso.isoform_id,
                    exon_chain_string(&iso.exon_chain, circ.strand)
                )
            })
            .collect::<Vec<_>>()
            .join(";");
        let isoform_len = circ
            .isoforms
            .iter()
            .map(|iso| format!("{}={}", iso.isoform_id, iso.seq.len()))
            .collect::<Vec<_>>()
            .join(";");
        let isoform_read_cnt = circ
            .isoforms
            .iter()
            .map(|iso| format!("{}={}", iso.isoform_id, iso.read_pair_count))
            .collect::<Vec<_>>()
            .join(";");
        let isoform_bsj_read_cnt = circ
            .isoforms
            .iter()
            .map(|iso| format!("{}={}", iso.isoform_id, iso.bsj_read_pair_count))
            .collect::<Vec<_>>()
            .join(";");
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            circ.circ_id,
            circ.chrom,
            circ.start,
            circ.end,
            circ.strand,
            circ.gene_id,
            circ.transcript_id,
            circ.sampled_circ_coverage,
            circ.read_pair_count,
            circ.bsj_read_pair_count,
            circ.isoforms.len(),
            isoform_exons,
            isoform_len,
            isoform_read_cnt,
            isoform_bsj_read_cnt
        )?;
    }
    Ok(())
}

/// Ensures the user arguments describe a feasible simulation.
fn validate_args(args: &SimulateArgs) -> Result<()> {
    if args.circ_count == 0 {
        bail!("--circ-count must be > 0");
    }
    if args.read_len == 0 {
        bail!("--read-len must be > 0");
    }
    if args.insert_len < args.read_len {
        bail!("--insert-len must be >= --read-len");
    }
    if args.insert_len_minor < args.read_len {
        bail!("--insert-len-minor must be >= --read-len");
    }
    if args.insert_sd < 0.0 || args.insert_sd_minor < 0.0 {
        bail!("insert standard deviations must be >= 0");
    }
    if !(0.0..=1.0).contains(&args.minor_insert_fraction) {
        bail!("--minor-insert-fraction must be in [0, 1]");
    }
    if args.circ_coverage < 0.0 || args.linear_coverage < 0.0 {
        bail!("coverage values must be >= 0");
    }
    if args.scale < 0.0 {
        bail!("--scale must be >= 0");
    }
    if !(0.0..=1.0).contains(&args.error_rate) {
        bail!("--error-rate must be in [0, 1]");
    }
    if !(0.0..=1.0).contains(&args.exon_exclusive_rate) {
        bail!("--exon-exclusive-rate must be in [0, 1]");
    }
    Ok(())
}

/// Loads annotation/reference data, simulates reads, and writes all formal outputs.
///
/// The simulator is run through the dedicated `ciri-simulator` binary in normal
/// use, but this function keeps the generator callable from module tests without
/// shelling out to a subprocess.
pub fn run(args: SimulateArgs) -> Result<SimulationSummary> {
    validate_args(&args)?;
    let fastq_compressor = select_fastq_compressor()?;
    if let Some(parent) = Path::new(&args.out_prefix).parent() {
        if !parent.as_os_str().is_empty() {
            create_dir_all(parent)?;
        }
    }

    let mut rng = Lcg64::new(args.seed);
    let mut fasta = FastaReader::new();
    fasta.read_fasta(&args.ref_fasta)?;
    let annotation = read_annotation_model(
        &args.gtf,
        args.chrom.as_deref(),
        args.exon_exclusive_rate,
        &mut rng,
    )?;
    if annotation.transcripts.is_empty() {
        bail!("no exon-bearing transcripts found for the requested filter");
    }
    if annotation.linear_transcripts.is_empty() {
        bail!("no usable linear transcripts remain after circRNA-exclusive exon filtering");
    }
    write_linear_annotation(&annotation, &args.out_prefix)?;

    let mut circs = build_circs(
        &annotation.transcripts,
        &fasta.chr_tcga_map,
        &args,
        &mut rng,
    )?;
    if circs.is_empty() {
        bail!("no circRNAs could be built; try lowering --read-len or --circ-count");
    }

    let r1_tmp_path = prefixed_path(&args.out_prefix, "_1.fq.tmp");
    let r2_tmp_path = prefixed_path(&args.out_prefix, "_2.fq.tmp");
    let r1_gz_path = prefixed_path(&args.out_prefix, "_1.fq.gz");
    let r2_gz_path = prefixed_path(&args.out_prefix, "_2.fq.gz");
    let mut r1_writer = BufWriter::new(File::create(&r1_tmp_path)?);
    let mut r2_writer = BufWriter::new(File::create(&r2_tmp_path)?);
    let mut read_truth =
        BufWriter::new(File::create(prefixed_path(&args.out_prefix, ".reads.tsv"))?);
    writeln!(
        read_truth,
        "read_id\tcirc_id\tchrom\tstart\tend\tstrand\tisoform_id\tis_circular\tis_bsj\tr1_segments\tr1_is_bsj\tr2_segments\tr2_is_bsj\ttype"
    )?;

    let mut next_read_index = 1usize;
    write_circular_reads(
        &mut circs,
        &args,
        &mut rng,
        &mut r1_writer,
        &mut r2_writer,
        &mut read_truth,
        &mut next_read_index,
    )?;
    let sampled_linear_coverage = sample_coverage(args.linear_coverage, args.scale, &mut rng);
    write_linear_reads(
        &annotation.linear_transcripts,
        &fasta.chr_tcga_map,
        &args,
        sampled_linear_coverage,
        &mut rng,
        &mut r1_writer,
        &mut r2_writer,
        &mut read_truth,
        &mut next_read_index,
    )?;
    r1_writer.flush()?;
    r2_writer.flush()?;
    read_truth.flush()?;
    drop(r1_writer);
    drop(r2_writer);
    drop(read_truth);
    compress_fastq_to_gz(&r1_tmp_path, &r1_gz_path, fastq_compressor)?;
    compress_fastq_to_gz(&r2_tmp_path, &r2_gz_path, fastq_compressor)?;
    write_isoform_truth(&circs, &args.out_prefix)?;

    let isoform_count: usize = circs.iter().map(|circ| circ.isoforms.len()).sum();
    let total_read_pairs = next_read_index - 1;
    let circ_read_pairs: usize = circs.iter().map(|circ| circ.read_pair_count).sum();
    let linear_read_pairs = total_read_pairs.saturating_sub(circ_read_pairs);
    let bsj_read_pairs: usize = circs.iter().map(|circ| circ.bsj_read_pair_count).sum();
    let bsj_reads: usize = circs.iter().map(|circ| circ.bsj_read_count).sum();

    Ok(SimulationSummary {
        circ_count: circs.len(),
        isoform_count,
        total_read_pairs,
        circ_read_pairs,
        linear_read_pairs,
        bsj_reads,
        bsj_read_pairs,
    })
}

#[cfg(test)]
mod tests;
