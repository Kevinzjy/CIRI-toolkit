//! Post-Summary circRNA read reconstruction entry points.
//!
//! This module remains intentionally separate from `Scan1 -> Scan2 -> Summary`.
//! The CIRI3 Rust path is parity-sensitive, so all circRNA internal-structure
//! work stays in a post-Summary phase that consumes the final circRNA table and
//! sidecar alignment evidence without feeding evidence back into BSJ detection.
//! The current deliverable is `<prefix>.segments`: a read-level, splice-aware
//! representation of confirmed BSJ reads that can be compared directly against
//! simulator truth before full-length path reconstruction is re-enabled on top
//! of it.

use anyhow::{anyhow, bail, Context, Result};
use noodles::bam;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::time::Instant;

use crate::annotation::Annotation;
use crate::sam_bam::{detect_format, InputFormat};
use crate::utils::{bsj_payload_start, is_bsj_mate_label, reverse_complement};

const MIN_INTRON: i32 = 70;
const MIN_EXON_LENGTH: i32 = 20;
const MAX_EXON_LENGTH: i32 = 2000;
const MAX_ISOFORM_PATHS_PER_CIRC: usize = 1024;
const Z_ALPHA: f64 = 1.6449;
const MIN_JUNCTION_SUPPORT_SEGMENT_LENGTH: i32 = 10;
const INTERNAL_SPLICE_CORRECTION_WINDOW: i32 = 4;
const PARTIAL_LOCAL_SPLICE_CORRECTION_WINDOW: i32 = 16;
const MAPQ_THRES: i32 = 5;
const MAPQ_UNI: i32 = 0;
const MAPQ_BOTH: i32 = 0;
const STRINGENCY: usize = 1;
/// Minimum selected XA anchor length accepted by the backward negative filter.
const XA_REJECT_MIN_ANCHOR_LEN: i32 = MIN_JUNCTION_SUPPORT_SEGMENT_LENGTH;
/// Minimum span reduction required before XA can reject a backward sidecar row.
const XA_REJECT_MIN_SPAN_REDUCTION: i32 = 1000;
/// Read-coordinate tolerance for treating an XA hit as the same aligned slice.
const XA_SEAMLESS_QUERY_TOLERANCE: i32 = 2;
/// Maximum exact XA alternatives evaluated per selected alignment record.
const XA_MAX_ALTERNATIVES_PER_RECORD: usize = 4;
/// CIRI3 default `-Max/--max_span`; sidecar backward rows should not exceed the
/// same circRNA spanning range when they are candidates for later graph work.
const BACKWARD_MAX_SPAN: i32 = 200000;

/// Runtime options for the post-Summary circRNA segments phase.
///
/// The fields intentionally reuse the existing CIRI-AS-style second-sweep inputs
/// so read recognition can keep its current parity-sensitive logic while the
/// user-facing output is simplified to `<prefix>.segments`.
pub struct AsConfig<'a> {
    /// Original queryname-sorted SAM/BAM used by the main CIRI run.
    pub input_path: &'a str,
    /// Final CIRI3-style circRNA result table used as CIRI-AS `-C` input.
    pub circ_path: &'a str,
    /// Final mate-level BSJ table written after Summary.
    ///
    /// This file contains both `priority=1` CIRI3-compatible evidence and
    /// `priority=0` mate-level evidence. The segments phase uses it only as
    /// post-Summary read/mate annotation; it never feeds these rows back into
    /// Summary counts or circRNA discovery.
    pub bsj_path: Option<&'a str>,
    /// Scan1/Scan2 segment-evidence sidecars produced before Summary.
    ///
    /// When present, the segments phase reads these compact alignment dumps and
    /// filters them by the final `.out` junction-read assignments. This avoids a
    /// second BAM/SAM scan for confirmed BSJ reads while keeping the sidecar out
    /// of Summary's parity-sensitive inputs.
    pub segment_evidence_paths: Vec<&'a str>,
    /// Output prefix; the current phase writes `<prefix>.segments`.
    pub out_prefix: &'a str,
    /// In-memory reference sequence map loaded by the main CLI.
    pub reference: &'a HashMap<String, String>,
    /// Optional exon-boundary annotation used to resolve ambiguous motif offsets.
    ///
    /// CIRI-AS v1.2 does not consult annotation during `index_compare`, so Perl
    /// falls back to hash iteration order when multiple motif offsets are valid.
    /// The Rust sidecar deliberately uses annotation as a deterministic
    /// biological tie-break before falling back to the lowest offset.
    pub annotation: Option<&'a Annotation>,
}

/// Compact alignment representation for CIRI-AS read-group matching.
///
/// CIRI-AS only needs a subset of SAM fields. `seq` stores the full read
/// sequence in the legacy SAM/BAM sweep and the compact sidecar clip payload for
/// auditability; the default sidecar path materializes validator-accepted local
/// clips as pseudo-alignment rows before the post-Summary chain builder runs.
#[derive(Debug, Clone)]
struct AsAlignment {
    flag: i32,
    chr: String,
    pos: i32,
    mapq: i32,
    cigar: String,
    seq: String,
    from_local_clip: bool,
    xa_alternatives: Vec<XaAlternative>,
}

/// One BWA `XA:Z` alternative alignment parsed for sidecar chain rejection.
///
/// CIRI3 parity code ignores optional tags. The segments sidecar uses XA only as
/// a conservative negative signal: if a short ambiguous supplementary block can
/// be placed by XA into a non-circular chain, the read should not be promoted to
/// `type=backward` solely because BWA materialized a different equivalent hit.
#[derive(Debug, Clone, PartialEq, Eq)]
struct XaAlternative {
    chr: String,
    strand: char,
    pos: i32,
    cigar: String,
    edit_distance: i32,
}

/// Read-level output row written to `<prefix>.segments`.
///
/// This is the formal handoff format between the CIRI-AS-style evidence capture
/// and future full-length reconstruction. The current default CLI path emits
/// confirmed `bsj` rows only; `backward` and `forward` remain supported by the
/// internal builders for later extra-scan phases.
struct SegmentRecord {
    read_id: String,
    type_name: &'static str,
    circ_id: String,
    chrom: String,
    start: String,
    end: String,
    strand: String,
    is_circular: usize,
    is_r1_bsj: usize,
    is_r2_bsj: usize,
    r1_cigar: String,
    r1_segments: String,
    r2_cigar: String,
    r2_segments: String,
}

/// One mate-level BSJ row parsed from the final `<prefix>.bsj` display file.
///
/// The final `.bsj` protocol can report both mates for one read. Keeping this
/// sidecar evidence separate from `junction_read_to_circ` lets segments mark
/// the R1/R2 mate that CIRI identified as a BSJ while preserving the Summary
/// rule that only `priority=1` rows define `.out` junction read counts.
#[derive(Debug, Clone)]
struct MateBsjEvidence {
    mate_bucket: usize,
    chr: String,
    start: i32,
    end: i32,
    strand: String,
    priority: u8,
    source_stage: String,
}

/// One genomic segment block in read order.
///
/// CIGAR parsing first walks the alignment in reference order, then flips the
/// query coordinates for reverse-strand records so downstream chain assembly can
/// reason in the same read-order protocol used by the simulator truth tables.
#[derive(Debug, Clone)]
struct SegmentBlock {
    read_start: i32,
    read_end: i32,
    ref_start: i32,
    ref_end: i32,
    from_local_clip: bool,
}

/// One alignment record with CIGAR-derived blocks in read order.
#[derive(Debug, Clone)]
struct ParsedAlignment {
    flag: i32,
    chrom: String,
    strand: char,
    mapq: i32,
    blocks: Vec<SegmentBlock>,
}

/// One candidate alignment chain for a single mate.
///
/// The chain is built from primary/supplementary alignments first, with
/// secondary alignments only allowed to replace that chain when they produce a
/// topology-compatible circRNA explanation. Keeping the chain explicit prevents
/// us from dumping raw chimeric records that do not reflect the CIRI-AS splice
/// interpretation.
#[derive(Debug, Clone)]
struct MateChain {
    chrom: String,
    order_strand: char,
    token_strand: char,
    blocks: Vec<SegmentBlock>,
    token_spans: Vec<(i32, i32)>,
    tokens: Vec<String>,
    cigar: String,
    is_bsj: bool,
    is_circular: bool,
    used_secondary: usize,
    used_supplementary: usize,
    query_coverage: i32,
}

/// One final circRNA row loaded from the CIRI3 Summary output.
///
/// The column parser is header-driven because modern CIRI3-style tables include
/// a trailing `Score` column, while CIRI-AS v1.2 expects junction read IDs to be
/// the final column. Header lookup lets the sidecar consume the current Rust
/// output directly and avoids temporary lossy conversion for normal use.
#[derive(Debug, Clone)]
struct CircRecord {
    id: String,
    chr: String,
    start: i32,
    end: i32,
    junction_reads: Vec<String>,
    junction_read_count: String,
    pcc: String,
    non_junction_reads: String,
    junction_reads_ratio: String,
    circ_type: String,
    gene_id: String,
    strand: String,
}

/// Genomic interval contributed by one alignment record of a junction read.
///
/// CIRI-AS uses these intervals twice: first to build junction-read-only
/// coverage inside a circRNA, and later to ask whether a junction read spans an
/// internal splice site by at least six bases on both sides. The read-coordinate
/// fields mirror Perl `MSID_start` so paired-end AS quantification can be added
/// without rescanning SAM/BAM.
#[derive(Debug, Clone)]
#[allow(dead_code)]
struct ReadMapping {
    chr: String,
    start: i32,
    end: i32,
    reverse: usize,
    mapq: i32,
    read_start: i32,
    read_end: i32,
}

/// Boundary support used while building candidate cirexons.
///
/// Perl stores the circRNA outer boundaries as strings (`start`/`end`) in the
/// same hashes as numeric splice support. A small enum keeps that output
/// contract explicit while avoiding string-to-number coercion in Rust.
#[derive(Debug, Clone, Copy)]
enum BoundarySupport {
    Count(usize),
    CircStart,
    CircEnd,
}

impl BoundarySupport {
    /// Numeric contribution used by CIRI-AS candidate-exon generation.
    fn numeric(self) -> usize {
        match self {
            Self::Count(count) => count,
            Self::CircStart | Self::CircEnd => 0,
        }
    }

    /// Legacy text written in the support columns of `<prefix>.list`.
    fn as_output(self) -> String {
        match self {
            Self::Count(count) => count.to_string(),
            Self::CircStart => "start".to_string(),
            Self::CircEnd => "end".to_string(),
        }
    }
}

/// One validated cirexon row ready for CIRI-AS `.list` output.
#[derive(Debug, Clone)]
struct CirexonRecord {
    circ_id: String,
    chr: String,
    circ_start: i32,
    circ_end: i32,
    strand: String,
    junction_read_count: String,
    pcc: String,
    non_junction_reads: String,
    junction_reads_ratio: String,
    circ_type: String,
    gene_id: String,
    cirexon_id: String,
    start: i32,
    end: i32,
    start_support: String,
    end_support: String,
    coverage_median: i32,
    is_icf: bool,
}

/// One reconstructed full-length circRNA isoform path.
///
/// The first implementation deliberately reports path structure rather than AS
/// event taxonomy. Each record is anchored to one Summary-confirmed BSJ and is
/// therefore allowed to become sequence output without changing the baseline
/// CIRI3 circRNA calls.
#[derive(Debug, Clone)]
struct IsoformRecord {
    circ_id: String,
    isoform_id: String,
    chr: String,
    start: i32,
    end: i32,
    strand: String,
    path_tier: String,
    bsj_supported: String,
    gene_id: String,
    exons: Vec<(i32, i32)>,
    junctions: Vec<(i32, i32)>,
    isoform_len: i32,
    bsj_seed_reads: usize,
    internal_split_reads: usize,
    boundary_clip_reads: usize,
    anomalous_pair_reads: usize,
    annotation_supported_junctions: usize,
    de_novo_supported_junctions: usize,
    score: i32,
    rank: usize,
}

/// Validation result from the CIRI-AS coverage filter.
#[derive(Debug, Clone, Copy)]
struct CoverageValidation {
    code: i32,
    median: i32,
}

/// CircRNA locus cluster used to decide whether non-BSJ reads should be checked.
///
/// Perl CIRI-AS materializes every covered base into a hash. Rust keeps compact
/// intervals and performs explicit containment tests, which preserves the same
/// boundary rule (`start - 70 .. end + 70`) without whole-genome per-base maps.
#[derive(Debug, Clone)]
struct CircCluster {
    chr: String,
    start: i32,
    end: i32,
}

/// Parsed result of the CIRI-AS `MSID` CIGAR classifier.
///
/// `kind`, `clip1`, `clip2`, and `ref_len` map directly to the four values
/// returned by the Perl routine. Invalid CIGARs keep a negative `ref_len`, which
/// lets the caller retain Perl's "skip by sentinel" behavior.
#[derive(Debug, Clone, Copy)]
struct Msid {
    kind: i32,
    clip1: i32,
    clip2: i32,
    ref_len: i32,
}

/// Candidate splice junction before motif validation and clustering.
#[derive(Debug, Clone)]
struct PositiveCandidate {
    read_id: String,
    index: usize,
    chr: String,
    site1: i32,
    site2: i32,
    adjust1: i32,
    adjust2: i32,
    strand_hint: i32,
    cigars: [Option<String>; 3],
}

/// Final non-redundant splice cluster written to `_splice.list`.
#[derive(Debug, Clone)]
struct SpliceCluster {
    chr: String,
    site1: i32,
    site2: i32,
    reads: Vec<usize>,
    cigar_counts: [usize; 3],
}

/// Summary counters reported in the CIRI-AS sidecar log.
#[derive(Debug, Default)]
struct AsStats {
    junction_reads_loaded: usize,
    junction_reads_seen: usize,
    known_candidates: usize,
    add_candidates: usize,
    motif_validated: usize,
    final_splice_clusters: usize,
    final_cirexons: usize,
    final_isoforms: usize,
    final_segments: usize,
}

/// Runs the post-Summary read-level circRNA segments phase.
///
/// The default CLI path consumes Scan1/Scan2 sidecar evidence and filters it by
/// Summary-confirmed junction reads, which avoids re-scanning BAM/SAM while
/// preserving CIRI3 `.out` parity. If no sidecar paths are supplied, the older
/// CIRI-AS-style second sweep remains available for focused development and
/// tests; full-length cirexon/path reconstruction stays dormant in this module.
pub fn run_ciri_as(config: AsConfig<'_>) -> Result<()> {
    let profile = segments_profile_enabled();
    let phase_started = profile.then(Instant::now);
    let (circ_records, junction_read_to_circ) = load_circ_records(config.circ_path)?;
    log_segments_profile(profile, "load_circ_records", phase_started);
    if circ_records.is_empty() {
        bail!("No circRNAs were loaded from {}", config.circ_path);
    }
    if junction_read_to_circ.is_empty() {
        bail!(
            "No circular junction read IDs were loaded from {}",
            config.circ_path
        );
    }

    let phase_started = profile.then(Instant::now);
    let mate_bsj_evidence = if let Some(path) = config.bsj_path {
        load_mate_bsj_evidence(path)?
    } else {
        HashMap::new()
    };
    log_segments_profile(profile, "load_mate_bsj_evidence", phase_started);
    let phase_started = profile.then(Instant::now);
    let clusters = build_circ_clusters(&circ_records);
    log_segments_profile(profile, "build_circ_clusters", phase_started);
    let phase_started = profile.then(Instant::now);
    let (segment_groups, read_len) = if config.segment_evidence_paths.is_empty() {
        let read_len = infer_read_length(config.input_path)?;
        if read_len < 40 {
            bail!(
                "CIRI-AS requires paired reads with inferred read length >= 40, got {}",
                read_len
            );
        }
        (HashMap::new(), read_len)
    } else {
        let (groups, read_len) =
            load_segment_evidence(&config.segment_evidence_paths, &junction_read_to_circ)?;
        if read_len < 40 {
            bail!(
                "segments sidecar read length must be >= 40 for CIRI-AS chain assembly, got {}",
                read_len
            );
        }
        (groups, read_len)
    };
    log_segments_profile(profile, "load_segment_evidence", phase_started);
    let phase_started = profile.then(Instant::now);
    let mut state = ScanState {
        circ_by_id: circ_records
            .iter()
            .map(|r| (r.id.clone(), r.clone()))
            .collect(),
        junction_read_to_circ,
        mate_bsj_evidence,
        reference: config.reference,
        annotation: config.annotation,
        clusters_by_chr: clusters_by_chr(clusters),
        read_len,
        candidates: Vec::new(),
        coverage: HashMap::new(),
        read_mappings: HashMap::new(),
        seen_junction_reads: HashSet::new(),
        segment_groups,
        stats: AsStats::default(),
    };
    state.stats.junction_reads_loaded = state.junction_read_to_circ.len();
    log_segments_profile(profile, "build_scan_state", phase_started);

    let phase_started = profile.then(Instant::now);
    let segments = if config.segment_evidence_paths.is_empty() {
        scan_alignment_groups(config.input_path, &mut state)?;
        validate_splice_motifs(&mut state.candidates, config.reference, config.annotation)?;
        state.stats.motif_validated = state.candidates.len();
        state.stats.junction_reads_seen = state.seen_junction_reads.len();
        build_segment_records(&state)?
    } else {
        let phase_started = profile.then(Instant::now);
        scan_backward_alignment_groups(config.input_path, &mut state)?;
        log_segments_profile(profile, "scan_backward_alignment_groups", phase_started);
        let phase_started = profile.then(Instant::now);
        validate_splice_motifs(&mut state.candidates, config.reference, config.annotation)?;
        log_segments_profile(profile, "validate_backward_splice_motifs", phase_started);
        state.stats.motif_validated = state.candidates.len();
        state.stats.junction_reads_seen = state.segment_groups.len();
        build_sidecar_segment_records(&state)?
    };
    log_segments_profile(profile, "build_segments", phase_started);
    state.stats.final_segments = segments.len();
    let phase_started = profile.then(Instant::now);
    write_segments(&format!("{}.segments", config.out_prefix), &segments)?;
    log_segments_profile(profile, "write_segments", phase_started);
    Ok(())
}

/// Returns whether segments-phase profiling is enabled.
///
/// This mirrors the existing `CIRI_PROFILE_SCAN1/2` convention: profiling is
/// opt-in, has no effect on normal output, and is scoped to post-Summary
/// sidecar reconstruction so CIRI3 parity logs remain stable by default.
fn segments_profile_enabled() -> bool {
    matches!(std::env::var("CIRI_PROFILE_SEGMENTS"), Ok(v) if !v.is_empty() && v != "0")
}

/// Prints one segments profiling phase when the opt-in switch is enabled.
fn log_segments_profile(enabled: bool, label: &str, started: Option<Instant>) {
    if let (true, Some(started)) = (enabled, started) {
        eprintln!(
            "[CIRI_PROFILE_SEGMENTS] {label}: {:.3}s",
            started.elapsed().as_secs_f64()
        );
    }
}

/// Mutable state accumulated while scanning the original alignment file.
///
/// Perl stores these as package globals. Rust keeps them in one struct so later
/// BAM parallelization can shard the scan and merge candidate/state fragments
/// without changing the CIRI-AS decision routines themselves.
struct ScanState<'a> {
    circ_by_id: HashMap<String, CircRecord>,
    junction_read_to_circ: HashMap<String, String>,
    mate_bsj_evidence: HashMap<String, Vec<MateBsjEvidence>>,
    reference: &'a HashMap<String, String>,
    annotation: Option<&'a Annotation>,
    clusters_by_chr: HashMap<String, Vec<CircCluster>>,
    read_len: i32,
    candidates: Vec<PositiveCandidate>,
    coverage: HashMap<String, HashMap<i32, u32>>,
    read_mappings: HashMap<String, [Vec<ReadMapping>; 2]>,
    seen_junction_reads: HashSet<String>,
    segment_groups: HashMap<String, Vec<AsAlignment>>,
    stats: AsStats,
}

/// Reference and annotation inputs used only when formatting read-level chains.
///
/// Keeping this context out of Scan1/Scan2 preserves CIRI3 parity: boundary
/// correction changes only the final `<prefix>.segments` representation after a
/// read has already been assigned to a confirmed circRNA.
struct SegmentCorrectionContext<'a> {
    reference: &'a HashMap<String, String>,
    annotation: Option<&'a Annotation>,
    junction_support: Option<&'a JunctionSupportMap>,
}

type JunctionSupportKey = (i32, i32, char);
type JunctionSupportMap = HashMap<String, HashMap<JunctionSupportKey, usize>>;
type JunctionSupportIndex = HashMap<String, HashMap<char, HashMap<i32, Vec<i32>>>>;

/// Loads the CIRI3 Summary table and maps junction reads back to circRNAs.
///
/// CIRI-AS v1.2 used positional columns. This parser is stricter and safer by
/// using header names, while preserving the Perl overwrite behavior if one read
/// ID appears under more than one circRNA: the later row wins.
fn load_circ_records(path: &str) -> Result<(Vec<CircRecord>, HashMap<String, String>)> {
    let file = File::open(path).with_context(|| format!("open circ table {}", path))?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .transpose()?
        .ok_or_else(|| anyhow!("empty circ table: {}", path))?;
    let headers: Vec<&str> = header.split('\t').collect();
    let idx = |name: &str| -> Result<usize> {
        headers
            .iter()
            .position(|h| *h == name)
            .ok_or_else(|| anyhow!("missing `{}` column in {}", name, path))
    };
    let id_idx = idx("circRNA_ID")?;
    let chr_idx = idx("chr")?;
    let start_idx = idx("circRNA_start")?;
    let end_idx = idx("circRNA_end")?;
    let junc_count_idx = idx("#junction_reads")?;
    let pcc_idx = idx("SM_MS_SMS")?;
    let non_junc_idx = idx("#non_junction_reads")?;
    let ratio_idx = idx("junction_reads_ratio")?;
    let type_idx = idx("circRNA_type")?;
    let gene_idx = idx("gene_id")?;
    let strand_idx = idx("strand")?;
    let read_ids_idx = idx("junction_reads_ID")?;

    let mut records = Vec::new();
    let mut read_to_circ = HashMap::new();
    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<String> = line.split('\t').map(str::to_string).collect();
        let get = |i: usize| -> Result<&str> {
            cols.get(i)
                .map(String::as_str)
                .ok_or_else(|| anyhow!("malformed circ row in {}: {}", path, line))
        };
        let id = get(id_idx)?.to_string();
        let chr = get(chr_idx)?.to_string();
        let start = get(start_idx)?.parse::<i32>()?;
        let end = get(end_idx)?.parse::<i32>()?;
        let junction_reads: Vec<String> = get(read_ids_idx)?
            .split(',')
            .filter(|s| !s.is_empty())
            .map(str::to_string)
            .collect();
        for read in &junction_reads {
            read_to_circ.insert(read.clone(), id.clone());
        }
        records.push(CircRecord {
            id,
            chr,
            start,
            end,
            junction_reads,
            junction_read_count: get(junc_count_idx)?.to_string(),
            pcc: get(pcc_idx)?.to_string(),
            non_junction_reads: get(non_junc_idx)?.to_string(),
            junction_reads_ratio: get(ratio_idx)?.to_string(),
            circ_type: get(type_idx)?.to_string(),
            gene_id: get(gene_idx)?.to_string(),
            strand: get(strand_idx)?.to_string(),
        });
    }
    Ok((records, read_to_circ))
}

/// Loads post-Summary mate-level BSJ evidence from the final `.bsj` display file.
///
/// New `.bsj` rows start with `read_id, mate_label, priority` before the legacy
/// CIRI payload and end with `source_stage`. Historical rows without an explicit
/// mate label are ignored here because `<prefix>.segments` needs R1/R2-specific
/// evidence; Summary compatibility remains handled by the core pipeline.
fn load_mate_bsj_evidence(path: &str) -> Result<HashMap<String, Vec<MateBsjEvidence>>> {
    let file = File::open(path).with_context(|| format!("open BSJ evidence {}", path))?;
    let reader = BufReader::new(file);
    let mut out: HashMap<String, Vec<MateBsjEvidence>> = HashMap::new();
    for line_res in reader.lines() {
        let line = line_res?;
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        let payload_start = bsj_payload_start(&parts);
        if payload_start != 3 || parts.len() < payload_start + 9 {
            continue;
        }
        let Some(read_id) = parts.first().copied() else {
            continue;
        };
        let mate_label = parts[1];
        if !is_bsj_mate_label(mate_label) {
            continue;
        }
        let priority = parts[2].parse::<u8>().unwrap_or(1);
        let start = parts[payload_start + 3]
            .parse::<i32>()
            .with_context(|| format!("parse BSJ start from {}", line))?;
        let end = parts[payload_start + 4]
            .parse::<i32>()
            .with_context(|| format!("parse BSJ end from {}", line))?;
        out.entry(read_id.to_string())
            .or_default()
            .push(MateBsjEvidence {
                mate_bucket: if mate_label == "R1" { 0 } else { 1 },
                chr: parts[payload_start + 2].to_string(),
                start,
                end,
                strand: parts[payload_start + 5].to_string(),
                priority,
                source_stage: parts.last().copied().unwrap_or("unknown").to_string(),
            });
    }
    for rows in out.values_mut() {
        rows.sort_by(|a, b| {
            a.mate_bucket
                .cmp(&b.mate_bucket)
                .then_with(|| a.priority.cmp(&b.priority))
                .then_with(|| a.source_stage.cmp(&b.source_stage))
                .then_with(|| a.chr.cmp(&b.chr))
                .then_with(|| a.start.cmp(&b.start))
                .then_with(|| a.end.cmp(&b.end))
        });
    }
    Ok(out)
}

/// Loads Scan1/Scan2 mapper-block evidence for Summary-confirmed BSJ reads.
///
/// The sidecar rows are deliberately compact:
/// `read_id, stage, mate, flag, chrom, pos, mapq, cigar, read_len, clips`. Only
/// `read_id` and the mapper fields are mandatory; `clips` stores compact
/// `L:<seq>,R:<seq>` soft-clip subsequences for traceability, while accepted
/// local clip placements arrive as normal-looking `scan*_local` pseudo rows.
/// Duplicate mapper rows can appear when Scan1 and Scan2 both capture the same
/// read group, so rows are de-duplicated before chain selection.
fn load_segment_evidence(
    paths: &[&str],
    confirmed_reads: &HashMap<String, String>,
) -> Result<(HashMap<String, Vec<AsAlignment>>, i32)> {
    let mut groups: HashMap<String, Vec<AsAlignment>> = HashMap::new();
    let mut seen = HashSet::new();
    let mut read_len = 0_i32;
    for path in paths {
        let file = File::open(path).with_context(|| format!("open segment evidence {}", path))?;
        let reader = BufReader::new(file);
        for line_res in reader.lines() {
            let line = line_res?;
            if line.trim().is_empty() {
                continue;
            }
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 9 {
                continue;
            }
            let read_id = parts[0];
            if !confirmed_reads.contains_key(read_id) {
                continue;
            }
            let flag = parts[3]
                .parse::<i32>()
                .with_context(|| format!("parse segment evidence flag from {}", line))?;
            let chrom = parts[4];
            let pos = parts[5]
                .parse::<i32>()
                .with_context(|| format!("parse segment evidence position from {}", line))?;
            let mapq = parts[6]
                .parse::<i32>()
                .with_context(|| format!("parse segment evidence MAPQ from {}", line))?;
            let cigar = parts[7];
            let row_read_len = parts[8]
                .parse::<i32>()
                .with_context(|| format!("parse segment evidence read length from {}", line))?;
            read_len = read_len.max(row_read_len);
            let seq = parts.get(9).copied().unwrap_or("*");
            let key = format!("{read_id}\t{flag}\t{chrom}\t{pos}\t{mapq}\t{cigar}\t{seq}");
            if !seen.insert(key) {
                continue;
            }
            groups
                .entry(read_id.to_string())
                .or_default()
                .push(AsAlignment {
                    flag,
                    chr: chrom.to_string(),
                    pos,
                    mapq,
                    cigar: cigar.to_string(),
                    seq: if seq == "*" {
                        String::new()
                    } else {
                        seq.to_string()
                    },
                    from_local_clip: parts[1].contains("_local"),
                    xa_alternatives: Vec::new(),
                });
        }
    }
    for records in groups.values_mut() {
        records.sort_by(|a, b| {
            mate_bucket(a.flag)
                .cmp(&mate_bucket(b.flag))
                .then_with(|| as_alignment_order_key(a).cmp(&as_alignment_order_key(b)))
                .then_with(|| a.chr.cmp(&b.chr))
                .then_with(|| a.pos.cmp(&b.pos))
                .then_with(|| a.cigar.cmp(&b.cigar))
        });
    }
    Ok((groups, read_len))
}

/// Builds CIRI-AS circRNA locus clusters.
///
/// The merge condition intentionally uses `< cluster_end + min_intron * 2 + 1`,
/// matching the Perl expression instead of a more conventional interval-overlap
/// predicate. That off-by-one-looking boundary is part of the legacy contract.
fn build_circ_clusters(records: &[CircRecord]) -> Vec<CircCluster> {
    let mut by_chr: BTreeMap<&str, Vec<&CircRecord>> = BTreeMap::new();
    for record in records {
        by_chr.entry(record.chr.as_str()).or_default().push(record);
    }

    let mut clusters = Vec::new();
    for (chr, mut chr_records) in by_chr {
        chr_records.sort_by_key(|r| (r.start, r.end));
        let Some(first) = chr_records.first() else {
            continue;
        };
        let mut cluster_start = first.start;
        let mut cluster_end = first.end;
        for record in chr_records.into_iter().skip(1) {
            if record.start < cluster_end + MIN_INTRON * 2 + 1 {
                cluster_end = cluster_end.max(record.end);
            } else {
                clusters.push(CircCluster {
                    chr: chr.to_string(),
                    start: cluster_start,
                    end: cluster_end,
                });
                cluster_start = record.start;
                cluster_end = record.end;
            }
        }
        clusters.push(CircCluster {
            chr: chr.to_string(),
            start: cluster_start,
            end: cluster_end,
        });
    }
    clusters
}

/// Indexes compact circ clusters by chromosome for fast overlap checks.
fn clusters_by_chr(clusters: Vec<CircCluster>) -> HashMap<String, Vec<CircCluster>> {
    let mut map: HashMap<String, Vec<CircCluster>> = HashMap::new();
    for cluster in clusters {
        map.entry(cluster.chr.clone()).or_default().push(cluster);
    }
    for clusters in map.values_mut() {
        clusters.sort_by_key(|c| (c.start, c.end));
    }
    map
}

/// Infers the global read length using the CIRI-AS sampling rule.
///
/// Perl samples the first 100 unique read IDs for each mate flag and requires
/// both mate categories to agree on one length. Counting unique IDs is important:
/// supplementary hard-clipped alignments can have short SEQ fields and must not
/// bias the global `read_length` used by `MSID`.
fn infer_read_length(path: &str) -> Result<i32> {
    let mut seen: [HashSet<String>; 2] = [HashSet::new(), HashSet::new()];
    let mut length_types: [HashMap<usize, usize>; 2] = [HashMap::new(), HashMap::new()];
    let mut fallback = None;
    for_group_records(path, |read_id, records| {
        for record in records {
            if record.seq.is_empty() || record.seq == "*" {
                continue;
            }
            fallback.get_or_insert(record.seq.len() as i32);
            let idx = if record.flag & 0x40 != 0 { 1 } else { 0 };
            if seen[idx].insert(read_id.to_string()) {
                if seen[idx].len() > 100 {
                    seen[idx].remove(read_id);
                    return Ok(true);
                }
                *length_types[idx].entry(record.seq.len()).or_insert(0) += 1;
            }
        }
        Ok(false)
    })?;
    if seen[0].len() >= 100
        && seen[1].len() >= 100
        && length_types[0].len() == 1
        && length_types[1].len() == 1
    {
        let len0 = *length_types[0].keys().next().unwrap();
        let len1 = *length_types[1].keys().next().unwrap();
        if len0 == len1 {
            return Ok(len0 as i32);
        }
    }
    fallback.ok_or_else(|| anyhow!("could not infer read length from {}", path))
}

/// Scans queryname-sorted SAM/BAM groups and records CIRI-AS splice candidates.
fn scan_alignment_groups(path: &str, state: &mut ScanState) -> Result<()> {
    for_group_records(path, |read_id, records| {
        process_group(read_id, records, state)?;
        Ok(false)
    })
}

/// Scans only non-BSJ read groups for backward/circular candidates.
///
/// The default segments path already receives Summary-confirmed BSJ mapper
/// blocks from Scan1/Scan2 sidecars. This pass deliberately skips those read IDs
/// so it can supplement `<prefix>.segments` with `type=backward` rows without
/// overwriting the richer BSJ sidecar evidence or changing Summary parity.
fn scan_backward_alignment_groups(path: &str, state: &mut ScanState) -> Result<()> {
    for_group_records(path, |read_id, records| {
        process_backward_group(read_id, records, state)?;
        Ok(false)
    })
}

/// Iterates queryname-sorted alignment groups from either SAM or BAM.
///
/// Both formats feed the same lightweight `AsAlignment` records so the CIRI-AS
/// decision code remains format-agnostic. BAM support is sequential in this first
/// parity stage; the state boundary above is deliberately ready for sharded
/// merging once splice output is aligned to the Perl reference.
fn for_group_records<F>(path: &str, mut on_group: F) -> Result<()>
where
    F: FnMut(&str, &[AsAlignment]) -> Result<bool>,
{
    match detect_format(path)? {
        InputFormat::Sam => scan_sam_groups(path, &mut on_group),
        InputFormat::Bam => scan_bam_groups(path, &mut on_group),
    }
}

/// Iterates queryname-sorted SAM groups.
fn scan_sam_groups<F>(path: &str, on_group: &mut F) -> Result<()>
where
    F: FnMut(&str, &[AsAlignment]) -> Result<bool>,
{
    let file = File::open(path).with_context(|| format!("open SAM {}", path))?;
    let mut current_id = String::new();
    let mut group = Vec::with_capacity(8);
    for line in BufReader::new(file).lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('@') {
            continue;
        }
        let mut cols = line.split('\t');
        let read_id = cols.next().unwrap_or_default();
        if !current_id.is_empty() && read_id != current_id {
            if on_group(&current_id, &group)? {
                return Ok(());
            }
            group.clear();
        }
        current_id.clear();
        current_id.push_str(read_id);
        let flag = cols.next().unwrap_or("0").parse::<i32>().unwrap_or(0);
        let chr = cols.next().unwrap_or("*").to_string();
        let pos = cols.next().unwrap_or("0").parse::<i32>().unwrap_or(0);
        let mapq = cols.next().unwrap_or("0").parse::<i32>().unwrap_or(0);
        let cigar = cols.next().unwrap_or("*").to_string();
        cols.next();
        cols.next();
        cols.next();
        let seq = cols.next().unwrap_or("*").to_string();
        cols.next();
        let xa_alternatives = parse_xa_from_sam_fields(cols);
        group.push(AsAlignment {
            flag,
            chr,
            pos,
            mapq,
            cigar,
            seq,
            from_local_clip: false,
            xa_alternatives,
        });
    }
    if !group.is_empty() {
        on_group(&current_id, &group)?;
    }
    Ok(())
}

/// Iterates queryname-sorted BAM groups through noodles.
fn scan_bam_groups<F>(path: &str, on_group: &mut F) -> Result<()>
where
    F: FnMut(&str, &[AsAlignment]) -> Result<bool>,
{
    use noodles::sam::alignment::Record as _;

    let file = File::open(path).with_context(|| format!("open BAM {}", path))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header()?;
    let mut record = bam::Record::default();
    let mut current_id: Vec<u8> = Vec::new();
    let mut group = Vec::with_capacity(8);
    let mut cigar_buf = String::with_capacity(64);
    let mut seq_buf = String::with_capacity(256);
    while reader.read_record(&mut record)? != 0 {
        let read_id = record
            .name()
            .ok_or_else(|| anyhow!("Missing read name in {}", path))?;
        if !current_id.is_empty() && read_id != current_id.as_slice() {
            let id = String::from_utf8_lossy(&current_id);
            if on_group(&id, &group)? {
                return Ok(());
            }
            group.clear();
        }
        current_id.clear();
        current_id.extend_from_slice(read_id);
        let chr = match record.reference_sequence(&header) {
            Some(Ok((name, _))) => String::from_utf8_lossy(name).to_string(),
            _ => "*".to_string(),
        };
        let flag = i32::from(u16::from(record.flags()));
        let pos = record
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
        let xa_alternatives = parse_xa_from_bam_record(&record)?;
        group.push(AsAlignment {
            flag,
            chr,
            pos,
            mapq,
            cigar: cigar_buf.clone(),
            seq: seq_buf.clone(),
            from_local_clip: false,
            xa_alternatives,
        });
    }
    if !group.is_empty() {
        let id = String::from_utf8_lossy(&current_id);
        on_group(&id, &group)?;
    }
    Ok(())
}

/// Parses `XA:Z` alternatives from raw SAM optional columns.
///
/// The parser intentionally accepts only BWA's four-field entries
/// `chr,[+-]pos,cigar,NM;`. Malformed entries are ignored because XA is a
/// sidecar-only negative signal and must never make alignment parsing fail.
fn parse_xa_from_sam_fields<'a>(fields: impl Iterator<Item = &'a str>) -> Vec<XaAlternative> {
    for field in fields {
        if let Some(raw) = field.strip_prefix("XA:Z:") {
            return parse_xa_tag(raw);
        }
    }
    Vec::new()
}

/// Parses `XA:Z` alternatives from one BAM record.
fn parse_xa_from_bam_record(record: &bam::Record) -> Result<Vec<XaAlternative>> {
    let tag = Tag::new(b'X', b'A');
    let data = record.data();
    let Some(value) = data.get(&tag).transpose()? else {
        return Ok(Vec::new());
    };
    let Value::String(raw) = value else {
        return Ok(Vec::new());
    };
    let Ok(raw) = std::str::from_utf8(raw.as_ref()) else {
        return Ok(Vec::new());
    };
    Ok(parse_xa_tag(raw))
}

/// Parses one BWA `XA:Z` payload into alternative alignment records.
fn parse_xa_tag(raw: &str) -> Vec<XaAlternative> {
    raw.split(';')
        .filter_map(|entry| {
            if entry.is_empty() {
                return None;
            }
            let mut parts = entry.split(',');
            let chr = parts.next()?.to_string();
            let signed_pos = parts.next()?;
            let mut chars = signed_pos.chars();
            let strand = chars.next()?;
            if !matches!(strand, '+' | '-') {
                return None;
            }
            let pos = chars.as_str().parse::<i32>().ok()?.abs();
            let cigar = parts.next()?.to_string();
            let edit_distance = parts.next()?.parse::<i32>().ok()?;
            Some(XaAlternative {
                chr,
                strand,
                pos,
                cigar,
                edit_distance,
            })
        })
        .collect()
}

/// Processes one queryname group through CIRI-AS known/non-known read paths.
fn process_group(read_id: &str, records: &[AsAlignment], state: &mut ScanState) -> Result<()> {
    if records.is_empty() {
        return Ok(());
    }
    record_cluster_coverage(records, state);
    if state.junction_read_to_circ.contains_key(read_id) {
        state.seen_junction_reads.insert(read_id.to_string());
        state
            .segment_groups
            .insert(read_id.to_string(), records.to_vec());
        record_mapping_detail(read_id, records, state);
        mapping_check1(read_id, records, true, state)?;
    } else if records.len() > 2 && overlaps_any_circ_cluster(records, state) {
        state
            .segment_groups
            .insert(read_id.to_string(), records.to_vec());
        mapping_check1(read_id, records, false, state)?;
    }
    Ok(())
}

/// Processes one non-BSJ read group for backward candidate discovery.
///
/// This is the sidecar-mode counterpart of the non-BSJ branch in
/// `process_group`. Confirmed BSJ reads are skipped because their final chains
/// must come from `.segments1/2`; non-BSJ groups are only retained when they
/// overlap a confirmed circRNA cluster and produce CIRI-AS-style circular
/// topology candidates.
fn process_backward_group(
    read_id: &str,
    records: &[AsAlignment],
    state: &mut ScanState,
) -> Result<()> {
    if records.is_empty() || state.junction_read_to_circ.contains_key(read_id) {
        return Ok(());
    }
    if records.len() > 2 && overlaps_any_circ_cluster(records, state) {
        state
            .segment_groups
            .entry(read_id.to_string())
            .or_insert_with(|| records.to_vec());
        mapping_check1(read_id, records, false, state)?;
    }
    Ok(())
}

/// Adds all alignment coverage that overlaps CIRI-AS circRNA locus clusters.
///
/// This mirrors the Perl scan-time `%coverage` population. Only the portion
/// anchored by a cluster boundary is counted, including the ±`MIN_INTRON` flank
/// that CIRI-AS uses when judging candidate exon edges against local background.
fn record_cluster_coverage(records: &[AsAlignment], state: &mut ScanState) {
    for record in records {
        if record.mapq < MAPQ_THRES {
            continue;
        }
        let msid = msid(&record.cigar, state.read_len);
        if msid.ref_len < 0 {
            continue;
        }
        let map_end = record.pos + msid.ref_len - 1;
        let Some((loci_start, loci_end)) = coverage_locus(record, map_end, state) else {
            continue;
        };
        if loci_start > loci_end {
            continue;
        }
        let chr_cov = state.coverage.entry(record.chr.clone()).or_default();
        for pos in loci_start..=loci_end {
            *chr_cov.entry(pos).or_insert(0) += 1;
        }
    }
}

/// Returns the clipped locus range counted for one alignment.
///
/// Perl tests the alignment start first and the alignment end second against a
/// per-base cluster range. Preserving that order matters for reads that touch
/// both sides of a cluster because the first branch chooses a different clipped
/// interval than the second branch.
fn coverage_locus(record: &AsAlignment, map_end: i32, state: &ScanState) -> Option<(i32, i32)> {
    let clusters = state.clusters_by_chr.get(&record.chr)?;
    if let Some(cluster) = clusters
        .iter()
        .find(|cluster| point_in_cluster_range(record.pos, cluster))
    {
        return Some((record.pos, map_end.min(cluster.end)));
    }
    if let Some(cluster) = clusters
        .iter()
        .find(|cluster| point_in_cluster_range(map_end, cluster))
    {
        return Some((cluster.start, map_end));
    }
    None
}

/// Tests membership in Perl's materialized `circ_cluster_range` hash.
fn point_in_cluster_range(pos: i32, cluster: &CircCluster) -> bool {
    pos >= cluster.start - MIN_INTRON && pos <= cluster.end + MIN_INTRON
}

/// Records junction-read mapping intervals for cirexon validation.
///
/// CIRI-AS only stores this detail for reads that were listed in the circRNA
/// Summary table. Non-BSJ reads may support `_splice.list`, but they are not
/// allowed to provide direct start/end BSJ support for `.list` cirexons.
fn record_mapping_detail(read_id: &str, records: &[AsAlignment], state: &mut ScanState) {
    let mut by_reverse: [Vec<ReadMapping>; 2] = [Vec::new(), Vec::new()];
    for record in records {
        let Some(msid) = msid_start(&record.cigar, state.read_len) else {
            continue;
        };
        if msid.ref_len <= 0 {
            continue;
        }
        let reverse = if record.flag & 0x10 != 0 { 1 } else { 0 };
        by_reverse[reverse].push(ReadMapping {
            chr: record.chr.clone(),
            start: record.pos,
            end: record.pos + msid.ref_len - 1,
            reverse,
            mapq: record.mapq,
            read_start: msid.clip1,
            read_end: msid.clip1 + msid.clip2 - 1,
        });
    }
    state.read_mappings.insert(read_id.to_string(), by_reverse);
}

/// True if any alignment start/end falls inside a clustered circRNA locus.
fn overlaps_any_circ_cluster(records: &[AsAlignment], state: &ScanState) -> bool {
    for record in records {
        if record.mapq < MAPQ_THRES {
            continue;
        }
        let msid = msid(&record.cigar, state.read_len);
        if msid.ref_len < 0 {
            continue;
        }
        let map_end = record.pos + msid.ref_len - 1;
        if let Some(clusters) = state.clusters_by_chr.get(&record.chr) {
            for cluster in clusters {
                let lo = cluster.start - MIN_INTRON;
                let hi = cluster.end + MIN_INTRON;
                if (record.pos >= lo && record.pos <= hi) || (map_end >= lo && map_end <= hi) {
                    return true;
                }
                if cluster.start > map_end + MIN_INTRON {
                    break;
                }
            }
        }
    }
    false
}

/// Splits one read group by reverse-strand flag and runs CIRI-AS pair matching.
///
/// The bucket is SAM flag bit 0x10, exactly Perl `ten2b(flag, 5)`. It is easy to
/// confuse with first/second mate bits because CIRI3's Scan1 grouping is mate
/// oriented; CIRI-AS splice detection is alignment-strand oriented.
fn mapping_check1(
    read_id: &str,
    records: &[AsAlignment],
    from_bsj: bool,
    state: &mut ScanState,
) -> Result<()> {
    let mut by_reverse: [Vec<&AsAlignment>; 2] = [Vec::new(), Vec::new()];
    for record in records {
        let idx = if record.flag & 0x10 != 0 { 1 } else { 0 };
        by_reverse[idx].push(record);
    }
    mapping_check2(read_id, &by_reverse[1], from_bsj, state)?;
    mapping_check2(read_id, &by_reverse[0], from_bsj, state)?;
    Ok(())
}

/// Mirrors CIRI-AS `mapping_check2` and `mapping_check2_add`.
///
/// The two Perl routines differ only by circRNA boundary checks and counters.
/// `from_bsj=true` enforces the original circRNA interval for known junction
/// reads; non-BSJ add-on reads are only pre-gated by clustered circ overlap.
fn mapping_check2(
    read_id: &str,
    records: &[&AsAlignment],
    from_bsj: bool,
    state: &mut ScanState,
) -> Result<()> {
    if records.len() < 2 {
        return Ok(());
    }
    let circ = if from_bsj {
        let circ_id = state
            .junction_read_to_circ
            .get(read_id)
            .ok_or_else(|| anyhow!("missing circ for junction read {}", read_id))?;
        Some(
            state
                .circ_by_id
                .get(circ_id)
                .ok_or_else(|| anyhow!("missing circ record {}", circ_id))?,
        )
    } else {
        None
    };

    let mut matches = Vec::with_capacity(records.len());
    for record in records {
        matches.push(msid(&record.cigar, state.read_len));
    }

    let initial_len = state.candidates.len();
    let mut candidate_chr: Option<String> = None;
    let mut na_tag = false;
    'read: for i in 0..records.len() - 1 {
        for j in i + 1..records.len() {
            let ri = records[i];
            let rj = records[j];
            if ri.chr != rj.chr {
                continue;
            }
            if let Some(circ) = circ {
                if ri.chr != circ.chr {
                    continue;
                }
            }
            if reverse_bit(ri.flag) != reverse_bit(rj.flag) {
                continue;
            }
            if let Some(existing_chr) = &candidate_chr {
                if existing_chr != &ri.chr {
                    na_tag = true;
                    break 'read;
                }
            }
            if let Some(candidate) = candidate_from_pair(
                read_id,
                ri,
                rj,
                matches[i],
                matches[j],
                i,
                j,
                circ,
                state.read_len,
            ) {
                candidate_chr.get_or_insert_with(|| candidate.chr.clone());
                if from_bsj {
                    state.stats.known_candidates += 1;
                } else {
                    state.stats.add_candidates += 1;
                }
                state.candidates.push(candidate);
            }
        }
    }

    if na_tag {
        state.candidates.truncate(initial_len);
    }
    Ok(())
}

/// Builds one candidate from a pair of CIRI-AS classified alignments.
fn candidate_from_pair(
    read_id: &str,
    ri: &AsAlignment,
    rj: &AsAlignment,
    mi: Msid,
    mj: Msid,
    i: usize,
    j: usize,
    circ: Option<&CircRecord>,
    read_len: i32,
) -> Option<PositiveCandidate> {
    let product = mi.kind * mj.kind;
    if product == -1 {
        let cir_scale = mi.kind * (ri.pos + mi.clip2) + mj.kind * (rj.pos + mj.clip2);
        if (mi.clip1 - mj.clip1).abs() <= 6
            && cir_scale < 0
            && boundary_ok(ri, mi, circ)
            && boundary_ok(rj, mj, circ)
            && ri.mapq >= MAPQ_UNI
            && rj.mapq >= MAPQ_UNI
            && ri.mapq + rj.mapq >= MAPQ_BOTH
        {
            let mut pos_com = [ri.pos + mi.clip2, rj.pos + mj.clip2];
            pos_com.sort_unstable();
            let end_adjustment1 = div_trunc(mi.clip1 * mi.kind + mj.clip1 * mj.kind, 2);
            let end_adjustment2 = mi.clip1 * mi.kind + mj.clip1 * mj.kind - end_adjustment1;
            let (x_record, y_record) = if mi.kind.cmp(&mj.kind) == Ordering::Greater {
                (rj, ri)
            } else {
                (ri, rj)
            };
            return Some(PositiveCandidate {
                read_id: read_id.to_string(),
                index: 0,
                chr: ri.chr.clone(),
                site2: pos_com[0] - end_adjustment2,
                site1: pos_com[1] + end_adjustment1,
                adjust1: end_adjustment1,
                adjust2: end_adjustment2,
                strand_hint: 0,
                cigars: [
                    Some(x_record.cigar.clone()),
                    Some(y_record.cigar.clone()),
                    None,
                ],
            });
        }
    } else if product.abs() == 10 {
        let (x, y, mx, my, rx, ry) = if mi.kind <= mj.kind {
            (i, j, mi, mj, ri, rj)
        } else {
            (j, i, mj, mi, rj, ri)
        };
        let _ = (x, y);
        if mx.kind == -1 {
            let cir_scale = ry.pos + my.ref_len - 1 - rx.pos;
            if (read_len - my.clip2 - mx.clip1).abs() <= 6
                && cir_scale < 0
                && boundary_start_ok(ry, circ)
                && boundary_end_ok(rx, mx, circ)
                && ri.mapq >= MAPQ_UNI
                && rj.mapq >= MAPQ_UNI
                && ri.mapq + rj.mapq >= MAPQ_BOTH
            {
                let end_adjustment1 = div_trunc(my.clip1 + my.ref_len - mx.clip1, 2);
                let end_adjustment2 = my.clip1 + my.ref_len - mx.clip1 - end_adjustment1;
                return Some(PositiveCandidate {
                    read_id: read_id.to_string(),
                    index: 0,
                    chr: rx.chr.clone(),
                    site2: ry.pos + my.ref_len - 1 - end_adjustment2,
                    site1: rx.pos + end_adjustment1,
                    adjust1: end_adjustment1,
                    adjust2: end_adjustment2,
                    strand_hint: 0,
                    cigars: [Some(rx.cigar.clone()), None, Some(ry.cigar.clone())],
                });
            }
        } else {
            let cir_scale = rx.pos + mx.ref_len - 1 - ry.pos;
            if (mx.clip1 - my.clip1).abs() <= 6
                && cir_scale < 0
                && boundary_start_ok(rx, circ)
                && boundary_end_ok(ry, my, circ)
                && ri.mapq >= MAPQ_UNI
                && rj.mapq >= MAPQ_UNI
                && ri.mapq + rj.mapq >= MAPQ_BOTH
            {
                let end_adjustment1 = div_trunc(mx.clip1 - my.clip1, 2);
                let end_adjustment2 = mx.clip1 - my.clip1 - end_adjustment1;
                return Some(PositiveCandidate {
                    read_id: read_id.to_string(),
                    index: 0,
                    chr: rx.chr.clone(),
                    site2: rx.pos + mx.ref_len - 1 - end_adjustment2,
                    site1: ry.pos + end_adjustment1,
                    adjust1: end_adjustment1,
                    adjust2: end_adjustment2,
                    strand_hint: 0,
                    cigars: [None, Some(rx.cigar.clone()), Some(ry.cigar.clone())],
                });
            }
        }
    }
    None
}

/// Perl-style integer division truncating toward zero.
fn div_trunc(value: i32, divisor: i32) -> i32 {
    value / divisor
}

/// Returns SAM flag bit 0x10 as CIRI-AS `ten2b(flag, 5)`.
fn reverse_bit(flag: i32) -> i32 {
    if flag & 0x10 != 0 {
        1
    } else {
        0
    }
}

/// Checks the full circ-boundary predicate used by the product `-1` branch.
fn boundary_ok(record: &AsAlignment, msid: Msid, circ: Option<&CircRecord>) -> bool {
    circ.is_none_or(|circ| {
        record.pos >= circ.start - 6 && record.pos + msid.ref_len - 1 <= circ.end + 6
    })
}

/// Checks start-side circ-boundary predicates from CIRI-AS branch `abs(product)==10`.
fn boundary_start_ok(record: &AsAlignment, circ: Option<&CircRecord>) -> bool {
    circ.is_none_or(|circ| record.pos >= circ.start - 6)
}

/// Checks end-side circ-boundary predicates from CIRI-AS branch `abs(product)==10`.
fn boundary_end_ok(record: &AsAlignment, msid: Msid, circ: Option<&CircRecord>) -> bool {
    circ.is_none_or(|circ| record.pos + msid.ref_len - 1 <= circ.end + 6)
}

/// Mirrors the CIRI-AS `MSID` CIGAR classifier.
///
/// Hard clips are normalized to soft clips because the Perl routine performs
/// `s/H/S/g` before classifying. The returned values intentionally retain
/// legacy sentinels `-1` and `-2` for unsupported full-length and malformed CIGAR
/// patterns.
fn msid(cigar: &str, read_len: i32) -> Msid {
    if cigar == "*" || cigar.is_empty() {
        return Msid {
            kind: 0,
            clip1: 0,
            clip2: 0,
            ref_len: -1,
        };
    }
    let mut counts = Vec::new();
    let mut styles = Vec::new();
    let mut number = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            number.push(ch);
        } else {
            counts.push(number.parse::<i32>().unwrap_or(0));
            number.clear();
            styles.push(if ch == 'H' { 'S' } else { ch });
        }
    }
    if styles.is_empty() {
        return Msid {
            kind: 0,
            clip1: 0,
            clip2: 0,
            ref_len: -2,
        };
    }
    if counts.len() == 1 {
        return if cigar == format!("{}M", read_len) {
            Msid {
                kind: 0,
                clip1: 0,
                clip2: 0,
                ref_len: read_len,
            }
        } else {
            Msid {
                kind: 0,
                clip1: 0,
                clip2: 0,
                ref_len: -1,
            }
        };
    }

    match counts.len() {
        2 => match (styles[0], styles[1]) {
            ('M', 'S') => Msid {
                kind: 1,
                clip1: counts[0],
                clip2: counts[0] - 1,
                ref_len: counts[0],
            },
            ('S', 'M') => Msid {
                kind: -1,
                clip1: counts[0],
                clip2: 0,
                ref_len: counts[1],
            },
            _ => invalid_msid(-2),
        },
        3 => match (styles[0], styles[1], styles[2]) {
            ('S', _, 'S') => Msid {
                kind: 10,
                clip1: counts[0],
                clip2: counts[2],
                ref_len: counts[1],
            },
            ('M', 'D', 'M') => Msid {
                kind: 0,
                clip1: 0,
                clip2: 0,
                ref_len: read_len + counts[1],
            },
            ('M', 'I', 'M') => Msid {
                kind: 0,
                clip1: 0,
                clip2: 0,
                ref_len: read_len - counts[1],
            },
            _ => invalid_msid(-2),
        },
        _ if styles[0] == 'M' && *styles.last().unwrap() == 'S' => {
            let (m_sum, d_sum) = sum_md(&styles[..styles.len() - 1], &counts[..counts.len() - 1]);
            Msid {
                kind: 1,
                clip1: read_len - counts[counts.len() - 1],
                clip2: m_sum + d_sum - 1,
                ref_len: m_sum + d_sum,
            }
        }
        _ if styles[0] == 'S' && *styles.last().unwrap() == 'M' => {
            let (m_sum, d_sum) = sum_md(&styles[1..], &counts[1..]);
            Msid {
                kind: -1,
                clip1: counts[0],
                clip2: 0,
                ref_len: m_sum + d_sum,
            }
        }
        _ if styles[0] == 'M' && *styles.last().unwrap() == 'M' => {
            let (m_sum, d_sum) = sum_md(&styles, &counts);
            Msid {
                kind: 0,
                clip1: 0,
                clip2: 0,
                ref_len: m_sum + d_sum,
            }
        }
        _ if styles[0] == 'S' && *styles.last().unwrap() == 'S' => {
            let (m_sum, d_sum) = sum_md(&styles[1..styles.len() - 1], &counts[1..counts.len() - 1]);
            Msid {
                kind: 10,
                clip1: counts[0],
                clip2: counts[counts.len() - 1],
                ref_len: m_sum + d_sum,
            }
        }
        _ => invalid_msid(-2),
    }
}

/// Mirrors CIRI-AS `MSID_start` for read-coordinate bookkeeping.
///
/// `MSID` and `MSID_start` look similar but their second and third fields have
/// different meanings. The splice detector needs clip offsets, while cirexon
/// path/coverage reporting needs read start/end coordinates. Keeping this as a
/// separate routine prevents accidental reuse of the wrong coordinate system.
fn msid_start(cigar: &str, read_len: i32) -> Option<Msid> {
    if cigar == "*" || cigar.is_empty() {
        return None;
    }
    let mut counts = Vec::new();
    let mut styles = Vec::new();
    let mut number = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            number.push(ch);
        } else {
            counts.push(number.parse::<i32>().unwrap_or(0));
            number.clear();
            styles.push(if ch == 'H' { 'S' } else { ch });
        }
    }
    if counts.is_empty() {
        return None;
    }
    if counts.len() == 1 {
        return Some(if cigar == format!("{}M", read_len) {
            Msid {
                kind: 0,
                clip1: 1,
                clip2: read_len,
                ref_len: read_len,
            }
        } else {
            Msid {
                kind: 0,
                clip1: 0,
                clip2: 0,
                ref_len: -1,
            }
        });
    }

    let result = match counts.len() {
        2 => match (styles[0], styles[1]) {
            ('M', 'S') => Msid {
                kind: 1,
                clip1: 1,
                clip2: counts[0],
                ref_len: counts[0],
            },
            ('S', 'M') => Msid {
                kind: -1,
                clip1: counts[0] + 1,
                clip2: counts[1],
                ref_len: counts[1],
            },
            _ => return None,
        },
        3 => match (styles[0], styles[1], styles[2]) {
            ('S', _, 'S') => Msid {
                kind: 10,
                clip1: counts[0] + 1,
                clip2: counts[1],
                ref_len: counts[1],
            },
            ('M', 'D', 'M') => Msid {
                kind: 0,
                clip1: 1,
                clip2: read_len,
                ref_len: read_len + counts[1],
            },
            ('M', 'I', 'M') => Msid {
                kind: 0,
                clip1: 1,
                clip2: read_len,
                ref_len: read_len - counts[1],
            },
            _ => return None,
        },
        _ if styles[0] == 'M' && *styles.last().unwrap() == 'S' => {
            let (m_sum, d_sum, i_sum) =
                sum_mdi(&styles[..styles.len() - 1], &counts[..counts.len() - 1]);
            Msid {
                kind: 1,
                clip1: 1,
                clip2: m_sum + i_sum,
                ref_len: m_sum + d_sum,
            }
        }
        _ if styles[0] == 'S' && *styles.last().unwrap() == 'M' => {
            let (m_sum, d_sum, i_sum) = sum_mdi(&styles[1..], &counts[1..]);
            Msid {
                kind: -1,
                clip1: counts[0] + 1,
                clip2: m_sum + i_sum,
                ref_len: m_sum + d_sum,
            }
        }
        _ if styles[0] == 'M' && *styles.last().unwrap() == 'M' => {
            let (m_sum, d_sum, _) = sum_mdi(&styles, &counts);
            Msid {
                kind: 0,
                clip1: 1,
                clip2: read_len,
                ref_len: m_sum + d_sum,
            }
        }
        _ if styles[0] == 'S' && *styles.last().unwrap() == 'S' => {
            let (m_sum, d_sum, i_sum) =
                sum_mdi(&styles[1..styles.len() - 1], &counts[1..counts.len() - 1]);
            Msid {
                kind: 10,
                clip1: counts[0] + 1,
                clip2: m_sum + i_sum,
                ref_len: m_sum + d_sum,
            }
        }
        _ => return None,
    };
    Some(result)
}

/// Returns an invalid `MSID` sentinel.
fn invalid_msid(ref_len: i32) -> Msid {
    Msid {
        kind: 0,
        clip1: 0,
        clip2: 0,
        ref_len,
    }
}

/// Sums reference-consuming `M` and `D` operations for `MSID`.
fn sum_md(styles: &[char], counts: &[i32]) -> (i32, i32) {
    let mut m_sum = 0;
    let mut d_sum = 0;
    for (style, count) in styles.iter().zip(counts) {
        match style {
            'M' => m_sum += count,
            'D' => d_sum += count,
            _ => {}
        }
    }
    (m_sum, d_sum)
}

/// Sums reference/read-consuming CIGAR operations for `MSID_start`.
fn sum_mdi(styles: &[char], counts: &[i32]) -> (i32, i32, i32) {
    let mut m_sum = 0;
    let mut d_sum = 0;
    let mut i_sum = 0;
    for (style, count) in styles.iter().zip(counts) {
        match style {
            'M' => m_sum += count,
            'D' => d_sum += count,
            'I' => i_sum += count,
            _ => {}
        }
    }
    (m_sum, d_sum, i_sum)
}

/// Applies CIRI-AS splice-signal motif checking and coordinate adjustment.
fn validate_splice_motifs(
    candidates: &mut Vec<PositiveCandidate>,
    reference: &HashMap<String, String>,
    annotation: Option<&Annotation>,
) -> Result<()> {
    let mut validated = Vec::with_capacity(candidates.len());
    for mut candidate in candidates.drain(..) {
        let Some(chr_seq) = reference.get(&candidate.chr) else {
            continue;
        };
        let total_adjustment = candidate.adjust1 + candidate.adjust2;
        let (end_string1, end_string2) = if candidate.adjust2 >= 0 {
            (
                perl_substr(
                    chr_seq,
                    candidate.site1 - candidate.adjust1 - 4,
                    4 + total_adjustment,
                ),
                perl_substr(
                    chr_seq,
                    candidate.site2 - candidate.adjust1 - 1,
                    4 + total_adjustment,
                ),
            )
        } else {
            (
                perl_substr(
                    chr_seq,
                    candidate.site1 + candidate.adjust2 - 4,
                    4 - total_adjustment,
                ),
                perl_substr(
                    chr_seq,
                    candidate.site2 + candidate.adjust2 - 1,
                    4 - total_adjustment,
                ),
            )
        };
        let (hint, index) = index_compare(&end_string1, &end_string2, &candidate, annotation);
        if hint != 0 {
            let diff_adjt = if candidate.adjust2 >= 0 {
                index - 1 - candidate.adjust1
            } else {
                index - 1 + candidate.adjust2
            };
            candidate.strand_hint = hint;
            candidate.adjust1 += diff_adjt;
            candidate.adjust2 = total_adjustment - candidate.adjust1;
            candidate.site1 += diff_adjt;
            candidate.site2 += diff_adjt;
            validated.push(candidate);
        }
    }
    for (idx, candidate) in validated.iter_mut().enumerate() {
        candidate.index = idx;
    }
    *candidates = validated;
    Ok(())
}

/// Perl-compatible positive-offset substring on ASCII FASTA sequences.
fn perl_substr(seq: &str, start: i32, len: i32) -> String {
    if len <= 0 || start >= seq.len() as i32 {
        return String::new();
    }
    let start = start.max(0) as usize;
    let end = (start + len as usize).min(seq.len());
    seq[start..end].to_string()
}

/// Matches splice motifs and resolves ambiguous strand/offset choices.
///
/// CIRI-AS v1.2 uses `if AC/CT ... elsif AG/GT`, so ambiguous windows are decided
/// by branch order and Perl hash iteration. Rust keeps that behavior only as the
/// fallback. When annotation is available, both negative-strand (`AC/CT`) and
/// positive-strand (`AG/GT`) motif explanations are scored against known exon
/// boundaries so strand and offset are chosen by biological support first.
fn index_compare(
    left: &str,
    right: &str,
    candidate: &PositiveCandidate,
    annotation: Option<&Annotation>,
) -> (i32, i32) {
    let left = left.to_ascii_uppercase();
    let right = right.to_ascii_uppercase();
    let has_negative_motifs = left.contains("AC") && right.contains("CT");
    let has_positive_motifs = left.contains("AG") && right.contains("GT");
    let negative = shared_motif_indices(&left, "AC", &right, "CT");
    let positive = shared_motif_indices(&left, "AG", &right, "GT");

    if annotation.is_some() {
        if let Some(choice) = choose_annotated_motif(&negative, &positive, candidate, annotation) {
            return choice;
        }
    }

    if has_negative_motifs {
        if negative.is_empty() {
            return (0, 0);
        }
        return (1, choose_motif_index(&negative, 1, candidate, annotation));
    }
    if has_positive_motifs && !positive.is_empty() {
        return (-1, choose_motif_index(&positive, -1, candidate, annotation));
    }
    (0, 0)
}

/// Finds all shared motif offsets between two splice windows.
fn shared_motif_indices(left: &str, motif_left: &str, right: &str, motif_right: &str) -> Vec<i32> {
    let left_indices: HashSet<usize> = left.match_indices(motif_left).map(|(i, _)| i).collect();
    let mut indices: Vec<i32> = right
        .match_indices(motif_right)
        .map(|(i, _)| i)
        .filter(|i| left_indices.contains(i))
        .map(|i| i as i32)
        .collect();
    indices.sort_unstable();
    indices
}

/// Chooses a motif offset using exon-boundary support before stable fallback.
fn choose_motif_index(
    indices: &[i32],
    hint: i32,
    candidate: &PositiveCandidate,
    annotation: Option<&Annotation>,
) -> i32 {
    indices
        .iter()
        .copied()
        .max_by_key(|&index| {
            (
                annotation_offset_score(candidate, index, hint, annotation),
                std::cmp::Reverse(index),
            )
        })
        .unwrap_or(0)
}

/// Chooses between negative- and positive-strand motif explanations.
///
/// Annotation can disambiguate both offset and strand. If both strands receive
/// the same support score, this falls back to CIRI-AS v1.2 branch order
/// (`AC/CT` before `AG/GT`) so uninformative annotation does not invent a new
/// strand convention.
fn choose_annotated_motif(
    negative: &[i32],
    positive: &[i32],
    candidate: &PositiveCandidate,
    annotation: Option<&Annotation>,
) -> Option<(i32, i32)> {
    let mut best: Option<(i32, i32, i32, i32)> = None;
    for &(hint, branch_rank, indices) in &[(1, 0, negative), (-1, 1, positive)] {
        for &index in indices {
            let score = annotation_offset_score(candidate, index, hint, annotation);
            if score == 0 {
                continue;
            }
            let candidate_key = (score, -branch_rank, -index, hint);
            if best.is_none_or(|current| candidate_key > current) {
                best = Some(candidate_key);
            }
        }
    }
    if let Some((_, _, neg_index, hint)) = best {
        return Some((hint, -neg_index));
    }
    None
}

/// Scores one possible motif offset by known exon-boundary support.
///
/// CIRI-AS later treats `site1` as a candidate exon start and `site2` as a
/// candidate exon end. The tie-break follows that downstream contract: exact
/// `site1 -> exon_start` plus `site2 -> exon_end` support is strongest, with a
/// smaller bonus when the motif-implied strand agrees with annotation.
fn annotation_offset_score(
    candidate: &PositiveCandidate,
    index: i32,
    hint: i32,
    annotation: Option<&Annotation>,
) -> i32 {
    let Some(annotation) = annotation else {
        return 0;
    };
    let diff_adjt = if candidate.adjust2 >= 0 {
        index - 1 - candidate.adjust1
    } else {
        index - 1 + candidate.adjust2
    };
    let site1 = candidate.site1 + diff_adjt;
    let site2 = candidate.site2 + diff_adjt;
    let start_key = format!("{}\t{}", candidate.chr, site1);
    let end_key = format!("{}\t{}", candidate.chr, site2);
    let start_hit = annotation.chr_exon_start_map.get(&start_key);
    let end_hit = annotation.chr_exon_end_map.get(&end_key);
    let mut score = match (start_hit, end_hit) {
        (Some(start), Some(end)) => {
            let start_gene = start.split('\t').next().unwrap_or("");
            let end_gene = end.split('\t').next().unwrap_or("");
            if start_gene == end_gene {
                100
            } else {
                80
            }
        }
        (Some(_), None) | (None, Some(_)) => 20,
        (None, None) => 0,
    };
    let motif_strand = if hint > 0 { "-" } else { "+" };
    if let Some(start) = start_hit {
        if start.split('\t').nth(1) == Some(motif_strand) {
            score += 5;
        }
    }
    if let Some(end) = end_hit {
        if end.split('\t').nth(1) == Some(motif_strand) {
            score += 5;
        }
    }
    score
}

/// Clusters motif-validated splice candidates using the CIRI-AS two-pass rule.
fn cluster_candidates(candidates: &[PositiveCandidate]) -> Vec<SpliceCluster> {
    if candidates.is_empty() {
        return Vec::new();
    }
    let mut sorted: Vec<usize> = (0..candidates.len()).collect();
    sorted.sort_by(|&a, &b| {
        candidates[a]
            .chr
            .cmp(&candidates[b].chr)
            .then(candidates[a].site2.cmp(&candidates[b].site2))
    });

    let mut groups: Vec<Vec<usize>> = Vec::new();
    for idx in sorted {
        let same_group = groups.last().and_then(|g| g.last()).is_some_and(|&prev| {
            candidates[idx].chr == candidates[prev].chr
                && candidates[idx].site2 - candidates[prev].site2 <= 3
        });
        if same_group {
            groups.last_mut().unwrap().push(idx);
        } else {
            groups.push(vec![idx]);
        }
    }

    let mut clusters = Vec::new();
    for mut group in groups {
        group.sort_by_key(|&idx| candidates[idx].site1);
        let mut current: Vec<usize> = Vec::new();
        for idx in group {
            let same_cluster = current.last().is_some_and(|&prev| {
                candidates[idx].chr == candidates[prev].chr
                    && candidates[idx].site1 - candidates[prev].site1 <= 3
            });
            if same_cluster {
                current.push(idx);
            } else {
                push_cluster(&mut clusters, candidates, &current);
                current = vec![idx];
            }
        }
        push_cluster(&mut clusters, candidates, &current);
    }
    clusters
}

/// Finalizes one candidate cluster if it passes CIRI-AS stringency.
fn push_cluster(
    clusters: &mut Vec<SpliceCluster>,
    candidates: &[PositiveCandidate],
    current: &[usize],
) {
    if current.is_empty() {
        return;
    }
    let median_idx = current[(current.len() - 1) / 2];
    let mut cigar_counts = [0usize; 3];
    for slot in 0..3 {
        let mut seen = HashSet::new();
        for &idx in current {
            if let Some(cigar) = &candidates[idx].cigars[slot] {
                seen.insert(cigar.clone());
            }
        }
        cigar_counts[slot] = seen.len();
    }
    if cigar_counts.iter().sum::<usize>() >= STRINGENCY {
        clusters.push(SpliceCluster {
            chr: candidates[median_idx].chr.clone(),
            site1: candidates[median_idx].site1,
            site2: candidates[median_idx].site2,
            reads: current.to_vec(),
            cigar_counts,
        });
    }
}

/// Writes the CIRI-AS `_splice.list` table.
fn write_splice_list(
    path: &str,
    clusters: &[SpliceCluster],
    candidates: &[PositiveCandidate],
    read_to_circ: &HashMap<String, String>,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "ID\tchr\tsplice_start\tsplice_end\t#supporting_reads\tSM_MS_SMS\tjunction_reads_ID"
    )?;
    for cluster in clusters {
        write!(
            writer,
            "{}:{}|{}\t{}\t{}\t{}\t{}\t{}_{}_{}\t",
            cluster.chr,
            cluster.site1,
            cluster.site2,
            cluster.chr,
            cluster.site1,
            cluster.site2,
            cluster.reads.len(),
            cluster.cigar_counts[0],
            cluster.cigar_counts[1],
            cluster.cigar_counts[2]
        )?;
        for &idx in &cluster.reads {
            let candidate = &candidates[idx];
            let circ = read_to_circ
                .get(&candidate.read_id)
                .map(String::as_str)
                .unwrap_or("");
            write!(writer, "{}({}),", candidate.read_id, circ)?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

/// Predicts CIRI-AS cirexons from validated splice clusters and coverage.
///
/// This is the first `.list` implementation stage. It ports the deterministic
/// exon-candidate and coverage-validation path from CIRI-AS v1.2, but leaves
/// intron-retention rescue and AS path classification for the next stage so the
/// cirexon table can be validated independently.
fn predict_cirexons(
    state: &ScanState,
    splice_clusters: &[SpliceCluster],
    reference: &HashMap<String, String>,
    annotation: Option<&Annotation>,
) -> (Vec<CirexonRecord>, Vec<IsoformRecord>) {
    let mut circ_by_chr: BTreeMap<&str, Vec<&CircRecord>> = BTreeMap::new();
    for circ in state.circ_by_id.values() {
        circ_by_chr.entry(circ.chr.as_str()).or_default().push(circ);
    }
    let mut clusters_by_chr: BTreeMap<&str, Vec<&SpliceCluster>> = BTreeMap::new();
    for cluster in splice_clusters {
        clusters_by_chr
            .entry(cluster.chr.as_str())
            .or_default()
            .push(cluster);
    }
    for clusters in clusters_by_chr.values_mut() {
        clusters.sort_by_key(|cluster| (cluster.site2, cluster.site1));
    }

    let mut cirexons = Vec::new();
    let mut isoforms = Vec::new();
    for (chr, mut circs) in circ_by_chr {
        circs.sort_by_key(|circ| (circ.start, circ.end));
        let chr_clusters = clusters_by_chr.get(chr).map(Vec::as_slice).unwrap_or(&[]);
        for circ in circs {
            let Some(strand) = infer_circ_strand(circ, reference, annotation) else {
                continue;
            };
            let in_circ: Vec<&SpliceCluster> = chr_clusters
                .iter()
                .copied()
                .filter(|cluster| {
                    cluster.site2 >= circ.start
                        && cluster.site2 <= circ.end
                        && cluster.site1 < circ.end
                        && cluster.site2 > circ.start
                })
                .collect();
            if in_circ.is_empty() {
                continue;
            }
            let junction_context = build_junction_context(circ, state);
            let Some(cirexon_context) =
                build_cirexon_context(circ, &in_circ, &junction_context, state)
            else {
                continue;
            };
            let known_exons = known_exons_for_circ(circ, annotation);
            let mut validated = validate_candidate_cirexons(
                circ,
                &strand,
                &cirexon_context,
                &junction_context,
                state,
                &known_exons,
                annotation.is_some(),
            );
            validated.sort_by_key(|record| (record.start, record.end));
            let mut circ_isoforms =
                build_full_length_isoforms(circ, &strand, &validated, &cirexon_context, annotation);
            isoforms.append(&mut circ_isoforms);
            append_ordered_cirexons(&mut cirexons, &mut validated, &strand);
        }
    }
    (cirexons, isoforms)
}

/// Junction-read intervals and per-base coverage for one circRNA.
struct JunctionContext {
    mappings: Vec<ReadMapping>,
    coverage: HashMap<i32, u32>,
}

/// CIRI-AS splice-derived boundary state for one circRNA.
struct CirexonContext {
    supporting_start: BTreeMap<i32, BoundarySupport>,
    supporting_end: BTreeMap<i32, BoundarySupport>,
    splice_read_start: HashMap<i32, usize>,
    splice_read_end: HashMap<i32, usize>,
    splice_links: HashMap<(i32, i32), usize>,
    splice_across_count: HashMap<i32, usize>,
}

/// Builds mapping intervals for BSJ reads assigned to one circRNA.
fn build_junction_context(circ: &CircRecord, state: &ScanState) -> JunctionContext {
    let mut mappings = Vec::new();
    let mut coverage = HashMap::new();
    for read in &circ.junction_reads {
        let Some(by_reverse) = state.read_mappings.get(read) else {
            continue;
        };
        for bucket in by_reverse {
            for mapping in bucket {
                if mapping.chr == circ.chr
                    && mapping.start >= circ.start - 6
                    && mapping.end <= circ.end + 6
                {
                    mappings.push(mapping.clone());
                    for pos in mapping.start..=mapping.end {
                        *coverage.entry(pos).or_insert(0) += 1;
                    }
                }
            }
        }
    }
    mappings.sort_by_key(|mapping| (mapping.start, mapping.end));
    JunctionContext { mappings, coverage }
}

/// Builds candidate exon boundary support from splice clusters inside a circRNA.
fn build_cirexon_context(
    circ: &CircRecord,
    splice_clusters: &[&SpliceCluster],
    junction_context: &JunctionContext,
    state: &ScanState,
) -> Option<CirexonContext> {
    let mut supporting_start: BTreeMap<i32, BoundarySupport> = BTreeMap::new();
    let mut supporting_end: BTreeMap<i32, BoundarySupport> = BTreeMap::new();
    let mut splice_read_start: HashMap<i32, usize> = HashMap::new();
    let mut splice_read_end: HashMap<i32, usize> = HashMap::new();
    let mut splice_links: HashMap<(i32, i32), usize> = HashMap::new();
    let mut splice_across_count: HashMap<i32, usize> = HashMap::new();
    let mut any_bsj_supported_splice = false;

    for cluster in splice_clusters {
        supporting_start
            .entry(cluster.site1)
            .or_insert(BoundarySupport::Count(0));
        supporting_end
            .entry(cluster.site2)
            .or_insert(BoundarySupport::Count(0));
        *splice_read_start.entry(cluster.site1).or_insert(0) += cluster.reads.len();
        *splice_read_end.entry(cluster.site2).or_insert(0) += cluster.reads.len();

        let mut bsj_support = 0usize;
        for &idx in &cluster.reads {
            let candidate = &state.candidates[idx];
            if state
                .junction_read_to_circ
                .get(&candidate.read_id)
                .is_some_and(|id| id == &circ.id)
            {
                bsj_support += 1;
            }
        }
        if bsj_support > 0 {
            any_bsj_supported_splice = true;
            supporting_start.insert(cluster.site1, BoundarySupport::Count(bsj_support));
            supporting_end.insert(cluster.site2, BoundarySupport::Count(bsj_support));
            splice_links.insert((cluster.site2, cluster.site1), bsj_support);
        }

        for site in [cluster.site2, cluster.site1] {
            splice_across_count
                .entry(site)
                .or_insert_with(|| count_mappings_across_site(&junction_context.mappings, site));
        }
    }

    supporting_start.insert(circ.start, BoundarySupport::CircStart);
    supporting_end.insert(circ.end, BoundarySupport::CircEnd);

    if any_bsj_supported_splice {
        Some(CirexonContext {
            supporting_start,
            supporting_end,
            splice_read_start,
            splice_read_end,
            splice_links,
            splice_across_count,
        })
    } else {
        None
    }
}

/// Counts junction-read alignments that span a splice site by six bases.
fn count_mappings_across_site(mappings: &[ReadMapping], site: i32) -> usize {
    mappings
        .iter()
        .filter(|mapping| mapping.start <= site - 6 && mapping.end >= site + 6)
        .count()
}

/// Validates candidate cirexons using CIRI-AS coverage and support rules.
fn validate_candidate_cirexons(
    circ: &CircRecord,
    strand: &str,
    context: &CirexonContext,
    junction_context: &JunctionContext,
    state: &ScanState,
    known_exons: &[(i32, i32)],
    has_annotation: bool,
) -> Vec<CirexonRecord> {
    let starts: Vec<i32> = context.supporting_start.keys().copied().collect();
    let ends: Vec<i32> = context.supporting_end.keys().copied().collect();
    let mut candidates: Vec<(i32, i32, i32)> = Vec::new();
    for (i, start) in starts.iter().enumerate() {
        for (j, end) in ends.iter().enumerate() {
            let start_support = context.supporting_start[start].numeric();
            let end_support = context.supporting_end[end].numeric();
            if *end >= *start + MIN_EXON_LENGTH - 1
                && *end <= *start + MAX_EXON_LENGTH - 1
                && (start_support + end_support > 0 || i == 0 || j + 1 == ends.len())
            {
                candidates.push((*start, *end, *end - *start + 1));
            }
        }
    }
    candidates.sort_by_key(|&(start, end, length)| (length, start, end));

    let mut validated = Vec::new();
    let mut validated_starts: HashSet<i32> = HashSet::new();
    let mut validated_ends: HashSet<i32> = HashSet::new();
    for (start, end, _) in candidates {
        if validated_starts.contains(&start) && validated_ends.contains(&end) {
            continue;
        }
        let start_support = context.supporting_start[&start];
        let end_support = context.supporting_end[&end];
        if validated_starts.contains(&start) && end_support.numeric() == 0 {
            continue;
        }
        if validated_ends.contains(&end) && start_support.numeric() == 0 {
            continue;
        }
        let validation = exon_coverage_validation_single(
            &circ.chr,
            start,
            end,
            start_support,
            end_support,
            *context.splice_read_start.get(&start).unwrap_or(&0),
            *context.splice_read_end.get(&end).unwrap_or(&0),
            *context.splice_across_count.get(&start).unwrap_or(&0),
            *context.splice_across_count.get(&end).unwrap_or(&0),
            &junction_context.coverage,
            state,
        );
        if validation.code >= 1 {
            validated_starts.insert(start);
            validated_ends.insert(end);
            let is_icf = has_annotation
                && !known_exons.iter().any(|&(known_start, known_end)| {
                    known_start <= start + 6 && known_end >= end - 6
                });
            validated.push(CirexonRecord {
                circ_id: circ.id.clone(),
                chr: circ.chr.clone(),
                circ_start: circ.start,
                circ_end: circ.end,
                strand: strand.to_string(),
                junction_read_count: circ.junction_read_count.clone(),
                pcc: circ.pcc.clone(),
                non_junction_reads: circ.non_junction_reads.clone(),
                junction_reads_ratio: circ.junction_reads_ratio.clone(),
                circ_type: circ.circ_type.clone(),
                gene_id: circ.gene_id.clone(),
                cirexon_id: String::new(),
                start,
                end,
                start_support: start_support.as_output(),
                end_support: end_support.as_output(),
                coverage_median: validation.median,
                is_icf,
            });
        }
    }
    validated
}

/// Builds full-length anchored isoform paths from validated cirexons.
///
/// The goal here is not CIRI-AS event classification. We use the same exon graph
/// ingredients, but stop at concrete exon-chain reconstruction: nodes are
/// validated cirexons, explicit edges are BSJ-read-supported internal splices,
/// and the only fallback edge is the nearest downstream neighbouring exon used
/// by CIRI-AS when no splice link exists for a node.
fn build_full_length_isoforms(
    circ: &CircRecord,
    strand: &str,
    validated: &[CirexonRecord],
    context: &CirexonContext,
    annotation: Option<&Annotation>,
) -> Vec<IsoformRecord> {
    if validated.is_empty() {
        return Vec::new();
    }
    let mut exons: Vec<CirexonRecord> = validated.to_vec();
    exons.sort_by_key(|record| (record.start, record.end));
    let links = build_isoform_links(&exons, context);
    let starts: Vec<usize> = exons
        .iter()
        .enumerate()
        .filter_map(|(idx, exon)| (exon.start == circ.start).then_some(idx))
        .collect();
    if starts.is_empty() {
        return Vec::new();
    }

    let mut paths = Vec::new();
    for start_idx in starts {
        let mut path = Vec::new();
        enumerate_isoform_paths(start_idx, &exons, &links, &mut path, &mut paths);
        if paths.len() >= MAX_ISOFORM_PATHS_PER_CIRC {
            break;
        }
    }
    if paths.is_empty() {
        return Vec::new();
    }

    let mut records = Vec::new();
    for path in paths {
        let path_exons: Vec<(i32, i32)> = path
            .iter()
            .map(|&idx| (exons[idx].start, exons[idx].end))
            .collect();
        let mut transcript_exons = path_exons.clone();
        if strand == "-" {
            transcript_exons.reverse();
        }
        let junctions: Vec<(i32, i32)> = path_exons
            .windows(2)
            .map(|pair| (pair[0].1, pair[1].0))
            .collect();
        let internal_split_reads = junctions
            .iter()
            .map(|junction| context.splice_links.get(junction).copied().unwrap_or(0))
            .sum::<usize>();
        let annotation_supported_junctions = junctions
            .iter()
            .filter(|&&(end, start)| junction_has_annotation(circ, end, start, annotation))
            .count();
        let de_novo_supported_junctions = junctions
            .len()
            .saturating_sub(annotation_supported_junctions);
        let isoform_len = transcript_exons
            .iter()
            .map(|(start, end)| end - start + 1)
            .sum::<i32>();
        let icf_penalty = path
            .iter()
            .filter(|&&idx| exons[idx].is_icf)
            .count()
            .saturating_mul(2) as i32;
        let bsj_seed_reads = circ.junction_read_count.parse::<usize>().unwrap_or(0);
        let score = bsj_seed_reads as i32
            + (internal_split_reads as i32 * 10)
            + (annotation_supported_junctions as i32 * 3)
            - (de_novo_supported_junctions as i32)
            - icf_penalty;
        records.push(IsoformRecord {
            circ_id: circ.id.clone(),
            isoform_id: String::new(),
            chr: circ.chr.clone(),
            start: circ.start,
            end: circ.end,
            strand: strand.to_string(),
            path_tier: "anchored".to_string(),
            bsj_supported: "yes".to_string(),
            gene_id: circ.gene_id.clone(),
            exons: transcript_exons,
            junctions,
            isoform_len,
            bsj_seed_reads,
            internal_split_reads,
            boundary_clip_reads: 0,
            anomalous_pair_reads: 0,
            annotation_supported_junctions,
            de_novo_supported_junctions,
            score,
            rank: 0,
        });
    }

    records.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.exons.len().cmp(&b.exons.len()))
            .then(a.exon_chain().cmp(&b.exon_chain()))
    });
    records.dedup_by(|a, b| a.exons == b.exons);
    for (idx, record) in records.iter_mut().enumerate() {
        record.rank = idx + 1;
        record.isoform_id = format!("{}.isoform{}", circ.id, idx + 1);
    }
    records
}

/// Builds graph edges between validated cirexons for full-length path search.
fn build_isoform_links(
    exons: &[CirexonRecord],
    context: &CirexonContext,
) -> HashMap<usize, Vec<(usize, usize)>> {
    let mut links: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
    for (idx, exon) in exons.iter().enumerate() {
        let mut explicit = Vec::new();
        for (next_idx, next) in exons.iter().enumerate().skip(idx + 1) {
            if next.start <= exon.end {
                continue;
            }
            if let Some(&support) = context.splice_links.get(&(exon.end, next.start)) {
                explicit.push((next_idx, support));
            }
        }
        if !explicit.is_empty() {
            links.insert(idx, explicit);
            continue;
        }
        if let Some(first_next_start) = exons
            .iter()
            .skip(idx + 1)
            .filter(|next| next.start > exon.end)
            .map(|next| next.start)
            .min()
        {
            let neighbours = exons
                .iter()
                .enumerate()
                .skip(idx + 1)
                .filter_map(|(next_idx, next)| {
                    (next.start == first_next_start).then_some((next_idx, 0))
                })
                .collect::<Vec<_>>();
            if !neighbours.is_empty() {
                links.insert(idx, neighbours);
            }
        }
    }
    links
}

/// Depth-first path enumeration with a hard cap per circRNA.
fn enumerate_isoform_paths(
    idx: usize,
    exons: &[CirexonRecord],
    links: &HashMap<usize, Vec<(usize, usize)>>,
    path: &mut Vec<usize>,
    paths: &mut Vec<Vec<usize>>,
) {
    if paths.len() >= MAX_ISOFORM_PATHS_PER_CIRC {
        return;
    }
    if path.contains(&idx) {
        return;
    }
    path.push(idx);
    if exons[idx].end == exons[idx].circ_end {
        paths.push(path.clone());
        path.pop();
        return;
    }
    if let Some(nexts) = links.get(&idx) {
        for &(next_idx, _) in nexts {
            enumerate_isoform_paths(next_idx, exons, links, path, paths);
            if paths.len() >= MAX_ISOFORM_PATHS_PER_CIRC {
                break;
            }
        }
    }
    path.pop();
}

/// True when a path junction is supported by annotated exon boundaries.
fn junction_has_annotation(
    circ: &CircRecord,
    exon_end: i32,
    next_start: i32,
    annotation: Option<&Annotation>,
) -> bool {
    let Some(annotation) = annotation else {
        return false;
    };
    let end_key = format!("{}\t{}", circ.chr, exon_end);
    let start_key = format!("{}\t{}", circ.chr, next_start);
    annotation.chr_exon_end_map.contains_key(&end_key)
        && annotation.chr_exon_start_map.contains_key(&start_key)
}

impl IsoformRecord {
    /// Formats exon coordinates in transcript order for `.isoforms`.
    fn exon_chain(&self) -> String {
        self.exons
            .iter()
            .map(|(start, end)| format!("{}:{}!{}", start, end, self.strand))
            .collect::<Vec<_>>()
            .join(",")
    }

    /// Formats internal splice junctions in genomic order for debugging.
    fn junction_chain(&self) -> String {
        if self.junctions.is_empty() {
            return "NA".to_string();
        }
        self.junctions
            .iter()
            .map(|(end, start)| format!("{}:{}", end, start))
            .collect::<Vec<_>>()
            .join(",")
    }
}

/// Ports CIRI-AS `exon_coverage_validation_single` return-code semantics.
fn exon_coverage_validation_single(
    chr: &str,
    start: i32,
    end: i32,
    tag_start: BoundarySupport,
    tag_end: BoundarySupport,
    tag_start2: usize,
    tag_end2: usize,
    start_across: usize,
    end_across: usize,
    junction_coverage: &HashMap<i32, u32>,
    state: &ScanState,
) -> CoverageValidation {
    let length = end - start + 1;
    let last_group = (end - start) / (MIN_INTRON - 5);
    let mut total_cov0 = 0usize;
    let mut max_junc_gap = 0usize;
    let mut cont_junc_cov0 = 0usize;
    let mut junc_cov0_count = 0usize;
    let mut u_pre = vec![0.0; last_group as usize + 1];
    let mut u_after = vec![0.0; last_group as usize + 1];
    let mut u_pre_total = 0.0;
    let mut u_after_total = 0.0;
    let mut cov_all = Vec::with_capacity(length.max(0) as usize);

    for pos in start..=end {
        let group = ((pos - start) / (MIN_INTRON - 5)) as usize;
        let cov = coverage_at(state, chr, pos);
        if cov == 0 {
            total_cov0 += 1;
        }
        if junction_coverage.get(&pos).copied().unwrap_or(0) > 0 {
            cont_junc_cov0 = 0;
        } else {
            junc_cov0_count += 1;
            cont_junc_cov0 += 1;
            max_junc_gap = max_junc_gap.max(cont_junc_cov0);
        }
        for flank in 6..=MIN_INTRON {
            let pre = coverage_at(state, chr, start - flank);
            if cov > pre {
                u_pre[group] += 1.0;
                u_pre_total += 1.0;
            } else if cov == pre {
                u_pre[group] += 0.5;
                u_pre_total += 0.5;
            }
            let after = coverage_at(state, chr, end + flank);
            if cov > after {
                u_after[group] += 1.0;
                u_after_total += 1.0;
            } else if cov == after {
                u_after[group] += 0.5;
                u_after_total += 0.5;
            }
        }
        cov_all.push(cov);
    }

    if last_group > 0 {
        let tail_start = end - MIN_INTRON + 6;
        let tail_end = end - (end - start) % (MIN_INTRON - 5) - 1;
        for pos in tail_start..=tail_end {
            let cov = coverage_at(state, chr, pos);
            for flank in 6..=MIN_INTRON {
                let pre = coverage_at(state, chr, start - flank);
                if cov > pre {
                    *u_pre.last_mut().unwrap() += 1.0;
                } else if cov == pre {
                    *u_pre.last_mut().unwrap() += 0.5;
                }
                let after = coverage_at(state, chr, end + flank);
                if cov > after {
                    *u_after.last_mut().unwrap() += 1.0;
                } else if cov == after {
                    *u_after.last_mut().unwrap() += 0.5;
                }
            }
        }
    }
    let sample_size = if last_group > 0 {
        MIN_INTRON - 5
    } else {
        length
    };
    let mut decline_pre = 0usize;
    let mut decline_after = 0usize;
    for idx in 0..=last_group as usize {
        if z_calculation(u_pre[idx], sample_size, MIN_INTRON - 5) <= Z_ALPHA {
            decline_pre += 1;
        }
        if z_calculation(u_after[idx], sample_size, MIN_INTRON - 5) <= Z_ALPHA {
            decline_after += 1;
        }
    }
    cov_all.sort_unstable();
    let median = cov_all[(cov_all.len() - 1) / 2] as i32;
    let z_pre_total = z_calculation(u_pre_total, length, MIN_INTRON - 5);
    let z_after_total = z_calculation(u_after_total, length, MIN_INTRON - 5);

    let code = if matches!(tag_start, BoundarySupport::Count(0)) && start_across > 0 {
        0
    } else if matches!(tag_end, BoundarySupport::Count(0)) && end_across > 0 {
        -1
    } else if total_cov0 > 0 {
        -2
    } else if tag_start2 < 2 && !matches!(tag_start, BoundarySupport::CircStart) && decline_pre > 0
    {
        -3
    } else if tag_end2 < 2 && !matches!(tag_end, BoundarySupport::CircEnd) && decline_after > 0 {
        -4
    } else if tag_start2 < 2
        && !matches!(tag_start, BoundarySupport::CircStart)
        && z_pre_total <= Z_ALPHA
    {
        -8
    } else if tag_end2 < 2
        && !matches!(tag_end, BoundarySupport::CircEnd)
        && z_after_total <= Z_ALPHA
    {
        -9
    } else if (junc_cov0_count as f64) / (length as f64) <= 0.2 {
        1
    } else if !matches!(tag_start, BoundarySupport::Count(0))
        && !matches!(tag_end, BoundarySupport::Count(0))
        && max_junc_gap < MIN_INTRON as usize
    {
        2
    } else {
        -10
    };
    CoverageValidation { code, median }
}

/// Mann-Whitney-style Z score used by CIRI-AS coverage validation.
fn z_calculation(wxy: f64, n: i32, m: i32) -> f64 {
    let n = n as f64;
    let m = m as f64;
    let wy = wxy + n * (n + 1.0) / 2.0;
    let total = n + m;
    (wy - n * (total + 1.0) / 2.0) / (m * n * (total + 1.0) / 12.0).sqrt()
}

/// Returns all known annotation exons that can explain a circRNA interval.
fn known_exons_for_circ(circ: &CircRecord, annotation: Option<&Annotation>) -> Vec<(i32, i32)> {
    let Some(annotation) = annotation else {
        return Vec::new();
    };
    circ.gene_id
        .split(',')
        .filter_map(|gene| annotation.gene_exon_map.get(gene))
        .flat_map(|exons| exons.iter().copied())
        .filter(|&(start, end)| start >= circ.start && end <= circ.end)
        .collect()
}

/// Infers circRNA strand from Summary, splice motifs, then annotation.
fn infer_circ_strand(
    circ: &CircRecord,
    reference: &HashMap<String, String>,
    annotation: Option<&Annotation>,
) -> Option<String> {
    if circ.strand == "+" || circ.strand == "-" {
        return Some(circ.strand.clone());
    }
    let seq = reference.get(&circ.chr)?;
    let start_2bp = perl_substr(seq, circ.start - 3, 2).to_ascii_uppercase();
    let end_2bp = perl_substr(seq, circ.end, 2).to_ascii_uppercase();
    if start_2bp.contains("AG") && end_2bp.contains("GT") {
        Some("+".to_string())
    } else if start_2bp.contains("AC") && end_2bp.contains("CT") {
        Some("-".to_string())
    } else if start_2bp.contains("AG") || end_2bp.contains("GT") {
        Some("+".to_string())
    } else if start_2bp.contains("AC") || end_2bp.contains("CT") {
        Some("-".to_string())
    } else {
        annotation.and_then(|anno| {
            let start_key = format!("{}\t{}", circ.chr, circ.start);
            let end_key = format!("{}\t{}", circ.chr, circ.end);
            anno.chr_exon_start_map
                .get(&start_key)
                .or_else(|| anno.chr_exon_end_map.get(&end_key))
                .and_then(|value| value.split('\t').nth(1))
                .map(str::to_string)
        })
    }
}

/// Appends validated cirexons in CIRI-AS strand-aware `cirexonN` order.
fn append_ordered_cirexons(
    out: &mut Vec<CirexonRecord>,
    validated: &mut [CirexonRecord],
    strand: &str,
) {
    if validated.is_empty() {
        return;
    }
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for idx in 0..validated.len() {
        if idx == 0 || validated[idx].start > validated[idx - 1].end {
            groups.push(vec![idx]);
        } else {
            groups.last_mut().unwrap().push(idx);
        }
    }
    let group_order: Vec<usize> = if strand == "-" {
        (0..groups.len()).rev().collect()
    } else {
        (0..groups.len()).collect()
    };
    for (display_idx, group_idx) in group_order.into_iter().enumerate() {
        for &record_idx in &groups[group_idx] {
            let mut record = validated[record_idx].clone();
            record.cirexon_id = format!("cirexon{}", display_idx + 1);
            out.push(record);
        }
    }
}

/// Returns global CIRI-AS coverage at one genomic position.
fn coverage_at(state: &ScanState, chr: &str, pos: i32) -> u32 {
    state
        .coverage
        .get(chr)
        .and_then(|chr_cov| chr_cov.get(&pos).copied())
        .unwrap_or(0)
}

/// Writes the `.list` table used by CIRI-AS cirexon output.
fn write_cirexon_list(path: &str, records: &[CirexonRecord]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "circRNA_id\tchr\tstart\tend\tstrand\t#junction_reads\tPCC(MS_SM_SMS)\t#non_junction_reads\tjunction_reads_ratio\tcircRNA_type\tgene_id\tcirexon_id\tcirexon_start\tcirexon_end\t#start_supporting_BSJ_read\t#end_supporting_BSJ_read\tsequencing_depth_median\tif_ICF"
    )?;
    for record in records {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.circ_id,
            record.chr,
            record.circ_start,
            record.circ_end,
            record.strand,
            record.junction_read_count,
            record.pcc,
            record.non_junction_reads,
            record.junction_reads_ratio,
            record.circ_type,
            record.gene_id,
            record.cirexon_id,
            record.start,
            record.end,
            record.start_support,
            record.end_support,
            record.coverage_median,
            if record.is_icf { "ICF" } else { "non_ICF" }
        )?;
    }
    Ok(())
}

/// Writes the `_AS.list` header used by CIRI-AS alternative-splicing output.
fn write_as_header(path: &str) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "circRNA_id\talternatively_spliced_exon\tAS_type\tpsi_estimation_without_correction\tpsi_estimation_after_correction"
    )?;
    Ok(())
}

/// Writes full-length isoform paths reconstructed from cirexon graph traversal.
fn write_isoforms(path: &str, records: &[IsoformRecord]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "circ_id\tisoform_id\tchr\tstart\tend\tstrand\tpath_tier\tbsj_supported\tgene_id\texon_chain\tjunction_chain\texon_count\tisoform_len\tbsj_seed_reads\tinternal_split_reads\tboundary_clip_reads\tanomalous_pair_reads\tannotation_supported_junctions\tde_novo_supported_junctions\tscore\trank"
    )?;
    for record in records {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.circ_id,
            record.isoform_id,
            record.chr,
            record.start,
            record.end,
            record.strand,
            record.path_tier,
            record.bsj_supported,
            record.gene_id,
            record.exon_chain(),
            record.junction_chain(),
            record.exons.len(),
            record.isoform_len,
            record.bsj_seed_reads,
            record.internal_split_reads,
            record.boundary_clip_reads,
            record.anomalous_pair_reads,
            record.annotation_supported_junctions,
            record.de_novo_supported_junctions,
            record.score,
            record.rank
        )?;
    }
    Ok(())
}

/// Writes per-circRNA full-length isoform counts.
///
/// This summary is the quickest output for the current development goal: it
/// reports how many anchored isoform paths were reconstructed for each
/// Summary-confirmed circRNA, including zeros so missing path cases are explicit.
fn write_isoform_summary(
    path: &str,
    state: &ScanState,
    cirexons: &[CirexonRecord],
    isoforms: &[IsoformRecord],
) -> Result<()> {
    let mut cirexon_count: HashMap<&str, usize> = HashMap::new();
    for cirexon in cirexons {
        *cirexon_count.entry(cirexon.circ_id.as_str()).or_insert(0) += 1;
    }
    let mut isoform_count: HashMap<&str, usize> = HashMap::new();
    for isoform in isoforms {
        *isoform_count.entry(isoform.circ_id.as_str()).or_insert(0) += 1;
    }

    let mut circ_records: Vec<&CircRecord> = state.circ_by_id.values().collect();
    circ_records.sort_by_key(|circ| (&circ.chr, circ.start, circ.end));
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "circ_id\tchr\tstart\tend\tstrand\tgene_id\tjunction_reads\tcirexon_count\tisoform_count"
    )?;
    for circ in circ_records {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            circ.id,
            circ.chr,
            circ.start,
            circ.end,
            circ.strand,
            circ.gene_id,
            circ.junction_read_count,
            cirexon_count.get(circ.id.as_str()).copied().unwrap_or(0),
            isoform_count.get(circ.id.as_str()).copied().unwrap_or(0)
        )?;
    }
    Ok(())
}

/// Writes FASTA sequences for reconstructed full-length isoforms.
///
/// Coordinates remain 1-based inclusive in the table, but FASTA extraction uses
/// Rust's 0-based end-exclusive slicing in this one helper so coordinate
/// conversion is not duplicated around the reconstruction code.
fn write_isoform_fasta(
    path: &str,
    records: &[IsoformRecord],
    reference: &HashMap<String, String>,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for record in records {
        let seq = isoform_sequence(record, reference)?;
        writeln!(
            writer,
            ">{}|{}|tier={}|bsj_supported={}|score={}|rank={}|exons={}",
            record.circ_id,
            record.isoform_id,
            record.path_tier,
            record.bsj_supported,
            record.score,
            record.rank,
            record.exon_chain()
        )?;
        for chunk in seq.as_bytes().chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }
    Ok(())
}

/// Extracts one isoform sequence from the reference FASTA.
fn isoform_sequence(record: &IsoformRecord, reference: &HashMap<String, String>) -> Result<String> {
    let chr_seq = reference
        .get(&record.chr)
        .ok_or_else(|| anyhow!("missing reference sequence for {}", record.chr))?;
    let mut seq = String::with_capacity(record.isoform_len.max(0) as usize);
    for &(start, end) in &record.exons {
        if start < 1 || end < start || end as usize > chr_seq.len() {
            bail!(
                "invalid isoform exon coordinate {}:{}-{} for {}",
                record.chr,
                start,
                end,
                record.isoform_id
            );
        }
        let fragment = &chr_seq[(start - 1) as usize..end as usize];
        if record.strand == "-" {
            seq.push_str(&reverse_complement(fragment));
        } else {
            seq.push_str(fragment);
        }
    }
    Ok(seq)
}

/// Writes a compact development log for CIRI-AS sidecar runs.
fn write_as_log(path: &str, state: &ScanState) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "{} circRNAs are loaded.", state.circ_by_id.len())?;
    writeln!(
        writer,
        "{} circular junction reads are loaded.",
        state.stats.junction_reads_loaded
    )?;
    writeln!(
        writer,
        "{} circular junction reads were found from the alignment file.",
        state.stats.junction_reads_seen
    )?;
    writeln!(
        writer,
        "{} candidate splice junctions are recognized, of which {} are from known BSJ reads.",
        state.stats.known_candidates + state.stats.add_candidates,
        state.stats.known_candidates
    )?;
    writeln!(
        writer,
        "{} candidate splice junctions have splicing signals.",
        state.stats.motif_validated
    )?;
    writeln!(
        writer,
        "In sum, {} non-redundant splice junctions within circRNAs are detected.",
        state.stats.final_splice_clusters
    )?;
    writeln!(
        writer,
        "{} cirexons are predicted from internal splice junctions and coverage.",
        state.stats.final_cirexons
    )?;
    writeln!(
        writer,
        "{} full-length isoform paths are reconstructed from cirexon graph traversal.",
        state.stats.final_isoforms
    )?;
    Ok(())
}

/// Converts validated read classes into `<prefix>.segments` rows.
///
/// The current phase emits only BSJ reads and non-BSJ backward reads. The
/// classification itself still comes from the existing CIRI-AS-style scan; this
/// helper is only responsible for selecting the best alignment chain and
/// materializing it as simulator-compatible segments.
fn build_segment_records(state: &ScanState) -> Result<Vec<SegmentRecord>> {
    let backward_ids: HashSet<&str> = state
        .candidates
        .iter()
        .map(|candidate| candidate.read_id.as_str())
        .filter(|read_id| !state.junction_read_to_circ.contains_key(*read_id))
        .collect();
    let read_strand_hints = build_read_strand_hints(state);
    let read_junction_hints = build_read_junction_hints(state);
    let correction = SegmentCorrectionContext {
        reference: state.reference,
        annotation: state.annotation,
        junction_support: None,
    };
    let mut read_ids: Vec<&str> = state.segment_groups.keys().map(String::as_str).collect();
    read_ids.sort_unstable();

    let mut out = Vec::new();
    for read_id in read_ids {
        let Some(records) = state.segment_groups.get(read_id) else {
            continue;
        };
        if let Some(circ_id) = state.junction_read_to_circ.get(read_id) {
            let circ = state
                .circ_by_id
                .get(circ_id)
                .ok_or_else(|| anyhow!("missing circ record {}", circ_id))?;
            if let Some(record) = build_bsj_segment_record(
                read_id,
                records,
                circ,
                state
                    .mate_bsj_evidence
                    .get(read_id)
                    .map(Vec::as_slice)
                    .unwrap_or(&[]),
                read_junction_hints
                    .get(read_id)
                    .map(Vec::as_slice)
                    .unwrap_or(&[]),
                Some(&correction),
                state.read_len,
            ) {
                out.push(record);
            }
        } else if backward_ids.contains(read_id) {
            if let Some(record) = build_backward_segment_record(
                read_id,
                records,
                read_strand_hints.get(read_id).copied().flatten(),
                read_junction_hints
                    .get(read_id)
                    .map(Vec::as_slice)
                    .unwrap_or(&[]),
                Some(&correction),
                state.read_len,
            ) {
                out.push(record);
            }
        }
    }
    sort_segment_records(&mut out);
    Ok(out)
}

/// Builds the default sidecar-mode `<prefix>.segments` rows.
///
/// Confirmed BSJ reads are materialized from Scan1/Scan2 sidecar alignments so
/// the output keeps local clip pseudo rows captured during validation. Backward
/// rows are appended from the supplemental non-BSJ scan and intentionally carry
/// no circRNA assignment because the current stage has not run circRNA-level
/// interval arbitration for them.
fn build_sidecar_segment_records(state: &ScanState) -> Result<Vec<SegmentRecord>> {
    let mut out = build_confirmed_bsj_segment_records(state)?;
    out.extend(build_backward_segment_records(state)?);
    sort_segment_records(&mut out);
    Ok(out)
}

/// Converts sidecar evidence into confirmed `type=bsj` segment rows only.
///
/// This is the default post-Summary path: Scan1/Scan2 have already captured
/// mapper blocks for BSJ-supporting read groups, and `.out` supplies the final
/// circRNA assignment. Backward/internal reads are deliberately excluded here
/// because they require a later extra scan and must not be inferred from the
/// CIRI3 Summary-only evidence set.
fn build_confirmed_bsj_segment_records(state: &ScanState) -> Result<Vec<SegmentRecord>> {
    let profile = segments_profile_enabled();
    let first_pass_correction = SegmentCorrectionContext {
        reference: state.reference,
        annotation: state.annotation,
        junction_support: None,
    };
    let phase_started = profile.then(Instant::now);
    let mut out =
        build_confirmed_bsj_segment_records_with_correction(state, &first_pass_correction)?;
    log_segments_profile(profile, "confirmed_bsj_first_pass", phase_started);
    let phase_started = profile.then(Instant::now);
    let junction_support = collect_junction_support(&out);
    log_segments_profile(profile, "collect_junction_support", phase_started);
    let phase_started = profile.then(Instant::now);
    let junction_support_index = build_junction_support_index(&junction_support);
    log_segments_profile(profile, "build_junction_support_index", phase_started);
    let correction = SegmentCorrectionContext {
        reference: state.reference,
        annotation: state.annotation,
        junction_support: Some(&junction_support),
    };

    let phase_started = profile.then(Instant::now);
    let replacements: Vec<(usize, SegmentRecord)> = out
        .par_iter()
        .enumerate()
        .filter_map(|(idx, record)| {
            if !segment_record_has_supported_alternative(record, &junction_support_index) {
                return None;
            }
            let read_id = record.read_id.as_str();
            let records = state.segment_groups.get(read_id)?;
            let circ_id = state.junction_read_to_circ.get(read_id)?;
            let circ = state.circ_by_id.get(circ_id)?;
            build_bsj_segment_record(
                read_id,
                records,
                circ,
                state
                    .mate_bsj_evidence
                    .get(read_id)
                    .map(Vec::as_slice)
                    .unwrap_or(&[]),
                &[],
                Some(&correction),
                state.read_len,
            )
            .map(|record| (idx, record))
        })
        .collect();
    let ambiguous_records = replacements.len();
    for (idx, record) in replacements {
        out[idx] = record;
    }
    log_segments_profile(profile, "selective_second_pass", phase_started);
    if profile {
        eprintln!(
            "[CIRI_PROFILE_SEGMENTS] selective_second_pass_records: {ambiguous_records}/{}",
            out.len()
        );
    }
    Ok(out)
}

/// Builds confirmed BSJ segment rows using one boundary-correction context.
fn build_confirmed_bsj_segment_records_with_correction(
    state: &ScanState,
    correction: &SegmentCorrectionContext<'_>,
) -> Result<Vec<SegmentRecord>> {
    let mut read_ids: Vec<&str> = state
        .junction_read_to_circ
        .keys()
        .map(String::as_str)
        .collect();
    read_ids.sort_unstable();
    let records = read_ids
        .par_iter()
        .map(|read_id| -> Result<Option<SegmentRecord>> {
            let Some(records) = state.segment_groups.get(*read_id) else {
                return Ok(None);
            };
            let Some(circ_id) = state.junction_read_to_circ.get(*read_id) else {
                return Ok(None);
            };
            let circ = state
                .circ_by_id
                .get(circ_id)
                .ok_or_else(|| anyhow!("missing circ record {}", circ_id))?;
            Ok(build_bsj_segment_record(
                read_id,
                records,
                circ,
                state
                    .mate_bsj_evidence
                    .get(*read_id)
                    .map(Vec::as_slice)
                    .unwrap_or(&[]),
                &[],
                Some(correction),
                state.read_len,
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(records.into_iter().flatten().collect::<Vec<_>>())
}

/// Converts supplemental non-BSJ candidates into `type=backward` segment rows.
///
/// These reads do not have a reliable circRNA ID at this stage. The row-level
/// topology and mate chains are still useful for downstream circRNA-level
/// assembly, but assignment to a specific circRNA remains deferred.
fn build_backward_segment_records(state: &ScanState) -> Result<Vec<SegmentRecord>> {
    let backward_ids: HashSet<&str> = state
        .candidates
        .iter()
        .map(|candidate| candidate.read_id.as_str())
        .filter(|read_id| !state.junction_read_to_circ.contains_key(*read_id))
        .collect();
    let read_strand_hints = build_read_strand_hints(state);
    let read_junction_hints = build_read_junction_hints(state);
    let correction = SegmentCorrectionContext {
        reference: state.reference,
        annotation: state.annotation,
        junction_support: None,
    };
    let mut read_ids: Vec<&str> = backward_ids.into_iter().collect();
    read_ids.sort_unstable();
    let records = read_ids
        .par_iter()
        .filter_map(|read_id| {
            let records = state.segment_groups.get(*read_id)?;
            build_backward_segment_record(
                *read_id,
                records,
                read_strand_hints.get(*read_id).copied().flatten(),
                read_junction_hints
                    .get(*read_id)
                    .map(Vec::as_slice)
                    .unwrap_or(&[]),
                Some(&correction),
                state.read_len,
            )
        })
        .collect();
    Ok(records)
}

/// Counts preliminary non-BSJ junctions for support-aware boundary ranking.
///
/// The support map is built only from post-Summary `<prefix>.segments` rows, so
/// it cannot change which reads are counted as BSJ evidence. It is used as a
/// second-pass tie breaker when multiple nearby splice boundaries have comparable
/// annotation or motif support.
fn collect_junction_support(records: &[SegmentRecord]) -> JunctionSupportMap {
    let mut support = HashMap::new();
    for record in records {
        collect_junction_support_from_segments(&record.chrom, &record.r1_segments, &mut support);
        collect_junction_support_from_segments(&record.chrom, &record.r2_segments, &mut support);
    }
    support
}

/// Adds retained read-chain `N` junctions from one mate's segment string.
fn collect_junction_support_from_segments(
    chrom: &str,
    segments: &str,
    support: &mut JunctionSupportMap,
) {
    let mut previous: Option<(i32, i32, char)> = None;
    let mut pending_bsj = false;
    for token in segments.split('|') {
        if token == "<bsj>" {
            pending_bsj = previous.is_some();
            continue;
        }
        let Some((start, end, strand)) = parse_segment_token(token) else {
            pending_bsj = false;
            continue;
        };
        if end - start + 1 < MIN_JUNCTION_SUPPORT_SEGMENT_LENGTH {
            continue;
        }
        if let Some((prev_start, prev_end, prev_strand)) = previous {
            if !pending_bsj && prev_strand == strand {
                let (site2, site1) = if prev_start <= start {
                    (prev_end, start)
                } else {
                    (end, prev_start)
                };
                *support
                    .entry(chrom.to_string())
                    .or_default()
                    .entry((site2, site1, strand))
                    .or_insert(0) += 1;
            }
        }
        previous = Some((start, end, strand));
        pending_bsj = false;
    }
}

/// Builds a sparse lookup for finding nearby supported junction alternatives.
///
/// The ranking path still needs exact `(end, start, strand)` counts, so the
/// support map above remains authoritative. This auxiliary index avoids scanning
/// every coordinate inside a correction window for every confirmed BSJ read.
fn build_junction_support_index(support: &JunctionSupportMap) -> JunctionSupportIndex {
    let mut index = HashMap::new();
    for (chrom, chrom_support) in support {
        let chrom_index = index.entry(chrom.clone()).or_insert_with(HashMap::new);
        for (&(end, start, strand), &count) in chrom_support {
            if count == 0 {
                continue;
            }
            chrom_index
                .entry(strand)
                .or_insert_with(HashMap::new)
                .entry(end)
                .or_insert_with(Vec::new)
                .push(start);
        }
    }
    index
}

/// Tests whether preliminary segments have a nearby supported splice alternative.
///
/// The first pass already applies annotation and motif correction. A support
/// aware rebuild is only useful when another preliminary BSJ read supports a
/// different internal `N` boundary near the current one. This keeps the second
/// pass focused on ambiguous internal splice sites instead of materializing all
/// confirmed BSJ reads twice.
fn segment_record_has_supported_alternative(
    record: &SegmentRecord,
    support: &JunctionSupportIndex,
) -> bool {
    let Some(chrom_support) = support.get(&record.chrom) else {
        return false;
    };
    segments_have_supported_alternative(&record.r1_segments, chrom_support)
        || segments_have_supported_alternative(&record.r2_segments, chrom_support)
}

/// Scans one mate's read-chain segments for nearby non-BSJ supported alternatives.
fn segments_have_supported_alternative(
    segments: &str,
    chrom_support: &HashMap<char, HashMap<i32, Vec<i32>>>,
) -> bool {
    let mut previous: Option<(i32, i32, char)> = None;
    let mut pending_bsj = false;
    for token in segments.split('|') {
        if token == "<bsj>" {
            pending_bsj = previous.is_some();
            continue;
        }
        let Some((start, end, strand)) = parse_segment_token(token) else {
            previous = None;
            pending_bsj = false;
            continue;
        };
        if let Some((prev_start, prev_end, prev_strand)) = previous {
            if !pending_bsj && prev_strand == strand {
                let (left_start, left_end, right_start, right_end) = if prev_start <= start {
                    (prev_start, prev_end, start, end)
                } else {
                    (start, end, prev_start, prev_end)
                };
                let left_len = left_end - left_start + 1;
                let right_len = right_end - right_start + 1;
                if supported_alternative_near_junction(
                    left_end,
                    right_start,
                    left_len,
                    right_len,
                    strand,
                    chrom_support,
                ) {
                    return true;
                }
            }
        }
        previous = Some((start, end, strand));
        pending_bsj = false;
    }
    false
}

/// Looks for a different supported junction in the correction search window.
fn supported_alternative_near_junction(
    current_end: i32,
    current_start: i32,
    left_len: i32,
    right_len: i32,
    strand: char,
    chrom_support: &HashMap<char, HashMap<i32, Vec<i32>>>,
) -> bool {
    let Some(strand_support) = chrom_support.get(&strand) else {
        return false;
    };
    let base_window = if left_len <= MIN_EXON_LENGTH || right_len <= MIN_EXON_LENGTH {
        PARTIAL_LOCAL_SPLICE_CORRECTION_WINDOW
    } else {
        INTERNAL_SPLICE_CORRECTION_WINDOW
    };
    // The preliminary pass may already have moved one side by `base_window`.
    // Doubling the lookup window keeps the selective second pass accurate while
    // materialization itself now runs in parallel, so this broader prefilter is
    // no longer the dominant cost.
    let search_window = base_window * 2;
    for end in current_end - search_window..=current_end + search_window {
        let Some(starts) = strand_support.get(&end) else {
            continue;
        };
        for &start in starts {
            if end == current_end && start == current_start {
                continue;
            }
            if start < current_start - search_window || start > current_start + search_window {
                continue;
            }
            if end < start {
                return true;
            }
        }
    }
    false
}

/// Parses one `start-end:strand` segment token.
fn parse_segment_token(token: &str) -> Option<(i32, i32, char)> {
    let (range, strand_text) = token.rsplit_once(':')?;
    let (start, end) = range.split_once('-')?;
    let strand = strand_text.chars().next()?;
    Some((start.parse().ok()?, end.parse().ok()?, strand))
}

/// Builds one `type=bsj` output row.
fn build_bsj_segment_record(
    read_id: &str,
    records: &[AsAlignment],
    circ: &CircRecord,
    mate_bsj_evidence: &[MateBsjEvidence],
    junction_hints: &[(i32, i32)],
    correction: Option<&SegmentCorrectionContext<'_>>,
    read_len: i32,
) -> Option<SegmentRecord> {
    let token_strand = circ
        .strand
        .chars()
        .next()
        .filter(|strand| matches!(strand, '+' | '-'))
        .unwrap_or('?');
    let mut chains = build_pair_chains(
        records,
        read_len,
        Some(circ),
        "bsj",
        token_strand,
        junction_hints,
        correction,
    );
    apply_mate_bsj_evidence(&mut chains, mate_bsj_evidence, circ);
    if chains[0].is_none() && chains[1].is_none() {
        return None;
    }
    let (r1_segments, r1_cigar, is_r1_bsj) = chain_text(chains[0].as_ref());
    let (r2_segments, r2_cigar, is_r2_bsj) = chain_text(chains[1].as_ref());
    Some(SegmentRecord {
        read_id: read_id.to_string(),
        type_name: "bsj",
        circ_id: circ.id.clone(),
        chrom: circ.chr.clone(),
        start: circ.start.to_string(),
        end: circ.end.to_string(),
        strand: circ.strand.clone(),
        is_circular: 1,
        is_r1_bsj,
        is_r2_bsj,
        r1_cigar,
        r1_segments,
        r2_cigar,
        r2_segments,
    })
}

/// Builds one `type=backward` output row.
fn build_backward_segment_record(
    read_id: &str,
    records: &[AsAlignment],
    inferred_strand: Option<char>,
    junction_hints: &[(i32, i32)],
    correction: Option<&SegmentCorrectionContext<'_>>,
    read_len: i32,
) -> Option<SegmentRecord> {
    let token_strand = inferred_strand.unwrap_or('?');
    let mut chains = build_pair_chains(
        records,
        read_len,
        None,
        "backward",
        token_strand,
        junction_hints,
        correction,
    );
    let is_initial_circular = usize::from(
        chains
            .iter()
            .flatten()
            .any(|chain| chain.is_circular && !chain.is_bsj),
    );
    if is_initial_circular == 0 {
        return None;
    }
    match resolve_backward_records_by_xa(
        records,
        &chains,
        read_len,
        token_strand,
        junction_hints,
        correction,
    ) {
        XaBackwardResolution::Keep => {}
        XaBackwardResolution::Reject => return None,
        XaBackwardResolution::Repair(repaired_chains) => chains = repaired_chains,
    }
    let is_circular = usize::from(
        chains
            .iter()
            .flatten()
            .any(|chain| chain.is_circular && !chain.is_bsj),
    );
    if is_circular == 0 {
        return None;
    }
    let chrom = selected_chain_chrom(&chains)?;
    let (span_start, span_end) = selected_chain_span(&chains)?;
    if span_end - span_start + 1 > BACKWARD_MAX_SPAN {
        return None;
    }
    let (r1_segments, r1_cigar, _) = chain_text(chains[0].as_ref());
    let (r2_segments, r2_cigar, _) = chain_text(chains[1].as_ref());
    Some(SegmentRecord {
        read_id: read_id.to_string(),
        type_name: "backward",
        circ_id: "NA".to_string(),
        chrom,
        start: span_start.to_string(),
        end: span_end.to_string(),
        strand: "NA".to_string(),
        is_circular,
        is_r1_bsj: 0,
        is_r2_bsj: 0,
        r1_cigar,
        r1_segments,
        r2_cigar,
        r2_segments,
    })
}

/// Result of topology-neutral XA selection for a candidate backward read.
enum XaBackwardResolution {
    Keep,
    Repair([Option<MateChain>; 2]),
    Reject,
}

/// Resolves weak backward chains by ranking seamless XA alternatives first.
///
/// BWA-MEM emits `XA` alternatives without splice-aware context. For segments
/// sidecar output we can safely compare read-slice-compatible XA hits before
/// applying the max-span guard. The ranking is topology-neutral and does not use
/// XA edit distance: linear and circular replacements compete by strand/query
/// compatibility and genomic spanning size first, then the selected best chain
/// decides whether the read is still a `type=backward` row. This keeps CIRI3
/// Scan1/Scan2 parity untouched while handling both `sim:3966179`-style
/// circular repair and `sim:432652`-style linear rejection.
fn resolve_backward_records_by_xa(
    records: &[AsAlignment],
    chains: &[Option<MateChain>; 2],
    read_len: i32,
    token_strand: char,
    junction_hints: &[(i32, i32)],
    correction: Option<&SegmentCorrectionContext<'_>>,
) -> XaBackwardResolution {
    let Some(current_rank) = xa_alignment_rank(chains) else {
        return XaBackwardResolution::Keep;
    };
    let Some((current_start, current_end)) = selected_chain_span(chains) else {
        return XaBackwardResolution::Keep;
    };
    let current_span = current_end - current_start + 1;
    let current_coverage = selected_chain_query_coverage(chains);
    let candidate_sets = xa_candidate_sets(records, chains, read_len);
    if candidate_sets.is_empty() {
        return XaBackwardResolution::Keep;
    }

    let mut best: Option<((i32, i32, i32, i32, i32), [Option<MateChain>; 2])> = None;
    let mut candidate_records = records.to_vec();
    enumerate_xa_candidate_sets(
        records,
        &candidate_sets,
        0,
        &mut candidate_records,
        read_len,
        token_strand,
        junction_hints,
        correction,
        current_span,
        current_coverage,
        current_rank,
        &mut best,
    );

    if let Some((_, best_chains)) = best {
        if chain_has_backward_topology(&best_chains) {
            XaBackwardResolution::Repair(best_chains)
        } else {
            XaBackwardResolution::Reject
        }
    } else {
        XaBackwardResolution::Keep
    }
}

/// One selected alignment plus its seamless XA replacements.
struct XaCandidateSet {
    record_index: usize,
    alternatives: Vec<AsAlignment>,
}

/// Builds XA choices from all selected MAPQ=0 records, not only supplementary ones.
///
/// Some BWA-MEM records keep the more local splice-aware placement as `XA` on a
/// primary MAPQ=0 record. Restricting the sidecar to supplementary records leaves
/// false backward rows such as `sim:4599148`; therefore every selected low-MAPQ
/// record used by the current chain can contribute exact same-slice alternatives.
fn xa_candidate_sets(
    records: &[AsAlignment],
    chains: &[Option<MateChain>; 2],
    read_len: i32,
) -> Vec<XaCandidateSet> {
    records
        .iter()
        .enumerate()
        .filter_map(|(record_index, record)| {
            if !ambiguous_alignment_with_xa(record, read_len)
                || !chains_use_record_as_backward(chains, record, read_len)
            {
                return None;
            }
            let alternatives: Vec<AsAlignment> = record
                .xa_alternatives
                .iter()
                .filter(|alternative| xa_replacement_is_seamless(record, alternative, read_len))
                .take(XA_MAX_ALTERNATIVES_PER_RECORD)
                .map(|alternative| alignment_from_xa(record, alternative))
                .collect();
            (!alternatives.is_empty()).then_some(XaCandidateSet {
                record_index,
                alternatives,
            })
        })
        .collect()
}

/// Enumerates joint XA replacements across selected records.
///
/// Greedy single-record replacement misses cases where R1 and R2 each carry one
/// remote equivalent hit: replacing either record alone does not reduce the
/// pair-level span, but replacing both reveals a linear local interpretation.
fn enumerate_xa_candidate_sets(
    original_records: &[AsAlignment],
    candidate_sets: &[XaCandidateSet],
    set_idx: usize,
    candidate_records: &mut Vec<AsAlignment>,
    read_len: i32,
    token_strand: char,
    junction_hints: &[(i32, i32)],
    correction: Option<&SegmentCorrectionContext<'_>>,
    current_span: i32,
    current_coverage: i32,
    current_rank: (i32, i32, i32, i32, i32),
    best: &mut Option<((i32, i32, i32, i32, i32), [Option<MateChain>; 2])>,
) {
    if set_idx == candidate_sets.len() {
        if candidate_records
            .iter()
            .zip(original_records.iter())
            .all(|(candidate, original)| {
                candidate.pos == original.pos
                    && candidate.flag == original.flag
                    && candidate.cigar == original.cigar
                    && candidate.chr == original.chr
            })
        {
            return;
        }
        let candidate_chains = build_pair_chains(
            candidate_records,
            read_len,
            None,
            "backward",
            token_strand,
            junction_hints,
            correction,
        );
        if selected_chain_query_coverage(&candidate_chains)
            < current_coverage - XA_SEAMLESS_QUERY_TOLERANCE
        {
            return;
        }
        let Some((candidate_start, candidate_end)) = selected_chain_span(&candidate_chains) else {
            return;
        };
        let candidate_span = candidate_end - candidate_start + 1;
        if current_span - candidate_span < XA_REJECT_MIN_SPAN_REDUCTION {
            return;
        }
        let Some(rank) = xa_alignment_rank(&candidate_chains) else {
            return;
        };
        if rank <= current_rank {
            return;
        }
        if best.as_ref().is_none_or(|(best_rank, _)| rank > *best_rank) {
            *best = Some((rank, candidate_chains));
        }
        return;
    }

    let candidate_set = &candidate_sets[set_idx];
    enumerate_xa_candidate_sets(
        original_records,
        candidate_sets,
        set_idx + 1,
        candidate_records,
        read_len,
        token_strand,
        junction_hints,
        correction,
        current_span,
        current_coverage,
        current_rank,
        best,
    );
    let original = candidate_records[candidate_set.record_index].clone();
    for alternative in &candidate_set.alternatives {
        candidate_records[candidate_set.record_index] = alternative.clone();
        enumerate_xa_candidate_sets(
            original_records,
            candidate_sets,
            set_idx + 1,
            candidate_records,
            read_len,
            token_strand,
            junction_hints,
            correction,
            current_span,
            current_coverage,
            current_rank,
            best,
        );
    }
    candidate_records[candidate_set.record_index] = original;
}

/// Ranks one XA-derived alignment interpretation without using topology.
///
/// Exact same-slice XA hits are alternative mapper placements, so the sidecar
/// first asks which placement is the most compatible with the read pair's
/// genomic span. Topology is intentionally excluded here; the caller handles a
/// linear best chain by dropping the `type=backward` row.
fn xa_alignment_rank(chains: &[Option<MateChain>; 2]) -> Option<(i32, i32, i32, i32, i32)> {
    selected_chain_chrom(chains)?;
    let (start, end) = selected_chain_span(chains)?;
    let span = end - start + 1;
    let used_supplementary: usize = chains
        .iter()
        .flatten()
        .map(|chain| chain.used_supplementary)
        .sum();
    let used_secondary: usize = chains
        .iter()
        .flatten()
        .map(|chain| chain.used_secondary)
        .sum();
    Some((
        i32::from(span <= BACKWARD_MAX_SPAN),
        -span,
        selected_chain_query_coverage(chains),
        -(used_supplementary as i32),
        -(used_secondary as i32),
    ))
}

/// Returns whether any selected mate chain still supports backward topology.
fn chain_has_backward_topology(chains: &[Option<MateChain>; 2]) -> bool {
    chains
        .iter()
        .flatten()
        .any(|chain| chain.is_circular && !chain.is_bsj)
}

/// Returns whether the selected backward chain currently uses this alignment.
fn chains_use_record_as_backward(
    chains: &[Option<MateChain>; 2],
    record: &AsAlignment,
    read_len: i32,
) -> bool {
    if !chains
        .iter()
        .flatten()
        .any(|chain| chain.is_circular && !chain.is_bsj)
    {
        return false;
    }
    chains
        .iter()
        .flatten()
        .any(|chain| chain_uses_record(chain, record, read_len))
}

/// Identifies ambiguous alignments that should not define backward alone.
///
/// The CIRI3 parity path does not inspect XA tags, but this sidecar resolver
/// uses them to choose a better read-level chain. We require a minimum mapped
/// anchor so tiny clips do not drive the decision, and deliberately avoid a
/// maximum anchor length because exact XA alternatives can also explain longer
/// repetitive supplementary blocks.
fn ambiguous_alignment_with_xa(record: &AsAlignment, read_len: i32) -> bool {
    if record.xa_alternatives.is_empty() || record.from_local_clip || record.mapq != 0 {
        return false;
    }
    alignment_mapped_len(record, read_len).is_some_and(|len| len >= XA_REJECT_MIN_ANCHOR_LEN)
}

/// Tests whether an XA hit can replace the same read slice as one SAM record.
///
/// The replacement is intentionally strict: it must stay on the same alignment
/// strand, parse into the same number of blocks, and cover nearly the same
/// query coordinates. BWA already marked these as alternatives, but this
/// read-slice check keeps the sidecar from swapping unrelated clipped pieces.
fn xa_replacement_is_seamless(
    record: &AsAlignment,
    alternative: &XaAlternative,
    read_len: i32,
) -> bool {
    if alternative.strand != strand_char(record.flag) {
        return false;
    }
    let alt_record = alignment_from_xa(record, alternative);
    let Some(source_blocks) = parse_alignment_blocks(record, read_len) else {
        return false;
    };
    let Some(alt_blocks) = parse_alignment_blocks(&alt_record, read_len) else {
        return false;
    };
    if source_blocks.len() != alt_blocks.len() {
        return false;
    }
    source_blocks
        .iter()
        .zip(alt_blocks.iter())
        .all(|(source, alt)| {
            (source.read_start - alt.read_start).abs() <= XA_SEAMLESS_QUERY_TOLERANCE
                && (source.read_end - alt.read_end).abs() <= XA_SEAMLESS_QUERY_TOLERANCE
                && (source.read_end - source.read_start + 1) >= XA_REJECT_MIN_ANCHOR_LEN
                && (alt.read_end - alt.read_start + 1) >= XA_REJECT_MIN_ANCHOR_LEN
        })
}

/// Returns whether a selected chain contains one block from the source record.
fn chain_uses_record(chain: &MateChain, record: &AsAlignment, read_len: i32) -> bool {
    let Some(blocks) = parse_alignment_blocks(record, read_len) else {
        return false;
    };
    blocks.iter().any(|source| {
        chain.blocks.iter().any(|selected| {
            source.read_start == selected.read_start
                && source.read_end == selected.read_end
                && source.ref_start <= selected.ref_end
                && selected.ref_start <= source.ref_end
        })
    })
}

/// Returns total selected query coverage across both mates.
fn selected_chain_query_coverage(chains: &[Option<MateChain>; 2]) -> i32 {
    chains
        .iter()
        .flatten()
        .map(|chain| chain.query_coverage)
        .sum()
}

/// Returns the aligned query length represented by one record.
fn alignment_mapped_len(record: &AsAlignment, read_len: i32) -> Option<i32> {
    parse_alignment_blocks(record, read_len).map(|blocks| {
        blocks
            .iter()
            .map(|block| block.read_end - block.read_start + 1)
            .sum()
    })
}

/// Converts a BWA XA hit into an alternative alignment for sidecar ranking.
fn alignment_from_xa(record: &AsAlignment, alternative: &XaAlternative) -> AsAlignment {
    let mut flag = record.flag;
    if alternative.strand == '-' {
        flag |= 0x10;
    } else {
        flag &= !0x10;
    }
    AsAlignment {
        flag,
        chr: alternative.chr.clone(),
        pos: alternative.pos,
        mapq: 0,
        cigar: alternative.cigar.clone(),
        seq: record.seq.clone(),
        from_local_clip: false,
        xa_alternatives: Vec::new(),
    }
}

/// Sorts final read-level segment rows by their graph-construction coordinates.
///
/// Each row already contains both mates, so row-level genomic sorting does not
/// break read-pair phasing. The ordering makes later circRNA-level graph
/// construction stream mostly contiguous BSJ/backward evidence instead of
/// holding every read in a read-ID keyed map.
fn sort_segment_records(records: &mut [SegmentRecord]) {
    records.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then_with(|| sort_position(&a.start).cmp(&sort_position(&b.start)))
            .then_with(|| sort_position(&a.end).cmp(&sort_position(&b.end)))
            .then_with(|| type_sort_rank(a.type_name).cmp(&type_sort_rank(b.type_name)))
            .then_with(|| a.circ_id.cmp(&b.circ_id))
            .then_with(|| a.read_id.cmp(&b.read_id))
    });
}

/// Parses a segment row coordinate for stable output ordering.
fn sort_position(value: &str) -> i32 {
    value.parse::<i32>().unwrap_or(i32::MAX)
}

/// Keeps confirmed BSJ rows ahead of supplemental rows at identical loci.
fn type_sort_rank(type_name: &str) -> u8 {
    match type_name {
        "bsj" => 0,
        "backward" => 1,
        "forward" => 2,
        _ => 3,
    }
}

/// Returns the single chromosome shared by selected mate chains.
///
/// Backward rows have no circRNA assignment yet, but downstream candidate
/// clustering still needs chromosome context. The row-level schema has one
/// `chrom` column, so mixed-chromosome mate chains are filtered out instead of
/// being reported as circular segments with ambiguous reference context.
fn selected_chain_chrom(chains: &[Option<MateChain>; 2]) -> Option<String> {
    let mut chrom: Option<&str> = None;
    for chain in chains.iter().flatten() {
        match chrom {
            Some(existing) if existing != chain.chrom => return None,
            Some(_) => {}
            None => chrom = Some(chain.chrom.as_str()),
        }
    }
    chrom.map(str::to_string)
}

/// Returns the genomic span covered by all selected mate-chain blocks.
///
/// Backward rows are not assigned to a candidate BSJ at read-level output time:
/// one backward read can be compatible with multiple later BSJ/circRNA groups.
/// The sortable coordinates therefore describe the retained alignment span
/// across the selected R1/R2 chains, not a putative circRNA boundary.
fn selected_chain_span(chains: &[Option<MateChain>; 2]) -> Option<(i32, i32)> {
    let mut start = i32::MAX;
    let mut end = i32::MIN;
    for block in chains.iter().flatten().flat_map(|chain| &chain.blocks) {
        start = start.min(block.ref_start);
        end = end.max(block.ref_end);
    }
    (start <= end).then_some((start, end))
}

/// Selects one best chain for each mate of a read pair.
fn build_pair_chains(
    records: &[AsAlignment],
    read_len: i32,
    circ: Option<&CircRecord>,
    type_name: &str,
    token_strand: char,
    junction_hints: &[(i32, i32)],
    correction: Option<&SegmentCorrectionContext<'_>>,
) -> [Option<MateChain>; 2] {
    let mut buckets: [Vec<&AsAlignment>; 2] = [Vec::new(), Vec::new()];
    for record in records {
        buckets[mate_bucket(record.flag)].push(record);
    }
    [
        select_best_chain_for_mate(
            &buckets[0],
            read_len,
            circ,
            type_name,
            token_strand,
            junction_hints,
            correction,
            false,
        ),
        select_best_chain_for_mate(
            &buckets[1],
            read_len,
            circ,
            type_name,
            token_strand,
            junction_hints,
            correction,
            true,
        ),
    ]
}

/// Chooses the best chain for one mate.
///
/// Primary plus supplementary alignments are always considered first. Secondary
/// alignments are only allowed to win if they produce a better circRNA-consistent
/// chain than the primary/supplementary pool.
fn select_best_chain_for_mate(
    records: &[&AsAlignment],
    read_len: i32,
    circ: Option<&CircRecord>,
    type_name: &str,
    token_strand: char,
    junction_hints: &[(i32, i32)],
    correction: Option<&SegmentCorrectionContext<'_>>,
    reverse_chain_order: bool,
) -> Option<MateChain> {
    if records.is_empty() {
        return None;
    }
    let non_secondary: Vec<&AsAlignment> = records
        .iter()
        .copied()
        .filter(|record| !is_secondary(record.flag))
        .collect();
    let primary_chain = build_chain_from_pool(
        &non_secondary,
        read_len,
        circ,
        token_strand,
        junction_hints,
        correction,
        reverse_chain_order,
    );
    let fallback_chain = build_chain_from_pool(
        records,
        read_len,
        circ,
        token_strand,
        junction_hints,
        correction,
        reverse_chain_order,
    );
    better_chain(primary_chain, fallback_chain, type_name, circ)
}

/// Picks the better of the primary/supplementary and fallback chain candidates.
fn better_chain(
    left: Option<MateChain>,
    right: Option<MateChain>,
    type_name: &str,
    circ: Option<&CircRecord>,
) -> Option<MateChain> {
    match (left, right) {
        (Some(left), Some(right)) => {
            if chain_rank(&right, type_name, circ) > chain_rank(&left, type_name, circ) {
                Some(right)
            } else {
                Some(left)
            }
        }
        (Some(left), None) => Some(left),
        (None, Some(right)) => Some(right),
        (None, None) => None,
    }
}

/// Builds one best-effort chain from a pool of candidate alignments.
fn build_chain_from_pool(
    records: &[&AsAlignment],
    read_len: i32,
    circ: Option<&CircRecord>,
    token_strand: char,
    junction_hints: &[(i32, i32)],
    correction: Option<&SegmentCorrectionContext<'_>>,
    reverse_chain_order: bool,
) -> Option<MateChain> {
    let mut by_key: BTreeMap<(String, char), Vec<ParsedAlignment>> = BTreeMap::new();
    for record in records {
        let Some(blocks) = parse_alignment_blocks(record, read_len) else {
            continue;
        };
        if blocks.is_empty() || record.chr == "*" {
            continue;
        }
        let strand = strand_char(record.flag);
        by_key
            .entry((record.chr.clone(), strand))
            .or_default()
            .push(ParsedAlignment {
                flag: record.flag,
                chrom: record.chr.clone(),
                strand,
                mapq: record.mapq,
                blocks,
            });
    }

    let mut best: Option<MateChain> = None;
    for ((_chrom, _strand), mut group) in by_key {
        group.sort_by_key(|record| {
            (
                record
                    .blocks
                    .first()
                    .map(|block| block.read_start)
                    .unwrap_or(i32::MAX),
                usize::from(is_secondary(record.flag)),
                usize::from(is_supplementary(record.flag)),
                -record.mapq,
            )
        });

        let mut selected: Vec<ParsedAlignment> = Vec::new();
        for record in group {
            if let Some(last) = selected.last_mut() {
                if query_overlap(last, &record) > 6 {
                    if parsed_alignment_rank(&record) > parsed_alignment_rank(last) {
                        *last = record;
                    }
                    continue;
                }
            }
            selected.push(record);
        }

        let chain = materialize_chain(
            &selected,
            read_len,
            circ,
            token_strand,
            junction_hints,
            correction,
            reverse_chain_order,
        )?;
        if best.as_ref().is_none_or(|current| {
            chain_generic_rank(&chain, circ) > chain_generic_rank(current, circ)
        }) {
            best = Some(chain);
        }
    }
    best
}

/// Parses one alignment into read-order segment blocks.
fn parse_alignment_blocks(record: &AsAlignment, read_len: i32) -> Option<Vec<SegmentBlock>> {
    if record.cigar == "*" || record.cigar.is_empty() || record.chr == "*" {
        return None;
    }
    let mut ref_pos = record.pos;
    let mut read_pos = 1;
    let mut blocks = Vec::new();
    let mut current: Option<SegmentBlock> = None;
    let mut count = 0i32;
    let mut has_count = false;
    for mut op in record.cigar.chars() {
        if op.is_ascii_digit() {
            count = count
                .checked_mul(10)?
                .checked_add(op.to_digit(10)? as i32)?;
            has_count = true;
            continue;
        }
        if !has_count {
            return None;
        }
        if op == 'H' {
            op = 'S';
        }
        match op {
            'M' | '=' | 'X' => {
                let block = current.get_or_insert(SegmentBlock {
                    read_start: read_pos,
                    read_end: read_pos + count - 1,
                    ref_start: ref_pos,
                    ref_end: ref_pos + count - 1,
                    from_local_clip: record.from_local_clip,
                });
                block.read_end = read_pos + count - 1;
                block.ref_end = ref_pos + count - 1;
                read_pos += count;
                ref_pos += count;
            }
            'I' => {
                if let Some(block) = current.as_mut() {
                    block.read_end += count;
                }
                read_pos += count;
            }
            'D' => {
                if let Some(block) = current.as_mut() {
                    block.ref_end += count;
                }
                ref_pos += count;
            }
            'N' => {
                if let Some(block) = current.take() {
                    blocks.push(block);
                }
                ref_pos += count;
            }
            'S' | 'H' => {
                if let Some(block) = current.take() {
                    blocks.push(block);
                }
                read_pos += count;
            }
            'P' => {}
            _ => return None,
        }
        count = 0;
        has_count = false;
    }
    if has_count {
        return None;
    }
    if let Some(block) = current.take() {
        blocks.push(block);
    }
    if blocks.is_empty() {
        return None;
    }
    if record.flag & 0x10 != 0 {
        for block in &mut blocks {
            let flipped_start = read_len - block.read_end + 1;
            let flipped_end = read_len - block.read_start + 1;
            block.read_start = flipped_start;
            block.read_end = flipped_end;
        }
        blocks.reverse();
    }
    Some(blocks)
}

/// Builds the final output chain and segment CIGAR from selected alignments.
///
/// The final `<prefix>.segments` fields are emitted in read-chain order, not
/// genomic order. For BSJ reads this is the only representation that preserves
/// circular topology: a plus-strand junction is `circ-end-side -> B ->
/// circ-start-side`, while a minus-strand junction is the reverse genomic
/// direction after the reverse-strand alignment has been normalized into read
/// order.
fn materialize_chain(
    records: &[ParsedAlignment],
    read_len: i32,
    circ: Option<&CircRecord>,
    token_strand: char,
    junction_hints: &[(i32, i32)],
    correction: Option<&SegmentCorrectionContext<'_>>,
    reverse_chain_order: bool,
) -> Option<MateChain> {
    if records.is_empty() {
        return None;
    }
    let chrom = records.first()?.chrom.clone();
    let mut order_strand = records.first()?.strand;
    let mut blocks = Vec::new();
    let mut used_secondary = 0usize;
    let mut used_supplementary = 0usize;
    let mut query_coverage = 0i32;
    for record in records {
        used_secondary += usize::from(is_secondary(record.flag));
        used_supplementary += usize::from(is_supplementary(record.flag));
        for block in &record.blocks {
            query_coverage += block.read_end - block.read_start + 1;
            blocks.push(block.clone());
        }
    }
    blocks.sort_by_key(|block| (block.read_start, block.read_end, block.ref_start));
    if reverse_chain_order {
        blocks.reverse();
        order_strand = opposite_strand(order_strand);
    }
    apply_circ_boundary_corrections(&mut blocks, circ);
    let mut is_bsj = false;
    let mut is_circular = false;
    for pair in blocks.windows(2) {
        let prev_block = &pair[0];
        let next_block = &pair[1];
        let wrapped = wraps_in_read_order(prev_block, next_block, order_strand);
        if let Some(circ) = circ {
            if is_bsj_transition(prev_block, next_block, circ, order_strand) {
                is_bsj = true;
                is_circular = true;
            } else if wrapped {
                is_circular = true;
            }
        } else if wrapped {
            is_circular = true;
        }
    }

    let mut output_blocks = blocks.clone();
    let bsj_gap_idx = if is_bsj {
        circ.and_then(|circ| read_order_bsj_gap_index(&output_blocks, circ, order_strand))
    } else {
        None
    };
    apply_segment_boundary_corrections(
        &mut output_blocks,
        &chrom,
        bsj_gap_idx,
        junction_hints,
        correction,
        token_strand,
    );
    let (tokens, token_spans, cigar) = materialize_read_chain_output(
        &output_blocks,
        token_strand,
        bsj_gap_idx,
        read_len,
        reverse_chain_order,
    );

    Some(MateChain {
        chrom,
        order_strand,
        token_strand,
        blocks: output_blocks,
        token_spans,
        tokens,
        cigar,
        is_bsj,
        is_circular,
        used_secondary,
        used_supplementary,
        query_coverage,
    })
}

/// Builds read-chain segment tokens and the matching CIRI-specific CIGAR.
///
/// `N` and `B` both encode the genomic interval skipped between adjacent
/// read-chain blocks. `B` is deliberately tied to the read-order circular wrap;
/// re-sorting these blocks by coordinate would turn `C|B|A` circRNA evidence
/// into a linear-looking `A|B|C` chain.
fn materialize_read_chain_output(
    blocks: &[SegmentBlock],
    token_strand: char,
    bsj_gap_idx: Option<usize>,
    read_len: i32,
    reverse_chain_order: bool,
) -> (Vec<String>, Vec<(i32, i32)>, String) {
    let mut token_spans = Vec::with_capacity(blocks.len());
    let mut tokens = Vec::with_capacity(blocks.len() + usize::from(bsj_gap_idx.is_some()));
    let mut cigar = String::new();
    if let Some(first) = blocks.first() {
        let leading_clip = if reverse_chain_order {
            read_len - first.read_end
        } else {
            first.read_start - 1
        };
        if leading_clip > 0 {
            let _ = write!(&mut cigar, "{}S", leading_clip);
        }
    }
    for (idx, block) in blocks.iter().enumerate() {
        if idx > 0 {
            let prev = &blocks[idx - 1];
            let gap = interval_gap(prev, block);
            let op = if bsj_gap_idx == Some(idx) { 'B' } else { 'N' };
            let _ = write!(&mut cigar, "{}{}", gap, op);
            if op == 'B' {
                tokens.push("<bsj>".to_string());
            }
        }
        let len = block.ref_end - block.ref_start + 1;
        let _ = write!(&mut cigar, "{}M", len.max(0));
        tokens.push(format!(
            "{}-{}:{}",
            block.ref_start, block.ref_end, token_strand
        ));
        token_spans.push((block.ref_start, block.ref_end));
    }
    if let Some(last) = blocks.last() {
        let trailing_clip = if reverse_chain_order {
            last.read_start - 1
        } else {
            read_len - last.read_end
        };
        if trailing_clip > 0 {
            let _ = write!(&mut cigar, "{}S", trailing_clip);
        }
    }
    (tokens, token_spans, cigar)
}

/// Returns the read-chain gap index where the selected blocks cross BSJ.
fn read_order_bsj_gap_index(
    read_blocks: &[SegmentBlock],
    circ: &CircRecord,
    order_strand: char,
) -> Option<usize> {
    read_blocks.windows(2).enumerate().find_map(|(idx, pair)| {
        is_bsj_transition(&pair[0], &pair[1], circ, order_strand).then_some(idx + 1)
    })
}

/// Returns the skipped genomic distance between two read-chain blocks.
fn interval_gap(left: &SegmentBlock, right: &SegmentBlock) -> i32 {
    if left.ref_end < right.ref_start {
        right.ref_start - left.ref_end - 1
    } else if right.ref_end < left.ref_start {
        left.ref_start - right.ref_end - 1
    } else {
        0
    }
}

/// Snaps read-order blocks onto confirmed circ outer boundaries.
///
/// This runs before topology detection because BSJ transitions are judged in
/// read order. Internal splice correction is intentionally deferred until the
/// final read-chain blocks are selected, then each adjacent non-BSJ pair is
/// interpreted in genomic low/high order only for choosing splice boundaries.
fn apply_circ_boundary_corrections(blocks: &mut [SegmentBlock], circ: Option<&CircRecord>) {
    if let Some(circ) = circ {
        for block in blocks.iter_mut() {
            if (block.ref_start - circ.start).abs() <= 2 {
                block.ref_start = circ.start;
            }
            if (block.ref_end - circ.end).abs() <= 2 {
                block.ref_end = circ.end;
            }
        }
    }
}

/// Applies corrected internal splice boundaries to read-chain output blocks.
///
/// The correction order mirrors the current CIRI-AS design: annotation-supported
/// exon end/start pairs win first, read-specific validated junction hints are
/// used next, and de novo splice motifs are a final fallback. Only the
/// `<prefix>.segments` representation is changed; CIRI3 Summary evidence has
/// already been finalized before this function runs.
fn apply_segment_boundary_corrections(
    blocks: &mut [SegmentBlock],
    chrom: &str,
    bsj_gap_idx: Option<usize>,
    junction_hints: &[(i32, i32)],
    correction: Option<&SegmentCorrectionContext<'_>>,
    token_strand: char,
) {
    if blocks.len() < 2 {
        return;
    }
    for idx in 0..blocks.len() - 1 {
        if bsj_gap_idx == Some(idx + 1) {
            continue;
        }
        let (left_idx, right_idx) = if blocks[idx].ref_start <= blocks[idx + 1].ref_start {
            (idx, idx + 1)
        } else {
            (idx + 1, idx)
        };
        let prev_end = blocks[left_idx].ref_end;
        let next_start = blocks[right_idx].ref_start;
        if prev_end >= next_start {
            continue;
        }
        let annotation_window = correction_window_for_pair(&blocks[left_idx], &blocks[right_idx]);
        let chosen = if let Some(ctx) = correction {
            choose_supported_splice_boundary(
                ctx,
                chrom,
                prev_end,
                next_start,
                &blocks[left_idx],
                &blocks[right_idx],
                junction_hints,
                token_strand,
                annotation_window,
            )
        } else {
            choose_hinted_splice_boundary(
                prev_end,
                next_start,
                &blocks[left_idx],
                &blocks[right_idx],
                junction_hints,
                annotation_window,
            )
        };
        if let Some((site2, site1)) = chosen {
            blocks[left_idx].ref_end = site2;
            blocks[right_idx].ref_start = site1;
        }
    }
}

/// Chooses the best nearby internal splice boundary using all available signals.
///
/// The ranking is deliberately sidecar-local: read-specific CIRI-AS validation
/// hints win first, then transcript-consistent annotation, then preliminary
/// read-level junction support, then boundary-level annotation and splice motif.
/// BSJ gaps have already been skipped by the caller, so this function cannot
/// rewrite circRNA outer-boundary topology.
fn choose_supported_splice_boundary(
    ctx: &SegmentCorrectionContext<'_>,
    chrom: &str,
    prev_end: i32,
    next_start: i32,
    left: &SegmentBlock,
    right: &SegmentBlock,
    junction_hints: &[(i32, i32)],
    token_strand: char,
    window: i32,
) -> Option<(i32, i32)> {
    let chr_seq = ctx.reference.get(chrom);
    let mut best: Option<((i32, i32, i32, i32, i32, i32, i32, i32), i32, i32)> = None;
    for end in left.ref_end - window..=left.ref_end + window {
        if end < left.ref_start || end >= right.ref_start {
            continue;
        }
        for start in right.ref_start - window..=right.ref_start + window {
            if start > right.ref_end || end >= start {
                continue;
            }
            let hint_score = junction_hints
                .iter()
                .filter(|&&(site2, site1)| site2 == end && site1 == start)
                .count() as i32;
            let (transcript_score, annotation_score) =
                annotation_splice_pair_score(ctx.annotation, chrom, end, start, token_strand)
                    .unwrap_or((0, 0));
            let support_score = ctx
                .junction_support
                .and_then(|support| support.get(chrom))
                .and_then(|support| support.get(&(end, start, token_strand)).copied())
                .unwrap_or(0)
                .min(i32::MAX as usize) as i32;
            let motif_score = if (end - left.ref_end).abs() <= INTERNAL_SPLICE_CORRECTION_WINDOW
                && (start - right.ref_start).abs() <= INTERNAL_SPLICE_CORRECTION_WINDOW
            {
                chr_seq
                    .and_then(|seq| splice_motif_score(seq, end, start, token_strand))
                    .unwrap_or(0)
            } else {
                0
            };
            if hint_score == 0
                && transcript_score == 0
                && annotation_score == 0
                && support_score == 0
                && motif_score == 0
            {
                continue;
            }
            let movement = (end - prev_end).abs() + (start - next_start).abs();
            let key = (
                hint_score,
                transcript_score,
                support_score,
                annotation_score,
                motif_score,
                -movement,
                -(end.abs_diff(prev_end) as i32),
                -(start.abs_diff(next_start) as i32),
            );
            if best.as_ref().is_none_or(|(current, _, _)| key > *current) {
                best = Some((key, end, start));
            }
        }
    }
    best.map(|(_, end, start)| (end, start))
}

/// Chooses a nearby annotated exon end/start pair for one internal junction.
#[allow(dead_code)]
fn choose_annotation_splice_boundary(
    annotation: Option<&Annotation>,
    chrom: &str,
    left: &SegmentBlock,
    right: &SegmentBlock,
    token_strand: char,
    window: i32,
) -> Option<(i32, i32)> {
    let mut best: Option<(i32, i32, i32, i32)> = None;
    for end in left.ref_end - window..=left.ref_end + window {
        if end < left.ref_start || end >= right.ref_start {
            continue;
        }
        for start in right.ref_start - window..=right.ref_start + window {
            if start > right.ref_end || end >= start {
                continue;
            }
            let Some((transcript_score, annotation_score)) =
                annotation_splice_pair_score(annotation, chrom, end, start, token_strand)
            else {
                continue;
            };
            let movement = (end - left.ref_end).abs() + (start - right.ref_start).abs();
            let score = transcript_score * 1000 + annotation_score;
            let key = (
                score,
                -movement,
                -(end.abs_diff(left.ref_end) as i32),
                -(start.abs_diff(right.ref_start) as i32),
            );
            if best.is_none_or(|current| key > current) {
                best = Some((key.0, key.1, end, start));
            }
        }
    }
    best.map(|(_, _, end, start)| (end, start))
}

/// Scores one annotation-supported internal splice boundary candidate.
fn annotation_splice_pair_score(
    annotation: Option<&Annotation>,
    chrom: &str,
    end: i32,
    start: i32,
    token_strand: char,
) -> Option<(i32, i32)> {
    let annotation = annotation?;
    let transcript_score =
        transcript_splice_pair_score(annotation, chrom, end, start, token_strand);
    let boundary_score = boundary_splice_pair_score(annotation, chrom, end, start, token_strand);
    match (transcript_score, boundary_score) {
        (0, None) => None,
        (score, Some(boundary)) => Some((score, boundary)),
        (score, None) => Some((score, 0)),
    }
}

/// Scores whether one intron is supported by adjacent exons from the same transcript.
fn transcript_splice_pair_score(
    annotation: &Annotation,
    chrom: &str,
    end: i32,
    start: i32,
    token_strand: char,
) -> i32 {
    let has_pair = |strand: char| {
        transcript_splice_pair_score_for_strand(annotation, chrom, end, start, strand)
    };
    match token_strand {
        '+' | '-' => i32::from(has_pair(token_strand)),
        '?' => i32::from(has_pair('+') || has_pair('-')),
        _ => 0,
    }
}

/// Scores one transcript splice-pair strand without allocating keys when possible.
fn transcript_splice_pair_score_for_strand(
    annotation: &Annotation,
    chrom: &str,
    end: i32,
    start: i32,
    strand: char,
) -> bool {
    if let Some(by_strand) = annotation.transcript_splice_index.get(chrom) {
        return by_strand
            .get(&strand)
            .is_some_and(|sites| sites.contains(&(end, start)));
    }
    annotation
        .transcript_splice_map
        .contains(&format!("{}\t{}\t{}\t{}", chrom, end, start, strand))
}

/// Scores boundary-level exon end/start support without requiring transcript identity.
fn boundary_splice_pair_score(
    annotation: &Annotation,
    chrom: &str,
    end: i32,
    start: i32,
    token_strand: char,
) -> Option<i32> {
    if annotation.chr_exon_end_index.contains_key(chrom)
        || annotation.chr_exon_start_index.contains_key(chrom)
    {
        return typed_boundary_splice_pair_score(annotation, chrom, end, start, token_strand);
    }
    let end_key = format!("{}\t{}", chrom, end);
    let start_key = format!("{}\t{}", chrom, start);
    let end_hit = annotation.chr_exon_end_map.get(&end_key)?;
    let start_hit = annotation.chr_exon_start_map.get(&start_key)?;
    let (end_gene, end_strand) = split_gene_strand(end_hit);
    let (start_gene, start_strand) = split_gene_strand(start_hit);
    let same_gene = end_gene == start_gene;
    if token_strand != '?' && (end_strand != token_strand || start_strand != token_strand) {
        return None;
    }
    Some(if same_gene { 100 } else { 80 })
}

/// Scores exon end/start support from numeric annotation indexes.
fn typed_boundary_splice_pair_score(
    annotation: &Annotation,
    chrom: &str,
    end: i32,
    start: i32,
    token_strand: char,
) -> Option<i32> {
    let end_hit = annotation.chr_exon_end_index.get(chrom)?.get(&end)?;
    let start_hit = annotation.chr_exon_start_index.get(chrom)?.get(&start)?;
    if token_strand != '?' && (end_hit.1 != token_strand || start_hit.1 != token_strand) {
        return None;
    }
    Some(if end_hit.0 == start_hit.0 { 100 } else { 80 })
}

/// Splits the legacy `"gene\tstrand"` annotation payload without allocating.
fn split_gene_strand(value: &str) -> (&str, char) {
    let (gene, strand) = value.split_once('\t').unwrap_or((value, ""));
    (gene, strand.chars().next().unwrap_or('?'))
}

/// Returns the annotation/hint correction window for one internal splice pair.
///
/// Normal mapper blocks keep the conservative CIRI-AS-sized window. Partial
/// local clip blocks and very short retained blocks can sit farther from the
/// true annotated splice boundary, so annotation and validated read-specific
/// hints may use a wider window. De novo motif correction intentionally remains
/// conservative to avoid overfitting microhomology around BSJ-adjacent gaps.
fn correction_window_for_pair(left: &SegmentBlock, right: &SegmentBlock) -> i32 {
    let has_partial_clip_like_block = [left, right].iter().any(|block| {
        block.from_local_clip || block.ref_end - block.ref_start + 1 <= MIN_EXON_LENGTH
    });
    if has_partial_clip_like_block {
        PARTIAL_LOCAL_SPLICE_CORRECTION_WINDOW
    } else {
        INTERNAL_SPLICE_CORRECTION_WINDOW
    }
}

/// Reuses read-specific splice hints produced by CIRI-AS candidate validation.
fn choose_hinted_splice_boundary(
    prev_end: i32,
    next_start: i32,
    left: &SegmentBlock,
    right: &SegmentBlock,
    junction_hints: &[(i32, i32)],
    window: i32,
) -> Option<(i32, i32)> {
    junction_hints
        .iter()
        .copied()
        .filter(|&(site2, site1)| {
            site2 >= left.ref_start
                && site2 < site1
                && site1 <= right.ref_end
                && (prev_end - site2).abs() <= window
                && (next_start - site1).abs() <= window
        })
        .min_by_key(|&(site2, site1)| (prev_end - site2).abs() + (next_start - site1).abs())
}

/// Chooses a nearby canonical splice signal when annotation gives no answer.
fn choose_motif_splice_boundary(
    reference: &HashMap<String, String>,
    chrom: &str,
    left: &SegmentBlock,
    right: &SegmentBlock,
    token_strand: char,
) -> Option<(i32, i32)> {
    let chr_seq = reference.get(chrom)?;
    let mut best: Option<(i32, i32, i32, i32)> = None;
    for end in left.ref_end - INTERNAL_SPLICE_CORRECTION_WINDOW
        ..=left.ref_end + INTERNAL_SPLICE_CORRECTION_WINDOW
    {
        if end < left.ref_start || end >= right.ref_start {
            continue;
        }
        for start in right.ref_start - INTERNAL_SPLICE_CORRECTION_WINDOW
            ..=right.ref_start + INTERNAL_SPLICE_CORRECTION_WINDOW
        {
            if start > right.ref_end || end >= start {
                continue;
            }
            let Some(score) = splice_motif_score(chr_seq, end, start, token_strand) else {
                continue;
            };
            let movement = (end - left.ref_end).abs() + (start - right.ref_start).abs();
            let key = (
                score,
                -movement,
                -(end.abs_diff(left.ref_end) as i32),
                -(start.abs_diff(right.ref_start) as i32),
            );
            if best.is_none_or(|current| key > current) {
                best = Some((key.0, key.1, end, start));
            }
        }
    }
    best.map(|(_, _, end, start)| (end, start))
}

/// Scores the intronic dinucleotides implied by one genomic-order junction.
fn splice_motif_score(chr_seq: &str, end: i32, start: i32, token_strand: char) -> Option<i32> {
    if start - end < 3 {
        return None;
    }
    let low = reference_dinucleotide(chr_seq, end + 1)?;
    let high = reference_dinucleotide(chr_seq, start - 2)?;
    match token_strand {
        '+' => match (low.as_str(), high.as_str()) {
            ("GT", "AG") => Some(100),
            ("GC", "AG") | ("AT", "AC") => Some(60),
            _ => None,
        },
        '-' => match (low.as_str(), high.as_str()) {
            ("CT", "AC") => Some(100),
            ("CT", "GC") | ("GT", "AT") => Some(60),
            _ => None,
        },
        _ => match (low.as_str(), high.as_str()) {
            ("GT", "AG") | ("CT", "AC") => Some(80),
            ("GC", "AG") | ("AT", "AC") | ("CT", "GC") | ("GT", "AT") => Some(50),
            _ => None,
        },
    }
}

/// Returns a 1-based two-base genomic sequence.
fn reference_dinucleotide(chr_seq: &str, start: i32) -> Option<String> {
    if start <= 0 {
        return None;
    }
    let start_idx = (start - 1) as usize;
    let end_idx = start_idx + 2;
    if end_idx > chr_seq.len() {
        return None;
    }
    Some(chr_seq[start_idx..end_idx].to_ascii_uppercase())
}

/// Returns the query-interval overlap between two parsed alignments.
fn query_overlap(left: &ParsedAlignment, right: &ParsedAlignment) -> i32 {
    let left_start = left
        .blocks
        .first()
        .map(|block| block.read_start)
        .unwrap_or(0);
    let left_end = left.blocks.last().map(|block| block.read_end).unwrap_or(0);
    let right_start = right
        .blocks
        .first()
        .map(|block| block.read_start)
        .unwrap_or(0);
    let right_end = right.blocks.last().map(|block| block.read_end).unwrap_or(0);
    (left_end.min(right_end) - left_start.max(right_start) + 1).max(0)
}

/// Ranking key for one parsed alignment when two candidates overlap the same read slice.
fn parsed_alignment_rank(record: &ParsedAlignment) -> (i32, i32, i32, i32) {
    (
        i32::from(!is_secondary(record.flag)),
        i32::from(!is_supplementary(record.flag)),
        record.mapq,
        record
            .blocks
            .iter()
            .map(|block| block.read_end - block.read_start + 1)
            .sum(),
    )
}

/// Generic chain ranking used before type-specific filtering.
fn chain_generic_rank(chain: &MateChain, circ: Option<&CircRecord>) -> (i32, i32, i32, i32, i32) {
    (
        i32::from(circ.is_none_or(|circ| chain.chrom == circ.chr)),
        i32::from(chain.is_circular),
        i32::from(chain.is_bsj),
        -(chain.used_secondary as i32),
        chain.query_coverage,
    )
}

/// Final chain ranking for one requested output type.
fn chain_rank(
    chain: &MateChain,
    type_name: &str,
    circ: Option<&CircRecord>,
) -> (i32, i32, i32, i32, i32, i32) {
    let valid_type = match type_name {
        "bsj" => i32::from(chain.is_bsj),
        "backward" => i32::from(chain.is_circular && !chain.is_bsj),
        _ => 0,
    };
    (
        valid_type,
        i32::from(circ.is_none_or(|circ| chain.chrom == circ.chr)),
        i32::from(chain.is_circular),
        -(chain.used_secondary as i32),
        -(chain.used_supplementary as i32),
        chain.query_coverage,
    )
}

/// Returns whether two consecutive blocks wrap against the strand-specific linear order.
fn wraps_in_read_order(prev: &SegmentBlock, next: &SegmentBlock, strand: char) -> bool {
    match strand {
        '+' => next.ref_start < prev.ref_start,
        '-' => next.ref_start > prev.ref_start,
        _ => false,
    }
}

/// Returns the opposite genomic traversal direction for mate2 chain reversal.
fn opposite_strand(strand: char) -> char {
    match strand {
        '+' => '-',
        '-' => '+',
        _ => strand,
    }
}

/// Detects the output-chain BSJ transition inside a known circRNA span.
///
/// The circRNA strand determines token strand, but the adjacent block order is
/// determined by the normalized alignment/read chain. Reverse-strand alignments
/// therefore use the opposite genomic wrap direction after query-coordinate
/// normalization.
fn is_bsj_transition(
    prev: &SegmentBlock,
    next: &SegmentBlock,
    circ: &CircRecord,
    order_strand: char,
) -> bool {
    let prev_near_start = prev.ref_start <= circ.start + 6;
    let prev_near_end = prev.ref_end >= circ.end - 6;
    let next_near_start = next.ref_start <= circ.start + 6;
    let next_near_end = next.ref_end >= circ.end - 6;
    ((prev_near_start && next_near_end) || (prev_near_end && next_near_start))
        && wraps_in_read_order(prev, next, order_strand)
}

/// Applies final `.bsj` mate-level evidence to selected mate chains.
///
/// The chain builder can miss a BSJ flag when the mapper reports a single
/// soft-clipped alignment and the CIRI validator recognizes the clipped mate as
/// BSJ evidence. The final `.bsj` row is a post-Summary annotation, so using it
/// here improves mate-level segment labels without changing which reads Summary
/// counted as CIRI3-compatible junction reads.
fn apply_mate_bsj_evidence(
    chains: &mut [Option<MateChain>; 2],
    mate_bsj_evidence: &[MateBsjEvidence],
    circ: &CircRecord,
) {
    for evidence in mate_bsj_evidence
        .iter()
        .filter(|evidence| evidence_matches_circ(evidence, circ))
    {
        if let Some(chain) = chains
            .get_mut(evidence.mate_bucket)
            .and_then(Option::as_mut)
        {
            chain.is_bsj = true;
            chain.is_circular = true;
            if chain.token_strand == '?' {
                chain.token_strand = evidence
                    .strand
                    .chars()
                    .next()
                    .filter(|strand| matches!(strand, '+' | '-'))
                    .unwrap_or('?');
            }
        }
    }
}

/// Returns whether a mate-level BSJ row belongs to the confirmed circRNA row.
fn evidence_matches_circ(evidence: &MateBsjEvidence, circ: &CircRecord) -> bool {
    evidence.chr == circ.chr
        && evidence.start == circ.start
        && evidence.end == circ.end
        && (evidence.strand == circ.strand || evidence.strand == "NA" || circ.strand == "NA")
}

/// Converts one optional chain to output segment text, CIGAR, and BSJ flag.
fn chain_text(chain: Option<&MateChain>) -> (String, String, usize) {
    chain.map_or_else(
        || ("NA".to_string(), "NA".to_string(), 0),
        |chain| {
            (
                chain.tokens.join("|"),
                chain.cigar.clone(),
                usize::from(chain.is_bsj),
            )
        },
    )
}

/// Infers the RNA strand to write into `<prefix>.segments` for one read.
///
/// Validated splice candidates already carry the motif/annotation strand decision
/// from `validate_splice_motifs`. If all candidates for the read agree, we emit
/// that RNA strand; otherwise the read stays `unknown` so alignment strand does
/// not get misreported as biological strand.
fn infer_read_token_strand(read_id: &str, state: &ScanState) -> Option<char> {
    build_read_strand_hints(state)
        .get(read_id)
        .copied()
        .flatten()
}

/// Precomputes one RNA-strand hint per read ID from validated splice candidates.
///
/// Re-scanning `state.candidates` for every backward read turns segments writing
/// into an accidental O(reads * candidates) pass on large fixtures. The strand
/// hints are cheap to summarize once after motif validation, so we cache them
/// into a read-keyed map before building the output rows.
fn build_read_strand_hints(state: &ScanState) -> HashMap<String, Option<char>> {
    let mut raw_hints: HashMap<&str, i32> = HashMap::new();
    let mut conflicted: HashSet<&str> = HashSet::new();
    for candidate in &state.candidates {
        if candidate.strand_hint == 0 {
            continue;
        }
        let read_id = candidate.read_id.as_str();
        match raw_hints.get(read_id).copied() {
            Some(existing) if existing != candidate.strand_hint => {
                conflicted.insert(read_id);
            }
            Some(_) => {}
            None => {
                raw_hints.insert(read_id, candidate.strand_hint);
            }
        }
    }

    let mut out = HashMap::with_capacity(raw_hints.len());
    for (read_id, hint) in raw_hints {
        let strand = if conflicted.contains(read_id) {
            None
        } else {
            Some(match hint {
                -1 => '+',
                1 => '-',
                _ => '?',
            })
        };
        out.insert(read_id.to_string(), strand);
    }
    out
}

/// Precomputes corrected splice-boundary hints for each read.
///
/// Each candidate contributes one `(site2, site1)` pair, matching the downstream
/// convention that segment `end -> site2` and following segment `start -> site1`
/// represent the validated splice boundary after motif/annotation adjustment.
fn build_read_junction_hints(state: &ScanState) -> HashMap<String, Vec<(i32, i32)>> {
    let mut out: HashMap<String, Vec<(i32, i32)>> = HashMap::new();
    for candidate in &state.candidates {
        out.entry(candidate.read_id.clone())
            .or_default()
            .push((candidate.site2, candidate.site1));
    }
    out
}

/// Writes the simplified `<prefix>.segments` protocol.
fn write_segments(path: &str, records: &[SegmentRecord]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "read_id\ttype\tcirc_id\tchrom\tstart\tend\tstrand\tis_circular\tis_r1_bsj\tis_r2_bsj\tr1_cigar\tr1_segments\tr2_cigar\tr2_segments"
    )?;
    for record in records {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.read_id,
            record.type_name,
            record.circ_id,
            record.chrom,
            record.start,
            record.end,
            record.strand,
            record.is_circular,
            record.is_r1_bsj,
            record.is_r2_bsj,
            record.r1_cigar,
            record.r1_segments,
            record.r2_cigar,
            record.r2_segments
        )?;
    }
    Ok(())
}

/// Returns the 0-based mate bucket used by `<prefix>.segments`.
fn mate_bucket(flag: i32) -> usize {
    if flag & 0x40 != 0 {
        0
    } else {
        1
    }
}

/// Stable sidecar-order key for raw mapper alignments.
///
/// The final chain builder still performs biological/topology ranking; this key
/// only makes de-duplicated sidecar input deterministic across Scan1/Scan2 shard
/// merge order without preferring secondary alignments ahead of primary ones.
fn as_alignment_order_key(record: &AsAlignment) -> (usize, usize, i32) {
    (
        usize::from(is_secondary(record.flag)),
        usize::from(is_supplementary(record.flag)),
        -record.mapq,
    )
}

/// Returns whether one SAM flag represents a secondary alignment.
fn is_secondary(flag: i32) -> bool {
    flag & 0x100 != 0
}

/// Returns whether one SAM flag represents a supplementary alignment.
fn is_supplementary(flag: i32) -> bool {
    flag & 0x800 != 0
}

/// Returns the alignment strand as a user-facing `+` / `-`.
fn strand_char(flag: i32) -> char {
    if flag & 0x10 != 0 {
        '-'
    } else {
        '+'
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn annotation_breaks_multi_motif_offset_ties() {
        let mut annotation = Annotation::new();
        annotation
            .chr_exon_start_map
            .insert("chr1\t149094535".to_string(), "GENE1\t-".to_string());
        annotation
            .chr_exon_end_map
            .insert("chr1\t16901187".to_string(), "GENE1\t-".to_string());
        let candidate = PositiveCandidate {
            read_id: "simulate:4875198".to_string(),
            index: 0,
            chr: "chr1".to_string(),
            site1: 149094532,
            site2: 16901184,
            adjust1: -2,
            adjust2: -2,
            strand_hint: 0,
            cigars: [None, None, None],
        };

        let (hint, index) = index_compare("GTACTCAC", "ATCTCCCT", &candidate, Some(&annotation));

        assert_eq!(hint, 1);
        assert_eq!(index, 6);
    }

    #[test]
    fn motif_tie_falls_back_to_lowest_offset_without_annotation() {
        let candidate = PositiveCandidate {
            read_id: "read".to_string(),
            index: 0,
            chr: "chr1".to_string(),
            site1: 149094532,
            site2: 16901184,
            adjust1: -2,
            adjust2: -2,
            strand_hint: 0,
            cigars: [None, None, None],
        };

        let (hint, index) = index_compare("GTACTCAC", "ATCTCCCT", &candidate, None);

        assert_eq!(hint, 1);
        assert_eq!(index, 2);
    }

    #[test]
    fn annotation_can_choose_positive_strand_over_negative_branch() {
        let mut annotation = Annotation::new();
        annotation
            .chr_exon_start_map
            .insert("chr1\t101".to_string(), "GENE1\t+".to_string());
        annotation
            .chr_exon_end_map
            .insert("chr1\t201".to_string(), "GENE1\t+".to_string());
        let candidate = PositiveCandidate {
            read_id: "read".to_string(),
            index: 0,
            chr: "chr1".to_string(),
            site1: 100,
            site2: 200,
            adjust1: 0,
            adjust2: 0,
            strand_hint: 0,
            cigars: [None, None, None],
        };

        let (hint, index) = index_compare("ACAG", "CTGT", &candidate, Some(&annotation));

        assert_eq!(hint, -1);
        assert_eq!(index, 2);
    }

    #[test]
    fn coverage_validation_accepts_supported_exon_body() {
        let mut coverage = HashMap::new();
        let mut chr_cov = HashMap::new();
        for pos in 100..=119 {
            chr_cov.insert(pos, 10);
        }
        coverage.insert("chr1".to_string(), chr_cov);
        let reference = HashMap::new();
        let state = ScanState {
            circ_by_id: HashMap::new(),
            junction_read_to_circ: HashMap::new(),
            mate_bsj_evidence: HashMap::new(),
            reference: &reference,
            annotation: None,
            clusters_by_chr: HashMap::new(),
            read_len: 100,
            candidates: Vec::new(),
            coverage,
            read_mappings: HashMap::new(),
            seen_junction_reads: HashSet::new(),
            segment_groups: HashMap::new(),
            stats: AsStats::default(),
        };
        let junction_coverage = (100..=119).map(|pos| (pos, 1)).collect();

        let result = exon_coverage_validation_single(
            "chr1",
            100,
            119,
            BoundarySupport::Count(2),
            BoundarySupport::Count(2),
            2,
            2,
            0,
            0,
            &junction_coverage,
            &state,
        );

        assert_eq!(result.code, 1);
        assert_eq!(result.median, 10);
    }

    #[test]
    fn negative_strand_cirexon_ids_are_reverse_ordered() {
        let mut records = vec![
            test_cirexon_record(100, 120),
            test_cirexon_record(200, 220),
            test_cirexon_record(300, 320),
        ];
        let mut out = Vec::new();

        append_ordered_cirexons(&mut out, &mut records, "-");

        assert_eq!(out[0].start, 300);
        assert_eq!(out[0].cirexon_id, "cirexon1");
        assert_eq!(out[2].start, 100);
        assert_eq!(out[2].cirexon_id, "cirexon3");
    }

    #[test]
    fn isoform_graph_reconstructs_anchored_path() {
        let circ = CircRecord {
            id: "chr1:100|320".to_string(),
            chr: "chr1".to_string(),
            start: 100,
            end: 320,
            junction_reads: Vec::new(),
            junction_read_count: "12".to_string(),
            pcc: "1_1_1".to_string(),
            non_junction_reads: "0".to_string(),
            junction_reads_ratio: "1.00".to_string(),
            circ_type: "exon".to_string(),
            gene_id: "GENE1".to_string(),
            strand: "+".to_string(),
        };
        let mut splice_links = HashMap::new();
        splice_links.insert((120, 200), 3);
        splice_links.insert((220, 300), 4);
        let context = CirexonContext {
            supporting_start: BTreeMap::new(),
            supporting_end: BTreeMap::new(),
            splice_read_start: HashMap::new(),
            splice_read_end: HashMap::new(),
            splice_links,
            splice_across_count: HashMap::new(),
        };
        let mut exons = vec![
            test_cirexon_record(100, 120),
            test_cirexon_record(200, 220),
            test_cirexon_record(300, 320),
        ];
        for exon in &mut exons {
            exon.strand = "+".to_string();
        }

        let isoforms = build_full_length_isoforms(&circ, "+", &exons, &context, None);

        assert_eq!(isoforms.len(), 1);
        assert_eq!(isoforms[0].exon_chain(), "100:120!+,200:220!+,300:320!+");
        assert_eq!(isoforms[0].junction_chain(), "120:200,220:300");
        assert_eq!(isoforms[0].internal_split_reads, 7);
    }

    #[test]
    fn negative_isoform_sequence_uses_transcript_order() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), "AACCGG".to_string());
        let record = IsoformRecord {
            circ_id: "chr1:1|6".to_string(),
            isoform_id: "iso1".to_string(),
            chr: "chr1".to_string(),
            start: 1,
            end: 6,
            strand: "-".to_string(),
            path_tier: "anchored".to_string(),
            bsj_supported: "yes".to_string(),
            gene_id: "GENE1".to_string(),
            exons: vec![(5, 6), (1, 2)],
            junctions: vec![(2, 5)],
            isoform_len: 4,
            bsj_seed_reads: 1,
            internal_split_reads: 1,
            boundary_clip_reads: 0,
            anomalous_pair_reads: 0,
            annotation_supported_junctions: 0,
            de_novo_supported_junctions: 1,
            score: 1,
            rank: 1,
        };

        let seq = isoform_sequence(&record, &reference).unwrap();

        assert_eq!(seq, "CCTT");
    }

    #[test]
    fn reverse_alignment_blocks_are_flipped_into_read_order() {
        let record = AsAlignment {
            flag: 0x10,
            chr: "chr1".to_string(),
            pos: 100,
            mapq: 60,
            cigar: "10S90M".to_string(),
            seq: "A".repeat(100),
            from_local_clip: false,
            xa_alternatives: Vec::new(),
        };

        let blocks = parse_alignment_blocks(&record, 100).unwrap();

        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].read_start, 1);
        assert_eq!(blocks[0].read_end, 90);
        assert_eq!(blocks[0].ref_start, 100);
        assert_eq!(blocks[0].ref_end, 189);
    }

    #[test]
    fn backward_segments_keep_non_bsj_circular_order() {
        let records = vec![
            AsAlignment {
                flag: 0x40,
                chr: "chr1".to_string(),
                pos: 200,
                mapq: 60,
                cigar: "50M50S".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x40 | 0x800,
                chr: "chr1".to_string(),
                pos: 100,
                mapq: 60,
                cigar: "50S50M".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
        ];
        let reference = HashMap::new();
        let state = ScanState {
            circ_by_id: HashMap::new(),
            junction_read_to_circ: HashMap::new(),
            mate_bsj_evidence: HashMap::new(),
            reference: &reference,
            annotation: None,
            clusters_by_chr: HashMap::new(),
            read_len: 100,
            candidates: Vec::new(),
            coverage: HashMap::new(),
            read_mappings: HashMap::new(),
            seen_junction_reads: HashSet::new(),
            segment_groups: HashMap::new(),
            stats: AsStats::default(),
        };

        let record = build_backward_segment_record(
            "read1",
            &records,
            infer_read_token_strand("read1", &state),
            &[],
            None,
            100,
        )
        .unwrap();

        assert_eq!(record.type_name, "backward");
        assert_eq!(record.circ_id, "NA");
        assert_eq!(record.chrom, "chr1");
        assert_eq!(record.start, "100");
        assert_eq!(record.end, "249");
        assert_eq!(record.strand, "NA");
        assert_eq!(record.is_circular, 1);
        assert_eq!(record.r1_segments, "200-249:?|100-149:?");
        assert_eq!(record.r1_cigar, "50M50N50M");
        assert_eq!(record.is_r1_bsj, 0);
    }

    #[test]
    fn backward_segments_reject_weak_supplementary_when_xa_makes_linear_chain() {
        let records = vec![
            AsAlignment {
                flag: 0x40 | 0x10,
                chr: "chr1".to_string(),
                pos: 92091,
                mapq: 0,
                cigar: "27S123M".to_string(),
                seq: "A".repeat(150),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x40 | 0x10 | 0x800,
                chr: "chr1".to_string(),
                pos: 521712,
                mapq: 0,
                cigar: "29M121H".to_string(),
                seq: "A".repeat(29),
                from_local_clip: false,
                xa_alternatives: vec![XaAlternative {
                    chr: "chr1".to_string(),
                    strand: '-',
                    pos: 91603,
                    cigar: "29M121S".to_string(),
                    edit_distance: 0,
                }],
            },
        ];

        let record =
            build_backward_segment_record("sim:433108", &records, Some('-'), &[], None, 150);

        assert!(record.is_none());
    }

    #[test]
    fn backward_segments_reject_long_anchor_when_xa_makes_linear_chain() {
        let records = vec![
            AsAlignment {
                flag: 0x40,
                chr: "chr1".to_string(),
                pos: 2000,
                mapq: 0,
                cigar: "50M50S".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x40 | 0x800,
                chr: "chr1".to_string(),
                pos: 100,
                mapq: 0,
                cigar: "50S50M".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: vec![XaAlternative {
                    chr: "chr1".to_string(),
                    strand: '+',
                    pos: 2050,
                    cigar: "50S50M".to_string(),
                    edit_distance: 0,
                }],
            },
        ];

        let record =
            build_backward_segment_record("long_anchor", &records, Some('+'), &[], None, 100);

        assert!(record.is_none());
    }

    #[test]
    fn backward_segments_reject_linear_xa_before_circular_repair() {
        let records = vec![
            AsAlignment {
                flag: 83,
                chr: "chr1".to_string(),
                pos: 92_091,
                mapq: 0,
                cigar: "42S108M".to_string(),
                seq: "A".repeat(150),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 2131,
                chr: "chr1".to_string(),
                pos: 521_697,
                mapq: 0,
                cigar: "44M106H".to_string(),
                seq: "A".repeat(44),
                from_local_clip: false,
                xa_alternatives: vec![
                    XaAlternative {
                        chr: "chr1".to_string(),
                        strand: '-',
                        pos: 237_914,
                        cigar: "44M106S".to_string(),
                        edit_distance: 0,
                    },
                    XaAlternative {
                        chr: "chr1".to_string(),
                        strand: '-',
                        pos: 91_588,
                        cigar: "44M106S".to_string(),
                        edit_distance: 0,
                    },
                ],
            },
        ];

        let record =
            build_backward_segment_record("sim:432652", &records, Some('-'), &[], None, 150);

        assert!(record.is_none());
    }

    #[test]
    fn backward_segments_joint_xa_replacement_can_reject_linear_pair() {
        let records = vec![
            AsAlignment {
                flag: 81,
                chr: "chr1".to_string(),
                pos: 120_775,
                mapq: 0,
                cigar: "60S90M".to_string(),
                seq: "A".repeat(150),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 2129,
                chr: "chr1".to_string(),
                pos: 238_508,
                mapq: 0,
                cigar: "62M88H".to_string(),
                seq: "A".repeat(62),
                from_local_clip: false,
                xa_alternatives: vec![XaAlternative {
                    chr: "chr1".to_string(),
                    strand: '-',
                    pos: 92_181,
                    cigar: "62M88S".to_string(),
                    edit_distance: 0,
                }],
            },
            AsAlignment {
                flag: 161,
                chr: "chr1".to_string(),
                pos: 91_547,
                mapq: 55,
                cigar: "85M65S".to_string(),
                seq: "A".repeat(150),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 2209,
                chr: "chr1".to_string(),
                pos: 238_418,
                mapq: 0,
                cigar: "83H67M".to_string(),
                seq: "A".repeat(67),
                from_local_clip: false,
                xa_alternatives: vec![XaAlternative {
                    chr: "chr1".to_string(),
                    strand: '+',
                    pos: 92_091,
                    cigar: "83S67M".to_string(),
                    edit_distance: 0,
                }],
            },
        ];

        let record =
            build_backward_segment_record("sim:432594", &records, Some('-'), &[], None, 150);

        assert!(record.is_none());
    }

    #[test]
    fn backward_segments_uses_primary_mapq0_xa_in_joint_ranking() {
        let records = vec![
            AsAlignment {
                flag: 81,
                chr: "chr1".to_string(),
                pos: 92_091,
                mapq: 0,
                cigar: "44S106M".to_string(),
                seq: "A".repeat(150),
                from_local_clip: false,
                xa_alternatives: vec![XaAlternative {
                    chr: "chr1".to_string(),
                    strand: '-',
                    pos: 238_418,
                    cigar: "44S106M".to_string(),
                    edit_distance: 0,
                }],
            },
            AsAlignment {
                flag: 2129,
                chr: "chr1".to_string(),
                pos: 521_695,
                mapq: 0,
                cigar: "46M104H".to_string(),
                seq: "A".repeat(46),
                from_local_clip: false,
                xa_alternatives: vec![XaAlternative {
                    chr: "chr1".to_string(),
                    strand: '-',
                    pos: 237_912,
                    cigar: "46M104S".to_string(),
                    edit_distance: 0,
                }],
            },
            AsAlignment {
                flag: 161,
                chr: "chr1".to_string(),
                pos: 237_835,
                mapq: 23,
                cigar: "123M27S".to_string(),
                seq: "A".repeat(150),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 2209,
                chr: "chr1".to_string(),
                pos: 238_418,
                mapq: 0,
                cigar: "121H29M".to_string(),
                seq: "A".repeat(29),
                from_local_clip: false,
                xa_alternatives: vec![XaAlternative {
                    chr: "chr1".to_string(),
                    strand: '+',
                    pos: 92_091,
                    cigar: "121S29M".to_string(),
                    edit_distance: 0,
                }],
            },
        ];

        let record =
            build_backward_segment_record("sim:4599148", &records, Some('-'), &[], None, 150);

        assert!(record.is_none());
    }

    #[test]
    fn backward_segments_accepts_nm1_xa_when_span_is_better() {
        let records = vec![
            AsAlignment {
                flag: 0x40,
                chr: "chr1".to_string(),
                pos: 2_000,
                mapq: 0,
                cigar: "50M50S".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x40 | 0x800,
                chr: "chr1".to_string(),
                pos: 100,
                mapq: 0,
                cigar: "50S50M".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: vec![XaAlternative {
                    chr: "chr1".to_string(),
                    strand: '+',
                    pos: 2_050,
                    cigar: "50S50M".to_string(),
                    edit_distance: 1,
                }],
            },
        ];

        let record =
            build_backward_segment_record("nm1_linear_xa", &records, Some('+'), &[], None, 100);

        assert!(record.is_none());
    }

    #[test]
    fn backward_segments_choose_best_spanning_xa_regardless_of_topology() {
        let records = vec![
            AsAlignment {
                flag: 0x40,
                chr: "chr1".to_string(),
                pos: 5_000,
                mapq: 60,
                cigar: "50M50S".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x40 | 0x800,
                chr: "chr1".to_string(),
                pos: 100,
                mapq: 0,
                cigar: "50S50M".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: vec![
                    XaAlternative {
                        chr: "chr1".to_string(),
                        strand: '+',
                        pos: 4_990,
                        cigar: "50S50M".to_string(),
                        edit_distance: 0,
                    },
                    XaAlternative {
                        chr: "chr1".to_string(),
                        strand: '+',
                        pos: 9_000,
                        cigar: "50S50M".to_string(),
                        edit_distance: 0,
                    },
                ],
            },
        ];

        let record =
            build_backward_segment_record("best_circular_xa", &records, Some('+'), &[], None, 100)
                .unwrap();

        assert_eq!(record.start, "4990");
        assert_eq!(record.end, "5049");
        assert_eq!(record.r1_segments, "5000-5049:+|4990-5039:+");
    }

    #[test]
    fn backward_segments_reject_span_above_ciri3_max_span() {
        let records = vec![
            AsAlignment {
                flag: 0x40,
                chr: "chr1".to_string(),
                pos: 300_000,
                mapq: 60,
                cigar: "50M50S".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x40 | 0x800,
                chr: "chr1".to_string(),
                pos: 100,
                mapq: 60,
                cigar: "50S50M".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
        ];

        let record =
            build_backward_segment_record("wide_span", &records, Some('+'), &[], None, 100);

        assert!(record.is_none());
    }

    #[test]
    fn backward_segments_repair_remote_supplementary_with_local_xa() {
        let records = vec![
            AsAlignment {
                flag: 99,
                chr: "chr1".to_string(),
                pos: 203_745_131,
                mapq: 60,
                cigar: "115M35S".to_string(),
                seq: "A".repeat(150),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 2147,
                chr: "chr1".to_string(),
                pos: 203_741_196,
                mapq: 60,
                cigar: "115H35M".to_string(),
                seq: "A".repeat(35),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 147,
                chr: "chr1".to_string(),
                pos: 203_741_196,
                mapq: 60,
                cigar: "22S80M48S".to_string(),
                seq: "A".repeat(150),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 2195,
                chr: "chr1".to_string(),
                pos: 203_743_002,
                mapq: 60,
                cigar: "101H49M".to_string(),
                seq: "A".repeat(49),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 2195,
                chr: "chr1".to_string(),
                pos: 90_982,
                mapq: 0,
                cigar: "23M127H".to_string(),
                seq: "A".repeat(23),
                from_local_clip: false,
                xa_alternatives: vec![XaAlternative {
                    chr: "chr1".to_string(),
                    strand: '-',
                    pos: 203_745_224,
                    cigar: "22M128S".to_string(),
                    edit_distance: 0,
                }],
            },
        ];

        let record =
            build_backward_segment_record("sim:3966179", &records, Some('+'), &[], None, 150)
                .unwrap();

        assert_eq!(record.start, "203741196");
        assert_eq!(record.end, "203745245");
        assert_eq!(
            record.r2_segments,
            "203745224-203745245:+|203741196-203741275:+|203743002-203743050:+"
        );
        assert_eq!(record.r2_cigar, "22M3948N80M1726N49M");
    }

    #[test]
    fn backward_segments_filter_mixed_chromosome_chains() {
        let records = vec![
            AsAlignment {
                flag: 0x40,
                chr: "chr1".to_string(),
                pos: 200,
                mapq: 60,
                cigar: "50M50S".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x40 | 0x800,
                chr: "chr1".to_string(),
                pos: 100,
                mapq: 60,
                cigar: "50S50M".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x80,
                chr: "chr2".to_string(),
                pos: 500,
                mapq: 60,
                cigar: "100M".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
        ];

        let record = build_backward_segment_record("read1", &records, None, &[], None, 100);

        assert!(record.is_none());
    }

    #[test]
    fn segment_records_sort_by_genomic_position() {
        let mut records = vec![
            SegmentRecord {
                read_id: "read3".to_string(),
                type_name: "backward",
                circ_id: "NA".to_string(),
                chrom: "chr2".to_string(),
                start: "20".to_string(),
                end: "40".to_string(),
                strand: "NA".to_string(),
                is_circular: 1,
                is_r1_bsj: 0,
                is_r2_bsj: 0,
                r1_cigar: "NA".to_string(),
                r1_segments: "NA".to_string(),
                r2_cigar: "NA".to_string(),
                r2_segments: "NA".to_string(),
            },
            SegmentRecord {
                read_id: "read2".to_string(),
                type_name: "backward",
                circ_id: "NA".to_string(),
                chrom: "chr1".to_string(),
                start: "100".to_string(),
                end: "250".to_string(),
                strand: "NA".to_string(),
                is_circular: 1,
                is_r1_bsj: 0,
                is_r2_bsj: 0,
                r1_cigar: "NA".to_string(),
                r1_segments: "NA".to_string(),
                r2_cigar: "NA".to_string(),
                r2_segments: "NA".to_string(),
            },
            SegmentRecord {
                read_id: "read1".to_string(),
                type_name: "bsj",
                circ_id: "chr1:100|200".to_string(),
                chrom: "chr1".to_string(),
                start: "100".to_string(),
                end: "200".to_string(),
                strand: "+".to_string(),
                is_circular: 1,
                is_r1_bsj: 1,
                is_r2_bsj: 0,
                r1_cigar: "NA".to_string(),
                r1_segments: "NA".to_string(),
                r2_cigar: "NA".to_string(),
                r2_segments: "NA".to_string(),
            },
        ];

        sort_segment_records(&mut records);

        assert_eq!(
            records
                .iter()
                .map(|record| record.read_id.as_str())
                .collect::<Vec<_>>(),
            vec!["read1", "read2", "read3"]
        );
    }

    #[test]
    fn bsj_segments_insert_boundary_marker() {
        let circ = CircRecord {
            id: "chr1:100|199".to_string(),
            chr: "chr1".to_string(),
            start: 100,
            end: 199,
            junction_reads: vec!["read1".to_string()],
            junction_read_count: "1".to_string(),
            pcc: "1_0_0".to_string(),
            non_junction_reads: "0".to_string(),
            junction_reads_ratio: "1.00".to_string(),
            circ_type: "exon".to_string(),
            gene_id: "GENE1".to_string(),
            strand: "+".to_string(),
        };
        let records = vec![
            AsAlignment {
                flag: 0x40,
                chr: "chr1".to_string(),
                pos: 150,
                mapq: 60,
                cigar: "50M50S".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x40 | 0x800,
                chr: "chr1".to_string(),
                pos: 100,
                mapq: 60,
                cigar: "50S50M".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
        ];

        let record =
            build_bsj_segment_record("read1", &records, &circ, &[], &[], None, 100).unwrap();

        assert_eq!(record.type_name, "bsj");
        assert_eq!(record.circ_id, circ.id);
        assert_eq!(record.r1_segments, "150-199:+|<bsj>|100-149:+");
        assert_eq!(record.r1_cigar, "50M0B50M");
        assert_eq!(record.is_r1_bsj, 1);
    }

    #[test]
    fn bsj_marker_follows_read_order_wrap_in_three_block_chain() {
        let circ = CircRecord {
            id: "chr1:100|349".to_string(),
            chr: "chr1".to_string(),
            start: 100,
            end: 349,
            junction_reads: vec!["read1".to_string()],
            junction_read_count: "1".to_string(),
            pcc: "1_0_0".to_string(),
            non_junction_reads: "0".to_string(),
            junction_reads_ratio: "1.00".to_string(),
            circ_type: "exon".to_string(),
            gene_id: "GENE1".to_string(),
            strand: "+".to_string(),
        };
        let parsed = vec![ParsedAlignment {
            flag: 0,
            chrom: "chr1".to_string(),
            strand: '+',
            mapq: 60,
            blocks: vec![
                SegmentBlock {
                    read_start: 1,
                    read_end: 10,
                    ref_start: 200,
                    ref_end: 209,
                    from_local_clip: false,
                },
                SegmentBlock {
                    read_start: 11,
                    read_end: 60,
                    ref_start: 300,
                    ref_end: 349,
                    from_local_clip: false,
                },
                SegmentBlock {
                    read_start: 61,
                    read_end: 110,
                    ref_start: 100,
                    ref_end: 149,
                    from_local_clip: false,
                },
            ],
        }];

        let chain = materialize_chain(&parsed, 110, Some(&circ), '+', &[], None, false).unwrap();

        assert_eq!(
            chain.tokens,
            vec![
                "200-209:+".to_string(),
                "300-349:+".to_string(),
                "<bsj>".to_string(),
                "100-149:+".to_string(),
            ]
        );
        assert_eq!(chain.cigar, "10M90N50M150B50M");
        assert!(chain.is_bsj);
    }

    #[test]
    fn bsj_marker_uses_start_side_block_for_reverse_wrap_order() {
        let circ = CircRecord {
            id: "chr1:100|349".to_string(),
            chr: "chr1".to_string(),
            start: 100,
            end: 349,
            junction_reads: vec!["read1".to_string()],
            junction_read_count: "1".to_string(),
            pcc: "1_0_0".to_string(),
            non_junction_reads: "0".to_string(),
            junction_reads_ratio: "1.00".to_string(),
            circ_type: "exon".to_string(),
            gene_id: "GENE1".to_string(),
            strand: "-".to_string(),
        };
        let parsed = vec![ParsedAlignment {
            flag: 0x10,
            chrom: "chr1".to_string(),
            strand: '-',
            mapq: 60,
            blocks: vec![
                SegmentBlock {
                    read_start: 1,
                    read_end: 50,
                    ref_start: 100,
                    ref_end: 149,
                    from_local_clip: false,
                },
                SegmentBlock {
                    read_start: 51,
                    read_end: 110,
                    ref_start: 300,
                    ref_end: 349,
                    from_local_clip: false,
                },
                SegmentBlock {
                    read_start: 111,
                    read_end: 120,
                    ref_start: 200,
                    ref_end: 209,
                    from_local_clip: false,
                },
            ],
        }];

        let chain = materialize_chain(&parsed, 120, Some(&circ), '-', &[], None, false).unwrap();

        assert_eq!(
            chain.tokens,
            vec![
                "100-149:-".to_string(),
                "<bsj>".to_string(),
                "300-349:-".to_string(),
                "200-209:-".to_string(),
            ]
        );
        assert_eq!(chain.cigar, "50M150B50M90N10M");
        assert!(chain.is_bsj);
    }

    #[test]
    fn mate_bsj_evidence_marks_soft_clipped_mate() {
        let circ = CircRecord {
            id: "chr1:100|199".to_string(),
            chr: "chr1".to_string(),
            start: 100,
            end: 199,
            junction_reads: vec!["read1".to_string()],
            junction_read_count: "1".to_string(),
            pcc: "1_0_0".to_string(),
            non_junction_reads: "0".to_string(),
            junction_reads_ratio: "1.00".to_string(),
            circ_type: "exon".to_string(),
            gene_id: "GENE1".to_string(),
            strand: "+".to_string(),
        };
        let records = vec![AsAlignment {
            flag: 0x80,
            chr: "chr1".to_string(),
            pos: 100,
            mapq: 60,
            cigar: "10S90M".to_string(),
            seq: "A".repeat(100),
            from_local_clip: false,
            xa_alternatives: Vec::new(),
        }];
        let evidence = vec![MateBsjEvidence {
            mate_bucket: 1,
            chr: "chr1".to_string(),
            start: 100,
            end: 199,
            strand: "+".to_string(),
            priority: 0,
            source_stage: "scan2".to_string(),
        }];

        let record =
            build_bsj_segment_record("read1", &records, &circ, &evidence, &[], None, 100).unwrap();

        assert_eq!(record.is_r1_bsj, 0);
        assert_eq!(record.r2_segments, "100-189:+");
        assert_eq!(record.r2_cigar, "90M10S");
        assert_eq!(record.is_r2_bsj, 1);
    }

    #[test]
    fn materialize_chain_preserves_terminal_soft_clips_in_output_cigar() {
        let parsed = vec![ParsedAlignment {
            flag: 0,
            chrom: "chr1".to_string(),
            strand: '+',
            mapq: 60,
            blocks: vec![SegmentBlock {
                read_start: 11,
                read_end: 90,
                ref_start: 100,
                ref_end: 179,
                from_local_clip: false,
            }],
        }];

        let chain = materialize_chain(&parsed, 100, None, '+', &[], None, false).unwrap();

        assert_eq!(chain.tokens, vec!["100-179:+".to_string()]);
        assert_eq!(chain.cigar, "10S80M10S");
    }

    #[test]
    fn bsj_transition_follows_output_chain_wrap() {
        let circ = CircRecord {
            id: "chr1:100|199".to_string(),
            chr: "chr1".to_string(),
            start: 100,
            end: 199,
            junction_reads: vec![],
            junction_read_count: "1".to_string(),
            pcc: "1_0_0".to_string(),
            non_junction_reads: "0".to_string(),
            junction_reads_ratio: "1.00".to_string(),
            circ_type: "exon".to_string(),
            gene_id: "GENE1".to_string(),
            strand: "+".to_string(),
        };

        let low_then_high = (
            SegmentBlock {
                read_start: 1,
                read_end: 90,
                ref_start: 100,
                ref_end: 149,
                from_local_clip: false,
            },
            SegmentBlock {
                read_start: 91,
                read_end: 150,
                ref_start: 150,
                ref_end: 199,
                from_local_clip: false,
            },
        );
        let high_then_low = (low_then_high.1.clone(), low_then_high.0.clone());

        assert!(is_bsj_transition(
            &low_then_high.0,
            &low_then_high.1,
            &circ,
            '-'
        ));
        assert!(!is_bsj_transition(
            &high_then_low.0,
            &high_then_low.1,
            &circ,
            '-'
        ));
        assert!(!is_bsj_transition(
            &low_then_high.0,
            &low_then_high.1,
            &circ,
            '+'
        ));
        assert!(is_bsj_transition(
            &high_then_low.0,
            &high_then_low.1,
            &circ,
            '+'
        ));
    }

    #[test]
    fn materialize_chain_keeps_distant_blocks_separate() {
        let parsed = vec![
            ParsedAlignment {
                flag: 0x10,
                chrom: "chr1".to_string(),
                strand: '-',
                mapq: 60,
                blocks: vec![SegmentBlock {
                    read_start: 1,
                    read_end: 29,
                    ref_start: 93811203,
                    ref_end: 93811231,
                    from_local_clip: false,
                }],
            },
            ParsedAlignment {
                flag: 0x10,
                chrom: "chr1".to_string(),
                strand: '-',
                mapq: 60,
                blocks: vec![SegmentBlock {
                    read_start: 29,
                    read_end: 83,
                    ref_start: 93806014,
                    ref_end: 93806068,
                    from_local_clip: false,
                }],
            },
            ParsedAlignment {
                flag: 0x10,
                chrom: "chr1".to_string(),
                strand: '-',
                mapq: 60,
                blocks: vec![SegmentBlock {
                    read_start: 83,
                    read_end: 150,
                    ref_start: 93811273,
                    ref_end: 93811340,
                    from_local_clip: false,
                }],
            },
        ];

        let chain = materialize_chain(&parsed, 150, None, '-', &[], None, false).unwrap();

        assert_eq!(
            chain.tokens,
            vec![
                "93811203-93811231:-".to_string(),
                "93806014-93806068:-".to_string(),
                "93811273-93811340:-".to_string(),
            ]
        );
        assert_eq!(chain.cigar, "29M5134N55M5204N68M");
    }

    #[test]
    fn materialize_chain_corrects_internal_splice_boundaries_with_annotation() {
        let mut annotation = Annotation::new();
        annotation
            .chr_exon_end_map
            .insert("chr1\t171292331".to_string(), "GENE1\t+".to_string());
        annotation
            .chr_exon_start_map
            .insert("chr1\t171300779".to_string(), "GENE1\t+".to_string());
        let reference = HashMap::new();
        let correction = SegmentCorrectionContext {
            reference: &reference,
            annotation: Some(&annotation),
            junction_support: None,
        };
        let parsed = vec![ParsedAlignment {
            flag: 0,
            chrom: "chr1".to_string(),
            strand: '+',
            mapq: 60,
            blocks: vec![
                SegmentBlock {
                    read_start: 1,
                    read_end: 88,
                    ref_start: 171292245,
                    ref_end: 171292332,
                    from_local_clip: false,
                },
                SegmentBlock {
                    read_start: 89,
                    read_end: 146,
                    ref_start: 171300777,
                    ref_end: 171300834,
                    from_local_clip: false,
                },
            ],
        }];

        let chain =
            materialize_chain(&parsed, 146, None, '+', &[], Some(&correction), false).unwrap();

        assert_eq!(
            chain.tokens,
            vec![
                "171292245-171292331:+".to_string(),
                "171300779-171300834:+".to_string(),
            ]
        );
        assert_eq!(chain.cigar, "87M8447N56M");
    }

    #[test]
    fn materialize_chain_prefers_transcript_consistent_boundary() {
        let mut annotation = Annotation::new();
        annotation
            .chr_exon_end_map
            .insert("chr1\t196".to_string(), "GENE1\t+".to_string());
        annotation
            .chr_exon_start_map
            .insert("chr1\t300".to_string(), "GENE1\t+".to_string());
        annotation
            .chr_exon_end_map
            .insert("chr1\t199".to_string(), "GENE1\t+".to_string());
        annotation
            .chr_exon_start_map
            .insert("chr1\t304".to_string(), "GENE1\t+".to_string());
        annotation
            .transcript_splice_map
            .insert("chr1\t199\t304\t+".to_string());
        let reference = HashMap::new();
        let correction = SegmentCorrectionContext {
            reference: &reference,
            annotation: Some(&annotation),
            junction_support: None,
        };
        let parsed = vec![ParsedAlignment {
            flag: 0,
            chrom: "chr1".to_string(),
            strand: '+',
            mapq: 60,
            blocks: vec![
                SegmentBlock {
                    read_start: 1,
                    read_end: 97,
                    ref_start: 100,
                    ref_end: 196,
                    from_local_clip: false,
                },
                SegmentBlock {
                    read_start: 98,
                    read_end: 148,
                    ref_start: 300,
                    ref_end: 350,
                    from_local_clip: false,
                },
            ],
        }];

        let chain =
            materialize_chain(&parsed, 148, None, '+', &[], Some(&correction), false).unwrap();

        assert_eq!(
            chain.tokens,
            vec!["100-199:+".to_string(), "304-350:+".to_string()]
        );
        assert_eq!(chain.cigar, "100M104N47M");
    }

    #[test]
    fn materialize_chain_prefers_supported_boundary_for_local_clip() {
        let reference = HashMap::new();
        let mut support: JunctionSupportMap = HashMap::new();
        support
            .entry("chr1".to_string())
            .or_default()
            .insert((199, 328, '+'), 12);
        let correction = SegmentCorrectionContext {
            reference: &reference,
            annotation: None,
            junction_support: Some(&support),
        };
        let parsed = vec![ParsedAlignment {
            flag: 0,
            chrom: "chr1".to_string(),
            strand: '+',
            mapq: 60,
            blocks: vec![
                SegmentBlock {
                    read_start: 1,
                    read_end: 100,
                    ref_start: 100,
                    ref_end: 199,
                    from_local_clip: false,
                },
                SegmentBlock {
                    read_start: 131,
                    read_end: 150,
                    ref_start: 340,
                    ref_end: 359,
                    from_local_clip: true,
                },
            ],
        }];

        let chain =
            materialize_chain(&parsed, 150, None, '+', &[], Some(&correction), false).unwrap();

        assert_eq!(
            chain.tokens,
            vec!["100-199:+".to_string(), "328-359:+".to_string()]
        );
        assert_eq!(chain.cigar, "100M128N32M");
    }

    #[test]
    fn supported_alternative_filter_marks_only_nearby_different_junctions() {
        let mut support_map: JunctionSupportMap = HashMap::new();
        support_map
            .entry("chr1".to_string())
            .or_default()
            .insert((199, 328, '+'), 12);
        let support = build_junction_support_index(&support_map);
        assert!(segments_have_supported_alternative(
            "100-199:+|340-359:+",
            support.get("chr1").unwrap()
        ));

        let mut same_support_map: JunctionSupportMap = HashMap::new();
        same_support_map
            .entry("chr1".to_string())
            .or_default()
            .insert((199, 340, '+'), 12);
        let same_support = build_junction_support_index(&same_support_map);
        assert!(!segments_have_supported_alternative(
            "100-199:+|340-359:+",
            same_support.get("chr1").unwrap()
        ));
    }

    #[test]
    fn materialize_chain_widens_annotation_window_for_local_clip_block() {
        let mut annotation = Annotation::new();
        annotation
            .chr_exon_end_map
            .insert("chr1\t199".to_string(), "GENE1\t+".to_string());
        annotation
            .chr_exon_start_map
            .insert("chr1\t328".to_string(), "GENE1\t+".to_string());
        let reference = HashMap::new();
        let correction = SegmentCorrectionContext {
            reference: &reference,
            annotation: Some(&annotation),
            junction_support: None,
        };
        let parsed = vec![ParsedAlignment {
            flag: 0,
            chrom: "chr1".to_string(),
            strand: '+',
            mapq: 60,
            blocks: vec![
                SegmentBlock {
                    read_start: 1,
                    read_end: 100,
                    ref_start: 100,
                    ref_end: 199,
                    from_local_clip: false,
                },
                SegmentBlock {
                    read_start: 131,
                    read_end: 150,
                    ref_start: 340,
                    ref_end: 359,
                    from_local_clip: true,
                },
            ],
        }];

        let chain =
            materialize_chain(&parsed, 150, None, '+', &[], Some(&correction), false).unwrap();

        assert_eq!(
            chain.tokens,
            vec!["100-199:+".to_string(), "328-359:+".to_string()]
        );
        assert_eq!(chain.cigar, "100M128N32M");
    }

    #[test]
    fn materialize_chain_skips_bsj_gap_but_corrects_adjacent_internal_junction() {
        let circ = CircRecord {
            id: "chr1:100|349".to_string(),
            chr: "chr1".to_string(),
            start: 100,
            end: 349,
            junction_reads: vec!["read1".to_string()],
            junction_read_count: "1".to_string(),
            pcc: "1_0_0".to_string(),
            non_junction_reads: "0".to_string(),
            junction_reads_ratio: "1.00".to_string(),
            circ_type: "exon".to_string(),
            gene_id: "GENE1".to_string(),
            strand: "+".to_string(),
        };
        let mut annotation = Annotation::new();
        annotation
            .chr_exon_end_map
            .insert("chr1\t149".to_string(), "GENE1\t+".to_string());
        annotation
            .chr_exon_start_map
            .insert("chr1\t200".to_string(), "GENE1\t+".to_string());
        annotation
            .chr_exon_end_map
            .insert("chr1\t345".to_string(), "GENE1\t+".to_string());
        annotation
            .chr_exon_start_map
            .insert("chr1\t104".to_string(), "GENE1\t+".to_string());
        let reference = HashMap::new();
        let correction = SegmentCorrectionContext {
            reference: &reference,
            annotation: Some(&annotation),
            junction_support: None,
        };
        let parsed = vec![ParsedAlignment {
            flag: 0,
            chrom: "chr1".to_string(),
            strand: '+',
            mapq: 60,
            blocks: vec![
                SegmentBlock {
                    read_start: 1,
                    read_end: 50,
                    ref_start: 300,
                    ref_end: 349,
                    from_local_clip: false,
                },
                SegmentBlock {
                    read_start: 51,
                    read_end: 101,
                    ref_start: 100,
                    ref_end: 150,
                    from_local_clip: false,
                },
                SegmentBlock {
                    read_start: 102,
                    read_end: 153,
                    ref_start: 198,
                    ref_end: 249,
                    from_local_clip: false,
                },
            ],
        }];

        let chain = materialize_chain(
            &parsed,
            153,
            Some(&circ),
            '+',
            &[],
            Some(&correction),
            false,
        )
        .unwrap();

        assert_eq!(
            chain.tokens,
            vec![
                "300-349:+".to_string(),
                "<bsj>".to_string(),
                "100-149:+".to_string(),
                "200-249:+".to_string(),
            ]
        );
        assert_eq!(chain.cigar, "50M150B50M50N50M");
    }

    #[test]
    fn materialize_chain_uses_splice_signal_without_annotation() {
        let mut chr = vec![b'N'; 300];
        chr[149] = b'G';
        chr[150] = b'T';
        chr[198] = b'A';
        chr[199] = b'G';
        let chr = String::from_utf8(chr).unwrap();
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), chr);
        let correction = SegmentCorrectionContext {
            reference: &reference,
            annotation: None,
            junction_support: None,
        };
        let parsed = vec![ParsedAlignment {
            flag: 0,
            chrom: "chr1".to_string(),
            strand: '+',
            mapq: 60,
            blocks: vec![
                SegmentBlock {
                    read_start: 1,
                    read_end: 52,
                    ref_start: 100,
                    ref_end: 151,
                    from_local_clip: false,
                },
                SegmentBlock {
                    read_start: 53,
                    read_end: 105,
                    ref_start: 198,
                    ref_end: 250,
                    from_local_clip: false,
                },
            ],
        }];

        let chain =
            materialize_chain(&parsed, 105, None, '+', &[], Some(&correction), false).unwrap();

        assert_eq!(
            chain.tokens,
            vec!["100-149:+".to_string(), "201-250:+".to_string()]
        );
        assert_eq!(chain.cigar, "50M51N50M");
    }

    fn test_cirexon_record(start: i32, end: i32) -> CirexonRecord {
        CirexonRecord {
            circ_id: "chr1:100|320".to_string(),
            chr: "chr1".to_string(),
            circ_start: 100,
            circ_end: 320,
            strand: "-".to_string(),
            junction_read_count: "1".to_string(),
            pcc: "1_0_0".to_string(),
            non_junction_reads: "0".to_string(),
            junction_reads_ratio: "1.00".to_string(),
            circ_type: "exon".to_string(),
            gene_id: "GENE1".to_string(),
            cirexon_id: String::new(),
            start,
            end,
            start_support: "1".to_string(),
            end_support: "1".to_string(),
            coverage_median: 10,
            is_icf: false,
        }
    }
}
