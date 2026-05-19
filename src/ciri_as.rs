//! Post-Summary circRNA read reconstruction entry points.
//!
//! This module remains intentionally separate from `Scan1 -> Scan2 -> Summary`.
//! The CIRI3 Rust path is parity-sensitive, so all circRNA internal-structure
//! work stays in a post-Summary phase that consumes the final circRNA table and
//! sidecar alignment evidence without feeding evidence back into BSJ detection.
//! The current deliverable is `<prefix>.segments`: a read-level, splice-aware
//! representation of confirmed BSJ reads that can be compared directly against
//! simulator truth before full-length path reconstruction is re-enabled on top
//! of it. Full-length reconstruction is likewise a sidecar: it consumes the
//! completed segment rows and writes one major isoform per Summary-confirmed
//! circRNA without changing CIRI3 `.out` or `.bsj` decisions.

use anyhow::{anyhow, bail, Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use memmap2::Mmap;
use noodles::bam;
use noodles::sam::{
    self,
    alignment::record::data::field::{Tag, Value},
};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::{hash_map::DefaultHasher, BTreeMap, HashMap, HashSet};
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufRead, BufReader, BufWriter, Seek, Write};
use std::time::{Duration, Instant};

use crate::annotation::Annotation;
use crate::sam_bam::{detect_format, InputFormat};
use crate::utils::{
    bam_shard_count, bsj_payload_start, clip_placement_cigar, clip_sequence_payload,
    exact_clip_match_positions, is_bsj_mate_label, parse_cigar_ops_basic, parse_clip_payload,
    part_path, reverse_complement,
};

const MIN_INTRON: i32 = 70;
const MIN_EXON_LENGTH: i32 = 20;
const MAX_EXON_LENGTH: i32 = 2000;
const MIN_JUNCTION_SUPPORT_SEGMENT_LENGTH: i32 = 10;
const INTERNAL_SPLICE_CORRECTION_WINDOW: i32 = 4;
const PARTIAL_LOCAL_SPLICE_CORRECTION_WINDOW: i32 = 16;
const MAPQ_THRES: i32 = 5;
/// Minimum terminal 3' clip required on both mates for exact-span outward-pair evidence.
///
/// BSJ detection relies on roughly 20 bp split anchors. Exact-span outward pairs
/// have no aligned outward length, so their terminal clip evidence should not be
/// weaker.
const OUTWARD_MIN_TERMINAL_CLIP: i32 = 19;
/// Minimum aligned outward length required on both sides of a non-identical
/// outward-facing pair. Terminal clips are not required when the aligned spans
/// already provide this much outward extension.
const OUTWARD_MIN_PAIR_OFFSET: i32 = 19;
const MAPQ_UNI: i32 = 0;
const MAPQ_BOTH: i32 = 0;
/// Maximum unspliced graph block kept as a single mature exon in isoform output.
///
/// This follows the CIRI-AS exon-length scale instead of a permissive genomic
/// span cutoff. Larger blocks usually mean the segment graph lacks internal
/// splice evidence; when annotation is available, those blocks are projected to
/// known exons so intronic genomic span is not reported as mature RNA sequence.
const MAJOR_MAX_UNSPLICED_EXON_LEN: i32 = MAX_EXON_LENGTH;
/// Minimum BSJ/backward read-chain support required for every adjacent selected
/// junction pair in a mature multi-exon isoform.
///
/// Outward reads are useful weak completion evidence during path selection, but
/// they do not localize a BSJ range precisely enough to certify mature phasing.
const MAJOR_MIN_MATURE_LINK_SUPPORT: f64 = 1.0;
/// Smallest non-BSJ assignment probability retained for isoform graph support.
///
/// Backward/outward rows have no exact Summary BSJ. The isoform stage therefore
/// treats their circRNA membership as probabilistic and leaves weak, internal
/// placements unassigned instead of hard-counting them for an enclosing locus.
const MAJOR_MIN_NON_BSJ_ASSIGNMENT_WEIGHT: f64 = 0.001;
/// Minimum distance scale used when comparing non-BSJ anchors with BSJ sites.
///
/// The dynamic scale starts from the read-level aligned bases. This lower bound
/// keeps normal read-pair jitter from erasing real boundary-adjacent support.
const MAJOR_NON_BSJ_MIN_DISTANCE_SCALE: f64 = 100.0;
/// Maximum distance scale used for non-BSJ circRNA assignment.
///
/// Long genomic spans often include introns or internal circular wraps; letting
/// that span define the scale would incorrectly make distant internal reads look
/// compatible with an outer BSJ.
const MAJOR_NON_BSJ_MAX_DISTANCE_SCALE: f64 = 1000.0;
/// Maximum genomic gap used to merge high-confidence read spans into one
/// estimate anchor.
///
/// Paired BSJ/backward mates can bracket a short unsequenced insert inside the
/// same exon even when the reads do not overlap. Keeping this window small lets
/// those mate-level anchors guide annotation projection without merging across
/// the kilobase-scale introns that should still be represented as exon gaps.
const MAJOR_ESTIMATE_ANCHOR_MERGE_GAP: i32 = 300;
/// Minimum segment-covered fraction required to trust unannotated long exons.
///
/// Annotation can legitimately contain rare mega-exons, but an unannotated
/// long exon with only tiny read anchors is usually an unresolved structure
/// estimate rather than a sequence-ready isoform. This threshold only gates
/// FASTA output for estimates; the GTF audit row is still emitted.
const MAJOR_MIN_FASTA_SEGMENT_COVERAGE_PCT: f64 = 10.0;
/// Minimum segment-covered fraction required for estimate FASTA output.
///
/// FASTA is intended to be the high-confidence sequence set. Estimates can
/// remain valuable audit records in GTF, but sequence output should require
/// broad direct segment support unless the isoform is classified as mature.
const MAJOR_MIN_TRUSTED_ESTIMATE_SEGMENT_COVERAGE_PCT: f64 = 50.0;
/// Minimum segment-covered fraction for high-confidence unphased candidates.
///
/// These records are not mature because adjacent junction-pair phasing is
/// incomplete, but a broadly segment-covered chain is still useful enough to
/// include in the sequence FASTA with its estimate reason preserved.
const MAJOR_MIN_CANDIDATE_SEGMENT_COVERAGE_PCT: f64 = 90.0;
/// Half-window around local anchors searched for non-BSJ clip placement.
const NON_BSJ_LOCAL_CLIP_ANCHOR_FLANK_MULTIPLIER: i32 = 2;
/// Maximum local clip pseudo-alignments retained per read group.
const NON_BSJ_LOCAL_CLIP_MAX_ROWS_PER_READ: usize = 16;
/// Maximum local clip placements retained for one clipped side before chain ranking.
const NON_BSJ_LOCAL_CLIP_MAX_PLACEMENTS_PER_SIDE: usize = 4;
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
    /// Scan2 non-BSJ topology sidecars produced before Summary.
    ///
    /// These rows are read-level candidates, not circRNA-level assignments. They
    /// can represent backward/outward evidence with no known BSJ, so the final
    /// segments phase skips only reads already assigned as BSJ and then rebuilds
    /// the supplemental rows from the read chain itself.
    pub non_bsj_segment_evidence_paths: Vec<&'a str>,
    /// Output prefix; the current phase writes `<prefix>.segments`.
    pub out_prefix: &'a str,
    /// Keep shard-local intermediate files after successful segments finalize.
    ///
    /// The CLI uses this for `--debug`. Most retained non-BSJ alignments live in
    /// temporary shard streams, so deleting them inside this module would make
    /// the user-facing debug flag incomplete.
    pub keep_temp_files: bool,
    /// In-memory reference sequence map loaded by the main CLI.
    pub reference: &'a HashMap<String, String>,
    /// Optional exon-boundary annotation used to resolve ambiguous motif offsets.
    ///
    /// CIRI-AS v1.2 does not consult annotation during `index_compare`, so Perl
    /// falls back to hash iteration order when multiple motif offsets are valid.
    /// The Rust sidecar deliberately uses annotation as a deterministic
    /// biological tie-break before falling back to the lowest offset.
    pub annotation: Option<&'a Annotation>,
    /// Minimum mapping quality shared with the main BSJ scan.
    ///
    /// Outward reads are weaker circRNA evidence than BSJ/backward reads, so the
    /// post-Summary detector uses the same user-facing MAPQ threshold instead
    /// of the permissive CIRI-AS sidecar constant used for local coverage.
    pub min_mapq: i32,
    /// Optional user-facing progress logger supplied by the CLI.
    ///
    /// The segments phase owns an internal BAM progress bar, but several
    /// correctness-preserving steps run immediately after that bar completes.
    /// Routing status through the caller keeps stdout and `<prefix>.log`
    /// consistent with Scan1/Scan2 without giving this module direct log-file
    /// ownership.
    pub progress_log: Option<&'a mut dyn FnMut(&str, &str) -> Result<()>>,
}

/// User-facing counts produced by the post-Summary segments phase.
///
/// These counts are computed after final segment rows have been selected, so
/// they reflect the actual `<prefix>.segments` output rather than intermediate
/// candidate pools.
#[derive(Debug, Clone, Copy, Default)]
pub struct SegmentRunSummary {
    /// Total rows written to `<prefix>.segments`.
    pub total_segments: usize,
    /// Confirmed Summary-assigned BSJ rows.
    pub bsj_segments: usize,
    /// Non-BSJ mate-chain wrap rows.
    pub backward_segments: usize,
    /// Non-BSJ pair-orientation rows.
    pub outward_segments: usize,
}

/// User-facing counts produced by the major-isoform reconstruction stage.
///
/// The GTF is the complete per-circRNA audit output, while FASTA is stricter
/// and contains only sequence-ready isoforms. Keeping both counts visible helps
/// users distinguish "reconstructed structure exists" from "trusted sequence
/// was emitted".
#[derive(Debug, Clone, Copy, Default)]
pub struct IsoformRunSummary {
    /// Isoform records written to `<prefix>.isoforms.gtf`.
    pub total_isoforms: usize,
    /// Distinct circRNAs represented by GTF isoform records.
    pub circ_rnas: usize,
    /// Sequence records written to `<prefix>.isoforms.fa`.
    pub fasta_isoforms: usize,
    /// Distinct circRNAs represented by FASTA sequence records.
    pub fasta_circ_rnas: usize,
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
/// confirmed `bsj` rows, mate-chain `backward` rows, and pair-orientation
/// `outward` rows; `forward` remains reserved for later circ-span extraction.
#[derive(Clone)]
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
    r1_align_strand: String,
    r2_align_strand: String,
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
    token_strand: char,
    blocks: Vec<SegmentBlock>,
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
    junction_read_count: String,
    gene_id: String,
    strand: String,
}

/// One internal splice edge used by the major-isoform graph.
///
/// The edge represents a transcript-compatible intron from `donor_end` to
/// `acceptor_start` in genomic coordinates. It is intentionally independent of
/// read ID because BSJ, backward, and outward reads are combined as support
/// tiers after the final segment rows are stable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct MajorEdge {
    donor_end: i32,
    acceptor_start: i32,
    strand: char,
}

/// Junction token used to phase the selected major-isoform path.
///
/// The implicit BSJ is part of the circular junction chain. Treating it as a
/// first-class token prevents two-exon circRNAs with an unsupported long block
/// from being labeled `mature` simply because they have only one internal
/// splice edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum MajorJunction {
    Bsj,
    Edge(MajorEdge),
}

/// Pair of adjacent junctions observed in one read-chain segment.
///
/// Individual junction support is not enough to prove a mature full-length
/// chain: two neighboring junctions can be observed by separate reads while the
/// block between them remains unphased. Link support records that one read
/// chain contains both junctions consecutively, including BSJ-to-internal
/// boundary links.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct MajorJunctionLink {
    left: MajorJunction,
    right: MajorJunction,
}

impl MajorJunctionLink {
    /// Builds a directional link key from two adjacent junctions in read-chain order.
    ///
    /// The two BSJ boundaries around a single internal junction must remain
    /// distinguishable. Collapsing `(BSJ, edge)` and `(edge, BSJ)` would let a
    /// read that only covers the short exon between one internal junction and
    /// one BSJ side certify both circular blocks as mature.
    fn new(a: MajorJunction, b: MajorJunction) -> Self {
        Self { left: a, right: b }
    }
}

/// Tiered read support for one major-isoform splice edge.
///
/// BSJ and backward reads are treated as high-confidence phasing evidence;
/// outward reads are retained only as weak completion evidence because they
/// support circular origin but do not localize a BSJ span by themselves.
#[derive(Debug, Clone, Copy, Default)]
struct MajorEdgeSupport {
    bsj: f64,
    backward: f64,
    outward: f64,
}

/// Tiered read support for a neighboring-junction link.
#[derive(Debug, Clone, Copy, Default)]
struct MajorLinkSupport {
    bsj: f64,
    backward: f64,
    outward: f64,
}

/// Continuous aligned segment used as junction-exclusive evidence.
///
/// If a read has one uninterrupted alignment block spanning a candidate intron,
/// that read is evidence that the candidate splice junction was not used in
/// this molecule. The signal is only used to keep annotation-only junctions from
/// over-splitting unphased blocks; read-supported junctions still take priority.
#[derive(Debug, Clone, Copy)]
struct MajorAlignedSpan {
    start: i32,
    end: i32,
    bsj: f64,
    backward: f64,
    outward: f64,
}

/// Weighted circRNA assignment for one `<prefix>.segments` row.
///
/// BSJ rows keep a single exact assignment with weight 1.0. Backward/outward
/// rows can be compatible with several confirmed circRNAs or with none; their
/// weights are probabilities used only by isoform reconstruction so the public
/// read-level segments remain unmodified.
#[derive(Debug, Clone, Copy)]
struct MajorCircAssignment {
    circ_index: usize,
    weight: f64,
}

/// Exact circRNA interval index entry for assigning non-BSJ segment rows.
///
/// `max_end_through` enables backward scans over a start-sorted vector to stop
/// once no earlier interval can contain the query span, which keeps assignment
/// bounded without retaining read-level state for every circRNA.
#[derive(Debug, Clone, Copy)]
struct MajorCircIndexEntry {
    start: i32,
    end: i32,
    circ_index: usize,
    max_end_through: i32,
}

/// Major isoform selected for one Summary-confirmed circRNA.
///
/// The output is deliberately single-isoform for now. `structure_hash` and
/// `sample_id` are emitted so later multi-sample code can compare major
/// isoform switching without changing the file contract again.
#[derive(Debug, Clone)]
struct MajorIsoformRecord {
    circ_id: String,
    isoform_id: String,
    sample_id: String,
    chr: String,
    start: i32,
    end: i32,
    strand: char,
    source_gene_id: String,
    exons: Vec<(i32, i32)>,
    cov: f64,
    segment_coverage_pct: f64,
    path_score: f64,
    bsj_reads: usize,
    structure_hash: u64,
    isoform_len: i32,
    isoform_origin: String,
    estimate_reason: String,
}

/// Result of converting graph blocks into mature exon intervals.
struct MajorExonBuild {
    exons: Vec<(i32, i32)>,
    projected_long_block: bool,
    unresolved_long_block: bool,
    unphased_junction_chain: bool,
    inferred_internal_block: bool,
    unphased_single_exon_block: bool,
}

impl MajorExonBuild {
    /// Returns whether this isoform is fully supported by segment graph blocks.
    fn origin(&self) -> &'static str {
        if self.projected_long_block || self.unresolved_long_block {
            "estimate"
        } else if self.unphased_junction_chain
            || self.inferred_internal_block
            || self.unphased_single_exon_block
        {
            "estimate"
        } else {
            "mature"
        }
    }

    /// Returns a compact reason string for estimate records.
    fn estimate_reason(&self) -> String {
        let mut reasons = Vec::new();
        if self.projected_long_block {
            reasons.push("gtf_long_block_projection");
        }
        if self.unresolved_long_block {
            reasons.push("unresolved_long_block");
        }
        if self.unphased_junction_chain {
            reasons.push("unphased_junction_chain");
        }
        if self.inferred_internal_block {
            reasons.push("inferred_internal_block");
        }
        if self.unphased_single_exon_block {
            reasons.push("unphased_single_exon_block");
        }
        if reasons.is_empty() {
            "none".to_string()
        } else {
            reasons.join(",")
        }
    }
}

/// Per-block result from annotation projection.
#[derive(Debug, Clone, Copy, Default)]
struct MajorBlockBuildStatus {
    projected: bool,
    unresolved: bool,
    inferred: bool,
    unphased_single_exon: bool,
}

/// Column indexes required to rebuild major isoforms from `<prefix>.segments`.
///
/// Header-driven lookup keeps the isoform stage coupled to the public segments
/// contract instead of to the in-memory `SegmentRecord` layout, which is the
/// interface needed for later standalone and multi-sample processing.
#[derive(Debug, Clone, Copy)]
struct MajorSegmentsColumns {
    type_name: usize,
    circ_id: usize,
    chrom: usize,
    start: usize,
    end: usize,
    r1_segments: usize,
    r2_segments: usize,
}

/// CircRNA locus cluster used to define local splice-repair search windows.
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

/// Exact Summary circRNA span used for local clip and splice-repair context.
///
/// Backward/outward rows are read-level topology evidence and do not require a
/// known BSJ. When a read is close to known circ loci, these compact spans still
/// provide bounded windows for annotation and splice-signal guided correction.
#[derive(Debug, Clone, Copy)]
struct CircSpan {
    start: i32,
    end: i32,
    max_end_through: i32,
}

/// Local reference window used to place non-BSJ soft clips.
///
/// Backward/outward reads are not assigned to one circRNA at read level, so clip
/// placement cannot be limited to a unique circ span. These windows stay local
/// to detected circ loci and the read's own mapped blocks, then the normal
/// chain ranking chooses whether any pseudo-alignment is biologically coherent.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct LocalClipWindow {
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

/// Runs the post-Summary read-level circRNA segments phase.
///
/// The default CLI path consumes Scan1/Scan2 sidecar evidence and filters it by
/// Summary-confirmed junction reads, which avoids re-scanning BAM/SAM while
/// preserving CIRI3 `.out` parity. If no sidecar paths are supplied, the older
/// CIRI-AS-style second sweep remains available for focused development and
/// tests; full-length cirexon/path reconstruction stays dormant in this module.
pub fn run_ciri_as(mut config: AsConfig<'_>) -> Result<SegmentRunSummary> {
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
    let (bsj_segment_shard_paths, read_len) = if config.segment_evidence_paths.is_empty() {
        let read_len = infer_read_length(config.input_path)?;
        if read_len < 40 {
            bail!(
                "CIRI-AS requires paired reads with inferred read length >= 40, got {}",
                read_len
            );
        }
        (Vec::new(), read_len)
    } else {
        let (shard_paths, read_len) = load_segment_evidence_spill(
            &config.segment_evidence_paths,
            &junction_read_to_circ,
            config.out_prefix,
            config.keep_temp_files,
        )?;
        if read_len < 40 {
            bail!(
                "segments sidecar read length must be >= 40 for CIRI-AS chain assembly, got {}",
                read_len
            );
        }
        (shard_paths, read_len)
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
        circ_spans_by_chr: circ_spans_by_chr(&circ_records),
        clusters_by_chr: clusters_by_chr(clusters),
        read_len,
        min_mapq: config.min_mapq,
        candidates: Vec::new(),
        outward_read_ids: HashSet::new(),
        prebuilt_segment_records: Vec::new(),
        segment_groups: HashMap::new(),
    };
    log_segments_profile(profile, "build_scan_state", phase_started);

    let phase_started = profile.then(Instant::now);
    let non_bsj_rescan_group_indexes = if config.non_bsj_segment_evidence_paths.is_empty() {
        None
    } else {
        log_segments_user_progress(
            &mut config,
            "Scanning exons 1/4",
            "Determining non-BSJ read topology...",
        )?;
        let non_bsj_paths = config.non_bsj_segment_evidence_paths.clone();
        let out_prefix = config.out_prefix;
        Some(load_non_bsj_segment_evidence(
            &non_bsj_paths,
            out_prefix,
            &mut state,
            || {
                log_segments_user_progress(
                    &mut config,
                    "Scanning exons 2/4",
                    "Preparing non-BSJ segments...",
                )
            },
        )?)
    };
    log_segments_profile(profile, "load_non_bsj_segment_evidence", phase_started);

    let phase_started = profile.then(Instant::now);
    let segments = if config.segment_evidence_paths.is_empty() {
        scan_alignment_groups(config.input_path, &mut state)?;
        validate_splice_motifs(&mut state.candidates, config.reference, config.annotation)?;
        build_segment_records(&state)?
    } else {
        let rescan_group_shard_paths = if let Some(paths) = non_bsj_rescan_group_indexes {
            paths
        } else {
            let phase_started = profile.then(Instant::now);
            let rescan_group_shards =
                scan_backward_alignment_groups(config.input_path, config.out_prefix, &mut state)?;
            log_segments_profile(profile, "scan_backward_alignment_groups", phase_started);
            rescan_group_shards.unwrap_or_default()
        };
        log_segments_user_progress(
            &mut config,
            "Scanning exons 3/4",
            "Collecting junction support from retained segments...",
        )?;
        let phase_started = profile.then(Instant::now);
        let keep_temp_files = config.keep_temp_files;
        let records = build_sidecar_segment_records(
            &mut state,
            bsj_segment_shard_paths,
            rescan_group_shard_paths,
            keep_temp_files,
            || {
                log_segments_user_progress(
                    &mut config,
                    "Scanning exons 4/4",
                    "Correcting ambiguous segments with junction support...",
                )
            },
        )?;
        log_segments_profile(profile, "build_sidecar_segment_records", phase_started);
        records
    };
    log_segments_profile(profile, "build_segments", phase_started);
    let summary = segment_run_summary(&segments);
    drop(state);
    let phase_started = profile.then(Instant::now);
    let segments_path = format!("{}.segments", config.out_prefix);
    write_segments(&segments_path, &segments)?;
    log_segments_profile(profile, "write_segments", phase_started);
    drop(segments);
    Ok(summary)
}

/// Rebuilds major circRNA isoform sidecars from completed `.out` and `.segments` files.
///
/// This is the resumable checkpoint used by `ciri --continue` when
/// `<prefix>.segments` already exists. It deliberately starts from the public
/// merged segment file instead of any shard-local temporary files, so isoform
/// debugging can rerun without depending on partially written parallel state.
pub fn rebuild_major_isoforms_from_segments(
    circ_path: &str,
    segments_path: &str,
    out_prefix: &str,
    reference: &HashMap<String, String>,
    annotation: Option<&Annotation>,
) -> Result<IsoformRunSummary> {
    let (circ_records, _junction_read_to_circ) = load_circ_records(circ_path)?;
    build_major_isoforms_from_segments_file(
        &circ_records,
        segments_path,
        out_prefix,
        reference,
        annotation,
    )
}

/// Summarizes final segment row counts by evidence type.
fn segment_run_summary(records: &[SegmentRecord]) -> SegmentRunSummary {
    let mut summary = SegmentRunSummary {
        total_segments: records.len(),
        ..SegmentRunSummary::default()
    };
    for record in records {
        match record.type_name {
            "bsj" => summary.bsj_segments += 1,
            "backward" => summary.backward_segments += 1,
            "outward" => summary.outward_segments += 1,
            _ => {}
        }
    }
    summary
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

/// Emits one user-facing segments progress line through the caller's logger.
fn log_segments_user_progress(config: &mut AsConfig<'_>, label: &str, message: &str) -> Result<()> {
    if let Some(logger) = config.progress_log.as_deref_mut() {
        logger(label, message)?;
    }
    Ok(())
}

/// Builds a byte-progress bar for retained segment shard passes.
///
/// The support and correction passes walk temporary `R/A/C/S` shard streams
/// instead of the original BAM. Summing shard byte sizes gives users the same
/// progress expectation as Scan1/Scan2 without retaining read-level state just
/// for reporting.
fn retained_segment_progress_bar(
    bsj_shards: &[SegmentScanShardPaths],
    non_bsj_shards: &[SegmentScanShardPaths],
    message: &'static str,
) -> Result<Option<ProgressBar>> {
    let total_bytes = retained_segment_shard_bytes(bsj_shards, non_bsj_shards)?;
    if total_bytes == 0 {
        return Ok(None);
    }
    let pb = ProgressBar::new(total_bytes);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {percent:>3}% ({eta}) {msg}")?
            .progress_chars("#>-"),
    );
    pb.set_message(message);
    pb.enable_steady_tick(Duration::from_millis(120));
    Ok(Some(pb))
}

/// Returns the total on-disk size of retained segment shards.
fn retained_segment_shard_bytes(
    bsj_shards: &[SegmentScanShardPaths],
    non_bsj_shards: &[SegmentScanShardPaths],
) -> Result<u64> {
    bsj_shards
        .iter()
        .chain(non_bsj_shards)
        .try_fold(0_u64, |acc, shard| {
            Ok(acc + std::fs::metadata(&shard.path)?.len())
        })
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
    circ_spans_by_chr: HashMap<String, Vec<CircSpan>>,
    clusters_by_chr: HashMap<String, Vec<CircCluster>>,
    read_len: i32,
    min_mapq: i32,
    candidates: Vec<PositiveCandidate>,
    outward_read_ids: HashSet<String>,
    prebuilt_segment_records: Vec<SegmentRecord>,
    segment_groups: HashMap<String, Vec<AsAlignment>>,
}

/// Read-only inputs shared by parallel segments BAM shards.
///
/// The third-pass segments scan needs the same circRNA indexes in every shard,
/// but cloning the full `ScanState` would multiply the already-loaded BSJ
/// sidecar evidence. This view keeps shard workers read-only and lets each
/// worker return only newly discovered backward/outward evidence.
struct SegmentScanContext<'a> {
    junction_read_to_circ: &'a HashMap<String, String>,
    reference: &'a HashMap<String, String>,
    annotation: Option<&'a Annotation>,
    circ_spans_by_chr: &'a HashMap<String, Vec<CircSpan>>,
    clusters_by_chr: &'a HashMap<String, Vec<CircCluster>>,
    read_len: i32,
    min_mapq: i32,
}

/// Shard-local mutable output from the segments BAM rescan.
///
/// Keeping this separate from `ScanState` makes BGZF sharding deterministic:
/// every read group is owned by exactly one shard, and the main thread performs
/// the only merge into the final state.
#[derive(Default)]
struct SegmentScanAccum {
    candidates: Vec<PositiveCandidate>,
    segment_groups: HashMap<String, Vec<AsAlignment>>,
    outward_read_ids: HashSet<String>,
    prebuilt_segment_records: Vec<SegmentRecord>,
}

/// Path for one shard-local segments rescan spill.
struct SegmentScanShardPaths {
    path: String,
}

/// One read-complete retained segment shard block.
///
/// The retained `R/A/C/S` protocol is the scalability boundary for segments
/// finalize: each block contains the alignments and local splice hints needed
/// to rebuild that read without a global read-id index.
struct RetainedSegmentGroup {
    read_id: String,
    records: Vec<AsAlignment>,
    candidates: Vec<PositiveCandidate>,
    segment_records: Vec<SegmentRecord>,
}

/// Parsed Scan1/Scan2 segment-evidence row before shard-local aggregation.
///
/// BSJ sidecars can be large enough that keeping all parsed rows in a global
/// read-keyed map dominates RSS. This owned row is written to hash partitions
/// first, then each partition is aggregated independently into the same retained
/// `R/A` stream used by non-BSJ finalize.
struct SegmentEvidenceRow {
    read_id: String,
    alignment: AsAlignment,
    read_len: i32,
}

/// Streaming writer for one segments rescan shard.
///
/// Shards can discover many non-BSJ groups on large BAMs. Writing those groups
/// immediately follows the Scan1/Scan2 pattern and prevents parallel workers
/// from retaining all candidate evidence in memory until the slowest shard
/// finishes.
struct SegmentScanShardWriter {
    writer: BufWriter<File>,
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
            junction_read_count: get(junc_count_idx)?.to_string(),
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

/// Spills Scan1/Scan2 mapper-block evidence for Summary-confirmed BSJ reads.
///
/// The sidecar rows are deliberately compact:
/// `read_id, stage, mate, flag, chrom, pos, mapq, cigar, read_len, clips`. Only
/// `read_id` and the mapper fields are mandatory; `clips` stores compact
/// `L:<seq>,R:<seq>` soft-clip subsequences for traceability, while accepted
/// local clip placements arrive as normal-looking `scan*_local` pseudo rows.
/// The important memory boundary is that rows are hash-partitioned by read ID
/// first, then each partition is aggregated independently. This mirrors the
/// Scan1/Scan2 spill pattern and avoids a whole-run
/// `HashMap<read_id, Vec<AsAlignment>>` during segments finalize.
fn load_segment_evidence_spill(
    paths: &[&str],
    confirmed_reads: &HashMap<String, String>,
    out_prefix: &str,
    keep_temp_files: bool,
) -> Result<(Vec<SegmentScanShardPaths>, i32)> {
    let partition_count = segment_spill_partition_count();
    let raw_paths: Vec<SegmentScanShardPaths> = (0..partition_count)
        .map(|idx| SegmentScanShardPaths::new_bsj_raw(out_prefix, idx))
        .collect();
    let mut writers = raw_paths
        .iter()
        .map(|paths| {
            Ok(BufWriter::with_capacity(
                512 * 1024,
                File::create(&paths.path)?,
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut read_len = 0_i32;
    for path in paths {
        let file = File::open(path).with_context(|| format!("open segment evidence {}", path))?;
        let reader = BufReader::new(file);
        for line_res in reader.lines() {
            let line = line_res?;
            if line.trim().is_empty() {
                continue;
            }
            let Some(row) = parse_segment_evidence_row(&line, confirmed_reads)? else {
                continue;
            };
            read_len = read_len.max(row.read_len);
            let idx = read_id_partition(&row.read_id, partition_count);
            write_segment_evidence_raw_row(&mut writers[idx], &row)?;
        }
    }
    for writer in &mut writers {
        writer.flush()?;
    }
    drop(writers);

    let shard_paths = raw_paths
        .par_iter()
        .enumerate()
        .map(|(idx, raw_path)| {
            let retained_path = SegmentScanShardPaths::new_bsj(out_prefix, idx);
            aggregate_segment_evidence_partition(&raw_path.path, &retained_path)?;
            if !keep_temp_files {
                let _ = std::fs::remove_file(&raw_path.path);
            }
            Ok(retained_path)
        })
        .collect::<Result<Vec<_>>>()?;
    Ok((shard_paths, read_len))
}

/// Chooses how many bounded BSJ sidecar partitions to use.
fn segment_spill_partition_count() -> usize {
    rayon::current_num_threads().saturating_mul(4).clamp(1, 256)
}

/// Returns the stable hash partition for one read ID.
fn read_id_partition(read_id: &str, partition_count: usize) -> usize {
    let mut hasher = DefaultHasher::new();
    read_id.hash(&mut hasher);
    (hasher.finish() as usize) % partition_count.max(1)
}

/// Parses one Scan1/Scan2 segment-evidence line.
fn parse_segment_evidence_row(
    line: &str,
    confirmed_reads: &HashMap<String, String>,
) -> Result<Option<SegmentEvidenceRow>> {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 9 {
        return Ok(None);
    }
    let read_id = parts[0];
    if !confirmed_reads.contains_key(read_id) {
        return Ok(None);
    }
    let flag = parts[3]
        .parse::<i32>()
        .with_context(|| format!("parse segment evidence flag from {}", line))?;
    let pos = parts[5]
        .parse::<i32>()
        .with_context(|| format!("parse segment evidence position from {}", line))?;
    let mapq = parts[6]
        .parse::<i32>()
        .with_context(|| format!("parse segment evidence MAPQ from {}", line))?;
    let read_len = parts[8]
        .parse::<i32>()
        .with_context(|| format!("parse segment evidence read length from {}", line))?;
    let seq = parts.get(9).copied().unwrap_or("*");
    Ok(Some(SegmentEvidenceRow {
        read_id: read_id.to_string(),
        alignment: AsAlignment {
            flag,
            chr: parts[4].to_string(),
            pos,
            mapq,
            cigar: parts[7].to_string(),
            seq: if seq == "*" {
                String::new()
            } else {
                seq.to_string()
            },
            from_local_clip: parts[1].contains("_local"),
            xa_alternatives: Vec::new(),
        },
        read_len,
    }))
}

/// Writes one parsed BSJ segment-evidence row to its raw hash partition.
fn write_segment_evidence_raw_row(
    writer: &mut BufWriter<File>,
    row: &SegmentEvidenceRow,
) -> Result<()> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        row.read_id,
        row.alignment.flag,
        row.alignment.chr,
        row.alignment.pos,
        row.alignment.mapq,
        row.alignment.cigar,
        if row.alignment.seq.is_empty() {
            "*"
        } else {
            &row.alignment.seq
        },
        u8::from(row.alignment.from_local_clip),
        row.read_len
    )?;
    Ok(())
}

/// Aggregates one raw BSJ hash partition into a retained `R/A` shard.
fn aggregate_segment_evidence_partition(
    raw_path: &str,
    retained_path: &SegmentScanShardPaths,
) -> Result<()> {
    let file = File::open(raw_path)
        .with_context(|| format!("open segment evidence partition {raw_path}"))?;
    let reader = BufReader::new(file);
    let mut groups: HashMap<String, Vec<AsAlignment>> = HashMap::new();
    for line_res in reader.lines() {
        let line = line_res?;
        if line.trim().is_empty() {
            continue;
        }
        let mut parts = line.split('\t');
        let Some(read_id) = parts.next() else {
            continue;
        };
        let Some(flag) = parts.next().and_then(|value| value.parse::<i32>().ok()) else {
            continue;
        };
        let Some(chr) = parts.next() else {
            continue;
        };
        let Some(pos) = parts.next().and_then(|value| value.parse::<i32>().ok()) else {
            continue;
        };
        let Some(mapq) = parts.next().and_then(|value| value.parse::<i32>().ok()) else {
            continue;
        };
        let Some(cigar) = parts.next() else {
            continue;
        };
        let Some(seq) = parts.next() else {
            continue;
        };
        let from_local_clip = parts.next().is_some_and(|value| value == "1");
        groups
            .entry(read_id.to_string())
            .or_default()
            .push(AsAlignment {
                flag,
                chr: chr.to_string(),
                pos,
                mapq,
                cigar: cigar.to_string(),
                seq: if seq == "*" {
                    String::new()
                } else {
                    seq.to_string()
                },
                from_local_clip,
                xa_alternatives: Vec::new(),
            });
    }
    let mut writer = SegmentScanShardWriter::new(retained_path)?;
    for (read_id, mut records) in groups {
        records.sort_by(|a, b| {
            mate_bucket(a.flag)
                .cmp(&mate_bucket(b.flag))
                .then_with(|| as_alignment_order_key(a).cmp(&as_alignment_order_key(b)))
                .then_with(|| a.chr.cmp(&b.chr))
                .then_with(|| a.pos.cmp(&b.pos))
                .then_with(|| a.cigar.cmp(&b.cigar))
        });
        // Duplicate mapper rows are local to a read group. Keeping one global
        // string key per sidecar row can dominate RSS on large BAMs, so dedupe
        // after grouping with a short-lived per-read key set instead.
        let mut seen = HashSet::with_capacity(records.len());
        records.retain(|record| {
            seen.insert((
                record.flag,
                record.chr.clone(),
                record.pos,
                record.mapq,
                record.cigar.clone(),
                record.seq.clone(),
            ))
        });
        writer.write_read_result(&read_id, &records, false, &[], &[])?;
    }
    writer.flush()?;
    Ok(())
}

/// Loads Scan2 read-level non-BSJ topology candidates.
///
/// These rows are deliberately re-evaluated after Summary instead of trusted as
/// final calls. Scan2 only knows read-level topology; this loader skips final
/// BSJ reads and then runs the same backward/outward materializer used by the
/// legacy third BAM scan. The sequence column is the compact `L:/R:` soft-clip
/// payload, matching the BSJ sidecar strategy to keep I/O bounded.
fn load_non_bsj_segment_evidence(
    paths: &[&str],
    out_prefix: &str,
    state: &mut ScanState,
    mut on_spill_done: impl FnMut() -> Result<()>,
) -> Result<Vec<SegmentScanShardPaths>> {
    let shard_paths = {
        let ctx = SegmentScanContext::from_state(state);
        let grouped_inputs = paths
            .iter()
            .map(|path| non_bsj_sidecar_is_grouped(path))
            .collect::<Result<Vec<_>>>()?;
        if paths.len() > 1 && grouped_inputs.iter().all(|is_grouped| *is_grouped) {
            spill_grouped_non_bsj_segment_evidence_shards(paths, out_prefix, 0, &ctx)?
        } else {
            let mut shard_paths = Vec::new();
            for (path, is_grouped) in paths.iter().zip(grouped_inputs) {
                if is_grouped {
                    let mut paths = spill_grouped_non_bsj_segment_evidence_mmap(
                        path,
                        out_prefix,
                        shard_paths.len(),
                        &ctx,
                    )?;
                    shard_paths.append(&mut paths);
                } else {
                    let mut paths = spill_legacy_non_bsj_segment_evidence(
                        path,
                        out_prefix,
                        shard_paths.len(),
                        &ctx,
                    )?;
                    shard_paths.append(&mut paths);
                }
            }
            shard_paths
        }
    };
    on_spill_done()?;
    Ok(shard_paths)
}

/// Tests whether a non-BSJ sidecar uses the current read-group line protocol.
fn non_bsj_sidecar_is_grouped(path: &str) -> Result<bool> {
    let file =
        File::open(path).with_context(|| format!("open non-BSJ segment evidence {}", path))?;
    let reader = BufReader::new(file);
    for line_res in reader.lines() {
        let line = line_res?;
        if line.trim().is_empty() {
            continue;
        }
        return Ok(line.split('\t').nth(1).is_some_and(is_non_bsj_group_stage));
    }
    Ok(true)
}

/// Returns true for current compact and older verbose non-BSJ group stages.
fn is_non_bsj_group_stage(stage: &str) -> bool {
    stage == "N" || stage == "scan2_non_bsj_group"
}

/// Loads the current one-line-per-read non-BSJ sidecar with sharded mmap parsing.
///
/// Whole-genome sidecars can be tens of GiB. Mapping the file keeps memory
/// bounded by page cache instead of copying the file into heap memory, while
/// newline-aligned shards let read groups be parsed and materialized in parallel.
fn spill_grouped_non_bsj_segment_evidence_mmap(
    path: &str,
    out_prefix: &str,
    shard_offset: usize,
    ctx: &SegmentScanContext<'_>,
) -> Result<Vec<SegmentScanShardPaths>> {
    let file =
        File::open(path).with_context(|| format!("open non-BSJ segment evidence {}", path))?;
    let mmap = unsafe { Mmap::map(&file)? };
    if mmap.is_empty() {
        return Ok(Vec::new());
    }
    unsafe {
        libc::madvise(
            mmap.as_ptr() as *mut libc::c_void,
            mmap.len(),
            libc::MADV_SEQUENTIAL,
        );
    }
    let pb = segment_scan_progress_bar(path)?;
    pb.set_message("");
    let ranges = line_aligned_shard_ranges(&mmap, rayon::current_num_threads().max(1));
    let shard_paths = ranges
        .into_par_iter()
        .enumerate()
        .map(|(idx, (start, end))| {
            let paths = SegmentScanShardPaths::new(out_prefix, shard_offset + idx);
            process_grouped_non_bsj_sidecar_shard(&mmap, start, end, ctx, &pb, &paths)?;
            Ok(paths)
        })
        .collect::<Result<Vec<_>>>()?;
    pb.set_position(mmap.len() as u64);

    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {msg}")?
            .progress_chars("#>-"),
    );
    pb.finish_with_message("Completed");
    Ok(shard_paths)
}

/// Spills already sharded grouped non-BSJ sidecars into retained read streams.
///
/// Scan2 writes one sidecar per BAM shard. Each source shard is already
/// read-group complete, so finalize can process one raw shard into one retained
/// `R/A/C/S` shard without first merging the raw 10+ GiB sidecar or keeping the
/// retained alignment groups in heap memory. The later merge stage loads only
/// preliminary `R/C/S` metadata and seeks back to `A` rows for the small
/// ambiguous subset that needs support-aware correction.
fn spill_grouped_non_bsj_segment_evidence_shards(
    paths: &[&str],
    out_prefix: &str,
    shard_offset: usize,
    ctx: &SegmentScanContext<'_>,
) -> Result<Vec<SegmentScanShardPaths>> {
    let total_bytes = paths.iter().try_fold(0_u64, |acc, path| {
        Ok::<u64, anyhow::Error>(acc + std::fs::metadata(path)?.len())
    })?;
    let pb = ProgressBar::new(total_bytes);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {percent:>3}% ({eta}) {msg}",
            )?
            .progress_chars("#>-"),
    );
    pb.set_message("");

    let shard_paths = paths
        .par_iter()
        .enumerate()
        .map(|(idx, path)| {
            let output_paths = SegmentScanShardPaths::new(out_prefix, shard_offset + idx);
            process_grouped_non_bsj_sidecar_file(path, ctx, &pb, &output_paths)?;
            Ok(output_paths)
        })
        .collect::<Result<Vec<_>>>()?;
    pb.set_position(total_bytes);

    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {msg}")?
            .progress_chars("#>-"),
    );
    pb.finish_with_message("Completed");
    Ok(shard_paths)
}

/// Streams one grouped non-BSJ Scan2 sidecar into a retained segment shard.
///
/// Scan2 already writes these inputs as read-complete shard files. Re-mapping
/// every 10+ GiB source shard in parallel can make RSS look like heap growth on
/// real BAMs because the kernel charges touched mapped pages to the process.
/// A buffered line reader keeps the same one-read-at-a-time semantics while
/// bounding resident memory to the reader and writer buffers for each worker.
fn process_grouped_non_bsj_sidecar_file(
    path: &str,
    ctx: &SegmentScanContext<'_>,
    pb: &ProgressBar,
    paths: &SegmentScanShardPaths,
) -> Result<()> {
    let file = File::open(path).with_context(|| format!("open non-BSJ segment evidence {path}"))?;
    if file.metadata()?.len() == 0 {
        SegmentScanShardWriter::new(paths)?.flush()?;
        return Ok(());
    }

    let mut reader = BufReader::with_capacity(8 * 1024 * 1024, file);
    let mut writer = SegmentScanShardWriter::new(paths)?;
    let mut line = String::new();
    let mut pending_progress = 0usize;
    loop {
        line.clear();
        let bytes = reader.read_line(&mut line)?;
        if bytes == 0 {
            break;
        }
        pending_progress += bytes;
        let line = line.trim_end_matches(['\n', '\r']);
        if !line.is_empty() {
            if let Some((read_id, records)) = parse_non_bsj_group_line(line, ctx)? {
                process_backward_group_to_writer(read_id, &records, ctx, &mut writer)?;
            }
        }
        if pending_progress >= 4 * 1024 * 1024 {
            pb.inc(pending_progress as u64);
            pending_progress = 0;
        }
    }
    if pending_progress > 0 {
        pb.inc(pending_progress as u64);
    }
    writer.flush()?;
    Ok(())
}

/// Splits a line-oriented mmap into non-overlapping newline-aligned ranges.
fn line_aligned_shard_ranges(mmap: &Mmap, wanted_shards: usize) -> Vec<(usize, usize)> {
    let shard_count = wanted_shards.max(1).min(mmap.len().max(1));
    let mut starts = Vec::with_capacity(shard_count + 1);
    starts.push(0);
    for idx in 1..shard_count {
        let raw = idx * mmap.len() / shard_count;
        starts.push(next_line_start(mmap, raw));
    }
    starts.push(mmap.len());
    starts.dedup();
    starts
        .windows(2)
        .filter_map(|pair| (pair[0] < pair[1]).then_some((pair[0], pair[1])))
        .collect()
}

/// Returns the first byte after the line containing `pos`.
fn next_line_start(mmap: &Mmap, pos: usize) -> usize {
    if pos == 0 || pos >= mmap.len() {
        return pos.min(mmap.len());
    }
    if mmap[pos - 1] == b'\n' {
        return pos;
    }
    let mut idx = pos;
    while idx < mmap.len() {
        if mmap[idx] == b'\n' {
            return idx + 1;
        }
        idx += 1;
    }
    mmap.len()
}

/// Processes one mmap shard of grouped non-BSJ read evidence.
fn process_grouped_non_bsj_sidecar_shard(
    mmap: &Mmap,
    start: usize,
    end: usize,
    ctx: &SegmentScanContext<'_>,
    pb: &ProgressBar,
    paths: &SegmentScanShardPaths,
) -> Result<()> {
    let mut writer = SegmentScanShardWriter::new(paths)?;
    let mut pos = start;
    let mut last_progress = start;
    while pos < end {
        let line_start = pos;
        while pos < end && mmap[pos] != b'\n' {
            pos += 1;
        }
        let mut line_end = pos;
        if line_end > line_start && mmap[line_end - 1] == b'\r' {
            line_end -= 1;
        }
        if pos < end {
            pos += 1;
        }
        if line_end > line_start {
            let line = std::str::from_utf8(&mmap[line_start..line_end]).with_context(|| {
                format!("parse UTF-8 non-BSJ sidecar line at byte {line_start}")
            })?;
            if let Some((read_id, records)) = parse_non_bsj_group_line(line, ctx)? {
                process_backward_group_to_writer(&read_id, &records, ctx, &mut writer)?;
            }
        }
        if pos.saturating_sub(last_progress) >= 4 * 1024 * 1024 {
            pb.inc((pos - last_progress) as u64);
            last_progress = pos;
        }
    }
    if end > last_progress {
        pb.inc((end - last_progress) as u64);
    }
    writer.flush()?;
    advise_segment_mmap_dontneed(mmap, start, end - start);
    Ok(())
}

/// Parses one current-format grouped non-BSJ sidecar line.
fn parse_non_bsj_group_line<'a>(
    line: &'a str,
    ctx: &SegmentScanContext<'_>,
) -> Result<Option<(&'a str, Vec<AsAlignment>)>> {
    let mut parts = line.splitn(3, '\t');
    let Some(read_id) = parts.next() else {
        return Ok(None);
    };
    let Some(stage) = parts.next() else {
        return Ok(None);
    };
    let Some(payload) = parts.next() else {
        return Ok(None);
    };
    if !is_non_bsj_group_stage(stage) || ctx.junction_read_to_circ.contains_key(read_id) {
        return Ok(None);
    }
    let records = parse_non_bsj_group_records(payload, line)?;
    Ok(Some((read_id, records)))
}

/// Loads legacy one-alignment-per-line non-BSJ sidecars with visible progress.
fn spill_legacy_non_bsj_segment_evidence(
    path: &str,
    out_prefix: &str,
    shard_idx: usize,
    ctx: &SegmentScanContext<'_>,
) -> Result<Vec<SegmentScanShardPaths>> {
    let file =
        File::open(path).with_context(|| format!("open non-BSJ segment evidence {}", path))?;
    let mut reader = BufReader::new(file);
    let pb = segment_scan_progress_bar(path)?;
    pb.set_message("");
    let paths = SegmentScanShardPaths::new(out_prefix, shard_idx);
    let mut writer = SegmentScanShardWriter::new(&paths)?;
    let mut last_progress_pos = 0_u64;
    let mut line = String::new();
    let mut bytes_since_progress = 0_u64;
    {
        let mut current_read_id = String::new();
        let mut current_records = Vec::with_capacity(8);
        let mut seen_in_group = HashSet::new();
        loop {
            line.clear();
            let bytes = reader.read_line(&mut line)?;
            if bytes == 0 {
                break;
            }
            bytes_since_progress += bytes as u64;
            if bytes_since_progress >= 4 * 1024 * 1024 {
                last_progress_pos += bytes_since_progress;
                pb.inc(bytes_since_progress);
                bytes_since_progress = 0;
            }
            let line = line.trim_end_matches(['\n', '\r']);
            if line.is_empty() {
                continue;
            }
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 && is_non_bsj_group_stage(parts[1]) {
                if !current_read_id.is_empty() {
                    process_backward_group_to_writer(
                        &current_read_id,
                        &current_records,
                        ctx,
                        &mut writer,
                    )?;
                    current_read_id.clear();
                    current_records.clear();
                    seen_in_group.clear();
                }
                let read_id = parts[0];
                if ctx.junction_read_to_circ.contains_key(read_id) {
                    continue;
                }
                let records = parse_non_bsj_group_records(parts[2], &line)?;
                process_backward_group_to_writer(read_id, &records, ctx, &mut writer)?;
                continue;
            }
            if parts.len() < 10 {
                continue;
            }
            let read_id = parts[0];
            if ctx.junction_read_to_circ.contains_key(read_id) {
                continue;
            }
            if !current_read_id.is_empty() && read_id != current_read_id {
                process_backward_group_to_writer(
                    &current_read_id,
                    &current_records,
                    ctx,
                    &mut writer,
                )?;
                current_records.clear();
                seen_in_group.clear();
            }
            if current_read_id != read_id {
                current_read_id.clear();
                current_read_id.push_str(read_id);
            }

            let flag = parts[3]
                .parse::<i32>()
                .with_context(|| format!("parse non-BSJ segment flag from {}", line))?;
            let chrom = parts[4];
            let pos = parts[5]
                .parse::<i32>()
                .with_context(|| format!("parse non-BSJ segment position from {}", line))?;
            let mapq = parts[6]
                .parse::<i32>()
                .with_context(|| format!("parse non-BSJ segment MAPQ from {}", line))?;
            let cigar = parts[7];
            let seq = parts[9];
            let key = format!("{flag}\t{chrom}\t{pos}\t{mapq}\t{cigar}\t{seq}");
            if !seen_in_group.insert(key) {
                continue;
            }
            current_records.push(AsAlignment {
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
        if !current_read_id.is_empty() {
            process_backward_group_to_writer(&current_read_id, &current_records, ctx, &mut writer)?;
        }
    }
    writer.flush()?;
    if bytes_since_progress > 0 {
        last_progress_pos += bytes_since_progress;
        pb.inc(bytes_since_progress);
    }
    let total = std::fs::metadata(path)?.len();
    if total > last_progress_pos {
        pb.inc(total - last_progress_pos);
    }
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {msg}")?
            .progress_chars("#>-"),
    );
    pb.finish_with_message("Completed");
    Ok(vec![paths])
}

/// Parses the compact one-line Scan2 non-BSJ read-group sidecar encoding.
fn parse_non_bsj_group_records(payload: &str, line: &str) -> Result<Vec<AsAlignment>> {
    let mut records = Vec::new();
    let mut seen = HashSet::new();
    for encoded in payload.split(';') {
        if encoded.is_empty() {
            continue;
        }
        // This parser is on the 10+ GiB non-BSJ sidecar path. The current
        // compact protocol has six fields, while the older verbose protocol
        // has eight fields; keeping the first eight slots on the stack avoids
        // one short Vec allocation for every alignment payload.
        let mut fields = [None; 8];
        let mut field_count = 0usize;
        for field in encoded.split('|') {
            if field_count < fields.len() {
                fields[field_count] = Some(field);
            }
            field_count += 1;
            if field_count >= fields.len() {
                break;
            }
        }
        let (flag_idx, chrom_idx, pos_idx, mapq_idx, cigar_idx, seq_idx) = match field_count {
            0..=5 => continue,
            6 => (0, 1, 2, 3, 4, 5),
            7 => continue,
            _ => (1, 2, 3, 4, 5, 7),
        };
        let field = |idx: usize| fields[idx].expect("validated non-BSJ sidecar field");
        let flag_text = field(flag_idx);
        let chrom = field(chrom_idx);
        let pos_text = field(pos_idx);
        let mapq_text = field(mapq_idx);
        let cigar = field(cigar_idx);
        let seq = field(seq_idx);
        let flag = flag_text
            .parse::<i32>()
            .with_context(|| format!("parse non-BSJ grouped flag from {}", line))?;
        let pos = pos_text
            .parse::<i32>()
            .with_context(|| format!("parse non-BSJ grouped position from {}", line))?;
        let mapq = mapq_text
            .parse::<i32>()
            .with_context(|| format!("parse non-BSJ grouped MAPQ from {}", line))?;
        let key = (flag, chrom, pos, mapq, cigar, seq);
        if !seen.insert(key) {
            continue;
        }
        records.push(AsAlignment {
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
            from_local_clip: false,
            xa_alternatives: Vec::new(),
        });
    }
    Ok(records)
}

/// Returns recoverable soft-clip payload for one segments alignment.
fn segment_clip_payload(record: &AsAlignment) -> String {
    if record.seq == "*" || record.seq.starts_with("L:") || record.seq.starts_with("R:") {
        record.seq.clone()
    } else {
        clip_sequence_payload(&record.cigar, &record.seq)
    }
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

/// Indexes exact Summary circRNA spans by chromosome for outward-pair gating.
fn circ_spans_by_chr(records: &[CircRecord]) -> HashMap<String, Vec<CircSpan>> {
    let mut map: HashMap<String, Vec<CircSpan>> = HashMap::new();
    for record in records {
        map.entry(record.chr.clone()).or_default().push(CircSpan {
            start: record.start,
            end: record.end,
            max_end_through: record.end,
        });
    }
    for spans in map.values_mut() {
        spans.sort_by_key(|span| (span.start, span.end));
        let mut max_end = i32::MIN;
        for span in spans {
            max_end = max_end.max(span.end);
            span.max_end_through = max_end;
        }
    }
    map
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
    for_group_records(path, false, |read_id, records| {
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
    for_group_records(path, true, |read_id, records| {
        process_group(read_id, records, state)?;
        Ok(false)
    })
}

/// Scans only non-BSJ read groups for backward/circular candidates.
///
/// The default segments path already receives Summary-confirmed BSJ mapper
/// blocks from Scan1/Scan2 sidecars. This pass deliberately skips those read IDs
/// so it can supplement `<prefix>.segments` with `type=backward` rows without
/// overwriting the richer BSJ sidecar evidence or changing Summary parity. BAM
/// shards are returned before candidate merging so the CLI can report the
/// post-progress finalize phase immediately after the scan progress bar ends.
fn scan_backward_alignment_groups(
    path: &str,
    out_prefix: &str,
    state: &mut ScanState,
) -> Result<Option<Vec<SegmentScanShardPaths>>> {
    let ctx = SegmentScanContext::from_state(state);
    match detect_format(path)? {
        InputFormat::Sam => {
            let mut accum = SegmentScanAccum::default();
            for_group_records(path, true, |read_id, records| {
                process_backward_group(read_id, records, &ctx, &mut accum)?;
                Ok(false)
            })?;
            merge_segment_scan_accum(state, accum);
            Ok(None)
        }
        InputFormat::Bam => {
            let shard_paths = scan_backward_bam_groups_parallel(path, out_prefix, &ctx)?;
            Ok(Some(shard_paths))
        }
    }
}

/// Iterates queryname-sorted alignment groups from either SAM or BAM.
///
/// Both formats feed the same lightweight `AsAlignment` records so the CIRI-AS
/// decision code remains format-agnostic. BAM support is sequential in this first
/// parity stage; the state boundary above is deliberately ready for sharded
/// merging once splice output is aligned to the Perl reference.
fn for_group_records<F>(path: &str, show_progress: bool, mut on_group: F) -> Result<()>
where
    F: FnMut(&str, &[AsAlignment]) -> Result<bool>,
{
    let pb = if show_progress {
        Some(segment_scan_progress_bar(path)?)
    } else {
        None
    };
    let result = match detect_format(path)? {
        InputFormat::Sam => scan_sam_groups(path, &mut on_group, pb.as_ref()),
        InputFormat::Bam => scan_bam_groups(path, &mut on_group, pb.as_ref()),
    };
    if let Some(pb) = pb {
        if result.is_ok() {
            if let Some(len) = pb.length() {
                pb.set_position(len);
            }
            pb.set_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {msg}")?
                    .progress_chars("#>-"),
            );
            pb.finish_with_message("Completed");
        } else {
            pb.abandon();
        }
    }
    result
}

/// Creates the byte-progress display used by the segments BAM/SAM sweep.
///
/// This sweep can dominate whole-genome runtime because it re-reads the original
/// queryname-sorted alignment file to collect non-BSJ backward/outward evidence.
/// Reusing the Scan1/Scan2 progress style makes that cost visible without
/// changing the parser or evidence semantics.
fn segment_scan_progress_bar(path: &str) -> Result<ProgressBar> {
    let file_size = std::fs::metadata(path)?.len();
    let pb = ProgressBar::new(file_size);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {percent:>3}% ({eta}) {msg}")?
            .progress_chars("#>-"),
    );
    pb.set_message("");
    pb.enable_steady_tick(Duration::from_millis(120));
    Ok(pb)
}

/// Updates an optional byte-progress bar from the reader position.
fn update_segment_progress<R: Seek>(
    reader: &mut R,
    pb: Option<&ProgressBar>,
    last_progress_pos: &mut u64,
) -> Result<()> {
    let pos = reader.stream_position()?;
    update_segment_progress_position(pos, pb, last_progress_pos);
    Ok(())
}

/// Applies a monotonic byte-position update to an optional progress bar.
fn update_segment_progress_position(
    pos: u64,
    pb: Option<&ProgressBar>,
    last_progress_pos: &mut u64,
) {
    if let Some(pb) = pb {
        if pos > *last_progress_pos {
            pb.inc(pos - *last_progress_pos);
            *last_progress_pos = pos;
        }
    }
}

/// Finishes any remaining byte-progress after the stream reaches EOF.
fn finish_segment_progress<R: Seek>(
    reader: &mut R,
    pb: Option<&ProgressBar>,
    last_progress_pos: &mut u64,
) -> Result<()> {
    update_segment_progress(reader, pb, last_progress_pos)
}

/// Iterates queryname-sorted SAM groups.
fn scan_sam_groups<F>(path: &str, on_group: &mut F, pb: Option<&ProgressBar>) -> Result<()>
where
    F: FnMut(&str, &[AsAlignment]) -> Result<bool>,
{
    let file = File::open(path).with_context(|| format!("open SAM {}", path))?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    let mut current_id = String::new();
    let mut group = Vec::with_capacity(8);
    let mut records_since_progress = 0usize;
    let mut last_progress_pos = reader.stream_position()?;
    while reader.read_line(&mut line)? != 0 {
        if line.is_empty() || line.starts_with('@') {
            line.clear();
            continue;
        }
        let record_line = line.trim_end_matches(&['\r', '\n'][..]);
        let mut cols = record_line.split('\t');
        let read_id = cols.next().unwrap_or_default();
        if !current_id.is_empty() && read_id != current_id {
            if on_group(&current_id, &group)? {
                return Ok(());
            }
            group.clear();
        }
        records_since_progress += 1;
        if records_since_progress >= 4096 {
            update_segment_progress(&mut reader, pb, &mut last_progress_pos)?;
            records_since_progress = 0;
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
        line.clear();
    }
    if !group.is_empty() {
        on_group(&current_id, &group)?;
    }
    finish_segment_progress(&mut reader, pb, &mut last_progress_pos)?;
    Ok(())
}

/// Iterates queryname-sorted BAM groups through noodles.
fn scan_bam_groups<F>(path: &str, on_group: &mut F, pb: Option<&ProgressBar>) -> Result<()>
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
    let mut records_since_progress = 0usize;
    let mut last_progress_pos = reader.get_mut().position();
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
        records_since_progress += 1;
        if records_since_progress >= 4096 {
            update_segment_progress_position(
                reader.get_mut().position(),
                pb,
                &mut last_progress_pos,
            );
            records_since_progress = 0;
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
    update_segment_progress_position(reader.get_mut().position(), pb, &mut last_progress_pos);
    Ok(())
}

impl<'a> SegmentScanContext<'a> {
    /// Borrows the read-only indexes needed by the non-BSJ segments rescan.
    fn from_state(state: &'a ScanState<'a>) -> Self {
        Self {
            junction_read_to_circ: &state.junction_read_to_circ,
            reference: state.reference,
            annotation: state.annotation,
            circ_spans_by_chr: &state.circ_spans_by_chr,
            clusters_by_chr: &state.clusters_by_chr,
            read_len: state.read_len,
            min_mapq: state.min_mapq,
        }
    }
}

/// Merges one completed segments scan accumulator into the main state.
fn merge_segment_scan_accum(state: &mut ScanState, accum: SegmentScanAccum) {
    merge_segment_scan_accum_fields(
        &mut state.candidates,
        &mut state.segment_groups,
        &mut state.outward_read_ids,
        &mut state.prebuilt_segment_records,
        accum,
    );
}

/// Merges one accumulator into disjoint mutable state fields.
///
/// This path is still used by the legacy BAM rescan fallback, where workers
/// return only retained non-BSJ groups. The Scan2 sidecar path now prefers
/// shard-local spill files so large whole-genome inputs do not keep alignment
/// groups in memory before the support-aware second pass knows which reads are
/// actually needed.
fn merge_segment_scan_accum_fields(
    candidates: &mut Vec<PositiveCandidate>,
    segment_groups: &mut HashMap<String, Vec<AsAlignment>>,
    outward_read_ids: &mut HashSet<String>,
    prebuilt_segment_records: &mut Vec<SegmentRecord>,
    accum: SegmentScanAccum,
) {
    candidates.extend(accum.candidates);
    outward_read_ids.extend(accum.outward_read_ids);
    prebuilt_segment_records.extend(accum.prebuilt_segment_records);
    for (read_id, records) in accum.segment_groups {
        segment_groups.entry(read_id).or_insert(records);
    }
}

impl SegmentScanShardPaths {
    /// Builds deterministic shard spill paths beside the final `.segments` file.
    fn new(out_prefix: &str, shard_idx: usize) -> Self {
        Self {
            path: part_path(&format!("{out_prefix}.segments"), shard_idx),
        }
    }

    /// Builds shard spill paths for Summary-confirmed BSJ read evidence.
    fn new_bsj(out_prefix: &str, shard_idx: usize) -> Self {
        Self {
            path: part_path(&format!("{out_prefix}.segments.bsj"), shard_idx),
        }
    }

    /// Builds raw hash-partition paths used while aggregating BSJ sidecars.
    fn new_bsj_raw(out_prefix: &str, shard_idx: usize) -> Self {
        Self {
            path: part_path(&format!("{out_prefix}.segments.bsj.raw"), shard_idx),
        }
    }
}

impl SegmentScanShardWriter {
    /// Opens the shard-local read-result stream for segments rescan evidence.
    fn new(paths: &SegmentScanShardPaths) -> Result<Self> {
        Ok(Self {
            writer: BufWriter::with_capacity(512 * 1024, File::create(&paths.path)?),
        })
    }

    /// Writes all retained results for one non-BSJ read group.
    ///
    /// The tagged stream keeps evidence read-local while still allowing a cheap
    /// two-pass merge: `R` and `C` rows are loaded before motif validation, and
    /// `A` rows are loaded afterward only for read IDs that survived validation
    /// or were classified as outward.
    fn write_read_result(
        &mut self,
        read_id: &str,
        records: &[AsAlignment],
        outward: bool,
        candidates: &[PositiveCandidate],
        segment_records: &[SegmentRecord],
    ) -> Result<()> {
        writeln!(self.writer, "R\t{read_id}\t{}", u8::from(outward))?;
        for record in records {
            writeln!(
                self.writer,
                "A\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                read_id,
                record.flag,
                record.chr,
                record.pos,
                record.mapq,
                record.cigar,
                if record.seq.is_empty() {
                    "*"
                } else {
                    &record.seq
                },
                u8::from(record.from_local_clip),
                format_xa_alternatives(&record.xa_alternatives)
            )?;
        }
        for candidate in candidates {
            writeln!(
                self.writer,
                "C\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                candidate.read_id,
                candidate.chr,
                candidate.site2,
                candidate.site1,
                candidate.adjust1,
                candidate.adjust2,
                candidate.strand_hint,
                option_field(candidate.cigars[0].as_deref()),
                option_field(candidate.cigars[1].as_deref()),
                option_field(candidate.cigars[2].as_deref())
            )?;
        }
        for record in segment_records {
            writeln!(
                self.writer,
                "S\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
                record.r1_align_strand,
                record.r2_align_strand,
                record.r1_cigar,
                record.r1_segments,
                record.r2_cigar,
                record.r2_segments
            )?;
        }
        Ok(())
    }

    /// Flushes all shard-local spill writers before the main thread merges them.
    fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}

/// Formats an optional TSV field without introducing empty-column ambiguity.
fn option_field(value: Option<&str>) -> &str {
    value.filter(|value| !value.is_empty()).unwrap_or("*")
}

/// Serializes XA alternatives into a shard-local sidecar field.
fn format_xa_alternatives(alternatives: &[XaAlternative]) -> String {
    let mut out = String::new();
    for alternative in alternatives {
        let _ = write!(
            out,
            "{},{},{},{},{};",
            alternative.chr,
            alternative.strand,
            alternative.pos,
            alternative.cigar,
            alternative.edit_distance
        );
    }
    out
}

/// Parses XA alternatives previously written by a segments rescan shard.
fn parse_xa_alternatives_field(raw: &str) -> Vec<XaAlternative> {
    raw.split(';')
        .filter_map(|entry| {
            if entry.is_empty() {
                return None;
            }
            let mut parts = entry.split(',');
            Some(XaAlternative {
                chr: parts.next()?.to_string(),
                strand: parts.next()?.chars().next()?,
                pos: parts.next()?.parse().ok()?,
                cigar: parts.next()?.to_string(),
                edit_distance: parts.next()?.parse().ok()?,
            })
        })
        .collect()
}

/// Streams retained read blocks from one segment shard.
///
/// This is the preferred large-data access pattern for finalize. It keeps only
/// one read group in memory, preserving the read-local `A/C/S` payload needed
/// for support-aware correction without building a global read-id offset map.
fn for_each_retained_segment_group<F>(
    path: &str,
    pb: Option<&ProgressBar>,
    mut on_group: F,
) -> Result<()>
where
    F: FnMut(RetainedSegmentGroup) -> Result<()>,
{
    let file = File::open(path).with_context(|| format!("open segments retained shard {path}"))?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    let mut current: Option<RetainedSegmentGroup> = None;
    let mut pending_progress = 0usize;
    loop {
        line.clear();
        let bytes = reader.read_line(&mut line)?;
        if bytes == 0 {
            break;
        }
        pending_progress += bytes;
        let line = line.trim_end_matches(['\n', '\r']);
        if line.is_empty() {
            if pending_progress >= 4 * 1024 * 1024 {
                if let Some(pb) = pb {
                    pb.inc(pending_progress as u64);
                }
                pending_progress = 0;
            }
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        match parts.first().copied() {
            Some("R") if parts.len() >= 2 => {
                if let Some(group) = current.take() {
                    on_group(group)?;
                }
                current = Some(RetainedSegmentGroup {
                    read_id: parts[1].to_string(),
                    records: Vec::new(),
                    candidates: Vec::new(),
                    segment_records: Vec::new(),
                });
            }
            Some("A") if parts.len() >= 10 => {
                if let Some(group) = current.as_mut() {
                    group.records.push(AsAlignment {
                        flag: parts[2].parse().unwrap_or(0),
                        chr: parts[3].to_string(),
                        pos: parts[4].parse().unwrap_or(0),
                        mapq: parts[5].parse().unwrap_or(0),
                        cigar: parts[6].to_string(),
                        seq: if parts[7] == "*" {
                            String::new()
                        } else {
                            parts[7].to_string()
                        },
                        from_local_clip: parts[8] == "1",
                        xa_alternatives: parse_xa_alternatives_field(parts[9]),
                    });
                }
            }
            Some("C") if parts.len() >= 11 => {
                if let Some(group) = current.as_mut() {
                    group.candidates.push(PositiveCandidate {
                        read_id: parts[1].to_string(),
                        index: 0,
                        chr: parts[2].to_string(),
                        site2: parts[3].parse().unwrap_or(0),
                        site1: parts[4].parse().unwrap_or(0),
                        adjust1: parts[5].parse().unwrap_or(0),
                        adjust2: parts[6].parse().unwrap_or(0),
                        strand_hint: parts[7].parse().unwrap_or(0),
                        cigars: [
                            parse_optional_field(parts[8]),
                            parse_optional_field(parts[9]),
                            parse_optional_field(parts[10]),
                        ],
                    });
                }
            }
            Some("S") if parts.len() >= 15 => {
                if let (Some(group), Some(record)) = (
                    current.as_mut(),
                    segment_record_from_shard_fields(&parts[1..]),
                ) {
                    group.segment_records.push(record);
                }
            }
            _ => {}
        }
        if pending_progress >= 4 * 1024 * 1024 {
            if let Some(pb) = pb {
                pb.inc(pending_progress as u64);
            }
            pending_progress = 0;
        }
    }
    if let Some(group) = current {
        on_group(group)?;
    }
    if pending_progress > 0 {
        if let Some(pb) = pb {
            pb.inc(pending_progress as u64);
        }
    }
    Ok(())
}

/// Parses one shard-local preliminary segment row.
fn segment_record_from_shard_fields(parts: &[&str]) -> Option<SegmentRecord> {
    let has_alignment_strands = parts.len() >= 16;
    let (
        r1_align_strand,
        r2_align_strand,
        r1_cigar_idx,
        r1_segments_idx,
        r2_cigar_idx,
        r2_segments_idx,
    ) = if has_alignment_strands {
        (
            parts.get(10)?.to_string(),
            parts.get(11)?.to_string(),
            12usize,
            13usize,
            14usize,
            15usize,
        )
    } else {
        (
            "NA".to_string(),
            "NA".to_string(),
            10usize,
            11usize,
            12usize,
            13usize,
        )
    };
    Some(SegmentRecord {
        read_id: parts.first()?.to_string(),
        type_name: match *parts.get(1)? {
            "backward" => "backward",
            "outward" => "outward",
            _ => return None,
        },
        circ_id: parts.get(2)?.to_string(),
        chrom: parts.get(3)?.to_string(),
        start: parts.get(4)?.to_string(),
        end: parts.get(5)?.to_string(),
        strand: parts.get(6)?.to_string(),
        is_circular: parts.get(7)?.parse().ok()?,
        is_r1_bsj: parts.get(8)?.parse().ok()?,
        is_r2_bsj: parts.get(9)?.parse().ok()?,
        r1_align_strand,
        r2_align_strand,
        r1_cigar: parts.get(r1_cigar_idx)?.to_string(),
        r1_segments: parts.get(r1_segments_idx)?.to_string(),
        r2_cigar: parts.get(r2_cigar_idx)?.to_string(),
        r2_segments: parts.get(r2_segments_idx)?.to_string(),
    })
}

/// Parses a shard-local optional text field.
fn parse_optional_field(value: &str) -> Option<String> {
    (value != "*").then(|| value.to_string())
}

/// Processes the segments BAM rescan in BGZF shards.
///
/// This mirrors the Scan1/Scan2 ownership rule: a non-zero shard starts
/// decoding from the previous BGZF block, skips the leading partial read group,
/// and then owns each subsequent queryname group until it crosses the shard end.
fn scan_backward_bam_groups_parallel(
    path: &str,
    out_prefix: &str,
    ctx: &SegmentScanContext<'_>,
) -> Result<Vec<SegmentScanShardPaths>> {
    let file = File::open(path).with_context(|| format!("open BAM {}", path))?;
    let file_size = std::fs::metadata(path)?.len();
    let mmap = unsafe { Mmap::map(&file)? };
    let num_threads = bam_shard_count(mmap.len(), rayon::current_num_threads());
    let shard_size = mmap.len() / num_threads;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header()?;
    let pb = segment_scan_progress_bar(path)?;

    unsafe {
        libc::madvise(
            mmap.as_ptr() as *mut libc::c_void,
            mmap.len(),
            libc::MADV_SEQUENTIAL,
        );
    }

    let header_ref = &header;
    let shard_paths: Vec<SegmentScanShardPaths> = (0..num_threads)
        .into_par_iter()
        .map(|i| {
            let start = i * shard_size;
            let end = if i == num_threads - 1 {
                mmap.len()
            } else {
                (i + 1) * shard_size
            };
            let paths = SegmentScanShardPaths::new(out_prefix, i);
            process_backward_bam_shard(i, &mmap, start, end, header_ref, ctx, &pb, &paths)?;
            Ok(paths)
        })
        .collect::<Result<Vec<_>>>()?;

    pb.set_position(file_size);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {msg}",
            )?
            .progress_chars("#>-"),
    );
    pb.finish_with_message("Completed");

    Ok(shard_paths)
}

/// Processes one compressed-byte shard for the segments BAM rescan.
fn process_backward_bam_shard(
    shard_idx: usize,
    mmap: &Mmap,
    start: usize,
    end: usize,
    header: &sam::Header,
    ctx: &SegmentScanContext<'_>,
    pb: &ProgressBar,
    paths: &SegmentScanShardPaths,
) -> Result<()> {
    let block_start = if start == 0 {
        0
    } else {
        find_next_bgzf_block(mmap, start).unwrap_or(mmap.len())
    };
    let pos = if start == 0 {
        0
    } else {
        find_previous_bgzf_block(mmap, block_start).unwrap_or(0)
    };
    if pos >= mmap.len() {
        let mut writer = SegmentScanShardWriter::new(paths)?;
        writer.flush()?;
        return Ok(());
    }

    let mut reader = bam::io::Reader::new(&mmap[pos..]);
    if start == 0 {
        let _ = reader.read_header()?;
    }
    let mut record = bam::Record::default();
    let mut current_id: Vec<u8> = Vec::new();
    let mut group = Vec::with_capacity(8);
    let mut crossed_start = start == 0;
    let mut leading_partial_id: Option<Vec<u8>> = None;
    let mut last_progress_pos = block_start.saturating_sub(pos);
    let mut last_evicted_pos = pos;
    let eviction_threshold = 64 * 1024 * 1024;
    let mut cigar_buf = String::with_capacity(64);
    let mut seq_buf = String::with_capacity(256);
    let mut writer = SegmentScanShardWriter::new(paths)?;

    while reader.read_record(&mut record)? != 0 {
        let curr_c_pos = reader.get_ref().virtual_position().compressed() as usize;
        if curr_c_pos > last_progress_pos {
            pb.inc((curr_c_pos - last_progress_pos) as u64);
            last_progress_pos = curr_c_pos;
        }
        if curr_c_pos > last_evicted_pos.saturating_sub(pos) + eviction_threshold {
            let evicted_rel = last_evicted_pos - pos;
            advise_segment_mmap_dontneed(mmap, last_evicted_pos, curr_c_pos - evicted_rel);
            last_evicted_pos = pos + curr_c_pos;
        }

        let read_id = record
            .name()
            .ok_or_else(|| anyhow!("Missing read name in segments BAM shard {shard_idx}"))?;
        let abs_c_pos = pos + curr_c_pos;
        if !crossed_start {
            if abs_c_pos < block_start {
                leading_partial_id = Some(read_id.to_vec());
                continue;
            }
            crossed_start = true;
        }
        if let Some(partial_id) = &leading_partial_id {
            if read_id == partial_id.as_slice() {
                continue;
            }
            leading_partial_id = None;
        }
        if read_id != current_id.as_slice() {
            if !current_id.is_empty() {
                let id = String::from_utf8_lossy(&current_id);
                process_backward_group_to_writer(&id, &group, ctx, &mut writer)?;
                if pos + curr_c_pos > end {
                    current_id.clear();
                    break;
                }
                group.clear();
            }
            current_id.clear();
            current_id.extend_from_slice(read_id);
        }
        group.push(bam_record_to_as_alignment(
            &record,
            header,
            &mut cigar_buf,
            &mut seq_buf,
        )?);
    }
    if !current_id.is_empty() {
        let id = String::from_utf8_lossy(&current_id);
        process_backward_group_to_writer(&id, &group, ctx, &mut writer)?;
    }
    writer.flush()?;
    let processed_end = (pos + last_progress_pos).min(mmap.len());
    advise_segment_mmap_dontneed(
        mmap,
        last_evicted_pos,
        processed_end.saturating_sub(last_evicted_pos),
    );
    Ok(())
}

/// Returns the first BGZF block header at or after `start`.
fn find_next_bgzf_block(mmap: &Mmap, start: usize) -> Option<usize> {
    (start..mmap.len().saturating_sub(3)).find(|&i| &mmap[i..i + 4] == b"\x1f\x8b\x08\x04")
}

/// Returns the nearest BGZF block header at or before `start`.
fn find_previous_bgzf_block(mmap: &Mmap, start: usize) -> Option<usize> {
    let mut i = start.saturating_sub(1).min(mmap.len().saturating_sub(4));
    loop {
        if i + 3 < mmap.len() && &mmap[i..i + 4] == b"\x1f\x8b\x08\x04" {
            return Some(i);
        }
        if i == 0 {
            return None;
        }
        i -= 1;
    }
}

/// Advises the kernel that an already-processed segments shard range is cold.
fn advise_segment_mmap_dontneed(mmap: &Mmap, offset: usize, len: usize) {
    if len == 0 {
        return;
    }
    let page_size = 4096usize;
    let aligned_offset = (offset / page_size) * page_size;
    let aligned_len = ((offset + len + page_size - 1) / page_size) * page_size - aligned_offset;
    unsafe {
        let ptr = mmap.as_ptr().add(aligned_offset);
        libc::madvise(ptr as *mut libc::c_void, aligned_len, libc::MADV_DONTNEED);
    }
}

/// Converts one BAM record into the lightweight segments alignment form.
fn bam_record_to_as_alignment(
    record: &bam::Record,
    header: &sam::Header,
    cigar_buf: &mut String,
    seq_buf: &mut String,
) -> Result<AsAlignment> {
    use noodles::sam::alignment::Record as _;

    let chr = match record.reference_sequence(header) {
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
        let _ = write!(cigar_buf, "{}{}", op.len(), op_char);
    }
    seq_buf.clear();
    for b in record.sequence().iter() {
        seq_buf.push(char::from(b));
    }
    Ok(AsAlignment {
        flag,
        chr,
        pos,
        mapq,
        cigar: cigar_buf.clone(),
        seq: seq_buf.clone(),
        from_local_clip: false,
        xa_alternatives: parse_xa_from_bam_record(record)?,
    })
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
    if state.junction_read_to_circ.contains_key(read_id) {
        state
            .segment_groups
            .insert(read_id.to_string(), records.to_vec());
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
/// show read-level backward or outward topology. These supplemental rows do not
/// require a known circRNA span or a matching BSJ; downstream graph construction
/// can later project them onto compatible circRNA hypotheses.
fn process_backward_group(
    read_id: &str,
    records: &[AsAlignment],
    ctx: &SegmentScanContext<'_>,
    accum: &mut SegmentScanAccum,
) -> Result<()> {
    if records.is_empty() || ctx.junction_read_to_circ.contains_key(read_id) {
        return Ok(());
    }
    let is_outward = is_outward_pair_group_ctx(records, ctx);
    let may_backward = may_have_backward_candidate_shape(records);
    if !is_outward && !may_backward {
        return Ok(());
    }

    let candidates = if may_backward {
        mapping_check1_backward_candidates(read_id, records, ctx)
    } else {
        Vec::new()
    };
    let mut candidates = validate_splice_candidates(candidates, ctx.reference, ctx.annotation);
    let (segment_records, enriched_records) =
        prebuild_non_bsj_segment_records(read_id, records, is_outward, &candidates, ctx);
    if is_outward {
        accum
            .segment_groups
            .entry(read_id.to_string())
            .or_insert_with(|| enriched_records.unwrap_or_else(|| records.to_vec()));
        accum.outward_read_ids.insert(read_id.to_string());
    }
    if !candidates.is_empty() {
        accum
            .segment_groups
            .entry(read_id.to_string())
            .or_insert_with(|| records.to_vec());
        accum.candidates.append(&mut candidates);
    }
    accum.prebuilt_segment_records.extend(segment_records);
    Ok(())
}

/// Streams one non-BSJ read group into shard-local spill files.
fn process_backward_group_to_writer(
    read_id: &str,
    records: &[AsAlignment],
    ctx: &SegmentScanContext<'_>,
    writer: &mut SegmentScanShardWriter,
) -> Result<()> {
    if records.is_empty() || ctx.junction_read_to_circ.contains_key(read_id) {
        return Ok(());
    }
    let is_outward = is_outward_pair_group_ctx(records, ctx);
    let may_backward = may_have_backward_candidate_shape(records);
    if !is_outward && !may_backward {
        return Ok(());
    }

    let candidates = if may_backward {
        mapping_check1_backward_candidates(read_id, records, ctx)
    } else {
        Vec::new()
    };
    let candidates = validate_splice_candidates(candidates, ctx.reference, ctx.annotation);
    let (segment_records, enriched_records) =
        prebuild_non_bsj_segment_records(read_id, records, is_outward, &candidates, ctx);
    if is_outward {
        let retained_records = enriched_records.as_deref().unwrap_or(records);
        writer.write_read_result(
            read_id,
            retained_records,
            true,
            &candidates,
            &segment_records,
        )?;
    } else if !candidates.is_empty() {
        writer.write_read_result(read_id, records, false, &candidates, &segment_records)?;
    }
    Ok(())
}

/// Cheap read-level prefilter for CIRI-AS-style backward split evidence.
///
/// This intentionally does not consult circRNA spans. It only asks whether the
/// group has at least two same-strand same-chromosome alignments with a clipped
/// CIGAR shape that could later produce a backward junction candidate.
fn may_have_backward_candidate_shape(records: &[AsAlignment]) -> bool {
    let mut buckets: HashMap<(&str, i32), (usize, bool)> = HashMap::new();
    for record in records {
        if record.chr == "*" || record.cigar == "*" {
            continue;
        }
        let key = (record.chr.as_str(), reverse_bit(record.flag));
        let entry = buckets.entry(key).or_insert((0, false));
        entry.0 += 1;
        entry.1 |= record.cigar.contains('S') || record.cigar.contains('H');
        if entry.0 >= 2 && entry.1 {
            return true;
        }
    }
    false
}

/// Validates one read group's splice candidates before shard-local materialization.
fn validate_splice_candidates(
    candidates: Vec<PositiveCandidate>,
    reference: &HashMap<String, String>,
    annotation: Option<&Annotation>,
) -> Vec<PositiveCandidate> {
    candidates
        .into_iter()
        .enumerate()
        .filter_map(|(idx, mut candidate)| {
            candidate.index = idx;
            validate_splice_candidate(candidate, reference, annotation)
        })
        .collect()
}

/// Builds preliminary non-BSJ segment rows inside the owning BAM shard.
///
/// This moves the read-local backward/outward first pass into the parallel scan
/// while leaving global junction-support correction for the later build phase.
fn prebuild_non_bsj_segment_records(
    read_id: &str,
    records: &[AsAlignment],
    is_outward: bool,
    candidates: &[PositiveCandidate],
    ctx: &SegmentScanContext<'_>,
) -> (Vec<SegmentRecord>, Option<Vec<AsAlignment>>) {
    let correction = SegmentCorrectionContext {
        reference: ctx.reference,
        annotation: ctx.annotation,
        junction_support: None,
    };
    let junction_hints = local_junction_hints(candidates);
    let strand_hint = local_strand_hint(candidates);
    let enriched_records = is_outward
        .then(|| try_add_non_bsj_local_clip_alignments_ctx(records, ctx))
        .flatten();
    let materialization_records = enriched_records.as_deref().unwrap_or(records);

    if !candidates.is_empty() {
        if let Some(record) = build_backward_segment_record_from_group_ctx(
            read_id,
            materialization_records,
            ctx,
            &correction,
            strand_hint,
            &junction_hints,
        ) {
            return (vec![record], enriched_records);
        }
    }
    if is_outward {
        if let Some(record) = build_outward_segment_record(
            read_id,
            materialization_records,
            Some(&correction),
            ctx.read_len,
            ctx.min_mapq,
        ) {
            return (vec![record], enriched_records);
        }
    }
    (Vec::new(), enriched_records)
}

/// Builds one backward row using only shard-local read context.
fn build_backward_segment_record_from_group_ctx(
    read_id: &str,
    records: &[AsAlignment],
    ctx: &SegmentScanContext<'_>,
    correction: &SegmentCorrectionContext<'_>,
    strand_hint: Option<char>,
    junction_hints: &[(i32, i32)],
) -> Option<SegmentRecord> {
    let base = build_backward_segment_record(
        read_id,
        records,
        strand_hint,
        junction_hints,
        Some(correction),
        ctx.read_len,
    )?;
    let Some(enriched_records) = try_add_non_bsj_local_clip_alignments_ctx(records, ctx) else {
        return Some(base);
    };
    build_backward_segment_record(
        read_id,
        &enriched_records,
        strand_hint,
        junction_hints,
        Some(correction),
        ctx.read_len,
    )
    .or(Some(base))
}

/// Returns a unique RNA-strand hint from one read's validated candidates.
fn local_strand_hint(candidates: &[PositiveCandidate]) -> Option<char> {
    let mut hint = 0;
    for candidate in candidates {
        if candidate.strand_hint == 0 {
            continue;
        }
        if hint != 0 && hint != candidate.strand_hint {
            return None;
        }
        hint = candidate.strand_hint;
    }
    match hint {
        -1 => Some('+'),
        1 => Some('-'),
        _ => None,
    }
}

/// Returns corrected internal junction hints for one read.
fn local_junction_hints(candidates: &[PositiveCandidate]) -> Vec<(i32, i32)> {
    candidates
        .iter()
        .map(|candidate| (candidate.site2, candidate.site1))
        .collect()
}

/// Adds local soft-clip pseudo-alignments for non-BSJ read groups.
///
/// This mirrors the BSJ sidecar idea without assuming a read-level circRNA
/// assignment. Candidate placements are searched in local windows around
/// detected circ loci and the read's mapped blocks; the final backward/outward
/// chain builder still decides whether a placement improves the selected
/// topology. The original mapper rows are always retained first so these local
/// rows cannot erase the primary evidence.
fn add_non_bsj_local_clip_alignments(
    records: &[AsAlignment],
    state: &ScanState,
) -> Vec<AsAlignment> {
    add_non_bsj_local_clip_alignments_ctx(records, &SegmentScanContext::from_state(state))
}

/// Adds local soft-clip pseudo-alignments using the shared scan context.
fn add_non_bsj_local_clip_alignments_ctx(
    records: &[AsAlignment],
    ctx: &SegmentScanContext<'_>,
) -> Vec<AsAlignment> {
    try_add_non_bsj_local_clip_alignments_ctx(records, ctx).unwrap_or_else(|| records.to_vec())
}

/// Adds local soft-clip pseudo-alignments only when the rescue creates rows.
///
/// Most non-BSJ groups have no usable soft clip placement. Returning `None` for
/// those cases avoids repeatedly cloning the read group while preserving the old
/// materialized output through `add_non_bsj_local_clip_alignments_ctx`.
fn try_add_non_bsj_local_clip_alignments_ctx(
    records: &[AsAlignment],
    ctx: &SegmentScanContext<'_>,
) -> Option<Vec<AsAlignment>> {
    if !records.iter().any(|record| {
        !record.from_local_clip
            && record.cigar.contains('S')
            && !record.seq.is_empty()
            && record.seq != "*"
    }) {
        return None;
    }
    let windows = non_bsj_local_clip_windows_ctx(records, ctx);
    if windows.is_empty() {
        return None;
    }
    let mut local = non_bsj_local_clip_alignments(records, &windows, ctx);
    if local.is_empty() {
        return None;
    }
    local.sort_by(|a, b| {
        mate_bucket(a.flag)
            .cmp(&mate_bucket(b.flag))
            .then_with(|| a.chr.cmp(&b.chr))
            .then_with(|| a.pos.cmp(&b.pos))
            .then_with(|| a.cigar.cmp(&b.cigar))
    });
    local.truncate(NON_BSJ_LOCAL_CLIP_MAX_ROWS_PER_READ);
    let mut out = records.to_vec();
    out.extend(local);
    Some(out)
}

/// Builds bounded local reference windows for non-BSJ clip placement.
///
/// The search space deliberately includes circ-cluster and read-neighborhood
/// anchors, not only exact Summary circ spans. Windows are centered on likely
/// junction anchors rather than whole loci, which keeps this read-local rescue
/// bounded while still allowing clips to place just outside a particular circ
/// span when that is the most coherent chain explanation.
/// Builds bounded local reference windows from the shared scan context.
fn non_bsj_local_clip_windows_ctx(
    records: &[AsAlignment],
    ctx: &SegmentScanContext<'_>,
) -> Vec<LocalClipWindow> {
    let mut by_chr: HashMap<&str, (i32, i32)> = HashMap::new();
    let mut anchors_by_chr: HashMap<&str, Vec<i32>> = HashMap::new();
    let anchor_flank = ctx.read_len * NON_BSJ_LOCAL_CLIP_ANCHOR_FLANK_MULTIPLIER;
    for record in records {
        if record.chr == "*" {
            continue;
        }
        let Some((start, end)) = alignment_ref_span(record, ctx.read_len) else {
            continue;
        };
        by_chr
            .entry(record.chr.as_str())
            .and_modify(|span| {
                span.0 = span.0.min(start);
                span.1 = span.1.max(end);
            })
            .or_insert((start, end));
        anchors_by_chr
            .entry(record.chr.as_str())
            .or_default()
            .extend([start, end]);
    }

    let mut windows = Vec::new();
    let mut seen = HashSet::new();
    for (chr, (group_start, group_end)) in by_chr {
        if let Some(spans) = ctx.circ_spans_by_chr.get(chr) {
            for span in
                overlapping_circ_spans(spans, group_start - MIN_INTRON, group_end + MIN_INTRON)
            {
                push_local_clip_anchor_window(
                    &mut windows,
                    &mut seen,
                    ctx.reference,
                    chr,
                    span.start,
                    anchor_flank,
                );
                push_local_clip_anchor_window(
                    &mut windows,
                    &mut seen,
                    ctx.reference,
                    chr,
                    span.end,
                    anchor_flank,
                );
            }
        }
        if let Some(clusters) = ctx.clusters_by_chr.get(chr) {
            for cluster in overlapping_circ_clusters(
                clusters,
                group_start - MIN_INTRON,
                group_end + MIN_INTRON,
            ) {
                push_local_clip_anchor_window(
                    &mut windows,
                    &mut seen,
                    ctx.reference,
                    chr,
                    cluster.start,
                    anchor_flank,
                );
                push_local_clip_anchor_window(
                    &mut windows,
                    &mut seen,
                    ctx.reference,
                    chr,
                    cluster.end,
                    anchor_flank,
                );
            }
        }
        push_local_clip_anchor_window(
            &mut windows,
            &mut seen,
            ctx.reference,
            chr,
            group_start,
            anchor_flank,
        );
        push_local_clip_anchor_window(
            &mut windows,
            &mut seen,
            ctx.reference,
            chr,
            group_end,
            anchor_flank,
        );
        if let Some(anchors) = anchors_by_chr.get(chr) {
            for anchor in anchors {
                push_local_clip_anchor_window(
                    &mut windows,
                    &mut seen,
                    ctx.reference,
                    chr,
                    *anchor,
                    anchor_flank,
                );
            }
        }
    }
    windows
}

/// Returns circ spans overlapping a query interval using the sorted span index.
fn overlapping_circ_spans(
    spans: &[CircSpan],
    query_start: i32,
    query_end: i32,
) -> impl Iterator<Item = &CircSpan> {
    let right = spans.partition_point(|span| span.start <= query_end);
    spans[..right]
        .iter()
        .rev()
        .take_while(move |span| span.max_end_through >= query_start)
        .filter(move |span| span.end >= query_start)
}

/// Returns circ clusters overlapping a query interval using sorted disjoint clusters.
fn overlapping_circ_clusters(
    clusters: &[CircCluster],
    query_start: i32,
    query_end: i32,
) -> impl Iterator<Item = &CircCluster> {
    let right = clusters.partition_point(|cluster| cluster.start <= query_end);
    clusters[..right]
        .iter()
        .rev()
        .take_while(move |cluster| cluster.end >= query_start)
}

/// Appends one clamped anchor-centered local clip window if it is non-duplicate.
fn push_local_clip_anchor_window(
    windows: &mut Vec<LocalClipWindow>,
    seen: &mut HashSet<(String, i32, i32)>,
    reference: &HashMap<String, String>,
    chr: &str,
    anchor: i32,
    flank: i32,
) {
    let Some(chr_seq) = reference.get(chr) else {
        return;
    };
    let chr_len = chr_seq.len() as i32;
    let start = (anchor - flank).max(1);
    let end = (anchor + flank).min(chr_len);
    if start > end {
        return;
    }
    let key = (chr.to_string(), start, end);
    if seen.insert(key.clone()) {
        windows.push(LocalClipWindow {
            chr: key.0,
            start,
            end,
        });
    }
}

/// Generates pseudo-alignment rows by full-exact soft-clip placement.
///
/// BSJ sidecars can afford longest-partial local search because they run only
/// after a BSJ validator accepts one read group. Backward/outward candidates are
/// more numerous, so this stage keeps only full-clip exact placements and leaves
/// partial/approximate clip rescue for a later indexed or circRNA-level pass.
fn non_bsj_local_clip_alignments(
    records: &[AsAlignment],
    windows: &[LocalClipWindow],
    ctx: &SegmentScanContext<'_>,
) -> Vec<AsAlignment> {
    let mut out = Vec::new();
    let mut seen = HashSet::new();
    let min_clip_len = MIN_JUNCTION_SUPPORT_SEGMENT_LENGTH as usize;
    for record in records {
        if record.from_local_clip
            || !record.cigar.contains('S')
            || record.seq.is_empty()
            || record.seq == "*"
            || record.chr == "*"
        {
            continue;
        }
        let payload = segment_clip_payload(record);
        if payload == "*" {
            continue;
        }
        for (side, clip_seq) in parse_clip_payload(&payload) {
            if clip_seq.len() < min_clip_len || clip_seq.contains('N') {
                continue;
            }
            let query = clip_seq.to_ascii_uppercase();
            if query.contains('N') {
                continue;
            }
            let mut candidates: Vec<(i32, AsAlignment)> = Vec::new();
            for window in windows.iter().filter(|window| window.chr == record.chr) {
                let Some(chr_seq) = ctx.reference.get(&window.chr) else {
                    continue;
                };
                let Some(local_seq) = chr_seq.get((window.start - 1) as usize..window.end as usize)
                else {
                    continue;
                };
                let mut positions =
                    exact_clip_match_positions(local_seq, &query, window.start, window.end);
                positions.truncate(NON_BSJ_LOCAL_CLIP_MAX_PLACEMENTS_PER_SIDE);
                for pos in positions {
                    let read_len = record.seq.len() as i32;
                    let Some(cigar) =
                        clip_placement_cigar(side, clip_seq.len(), 0, clip_seq.len(), read_len)
                    else {
                        continue;
                    };
                    let flag = record.flag | 0x800;
                    let key = (flag, window.chr.clone(), pos, cigar.clone(), side);
                    if !seen.insert(key) {
                        continue;
                    }
                    let distance =
                        local_clip_distance_to_group(records, &window.chr, pos, ctx.read_len);
                    candidates.push((
                        -((clip_seq.len() as i32) * 1000) + distance,
                        AsAlignment {
                            flag,
                            chr: window.chr.clone(),
                            pos,
                            mapq: record.mapq.saturating_sub(1),
                            cigar,
                            seq: String::new(),
                            from_local_clip: true,
                            xa_alternatives: Vec::new(),
                        },
                    ));
                }
            }
            candidates.sort_by(|a, b| a.0.cmp(&b.0));
            out.extend(
                candidates
                    .into_iter()
                    .take(NON_BSJ_LOCAL_CLIP_MAX_PLACEMENTS_PER_SIDE)
                    .map(|(_, alignment)| alignment),
            );
        }
    }
    out
}

/// Scores how close a local clip placement is to the read group's mapped blocks.
fn local_clip_distance_to_group(
    records: &[AsAlignment],
    chr: &str,
    pos: i32,
    read_len: i32,
) -> i32 {
    records
        .iter()
        .filter(|record| record.chr == chr)
        .filter_map(|record| alignment_ref_span(record, read_len))
        .map(|(start, end)| (pos - start).abs().min((pos - end).abs()))
        .min()
        .unwrap_or(i32::MAX / 4)
}

/// Tests membership in Perl's materialized `circ_cluster_range` hash.
fn point_in_cluster_range(pos: i32, cluster: &CircCluster) -> bool {
    pos >= cluster.start - MIN_INTRON && pos <= cluster.end + MIN_INTRON
}

/// True if any alignment start/end falls inside a clustered circRNA locus.
fn overlaps_any_circ_cluster(records: &[AsAlignment], state: &ScanState) -> bool {
    overlaps_any_circ_cluster_ctx(records, &SegmentScanContext::from_state(state))
}

/// True if any alignment start/end falls inside a clustered circRNA locus.
fn overlaps_any_circ_cluster_ctx(records: &[AsAlignment], ctx: &SegmentScanContext<'_>) -> bool {
    for record in records {
        if record.mapq < MAPQ_THRES {
            continue;
        }
        let msid = msid(&record.cigar, ctx.read_len);
        if msid.ref_len < 0 {
            continue;
        }
        let map_end = record.pos + msid.ref_len - 1;
        if let Some(clusters) = ctx.clusters_by_chr.get(&record.chr) {
            for cluster in overlapping_circ_clusters(
                clusters,
                record.pos.min(map_end) - MIN_INTRON,
                record.pos.max(map_end) + MIN_INTRON,
            ) {
                if point_in_cluster_range(record.pos, cluster)
                    || point_in_cluster_range(map_end, cluster)
                {
                    return true;
                }
            }
        }
    }
    false
}

/// Returns whether one non-BSJ read group is pair-level outward circ evidence.
///
/// Unlike `type=backward`, this detector does not require a mate-internal read
/// chain wrap. It keeps primary-pair geometry where the reverse mate and
/// forward mate face outward by coordinate order. Gap length is intentionally
/// not capped because a gap-facing pair can still become useful once projected
/// onto confirmed BSJ loci; exact same-span pairs require 3' terminal clipping
/// on both mates so fully overlapping artifacts are not promoted. The row
/// remains `circ_id=NA`; later graph construction may project it onto every
/// compatible circRNA instead of forcing a read-level unique assignment.
#[cfg(test)]
fn is_outward_pair_group(records: &[AsAlignment], state: &ScanState) -> bool {
    is_outward_pair_group_ctx(records, &SegmentScanContext::from_state(state))
}

/// Shared-context form of outward pair detection for parallel segments shards.
fn is_outward_pair_group_ctx(records: &[AsAlignment], ctx: &SegmentScanContext<'_>) -> bool {
    let Some((r1, r2)) = primary_mate_pair(records) else {
        return false;
    };
    if r1.chr != r2.chr || r1.chr == "*" || r1.mapq < ctx.min_mapq || r2.mapq < ctx.min_mapq {
        return false;
    }
    if has_linear_mate_pair_mapping(records, ctx) {
        return false;
    }
    let Some((r1_start, r1_end)) = alignment_ref_span(r1, ctx.read_len) else {
        return false;
    };
    let Some((r2_start, r2_end)) = alignment_ref_span(r2, ctx.read_len) else {
        return false;
    };
    if !has_3p_outward_pair_geometry(r1, (r1_start, r1_end), r2, (r2_start, r2_end)) {
        return false;
    }
    let span_start = r1_start.min(r2_start);
    let span_end = r1_end.max(r2_end);
    span_end - span_start + 1 <= BACKWARD_MAX_SPAN
}

/// Returns whether any R1/R2 alignment pair can explain the read as linear.
///
/// Outward is weaker than backward because it is inferred from mate orientation
/// rather than a read-internal wrap. Before accepting the primary outward pair,
/// inspect all retained mapper alignments, including secondary/supplementary
/// alternatives, and reject the read if any mate pair forms a conventional
/// forward-left / reverse-right linear mapping within the allowed span.
fn has_linear_mate_pair_mapping(records: &[AsAlignment], ctx: &SegmentScanContext<'_>) -> bool {
    for i in 0..records.len() {
        let left = &records[i];
        if left.flag & 0x4 != 0 || left.chr == "*" || left.mapq < ctx.min_mapq {
            continue;
        }
        for right in &records[i + 1..] {
            if right.flag & 0x4 != 0
                || right.chr == "*"
                || right.mapq < ctx.min_mapq
                || left.chr != right.chr
                || mate_bucket(left.flag) == mate_bucket(right.flag)
            {
                continue;
            }
            let Some(left_span) = alignment_ref_span(left, ctx.read_len) else {
                continue;
            };
            let Some(right_span) = alignment_ref_span(right, ctx.read_len) else {
                continue;
            };
            let span = left_span.0.min(right_span.0)..=left_span.1.max(right_span.1);
            if *span.end() - *span.start() + 1 <= BACKWARD_MAX_SPAN
                && has_linear_pair_geometry(left, left_span, right, right_span)
            {
                return true;
            }
        }
    }
    false
}

/// Tests whether two mate alignments look like an ordinary linear pair.
fn has_linear_pair_geometry(
    left: &AsAlignment,
    left_span: (i32, i32),
    right: &AsAlignment,
    right_span: (i32, i32),
) -> bool {
    let left_reverse = is_reverse_strand(left.flag);
    let right_reverse = is_reverse_strand(right.flag);
    if left_reverse == right_reverse {
        return false;
    }
    let (forward_span, reverse_span) = if left_reverse {
        (right_span, left_span)
    } else {
        (left_span, right_span)
    };
    forward_span != reverse_span
        && forward_span.0 <= reverse_span.0
        && forward_span.1 <= reverse_span.1
}

/// Tests broad outward-facing geometry for a primary R1/R2 pair.
///
/// A valid pair has one reverse and one forward primary alignment. The reverse
/// mate must sit to the left of the forward mate by both start and end
/// coordinates. Non-identical spans are accepted when both aligned outward
/// extensions are at least 19 bp; terminal clips are not required in that case.
/// Exact same-span pairs have zero aligned outward length, so they still need
/// paired 3' terminal clip evidence, with the opposite mate's 5' clip
/// subtracted, to avoid promoting fully overlapping artifacts.
fn has_3p_outward_pair_geometry(
    r1: &AsAlignment,
    r1_span: (i32, i32),
    r2: &AsAlignment,
    r2_span: (i32, i32),
) -> bool {
    let r1_reverse = is_reverse_strand(r1.flag);
    let r2_reverse = is_reverse_strand(r2.flag);
    if r1_reverse == r2_reverse {
        return false;
    }
    let (reverse_span, forward_span) = if r1_reverse {
        (r1_span, r2_span)
    } else {
        (r2_span, r1_span)
    };
    let (reverse_record, forward_record) = if r1_reverse { (r1, r2) } else { (r2, r1) };
    if reverse_span == forward_span {
        return has_paired_outward_terminal_clip(reverse_record, forward_record);
    }
    if forward_span.0 - reverse_span.0 >= OUTWARD_MIN_PAIR_OFFSET
        && forward_span.1 - reverse_span.1 >= OUTWARD_MIN_PAIR_OFFSET
    {
        return true;
    }
    false
}

/// Returns whether both mates carry outward-specific terminal clip evidence.
///
/// SAM CIGAR is query-order regardless of mapper strand, so the first and last
/// operations represent the sequenced 5' and 3' read ends. A real outward pair
/// should have strong 3' tails on both mates; if the opposite mate has a similar
/// 5' clip, that tail can be explained as unaligned mate overlap and is not
/// strong outward evidence.
fn has_paired_outward_terminal_clip(
    reverse_record: &AsAlignment,
    forward_record: &AsAlignment,
) -> bool {
    let Some((reverse_5p, reverse_3p)) = terminal_query_clips(reverse_record) else {
        return false;
    };
    let Some((forward_5p, forward_3p)) = terminal_query_clips(forward_record) else {
        return false;
    };
    reverse_3p >= OUTWARD_MIN_TERMINAL_CLIP
        && forward_3p >= OUTWARD_MIN_TERMINAL_CLIP
        && reverse_3p - forward_5p >= OUTWARD_MIN_TERMINAL_CLIP
        && forward_3p - reverse_5p >= OUTWARD_MIN_TERMINAL_CLIP
}

/// Returns 5' and 3' terminal query clip lengths from a CIGAR string.
fn terminal_query_clips(record: &AsAlignment) -> Option<(i32, i32)> {
    let ops = parse_cigar_ops_basic(&record.cigar)?;
    let five_prime = ops
        .first()
        .copied()
        .filter(|(_, op)| matches!(op, 'S' | 'H'))
        .map_or(0, |(len, _)| len);
    let three_prime = ops
        .last()
        .copied()
        .filter(|(_, op)| matches!(op, 'S' | 'H'))
        .map_or(0, |(len, _)| len);
    Some((five_prime, three_prime))
}

/// Selects the single primary R1/R2 pair used for outward orientation checks.
fn primary_mate_pair(records: &[AsAlignment]) -> Option<(&AsAlignment, &AsAlignment)> {
    let mut r1: Option<&AsAlignment> = None;
    let mut r2: Option<&AsAlignment> = None;
    for record in records {
        if is_secondary(record.flag) || is_supplementary(record.flag) || record.flag & 0x4 != 0 {
            continue;
        }
        match mate_bucket(record.flag) {
            0 if r1.is_none() => r1 = Some(record),
            0 => return None,
            1 if r2.is_none() => r2 = Some(record),
            1 => return None,
            _ => return None,
        }
    }
    Some((r1?, r2?))
}

/// Returns the primary mapper alignment strand for R1 and R2.
///
/// Segment token strands describe RNA/circRNA interpretation and can remain
/// unknown for pair-orientation-only `outward` rows. These fields expose the
/// raw primary alignment orientation used to audit outward geometry without
/// changing the existing read-chain segment semantics.
fn primary_mate_alignment_strands(records: &[AsAlignment]) -> (String, String) {
    let mut strands: [Option<char>; 2] = [None, None];
    for record in records {
        if record.flag & 0x4 != 0 || is_secondary(record.flag) || is_supplementary(record.flag) {
            continue;
        }
        let bucket = mate_bucket(record.flag);
        if strands[bucket].is_none() {
            strands[bucket] = Some(strand_char(record.flag));
        }
    }
    (
        strands[0].map_or_else(|| "NA".to_string(), |strand| strand.to_string()),
        strands[1].map_or_else(|| "NA".to_string(), |strand| strand.to_string()),
    )
}

/// Returns the reference span covered by one alignment's retained blocks.
fn alignment_ref_span(record: &AsAlignment, read_len: i32) -> Option<(i32, i32)> {
    let blocks = parse_alignment_blocks(record, read_len)?;
    let start = blocks.iter().map(|block| block.ref_start).min()?;
    let end = blocks.iter().map(|block| block.ref_end).max()?;
    Some((start, end))
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

/// Runs the non-BSJ mapping check against shard-local segments state.
///
/// In sidecar mode BSJ reads already came from Scan1/Scan2 segment evidence, so
/// the parallel rescan only needs the `from_bsj=false` branch. Keeping this
/// accumulator-specific copy avoids cloning the full `ScanState` per BAM shard.
fn mapping_check1_backward_candidates(
    read_id: &str,
    records: &[AsAlignment],
    ctx: &SegmentScanContext<'_>,
) -> Vec<PositiveCandidate> {
    let mut by_reverse: [Vec<&AsAlignment>; 2] = [Vec::new(), Vec::new()];
    for record in records {
        let idx = if record.flag & 0x10 != 0 { 1 } else { 0 };
        by_reverse[idx].push(record);
    }
    let mut candidates = Vec::new();
    mapping_check2_backward(read_id, &by_reverse[1], ctx, &mut candidates);
    mapping_check2_backward(read_id, &by_reverse[0], ctx, &mut candidates);
    candidates
}

/// Accumulator-backed form of CIRI-AS `mapping_check2_add`.
fn mapping_check2_backward(
    read_id: &str,
    records: &[&AsAlignment],
    ctx: &SegmentScanContext<'_>,
    candidates: &mut Vec<PositiveCandidate>,
) {
    if records.len() < 2 {
        return;
    }

    let mut matches = Vec::with_capacity(records.len());
    for record in records {
        matches.push(msid(&record.cigar, ctx.read_len));
    }

    let initial_len = candidates.len();
    let mut candidate_chr: Option<String> = None;
    let mut na_tag = false;
    'read: for i in 0..records.len() - 1 {
        for j in i + 1..records.len() {
            let ri = records[i];
            let rj = records[j];
            if ri.chr != rj.chr || reverse_bit(ri.flag) != reverse_bit(rj.flag) {
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
                None,
                ctx.read_len,
            ) {
                candidate_chr.get_or_insert_with(|| candidate.chr.clone());
                candidates.push(candidate);
            }
        }
    }

    if na_tag {
        candidates.truncate(initial_len);
    }
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

/// Returns whether the SAM reverse-complement flag is set.
fn is_reverse_strand(flag: i32) -> bool {
    flag & 0x10 != 0
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

/// Applies CIRI-AS splice-signal motif checking and coordinate adjustment.
fn validate_splice_motifs(
    candidates: &mut Vec<PositiveCandidate>,
    reference: &HashMap<String, String>,
    annotation: Option<&Annotation>,
) -> Result<()> {
    // Motif validation is read-local and reference/annotation are immutable, so
    // this can run on Rayon workers without changing the final CIRI-AS ordering.
    let mut validated: Vec<(usize, PositiveCandidate)> = std::mem::take(candidates)
        .into_par_iter()
        .enumerate()
        .filter_map(|(idx, candidate)| {
            validate_splice_candidate(candidate, reference, annotation)
                .map(|candidate| (idx, candidate))
        })
        .collect();
    validated.sort_by_key(|(idx, _)| *idx);
    let mut validated: Vec<PositiveCandidate> = validated
        .into_iter()
        .map(|(_, candidate)| candidate)
        .collect();
    for (idx, candidate) in validated.iter_mut().enumerate() {
        candidate.index = idx;
    }
    *candidates = validated;
    Ok(())
}

/// Validates and motif-adjusts one splice candidate independently.
///
/// This helper keeps the per-candidate work separate from collection ordering so
/// callers can run it in parallel while retaining the serial output contract.
fn validate_splice_candidate(
    mut candidate: PositiveCandidate,
    reference: &HashMap<String, String>,
    annotation: Option<&Annotation>,
) -> Option<PositiveCandidate> {
    let chr_seq = reference.get(&candidate.chr)?;
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
    if hint == 0 {
        return None;
    }
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
    Some(candidate)
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

/// Converts validated read classes into `<prefix>.segments` rows.
///
/// The current phase emits only BSJ reads and non-BSJ backward reads. The
/// classification itself still comes from the existing CIRI-AS-style scan; this
/// helper is only responsible for selecting the best alignment chain and
/// materializing it as simulator-compatible segments.
fn build_segment_records(state: &ScanState) -> Result<Vec<SegmentRecord>> {
    let profile = segments_profile_enabled();
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
    let phase_started = profile.then(Instant::now);
    sort_segment_records(&mut out);
    log_segments_profile(profile, "sort_segment_records", phase_started);
    Ok(out)
}

/// Builds the default sidecar-mode `<prefix>.segments` rows.
///
/// Confirmed BSJ reads are materialized from Scan1/Scan2 sidecar alignments so
/// the output keeps local clip pseudo rows captured during validation. Backward
/// and outward rows come from the supplemental non-BSJ scan and intentionally
/// carry no circRNA assignment. All three evidence classes share the same
/// internal-junction support pass; only the confirmed BSJ gap itself remains
/// type-specific.
fn build_sidecar_segment_records(
    state: &mut ScanState,
    bsj_segment_shard_paths: Vec<SegmentScanShardPaths>,
    rescan_group_shard_paths: Vec<SegmentScanShardPaths>,
    keep_temp_files: bool,
    mut on_correction_start: impl FnMut() -> Result<()>,
) -> Result<Vec<SegmentRecord>> {
    let profile = segments_profile_enabled();
    let first_pass_correction = SegmentCorrectionContext {
        reference: state.reference,
        annotation: state.annotation,
        junction_support: None,
    };

    let phase_started = profile.then(Instant::now);
    let junction_support = collect_sidecar_junction_support_from_shards(
        state,
        &bsj_segment_shard_paths,
        &rescan_group_shard_paths,
        &first_pass_correction,
    )?;
    log_segments_profile(profile, "unified_collect_junction_support", phase_started);
    let phase_started = profile.then(Instant::now);
    let junction_support_index = build_junction_support_index(&junction_support);
    log_segments_profile(
        profile,
        "unified_build_junction_support_index",
        phase_started,
    );
    on_correction_start()?;
    let phase_started = profile.then(Instant::now);
    let correction = SegmentCorrectionContext {
        reference: state.reference,
        annotation: state.annotation,
        junction_support: Some(&junction_support),
    };
    let mut out = build_corrected_sidecar_records_from_shards(
        state,
        &bsj_segment_shard_paths,
        &rescan_group_shard_paths,
        &junction_support_index,
        &first_pass_correction,
        &correction,
    )?;
    let ambiguous_records = out
        .iter()
        .filter(|record| segment_record_has_supported_alternative(record, &junction_support_index))
        .count();
    log_segments_profile(profile, "unified_selective_second_pass", phase_started);
    if profile {
        eprintln!(
            "[CIRI_PROFILE_SEGMENTS] unified_second_pass_records: {ambiguous_records}/{}",
            out.len()
        );
    }
    let phase_started = profile.then(Instant::now);
    sort_segment_records(&mut out);
    log_segments_profile(profile, "sort_segment_records", phase_started);
    if !keep_temp_files {
        cleanup_segment_shards(&bsj_segment_shard_paths);
        cleanup_segment_shards(&rescan_group_shard_paths);
    }
    Ok(out)
}

/// Collects preliminary junction support by streaming retained segment shards.
///
/// Only the support count is global. Preliminary rows and read-local alignments
/// are rebuilt from each shard block and dropped immediately, avoiding the old
/// `Vec<SegmentRecord>` plus read-id offset index during the first pass.
fn collect_sidecar_junction_support_from_shards(
    state: &ScanState,
    bsj_shards: &[SegmentScanShardPaths],
    non_bsj_shards: &[SegmentScanShardPaths],
    correction: &SegmentCorrectionContext<'_>,
) -> Result<JunctionSupportMap> {
    let mut support = HashMap::new();
    let pb = retained_segment_progress_bar(bsj_shards, non_bsj_shards, "")?;
    let pb_ref = pb.as_ref();
    for shard in bsj_shards {
        for_each_retained_segment_group(&shard.path, pb_ref, |group| {
            if let Some(record) =
                build_bsj_segment_record_from_retained_group(state, &group, correction, &[])?
            {
                collect_junction_support_from_record(&record, &mut support);
            }
            Ok(())
        })?;
    }
    for shard in non_bsj_shards {
        for_each_retained_segment_group(&shard.path, pb_ref, |group| {
            for record in &group.segment_records {
                collect_junction_support_from_record(record, &mut support);
            }
            Ok(())
        })?;
    }
    if let Some(pb) = pb {
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {msg}")?
                .progress_chars("#>-"),
        );
        pb.finish_with_message("Completed");
    }
    Ok(support)
}

/// Builds final segment rows by streaming retained shards after support is known.
fn build_corrected_sidecar_records_from_shards(
    state: &ScanState,
    bsj_shards: &[SegmentScanShardPaths],
    non_bsj_shards: &[SegmentScanShardPaths],
    support_index: &JunctionSupportIndex,
    first_pass_correction: &SegmentCorrectionContext<'_>,
    correction: &SegmentCorrectionContext<'_>,
) -> Result<Vec<SegmentRecord>> {
    let mut out = Vec::new();
    let pb = retained_segment_progress_bar(bsj_shards, non_bsj_shards, "")?;
    let pb_ref = pb.as_ref();
    for shard in bsj_shards {
        for_each_retained_segment_group(&shard.path, pb_ref, |group| {
            let Some(preliminary) = build_bsj_segment_record_from_retained_group(
                state,
                &group,
                first_pass_correction,
                &[],
            )?
            else {
                return Ok(());
            };
            let record = if segment_record_has_supported_alternative(&preliminary, support_index) {
                build_bsj_segment_record_from_retained_group(state, &group, correction, &[])?
                    .unwrap_or(preliminary)
            } else {
                preliminary
            };
            out.push(record);
            Ok(())
        })?;
    }
    for shard in non_bsj_shards {
        for_each_retained_segment_group(&shard.path, pb_ref, |group| {
            let strand_hint = read_strand_hint_from_candidates(&group.candidates);
            let junction_hints = read_junction_hints_from_candidates(&group.candidates);
            for preliminary in &group.segment_records {
                let record = if segment_record_has_supported_alternative(preliminary, support_index)
                {
                    rebuild_segment_record_from_retained_group(
                        preliminary,
                        &group,
                        state,
                        correction,
                        strand_hint,
                        &junction_hints,
                    )
                    .unwrap_or_else(|| preliminary.clone())
                } else {
                    preliminary.clone()
                };
                out.push(record);
            }
            Ok(())
        })?;
    }
    if let Some(pb) = pb {
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {msg}")?
                .progress_chars("#>-"),
        );
        pb.finish_with_message("Completed");
    }
    Ok(out)
}

/// Builds one confirmed BSJ row from a retained read block.
fn build_bsj_segment_record_from_retained_group(
    state: &ScanState,
    group: &RetainedSegmentGroup,
    correction: &SegmentCorrectionContext<'_>,
    junction_hints: &[(i32, i32)],
) -> Result<Option<SegmentRecord>> {
    let read_id = group.read_id.as_str();
    let Some(circ_id) = state.junction_read_to_circ.get(read_id) else {
        return Ok(None);
    };
    let circ = state
        .circ_by_id
        .get(circ_id)
        .ok_or_else(|| anyhow!("missing circ record {}", circ_id))?;
    Ok(build_bsj_segment_record(
        read_id,
        &group.records,
        circ,
        state
            .mate_bsj_evidence
            .get(read_id)
            .map(Vec::as_slice)
            .unwrap_or(&[]),
        junction_hints,
        Some(correction),
        state.read_len,
    ))
}

/// Rebuilds one retained non-BSJ preliminary row with support-aware correction.
fn rebuild_segment_record_from_retained_group(
    record: &SegmentRecord,
    group: &RetainedSegmentGroup,
    state: &ScanState,
    correction: &SegmentCorrectionContext<'_>,
    strand_hint: Option<char>,
    junction_hints: &[(i32, i32)],
) -> Option<SegmentRecord> {
    match record.type_name {
        "backward" => build_backward_segment_record_from_records(
            group.read_id.as_str(),
            &group.records,
            state,
            correction,
            strand_hint,
            junction_hints,
        ),
        "outward" => build_outward_segment_record(
            group.read_id.as_str(),
            &group.records,
            Some(correction),
            state.read_len,
            state.min_mapq,
        ),
        _ => None,
    }
}

/// Deletes retained segment shard files after streamed finalize is complete.
fn cleanup_segment_shards(shards: &[SegmentScanShardPaths]) {
    for shard in shards {
        let _ = std::fs::remove_file(&shard.path);
    }
}

/// Builds one backward row from a streamed retained group and local hints.
fn build_backward_segment_record_from_records(
    read_id: &str,
    records: &[AsAlignment],
    state: &ScanState,
    correction: &SegmentCorrectionContext<'_>,
    strand_hint: Option<char>,
    junction_hints: &[(i32, i32)],
) -> Option<SegmentRecord> {
    let base = build_backward_segment_record(
        read_id,
        records,
        strand_hint,
        junction_hints,
        Some(correction),
        state.read_len,
    )?;
    let enriched_records = add_non_bsj_local_clip_alignments(records, state);
    if enriched_records.len() == records.len() {
        return Some(base);
    }
    build_backward_segment_record(
        read_id,
        &enriched_records,
        strand_hint,
        junction_hints,
        Some(correction),
        state.read_len,
    )
    .or(Some(base))
}

/// Adds one preliminary segment row to the global junction-support map.
fn collect_junction_support_from_record(record: &SegmentRecord, support: &mut JunctionSupportMap) {
    let chrom_support = support.entry(record.chrom.clone()).or_default();
    collect_junction_support_from_segments(&record.r1_segments, chrom_support);
    collect_junction_support_from_segments(&record.r2_segments, chrom_support);
}

/// Adds retained read-chain `N` junctions from one mate's segment string.
fn collect_junction_support_from_segments(
    segments: &str,
    chrom_support: &mut HashMap<JunctionSupportKey, usize>,
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
                *chrom_support.entry((site2, site1, strand)).or_insert(0) += 1;
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
    let (r1_align_strand, r2_align_strand) = primary_mate_alignment_strands(records);
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
        r1_align_strand,
        r2_align_strand,
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
    let (r1_align_strand, r2_align_strand) = primary_mate_alignment_strands(records);
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
        r1_align_strand,
        r2_align_strand,
        r1_cigar,
        r1_segments,
        r2_cigar,
        r2_segments,
    })
}

/// Builds one `type=outward` row from a primary outward R1/R2 pair.
///
/// This materializes only ordinary linear mate chains. Pair orientation is the
/// circular signal, so the row must not write `<bsj>` or `B`; any internal `N`
/// operators come solely from the original mate CIGAR and can support graph
/// exon-exon edges independently of the outward pair-level evidence.
fn build_outward_segment_record(
    read_id: &str,
    records: &[AsAlignment],
    correction: Option<&SegmentCorrectionContext<'_>>,
    read_len: i32,
    min_mapq: i32,
) -> Option<SegmentRecord> {
    let (r1, r2) = primary_mate_pair(records)?;
    if r1.mapq < min_mapq || r2.mapq < min_mapq {
        return None;
    }
    let mut chains = build_pair_chains(records, read_len, None, "outward", '?', &[], correction);
    if chains.iter().flatten().any(|chain| chain.is_circular) {
        let selected = [r1.clone(), r2.clone()];
        chains = build_pair_chains(&selected, read_len, None, "outward", '?', &[], correction);
    }
    if chains.iter().flatten().any(|chain| chain.is_circular) {
        return None;
    }
    let token_strand = infer_outward_token_strand(&chains, correction).unwrap_or('?');
    if token_strand != '?' {
        chains = build_pair_chains(
            records,
            read_len,
            None,
            "outward",
            token_strand,
            &[],
            correction,
        );
        if chains.iter().flatten().any(|chain| chain.is_circular) {
            let selected = [r1.clone(), r2.clone()];
            chains = build_pair_chains(
                &selected,
                read_len,
                None,
                "outward",
                token_strand,
                &[],
                correction,
            );
        }
        if chains.iter().flatten().any(|chain| chain.is_circular) {
            return None;
        }
    }
    let chrom = selected_chain_chrom(&chains)?;
    let (span_start, span_end) = selected_chain_span(&chains)?;
    if span_end - span_start + 1 > BACKWARD_MAX_SPAN {
        return None;
    }
    let r1_span = chain_ref_span(chains[0].as_ref())?;
    let r2_span = chain_ref_span(chains[1].as_ref())?;
    if !has_3p_outward_pair_geometry(r1, r1_span, r2, r2_span) {
        return None;
    }
    let (r1_align_strand, r2_align_strand) = primary_mate_alignment_strands(records);
    let (r1_segments, r1_cigar, _) = chain_text(chains[0].as_ref());
    let (r2_segments, r2_cigar, _) = chain_text(chains[1].as_ref());
    Some(SegmentRecord {
        read_id: read_id.to_string(),
        type_name: "outward",
        circ_id: "NA".to_string(),
        chrom,
        start: span_start.to_string(),
        end: span_end.to_string(),
        strand: token_strand.to_string().replace('?', "NA"),
        is_circular: 1,
        is_r1_bsj: 0,
        is_r2_bsj: 0,
        r1_align_strand,
        r2_align_strand,
        r1_cigar,
        r1_segments,
        r2_cigar,
        r2_segments,
    })
}

/// Infers an outward row's RNA strand from mate-internal splice junctions.
///
/// Outward evidence itself is pair-orientation based and therefore does not
/// imply a transcript strand. When one mate already contains ordinary `N`
/// junctions, however, annotation or splice motifs can provide the same strand
/// evidence used by internal-boundary correction. The helper accepts only a
/// unique non-conflicting strand so alignment orientation is never reported as
/// biological strand by fallback.
fn infer_outward_token_strand(
    chains: &[Option<MateChain>; 2],
    correction: Option<&SegmentCorrectionContext<'_>>,
) -> Option<char> {
    let correction = correction?;
    let mut inferred: Option<char> = None;
    for chain in chains.iter().flatten() {
        for pair in chain.blocks.windows(2) {
            let (left, right) = if pair[0].ref_start <= pair[1].ref_start {
                (&pair[0], &pair[1])
            } else {
                (&pair[1], &pair[0])
            };
            let end = left.ref_end;
            let start = right.ref_start;
            if end >= start {
                continue;
            }
            let plus = splice_strand_evidence_score(correction, &chain.chrom, end, start, '+');
            let minus = splice_strand_evidence_score(correction, &chain.chrom, end, start, '-');
            let strand = match plus.cmp(&minus) {
                std::cmp::Ordering::Greater if plus > 0 => '+',
                std::cmp::Ordering::Less if minus > 0 => '-',
                _ => continue,
            };
            match inferred {
                Some(existing) if existing != strand => return None,
                Some(_) => {}
                None => inferred = Some(strand),
            }
        }
    }
    inferred
}

/// Scores one exact internal junction as strand evidence.
///
/// Transcript-consistent annotation is strongest, boundary-level annotation is
/// next, and de novo splice motifs are used only when annotation is absent.
/// Scores are comparable only between `+` and `-` for the same junction.
fn splice_strand_evidence_score(
    ctx: &SegmentCorrectionContext<'_>,
    chrom: &str,
    end: i32,
    start: i32,
    strand: char,
) -> i32 {
    let (transcript_score, annotation_score) =
        annotation_splice_pair_score(ctx.annotation, chrom, end, start, strand).unwrap_or((0, 0));
    let motif_score = ctx
        .reference
        .get(chrom)
        .and_then(|seq| splice_motif_score(seq, end, start, strand))
        .unwrap_or(0);
    transcript_score * 10_000 + annotation_score * 100 + motif_score
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
/// anchor so tiny clips do not drive the decision, and deliberately avoid both
/// MAPQ and maximum-anchor limits because exact XA alternatives can also
/// explain longer repetitive supplementary blocks that BWA scored confidently
/// in a non-splice-aware placement.
fn ambiguous_alignment_with_xa(record: &AsAlignment, read_len: i32) -> bool {
    if record.xa_alternatives.is_empty() || record.from_local_clip {
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
    records.sort_by_cached_key(|record| {
        (
            record.chrom.clone(),
            sort_position(&record.start),
            sort_position(&record.end),
            type_sort_rank(record.type_name),
            record.circ_id.clone(),
            record.read_id.clone(),
        )
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
        "outward" => 2,
        "forward" => 3,
        _ => 4,
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

/// Returns the genomic span covered by one selected mate chain.
fn chain_ref_span(chain: Option<&MateChain>) -> Option<(i32, i32)> {
    let chain = chain?;
    let mut start = i32::MAX;
    let mut end = i32::MIN;
    for block in &chain.blocks {
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
        type_name,
        token_strand,
        junction_hints,
        correction,
        reverse_chain_order,
    );
    let fallback_chain = build_chain_from_pool(
        records,
        read_len,
        circ,
        type_name,
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
    type_name: &str,
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
            chain_generic_rank(&chain, type_name, circ)
                > chain_generic_rank(current, type_name, circ)
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
    let boundary_gap_idxs = if is_bsj {
        circ.and_then(|circ| read_order_bsj_gap_index(&output_blocks, circ, order_strand))
            .into_iter()
            .collect::<Vec<_>>()
    } else if is_circular {
        read_order_wrap_gap_indexes(&output_blocks, order_strand)
    } else {
        Vec::new()
    };
    apply_segment_boundary_corrections(
        &mut output_blocks,
        &chrom,
        &boundary_gap_idxs,
        junction_hints,
        correction,
        token_strand,
    );
    let (tokens, _token_spans, cigar) = materialize_read_chain_output(
        &output_blocks,
        token_strand,
        &boundary_gap_idxs,
        read_len,
        reverse_chain_order,
    );

    Some(MateChain {
        chrom,
        token_strand,
        blocks: output_blocks,
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
/// read-chain blocks. `B` is deliberately tied to read-order circular wraps:
/// confirmed BSJ reads mark the annotated circ boundary, while `type=backward`
/// rows mark every selected wrap even when the boundary is not a confirmed BSJ.
/// Re-sorting these blocks by coordinate would turn `C|B|A` circRNA evidence
/// into a linear-looking `A|B|C` chain.
fn materialize_read_chain_output(
    blocks: &[SegmentBlock],
    token_strand: char,
    boundary_gap_idxs: &[usize],
    read_len: i32,
    reverse_chain_order: bool,
) -> (Vec<String>, Vec<(i32, i32)>, String) {
    let mut token_spans = Vec::with_capacity(blocks.len());
    let mut tokens = Vec::with_capacity(blocks.len() + boundary_gap_idxs.len());
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
            let op = if boundary_gap_idxs.contains(&idx) {
                'B'
            } else {
                'N'
            };
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

/// Returns read-chain gap indexes that wrap against strand-specific linear order.
///
/// These indexes are used for non-BSJ backward rows. They are not treated as
/// strong BSJ evidence, but the `.segments` CIGAR still needs `B` operators at
/// the actual read-order wrap positions so downstream graph code can distinguish
/// backward topology from ordinary splice `N` gaps.
fn read_order_wrap_gap_indexes(read_blocks: &[SegmentBlock], order_strand: char) -> Vec<usize> {
    read_blocks
        .windows(2)
        .enumerate()
        .filter_map(|(idx, pair)| {
            wraps_in_read_order(&pair[0], &pair[1], order_strand).then_some(idx + 1)
        })
        .collect()
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
    boundary_gap_idxs: &[usize],
    junction_hints: &[(i32, i32)],
    correction: Option<&SegmentCorrectionContext<'_>>,
    token_strand: char,
) {
    if blocks.len() < 2 {
        return;
    }
    for idx in 0..blocks.len() - 1 {
        if boundary_gap_idxs.contains(&(idx + 1)) {
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
fn chain_generic_rank(
    chain: &MateChain,
    type_name: &str,
    circ: Option<&CircRecord>,
) -> (i32, i32, i32, i32, i32) {
    let topology_rank = match type_name {
        "outward" => i32::from(!chain.is_circular),
        _ => i32::from(chain.is_circular),
    };
    (
        i32::from(circ.is_none_or(|circ| chain.chrom == circ.chr)),
        topology_rank,
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
        "outward" => i32::from(!chain.is_circular),
        _ => 0,
    };
    let topology_rank = match type_name {
        "outward" => i32::from(!chain.is_circular),
        _ => i32::from(chain.is_circular),
    };
    (
        valid_type,
        i32::from(circ.is_none_or(|circ| chain.chrom == circ.chr)),
        topology_rank,
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

/// Summarizes a streamed read group's validated candidate strand hint.
fn read_strand_hint_from_candidates(candidates: &[PositiveCandidate]) -> Option<char> {
    let mut hint = 0;
    for candidate in candidates {
        if candidate.strand_hint == 0 {
            continue;
        }
        if hint != 0 && hint != candidate.strand_hint {
            return None;
        }
        hint = candidate.strand_hint;
    }
    match hint {
        -1 => Some('+'),
        1 => Some('-'),
        _ => None,
    }
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

/// Returns local `(site2, site1)` hints for one streamed retained read group.
fn read_junction_hints_from_candidates(candidates: &[PositiveCandidate]) -> Vec<(i32, i32)> {
    candidates
        .iter()
        .map(|candidate| (candidate.site2, candidate.site1))
        .collect()
}

/// Writes the simplified `<prefix>.segments` protocol.
fn write_segments(path: &str, records: &[SegmentRecord]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "read_id\ttype\tcirc_id\tchrom\tstart\tend\tstrand\tis_circular\tis_r1_bsj\tis_r2_bsj\tr1_align_strand\tr2_align_strand\tr1_cigar\tr1_segments\tr2_cigar\tr2_segments"
    )?;
    for record in records {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
            record.r1_align_strand,
            record.r2_align_strand,
            record.r1_cigar,
            record.r1_segments,
            record.r2_cigar,
            record.r2_segments
        )?;
    }
    Ok(())
}

/// Builds and writes one major full-length isoform for every confirmed circRNA.
///
/// This is the Rust integration of the validated phase-seed-union prototype:
/// BSJ/backward reads define high-confidence phasing edges, outward reads can
/// fill compatible gaps with low weight, and ambiguous non-BSJ rows are ignored
/// unless their segment span is uniquely contained by one final circRNA.
fn build_major_isoforms_from_segments_file(
    circ_records: &[CircRecord],
    segments_path: &str,
    out_prefix: &str,
    reference: &HashMap<String, String>,
    annotation: Option<&Annotation>,
) -> Result<IsoformRunSummary> {
    let circ_by_id: HashMap<&str, usize> = circ_records
        .iter()
        .enumerate()
        .map(|(idx, circ)| (circ.id.as_str(), idx))
        .collect();
    let circ_index_by_chr = build_major_circ_index(circ_records);
    let mut edge_support_by_circ: HashMap<usize, HashMap<MajorEdge, MajorEdgeSupport>> =
        HashMap::new();
    let mut link_support_by_circ: HashMap<usize, HashMap<MajorJunctionLink, MajorLinkSupport>> =
        HashMap::new();
    let link_exclusion_by_circ: HashMap<usize, HashMap<MajorJunctionLink, MajorLinkSupport>> =
        HashMap::new();
    let mut span_support_by_circ: HashMap<usize, Vec<MajorAlignedSpan>> = HashMap::new();

    let file =
        File::open(segments_path).with_context(|| format!("open segments {}", segments_path))?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .transpose()?
        .ok_or_else(|| anyhow!("empty segments file: {}", segments_path))?;
    let columns = major_segments_columns(&header, segments_path)?;
    for line_result in lines {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        let Some(row) = major_segment_row(&cols, columns) else {
            continue;
        };
        let assignments = assign_major_circs_from_segments(
            row.type_name,
            row.circ_id,
            row.chrom,
            row.start,
            row.end,
            row.r1_segments,
            row.r2_segments,
            &circ_by_id,
            &circ_index_by_chr,
            circ_records,
        );
        if assignments.is_empty() {
            continue;
        }
        for assignment in assignments {
            let circ_index = assignment.circ_index;
            let circ = &circ_records[circ_index];
            let r1_is_chimeric = major_segment_has_reused_overlap(row.r1_segments, circ);
            let r2_is_chimeric = major_segment_has_reused_overlap(row.r2_segments, circ);
            let r1_junctions = if r1_is_chimeric {
                Vec::new()
            } else {
                major_segment_junction_chain(row.r1_segments, circ)
            };
            let r2_junctions = if r2_is_chimeric {
                Vec::new()
            } else {
                major_segment_junction_chain(row.r2_segments, circ)
            };
            let record_spans: HashSet<(i32, i32)> = (!r1_is_chimeric)
                .then(|| major_contiguous_segment_spans(row.r1_segments, circ))
                .into_iter()
                .flatten()
                .chain(
                    (!r2_is_chimeric)
                        .then(|| major_contiguous_segment_spans(row.r2_segments, circ))
                        .into_iter()
                        .flatten(),
                )
                .collect();
            if !record_spans.is_empty() {
                let circ_spans = span_support_by_circ.entry(circ_index).or_default();
                for (start, end) in record_spans {
                    let mut span = MajorAlignedSpan {
                        start,
                        end,
                        bsj: 0.0,
                        backward: 0.0,
                        outward: 0.0,
                    };
                    match row.type_name {
                        "bsj" => span.bsj = assignment.weight,
                        "backward" => span.backward = assignment.weight,
                        "outward" => span.outward = assignment.weight,
                        _ => {}
                    }
                    circ_spans.push(span);
                }
            }
            let record_edges: HashSet<MajorEdge> = r1_junctions
                .iter()
                .chain(r2_junctions.iter())
                .filter_map(|junction| match junction {
                    MajorJunction::Edge(edge) => Some(*edge),
                    MajorJunction::Bsj => None,
                })
                .collect();
            if !record_edges.is_empty() {
                let circ_support = edge_support_by_circ.entry(circ_index).or_default();
                for edge in record_edges {
                    let support = circ_support.entry(edge).or_default();
                    match row.type_name {
                        "bsj" => support.bsj += assignment.weight,
                        "backward" => support.backward += assignment.weight,
                        "outward" => support.outward += assignment.weight,
                        _ => {}
                    }
                }
            }
            let record_links: HashSet<MajorJunctionLink> =
                major_adjacent_junction_links(&r1_junctions)
                    .into_iter()
                    .chain(major_adjacent_junction_links(&r2_junctions))
                    .collect();
            if !record_links.is_empty() {
                let circ_links = link_support_by_circ.entry(circ_index).or_default();
                for &link in &record_links {
                    let support = circ_links.entry(link).or_default();
                    match row.type_name {
                        "bsj" => support.bsj += assignment.weight,
                        "backward" => support.backward += assignment.weight,
                        "outward" => support.outward += assignment.weight,
                        _ => {}
                    }
                }
            }
        }
    }

    let sample_id = major_sample_id(out_prefix);
    let mut isoforms = Vec::with_capacity(circ_records.len());
    for (idx, circ) in circ_records.iter().enumerate() {
        let empty = HashMap::new();
        let empty_links = HashMap::new();
        let empty_exclusions = HashMap::new();
        let empty_spans = Vec::new();
        let support = edge_support_by_circ.get(&idx).unwrap_or(&empty);
        let link_support = link_support_by_circ.get(&idx).unwrap_or(&empty_links);
        let link_exclusion = link_exclusion_by_circ
            .get(&idx)
            .unwrap_or(&empty_exclusions);
        let span_support = span_support_by_circ.get(&idx).unwrap_or(&empty_spans);
        isoforms.push(select_major_isoform(
            circ,
            support,
            link_support,
            link_exclusion,
            span_support,
            &sample_id,
            annotation,
        ));
    }
    isoforms.sort_by(|a, b| {
        a.chr
            .cmp(&b.chr)
            .then_with(|| a.start.cmp(&b.start))
            .then_with(|| a.end.cmp(&b.end))
            .then_with(|| a.circ_id.cmp(&b.circ_id))
    });

    let total_isoforms = isoforms.len();
    let circ_rnas = isoforms
        .iter()
        .map(|record| record.circ_id.as_str())
        .collect::<HashSet<_>>()
        .len();
    write_major_isoform_gtf(&format!("{}.isoforms.gtf", out_prefix), &isoforms)?;
    let fasta_summary =
        write_major_isoform_fasta(&format!("{}.isoforms.fa", out_prefix), &isoforms, reference)?;
    Ok(IsoformRunSummary {
        total_isoforms,
        circ_rnas,
        fasta_isoforms: fasta_summary.isoforms,
        fasta_circ_rnas: fasta_summary.circ_rnas,
    })
}

/// Borrowed view of one `<prefix>.segments` row used by the isoform stage.
struct MajorSegmentRow<'a> {
    type_name: &'a str,
    circ_id: &'a str,
    chrom: &'a str,
    start: &'a str,
    end: &'a str,
    r1_segments: &'a str,
    r2_segments: &'a str,
}

/// Resolves required `<prefix>.segments` columns by header name.
fn major_segments_columns(header: &str, path: &str) -> Result<MajorSegmentsColumns> {
    let headers: Vec<&str> = header.split('\t').collect();
    let idx = |name: &str| -> Result<usize> {
        headers
            .iter()
            .position(|field| *field == name)
            .ok_or_else(|| anyhow!("missing `{}` column in {}", name, path))
    };
    Ok(MajorSegmentsColumns {
        type_name: idx("type")?,
        circ_id: idx("circ_id")?,
        chrom: idx("chrom")?,
        start: idx("start")?,
        end: idx("end")?,
        r1_segments: idx("r1_segments")?,
        r2_segments: idx("r2_segments")?,
    })
}

/// Projects a parsed column slice into the fields needed for graph building.
fn major_segment_row<'a>(
    cols: &'a [&'a str],
    indexes: MajorSegmentsColumns,
) -> Option<MajorSegmentRow<'a>> {
    Some(MajorSegmentRow {
        type_name: *cols.get(indexes.type_name)?,
        circ_id: *cols.get(indexes.circ_id)?,
        chrom: *cols.get(indexes.chrom)?,
        start: *cols.get(indexes.start)?,
        end: *cols.get(indexes.end)?,
        r1_segments: *cols.get(indexes.r1_segments)?,
        r2_segments: *cols.get(indexes.r2_segments)?,
    })
}

/// Builds a start-sorted exact circRNA span index for non-BSJ row assignment.
fn build_major_circ_index(
    circ_records: &[CircRecord],
) -> HashMap<String, Vec<MajorCircIndexEntry>> {
    let mut by_chr: HashMap<String, Vec<MajorCircIndexEntry>> = HashMap::new();
    for (idx, circ) in circ_records.iter().enumerate() {
        by_chr
            .entry(circ.chr.clone())
            .or_default()
            .push(MajorCircIndexEntry {
                start: circ.start,
                end: circ.end,
                circ_index: idx,
                max_end_through: circ.end,
            });
    }
    for entries in by_chr.values_mut() {
        entries.sort_by_key(|entry| (entry.start, entry.end, entry.circ_index));
        let mut max_end = i32::MIN;
        for entry in entries {
            max_end = max_end.max(entry.end);
            entry.max_end_through = max_end;
        }
    }
    by_chr
}

/// Assigns a segment row to final circRNAs for graph reconstruction.
///
/// BSJ rows use their Summary circRNA ID directly. Backward/outward rows do not
/// have a trustworthy BSJ assignment, so every containing confirmed circRNA is
/// scored by how close the read's observed segment endpoints are to that
/// circRNA's BSJ boundaries. The normalized weights sum to at most one, with
/// the remaining probability intentionally left unassigned when the read looks
/// internal to an enclosing locus.
fn assign_major_circs_from_segments(
    type_name: &str,
    circ_id: &str,
    chrom: &str,
    start_text: &str,
    end_text: &str,
    r1_segments: &str,
    r2_segments: &str,
    circ_by_id: &HashMap<&str, usize>,
    circ_index_by_chr: &HashMap<String, Vec<MajorCircIndexEntry>>,
    circ_records: &[CircRecord],
) -> Vec<MajorCircAssignment> {
    if type_name == "bsj" {
        return circ_by_id
            .get(circ_id)
            .copied()
            .map(|circ_index| {
                vec![MajorCircAssignment {
                    circ_index,
                    weight: 1.0,
                }]
            })
            .unwrap_or_default();
    }
    let Ok(start) = start_text.parse::<i32>() else {
        return Vec::new();
    };
    let Ok(end) = end_text.parse::<i32>() else {
        return Vec::new();
    };
    let Some(entries) = circ_index_by_chr.get(chrom) else {
        return Vec::new();
    };
    let endpoints = major_segment_endpoints_for_assignment(r1_segments, r2_segments);
    if endpoints.is_empty() {
        return Vec::new();
    }
    let scale = major_non_bsj_assignment_distance_scale(&endpoints);
    let right = entries.partition_point(|entry| entry.start <= start);
    let mut raw = Vec::new();
    for entry in entries[..right]
        .iter()
        .rev()
        .take_while(|entry| entry.max_end_through >= end)
        .filter(|entry| entry.end >= end)
    {
        let circ = &circ_records[entry.circ_index];
        if circ.chr == chrom && circ.start <= start && circ.end >= end {
            let weight = major_non_bsj_assignment_weight(circ, &endpoints, scale);
            if weight >= MAJOR_MIN_NON_BSJ_ASSIGNMENT_WEIGHT {
                raw.push((entry.circ_index, weight));
            }
        }
    }
    let total: f64 = raw.iter().map(|(_, weight)| *weight).sum();
    if total <= 0.0 {
        return Vec::new();
    }
    let normalizer = total.max(1.0);
    raw.into_iter()
        .filter_map(|(circ_index, weight)| {
            let weight = weight / normalizer;
            (weight >= MAJOR_MIN_NON_BSJ_ASSIGNMENT_WEIGHT)
                .then_some(MajorCircAssignment { circ_index, weight })
        })
        .collect()
}

/// Collects read-level segment endpoints used to score non-BSJ circ assignment.
///
/// Only parsed `start-end:strand` tokens are considered; `<bsj>` markers remain
/// topology labels and should not create artificial coordinates. Both block
/// boundaries are retained because a circular fragment can expose either side
/// of the BSJ depending on mate orientation and clipping.
fn major_segment_endpoints_for_assignment(r1_segments: &str, r2_segments: &str) -> Vec<(i32, i32)> {
    r1_segments
        .split('|')
        .chain(r2_segments.split('|'))
        .filter_map(parse_segment_token)
        .map(|(start, end, _)| (start, end))
        .collect()
}

/// Estimates the fragment-scale distance window for non-BSJ assignment.
///
/// The scale follows aligned read bases rather than genomic span so introns and
/// internal circular wraps cannot make a far-away read look boundary-compatible
/// with a large outer circRNA.
fn major_non_bsj_assignment_distance_scale(endpoints: &[(i32, i32)]) -> f64 {
    let aligned_bases: i32 = endpoints
        .iter()
        .map(|(start, end)| end - start + 1)
        .filter(|len| *len > 0)
        .sum();
    (aligned_bases as f64 * 2.0).clamp(
        MAJOR_NON_BSJ_MIN_DISTANCE_SCALE,
        MAJOR_NON_BSJ_MAX_DISTANCE_SCALE,
    )
}

/// Scores how likely one non-BSJ row belongs to a candidate circRNA.
///
/// A true boundary-compatible backward/outward read should place observed
/// segment ends close to both BSJ sides within the read-pair fragment scale.
/// Internal reads from a large enclosing locus receive an exponentially small
/// probability and are mostly left unassigned.
fn major_non_bsj_assignment_weight(circ: &CircRecord, endpoints: &[(i32, i32)], scale: f64) -> f64 {
    let left_dist = major_min_endpoint_distance(circ.start, endpoints);
    let right_dist = major_min_endpoint_distance(circ.end, endpoints);
    let combined_dist = left_dist.saturating_add(right_dist) as f64;
    (-combined_dist / scale).exp()
}

/// Returns the nearest distance from a coordinate to any observed segment end.
fn major_min_endpoint_distance(position: i32, endpoints: &[(i32, i32)]) -> i32 {
    endpoints
        .iter()
        .flat_map(|(start, end)| [*start, *end])
        .map(|edge| (edge - position).abs())
        .min()
        .unwrap_or(i32::MAX)
}

/// Returns BSJ and internal splice junctions from one mate segment string.
///
/// Edges crossing an explicit `<bsj>` marker are represented by a `Bsj` token
/// instead of an internal edge. Keeping that boundary token in read order lets
/// mature classification require support for BSJ-to-internal neighboring
/// junctions, which is the only way to phase two-exon circRNA blocks from the
/// segment graph itself.
fn major_segment_junction_chain(segments: &str, circ: &CircRecord) -> Vec<MajorJunction> {
    let mut out = Vec::new();
    let mut previous: Option<(i32, i32, char)> = None;
    let mut pending_bsj = false;
    for token in segments.split('|') {
        if token == "<bsj>" {
            if previous.is_some() {
                out.push(MajorJunction::Bsj);
            }
            pending_bsj = previous.is_some();
            continue;
        }
        let Some((start, end, strand)) = parse_segment_token(token) else {
            previous = None;
            pending_bsj = false;
            continue;
        };
        if end - start + 1 < MIN_JUNCTION_SUPPORT_SEGMENT_LENGTH {
            continue;
        }
        if let Some((prev_start, prev_end, prev_strand)) = previous {
            if !pending_bsj && prev_strand == strand {
                let (donor_end, acceptor_start) = if prev_start <= start {
                    (prev_end, start)
                } else {
                    (end, prev_start)
                };
                if donor_end < acceptor_start
                    && donor_end >= circ.start
                    && acceptor_start <= circ.end
                {
                    out.push(MajorJunction::Edge(MajorEdge {
                        donor_end,
                        acceptor_start,
                        strand,
                    }));
                }
            }
        }
        previous = Some((start, end, strand));
        pending_bsj = false;
    }
    out
}

/// Returns continuous aligned blocks from one mate segment string.
///
/// These spans are not positive splice evidence. They are used later as
/// junction-exclusive evidence when an annotation-only candidate junction falls
/// entirely inside one uninterrupted alignment block.
fn major_contiguous_segment_spans(segments: &str, circ: &CircRecord) -> Vec<(i32, i32)> {
    let mut out = Vec::new();
    for token in segments.split('|') {
        if token == "<bsj>" {
            continue;
        }
        let Some((start, end, _)) = parse_segment_token(token) else {
            continue;
        };
        if end - start + 1 < MIN_JUNCTION_SUPPORT_SEGMENT_LENGTH {
            continue;
        }
        let clipped_start = start.max(circ.start);
        let clipped_end = end.min(circ.end);
        if clipped_start <= clipped_end {
            out.push((clipped_start, clipped_end));
        }
    }
    out
}

/// Detects cyclic reuse of aligned bases within one mate segment chain.
///
/// A mate that crosses multiple circular/splice junctions can contain positive
/// neighboring-junction tokens while also reusing the same genomic bases after a
/// later junction. In isoform reconstruction this is treated as a chimeric
/// reverse-transcription artifact and the affected mate is filtered before it
/// can add spans, splice edges, or phasing links. The check is intentionally
/// limited to overlaps between different junction-separated groups so adjacent
/// blocks from the same uninterrupted mapping are not penalized.
fn major_segment_has_reused_overlap(segments: &str, circ: &CircRecord) -> bool {
    let mut groups: Vec<(i32, i32, usize)> = Vec::new();
    let mut previous: Option<(i32, i32, char)> = None;
    let mut pending_bsj = false;
    let mut group_id = 0usize;
    for token in segments.split('|') {
        if token == "<bsj>" {
            if previous.is_some() {
                group_id += 1;
            }
            pending_bsj = previous.is_some();
            continue;
        }
        let Some((start, end, strand)) = parse_segment_token(token) else {
            previous = None;
            pending_bsj = false;
            continue;
        };
        if end - start + 1 < MIN_JUNCTION_SUPPORT_SEGMENT_LENGTH {
            continue;
        }
        if let Some((prev_start, prev_end, prev_strand)) = previous {
            if !pending_bsj && prev_strand == strand {
                let (donor_end, acceptor_start) = if prev_start <= start {
                    (prev_end, start)
                } else {
                    (end, prev_start)
                };
                if donor_end < acceptor_start
                    && donor_end >= circ.start
                    && acceptor_start <= circ.end
                {
                    group_id += 1;
                }
            }
        }
        let clipped_start = start.max(circ.start);
        let clipped_end = end.min(circ.end);
        if clipped_start <= clipped_end {
            for &(prior_start, prior_end, prior_group) in &groups {
                if prior_group != group_id {
                    let overlap = prior_end.min(clipped_end) - prior_start.max(clipped_start) + 1;
                    if overlap >= MIN_JUNCTION_SUPPORT_SEGMENT_LENGTH {
                        return true;
                    }
                }
            }
            groups.push((clipped_start, clipped_end, group_id));
        }
        previous = Some((start, end, strand));
        pending_bsj = false;
    }
    false
}

/// Builds unique neighboring-junction link keys from one read-chain list.
fn major_adjacent_junction_links(junctions: &[MajorJunction]) -> Vec<MajorJunctionLink> {
    junctions
        .windows(2)
        .map(|pair| MajorJunctionLink::new(pair[0], pair[1]))
        .collect()
}

/// Returns whether every neighboring junction pair in the selected path is phased.
///
/// The selected path is circular, so the implicit BSJ is adjacent to the first
/// and last internal splice edges. A single-exon circRNA has no internal
/// junction, so this function leaves it to block-level evidence checks; any
/// multi-exon path must have high-confidence BSJ/backward read-chain support
/// for all BSJ/internal boundary links and all internal neighboring-junction
/// links.
fn major_chain_has_adjacent_link_support(
    edges: &[MajorEdge],
    link_support: &HashMap<MajorJunctionLink, MajorLinkSupport>,
    link_exclusion: &HashMap<MajorJunctionLink, MajorLinkSupport>,
) -> bool {
    if edges.is_empty() {
        return true;
    }
    let mut junctions = Vec::with_capacity(edges.len() + 2);
    junctions.push(MajorJunction::Bsj);
    junctions.extend(edges.iter().copied().map(MajorJunction::Edge));
    junctions.push(MajorJunction::Bsj);
    junctions.windows(2).all(|pair| {
        major_link_is_mature_supported(
            MajorJunctionLink::new(pair[0], pair[1]),
            link_support,
            link_exclusion,
        )
    })
}

/// Returns whether one neighboring-junction link is supported and not excluded.
///
/// Positive phasing and junction-exclusive evidence are treated symmetrically
/// for mature classification: one high-confidence BSJ/backward read can support
/// a link, but one high-confidence BSJ/backward read showing contradictory
/// cyclic reuse prevents that link from certifying a mature path.
fn major_link_is_mature_supported(
    link: MajorJunctionLink,
    link_support: &HashMap<MajorJunctionLink, MajorLinkSupport>,
    link_exclusion: &HashMap<MajorJunctionLink, MajorLinkSupport>,
) -> bool {
    let supported = link_support
        .get(&link)
        .is_some_and(|support| support.bsj + support.backward >= MAJOR_MIN_MATURE_LINK_SUPPORT);
    if !supported {
        return false;
    }
    !link_exclusion
        .get(&link)
        .is_some_and(|support| support.bsj + support.backward >= MAJOR_MIN_MATURE_LINK_SUPPORT)
}

/// Selects the major isoform path for one circRNA with phase-seed-union.
fn select_major_isoform(
    circ: &CircRecord,
    support: &HashMap<MajorEdge, MajorEdgeSupport>,
    link_support: &HashMap<MajorJunctionLink, MajorLinkSupport>,
    link_exclusion: &HashMap<MajorJunctionLink, MajorLinkSupport>,
    span_support: &[MajorAlignedSpan],
    sample_id: &str,
    annotation: Option<&Annotation>,
) -> MajorIsoformRecord {
    let phase_edges = major_phase_chain(circ, support);
    let seed_edges = major_seed_chain(support);
    let union_edges = major_union_chain(seed_edges, phase_edges);
    let strand = major_isoform_strand(circ, &union_edges);
    let annotation_exons =
        major_projection_annotation_exons_for_circ(circ, annotation, span_support);
    let chain_is_phased =
        major_chain_has_adjacent_link_support(&union_edges, link_support, link_exclusion);
    let exon_build = major_exons_from_edges(
        circ.start,
        circ.end,
        strand,
        &union_edges,
        &annotation_exons,
        support,
        link_support,
        link_exclusion,
        span_support,
        chain_is_phased,
    );
    let isoform_origin = exon_build.origin().to_string();
    let mut estimate_reason = exon_build.estimate_reason();
    let exons = exon_build.exons;
    let isoform_len = exons.iter().map(|(start, end)| end - start + 1).sum();
    let path_score: f64 = union_edges
        .iter()
        .filter_map(|edge| support.get(edge))
        .map(major_edge_weight)
        .sum();
    let cov = union_edges
        .iter()
        .filter_map(|edge| support.get(edge))
        .map(major_edge_weight)
        .fold(None, |best: Option<f64>, value| {
            Some(best.map_or(value, |current| current.min(value)))
        })
        .unwrap_or_else(|| circ.junction_read_count.parse::<usize>().unwrap_or(0) as f64);
    let segment_coverage_pct = major_segment_coverage_pct(&exons, span_support);
    if isoform_origin == "estimate"
        && segment_coverage_pct < MAJOR_MIN_FASTA_SEGMENT_COVERAGE_PCT
        && major_has_unannotated_long_exon(&exons, &annotation_exons)
    {
        major_append_estimate_reason(
            &mut estimate_reason,
            "low_segment_coverage_unannotated_long_exon",
        );
    }
    let structure_hash = major_structure_hash(&circ.chr, strand, circ.start, circ.end, &exons);
    MajorIsoformRecord {
        circ_id: circ.id.clone(),
        isoform_id: format!("{}.major", circ.id),
        sample_id: sample_id.to_string(),
        chr: circ.chr.clone(),
        start: circ.start,
        end: circ.end,
        strand,
        source_gene_id: circ.gene_id.clone(),
        exons,
        cov,
        segment_coverage_pct,
        path_score,
        bsj_reads: circ.junction_read_count.parse::<usize>().unwrap_or(0),
        structure_hash,
        isoform_len,
        isoform_origin,
        estimate_reason,
    }
}

/// Builds the BSJ-phased chain and fills gaps with backward/outward evidence.
fn major_phase_chain(
    circ: &CircRecord,
    support: &HashMap<MajorEdge, MajorEdgeSupport>,
) -> Vec<MajorEdge> {
    let skeleton = choose_major_edge_chain(
        &support
            .iter()
            .filter_map(|(&edge, counts)| (counts.bsj > 0.0).then_some((edge, counts.bsj)))
            .collect::<Vec<_>>(),
    );
    let mut phase = skeleton.clone();
    for (start, end) in major_completion_intervals(circ.start, circ.end, &skeleton) {
        let completion = choose_major_edge_chain(
            &support
                .iter()
                .filter_map(|(&edge, counts)| {
                    if edge.donor_end >= start && edge.acceptor_start <= end {
                        let weight = counts.backward as f64 + counts.outward as f64 * 0.05;
                        (weight > 0.0).then_some((edge, weight))
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>(),
        );
        phase.extend(completion);
    }
    major_sorted_unique_chain(phase)
}

/// Returns circRNA intervals between BSJ-phased skeleton edges.
fn major_completion_intervals(
    circ_start: i32,
    circ_end: i32,
    skeleton: &[MajorEdge],
) -> Vec<(i32, i32)> {
    let mut intervals = Vec::with_capacity(skeleton.len() + 1);
    let mut cursor = circ_start;
    for edge in skeleton {
        if cursor <= edge.donor_end {
            intervals.push((cursor, edge.donor_end));
        }
        cursor = edge.acceptor_start;
    }
    if cursor <= circ_end {
        intervals.push((cursor, circ_end));
    }
    intervals
}

/// Builds a high-confidence seed chain from strongest compatible edges.
fn major_seed_chain(support: &HashMap<MajorEdge, MajorEdgeSupport>) -> Vec<MajorEdge> {
    let mut ranked: Vec<(MajorEdge, f64)> = support
        .iter()
        .filter_map(|(&edge, counts)| {
            let seed = counts.bsj as f64 * 3.0 + counts.backward as f64 * 2.0;
            (seed > 0.0).then_some((edge, seed))
        })
        .collect();
    ranked.sort_by(|(edge_a, score_a), (edge_b, score_b)| {
        score_b
            .partial_cmp(score_a)
            .unwrap_or(Ordering::Equal)
            .then_with(|| edge_a.cmp(edge_b))
    });
    let Some((_, max_score)) = ranked.first().copied() else {
        return Vec::new();
    };
    let min_score = max_score * 0.05;
    let mut chain = Vec::new();
    for (edge, score) in ranked {
        if score < min_score {
            continue;
        }
        if chain
            .iter()
            .all(|selected| major_edges_compatible(*selected, edge))
        {
            chain.push(edge);
        }
    }
    major_sorted_unique_chain(chain)
}

/// Merges the seed and phase paths, keeping the seed path as the primary path.
fn major_union_chain(mut seed: Vec<MajorEdge>, phase: Vec<MajorEdge>) -> Vec<MajorEdge> {
    for edge in phase {
        if seed
            .iter()
            .all(|selected| major_edges_compatible(*selected, edge))
        {
            seed.push(edge);
        }
    }
    major_sorted_unique_chain(seed)
}

/// Chooses the highest-weight compatible edge chain by dynamic programming.
fn choose_major_edge_chain(weighted_edges: &[(MajorEdge, f64)]) -> Vec<MajorEdge> {
    if weighted_edges.is_empty() {
        return Vec::new();
    }
    let mut edges = weighted_edges.to_vec();
    edges.sort_by(|(edge_a, weight_a), (edge_b, weight_b)| {
        edge_a
            .donor_end
            .cmp(&edge_b.donor_end)
            .then_with(|| edge_a.acceptor_start.cmp(&edge_b.acceptor_start))
            .then_with(|| edge_a.strand.cmp(&edge_b.strand))
            .then_with(|| weight_b.partial_cmp(weight_a).unwrap_or(Ordering::Equal))
    });
    let n = edges.len();
    let mut best = vec![0.0; n];
    let mut prev = vec![None; n];
    for i in 0..n {
        best[i] = edges[i].1;
        for j in 0..i {
            if major_edges_compatible(edges[j].0, edges[i].0) && best[j] + edges[i].1 > best[i] {
                best[i] = best[j] + edges[i].1;
                prev[i] = Some(j);
            }
        }
    }
    let mut best_idx = 0;
    for i in 1..n {
        if best[i] > best[best_idx] {
            best_idx = i;
        }
    }
    let mut chain = Vec::new();
    let mut current = Some(best_idx);
    while let Some(idx) = current {
        chain.push(edges[idx].0);
        current = prev[idx];
    }
    chain.reverse();
    major_sorted_unique_chain(chain)
}

/// Tests whether two internal splice edges can coexist in one transcript path.
fn major_edges_compatible(a: MajorEdge, b: MajorEdge) -> bool {
    if a == b {
        return true;
    }
    if a.strand != b.strand {
        return false;
    }
    let (left, right) = if a.donor_end <= b.donor_end {
        (a, b)
    } else {
        (b, a)
    };
    left.acceptor_start <= right.donor_end
}

/// Sorts a chain in genomic order and removes duplicate edges.
fn major_sorted_unique_chain(mut chain: Vec<MajorEdge>) -> Vec<MajorEdge> {
    chain.sort_by_key(|edge| (edge.donor_end, edge.acceptor_start, edge.strand));
    chain.dedup();
    chain
}

/// Returns annotation exon intervals overlapping a circRNA span.
///
/// The segment graph defines which introns are skipped, but mature RNA length
/// must not count every genomic base between two supported splice edges. When
/// annotation is available, each graph block is projected through one
/// transcript-consistent exon chain from the Summary-assigned gene. Using the
/// gene-level exon union here can combine mutually exclusive transcript exons
/// and create unsupported isoform structures, so the union is only a fallback
/// for annotations that lack transcript IDs.
fn major_annotation_exons_for_circ(
    circ: &CircRecord,
    annotation: Option<&Annotation>,
) -> Vec<(i32, i32)> {
    let Some(annotation) = annotation else {
        return Vec::new();
    };
    let mut best_transcript: Option<((i32, i32, i32, i32), Vec<(i32, i32)>)> = None;
    for gene in circ.gene_id.split(',') {
        let Some(transcripts) = annotation.gene_transcript_exon_map.get(gene) else {
            continue;
        };
        for transcript_exons in transcripts {
            let Some(clipped) = major_clip_annotation_chain(transcript_exons, circ.start, circ.end)
            else {
                continue;
            };
            let rank = major_annotation_transcript_rank(transcript_exons, &clipped, circ);
            if best_transcript
                .as_ref()
                .is_none_or(|(best_rank, _)| rank > *best_rank)
            {
                best_transcript = Some((rank, clipped));
            }
        }
    }
    if let Some((_rank, exons)) = best_transcript {
        return exons;
    }
    let mut exons: Vec<(i32, i32)> = circ
        .gene_id
        .split(',')
        .filter_map(|gene| annotation.gene_exon_map.get(gene))
        .flat_map(|items| items.iter().copied())
        .filter_map(|(start, end)| {
            let clipped_start = start.max(circ.start);
            let clipped_end = end.min(circ.end);
            (clipped_start <= clipped_end).then_some((clipped_start, clipped_end))
        })
        .collect();
    exons.sort_unstable();
    exons.dedup();
    exons
}

/// Chooses annotation exons used to project unphased major-isoform blocks.
///
/// Transcript-consistent annotation remains the default because it avoids
/// impossible exon combinations. When high-confidence BSJ/backward spans show
/// that one transcript chain fails to cover observed mate-level anchors, the
/// projection can switch to a gene-level hybrid chain that better explains the
/// anchored blocks. This keeps the result as an estimate while avoiding the
/// worse fallback of filling the whole circRNA span as one exon.
fn major_projection_annotation_exons_for_circ(
    circ: &CircRecord,
    annotation: Option<&Annotation>,
    span_support: &[MajorAlignedSpan],
) -> Vec<(i32, i32)> {
    let primary = major_annotation_exons_for_circ(circ, annotation);
    let Some(annotation) = annotation else {
        return primary;
    };
    let anchors = major_estimate_anchor_blocks(circ.start, circ.end, span_support);
    if anchors.is_empty() {
        return primary;
    }
    let gene_exons = major_gene_annotation_exons_for_circ(circ, annotation);
    if gene_exons.is_empty() {
        return primary;
    }
    let hybrid = major_anchor_guided_annotation_exons(circ, &gene_exons, &anchors);
    if hybrid.is_empty() {
        return primary;
    }
    let primary_score = major_annotation_anchor_score(&primary, &anchors, circ);
    let hybrid_score = major_annotation_anchor_score(&hybrid, &anchors, circ);
    if hybrid_score > primary_score {
        hybrid
    } else {
        primary
    }
}

/// Returns the gene-level exon union clipped to one circRNA span.
///
/// This is intentionally separate from `major_annotation_exons_for_circ`:
/// full-length mature calls prefer transcript consistency, while estimate
/// projection sometimes needs a retained-intron block from one transcript and a
/// terminal exon from another to explain observed read anchors.
fn major_gene_annotation_exons_for_circ(
    circ: &CircRecord,
    annotation: &Annotation,
) -> Vec<(i32, i32)> {
    let mut exons: Vec<(i32, i32)> = circ
        .gene_id
        .split(',')
        .filter_map(|gene| annotation.gene_exon_map.get(gene))
        .flat_map(|items| items.iter().copied())
        .filter_map(|(start, end)| {
            let clipped_start = start.max(circ.start);
            let clipped_end = end.min(circ.end);
            (clipped_start <= clipped_end).then_some((clipped_start, clipped_end))
        })
        .collect();
    exons.sort_unstable();
    exons.dedup();
    exons
}

/// Merges BSJ/backward continuous spans into estimate anchors.
///
/// These anchors are not used to certify maturity. They only tell annotation
/// projection which parts of a long unphased block have direct read support, so
/// the fallback can choose retained-intron or terminal exons that actually cover
/// the observed alignments.
fn major_estimate_anchor_blocks(
    circ_start: i32,
    circ_end: i32,
    span_support: &[MajorAlignedSpan],
) -> Vec<(i32, i32)> {
    let mut spans: Vec<(i32, i32)> = span_support
        .iter()
        .filter(|span| span.bsj + span.backward >= MAJOR_MIN_MATURE_LINK_SUPPORT)
        .filter_map(|span| {
            let start = span.start.max(circ_start);
            let end = span.end.min(circ_end);
            (start <= end).then_some((start, end))
        })
        .collect();
    spans.sort_unstable();
    let mut anchors: Vec<(i32, i32)> = Vec::new();
    for (start, end) in spans {
        if let Some(last) = anchors.last_mut() {
            if start <= last.1 + MAJOR_ESTIMATE_ANCHOR_MERGE_GAP + 1 {
                last.1 = last.1.max(end);
                continue;
            }
        }
        anchors.push((start, end));
    }
    anchors
}

/// Selects the best gene-level annotation exon for every read-supported anchor.
///
/// The rank strongly prefers exons that contain the whole anchor, then the
/// largest overlap, then terminal-boundary agreement and compact exon length.
/// This favors retained-intron exons when mates extend past a canonical exon
/// boundary, while still choosing the shortest common terminal exon for a BSJ
/// side anchor.
fn major_anchor_guided_annotation_exons(
    circ: &CircRecord,
    gene_exons: &[(i32, i32)],
    anchors: &[(i32, i32)],
) -> Vec<(i32, i32)> {
    let mut selected = Vec::new();
    for &(anchor_start, anchor_end) in anchors {
        let best = gene_exons
            .iter()
            .copied()
            .filter(|&(exon_start, exon_end)| exon_start <= anchor_end && exon_end >= anchor_start)
            .max_by_key(|&(exon_start, exon_end)| {
                let overlap_start = exon_start.max(anchor_start);
                let overlap_end = exon_end.min(anchor_end);
                let overlap = overlap_end - overlap_start + 1;
                let contains = (exon_start <= anchor_start && exon_end >= anchor_end) as i32;
                let terminal_match =
                    (exon_start == circ.start) as i32 + (exon_end == circ.end) as i32;
                let exon_len = exon_end - exon_start + 1;
                let anchor_delta =
                    (exon_start - anchor_start).abs() + (exon_end - anchor_end).abs();
                (contains, overlap, terminal_match, -exon_len, -anchor_delta)
            });
        if let Some(exon) = best {
            selected.push(exon);
        }
    }
    major_normalize_annotation_blocks(selected)
}

/// Scores how well an annotation projection explains read-supported anchors.
fn major_annotation_anchor_score(
    exons: &[(i32, i32)],
    anchors: &[(i32, i32)],
    circ: &CircRecord,
) -> (i32, i32, i32, i32, i32) {
    let mut contained = 0;
    let mut covered_bases = 0;
    for &(anchor_start, anchor_end) in anchors {
        if exons
            .iter()
            .any(|&(exon_start, exon_end)| exon_start <= anchor_start && exon_end >= anchor_end)
        {
            contained += 1;
        }
        covered_bases += major_interval_union_overlap(exons, anchor_start, anchor_end);
    }
    let terminal_matches = exons.iter().any(|&(start, _end)| start == circ.start) as i32
        + exons.iter().any(|&(_start, end)| end == circ.end) as i32;
    let total_len: i32 = exons.iter().map(|(start, end)| end - start + 1).sum();
    (
        contained,
        covered_bases,
        terminal_matches,
        -(exons.len() as i32),
        -total_len,
    )
}

/// Returns covered bases between an exon set and one anchor interval.
fn major_interval_union_overlap(exons: &[(i32, i32)], anchor_start: i32, anchor_end: i32) -> i32 {
    let mut pieces: Vec<(i32, i32)> = exons
        .iter()
        .filter_map(|&(start, end)| {
            let clipped_start = start.max(anchor_start);
            let clipped_end = end.min(anchor_end);
            (clipped_start <= clipped_end).then_some((clipped_start, clipped_end))
        })
        .collect();
    pieces.sort_unstable();
    let mut total = 0;
    let mut current: Option<(i32, i32)> = None;
    for (start, end) in pieces {
        if let Some((current_start, current_end)) = current {
            if start <= current_end + 1 {
                current = Some((current_start, current_end.max(end)));
            } else {
                total += current_end - current_start + 1;
                current = Some((start, end));
            }
        } else {
            current = Some((start, end));
        }
    }
    if let Some((start, end)) = current {
        total += end - start + 1;
    }
    total
}

/// Sorts, deduplicates and removes overlaps from hybrid annotation blocks.
///
/// Hybrid estimate projection may pick exons from different transcripts. The
/// GTF sidecar must still emit a non-overlapping exon chain, so overlapping
/// selected blocks are conservatively merged instead of emitted as conflicting
/// features.
fn major_normalize_annotation_blocks(mut exons: Vec<(i32, i32)>) -> Vec<(i32, i32)> {
    exons.sort_unstable();
    exons.dedup();
    let mut normalized: Vec<(i32, i32)> = Vec::with_capacity(exons.len());
    for (start, end) in exons {
        if let Some(last) = normalized.last_mut() {
            if start <= last.1 {
                last.1 = last.1.max(end);
                continue;
            }
        }
        normalized.push((start, end));
    }
    normalized
}

/// Clips one transcript exon chain to a circRNA span.
fn major_clip_annotation_chain(
    exons: &[(i32, i32)],
    circ_start: i32,
    circ_end: i32,
) -> Option<Vec<(i32, i32)>> {
    let mut clipped = Vec::new();
    for &(start, end) in exons {
        let clipped_start = start.max(circ_start);
        let clipped_end = end.min(circ_end);
        if clipped_start <= clipped_end {
            clipped.push((clipped_start, clipped_end));
        }
    }
    if clipped.is_empty() {
        None
    } else {
        clipped.sort_unstable();
        clipped.dedup();
        Some(clipped)
    }
}

/// Ranks transcript candidates for annotation fallback projection.
///
/// Exact terminal exon-boundary agreement with the circRNA span is preferred.
/// This avoids selecting a transcript merely because a long terminal exon
/// overlaps the circ boundary after clipping, which would erase the intended
/// transcript-specific chain.
fn major_annotation_transcript_rank(
    transcript_exons: &[(i32, i32)],
    clipped: &[(i32, i32)],
    circ: &CircRecord,
) -> (i32, i32, i32, i32) {
    let boundary_matches = transcript_exons
        .iter()
        .any(|(start, _end)| *start == circ.start) as i32
        + transcript_exons
            .iter()
            .any(|(_start, end)| *end == circ.end) as i32;
    let terminal_delta = transcript_exons
        .iter()
        .map(|(start, _end)| (*start - circ.start).abs())
        .min()
        .unwrap_or(i32::MAX / 4)
        + transcript_exons
            .iter()
            .map(|(_start, end)| (*end - circ.end).abs())
            .min()
            .unwrap_or(i32::MAX / 4);
    let clipped_len = clipped.iter().map(|(start, end)| end - start + 1).sum();
    (
        boundary_matches,
        -terminal_delta,
        clipped_len,
        -(clipped.len() as i32),
    )
}

/// Converts an internal edge chain into 1-based inclusive exon intervals.
fn major_exons_from_edges(
    circ_start: i32,
    circ_end: i32,
    strand: char,
    edges: &[MajorEdge],
    annotation_exons: &[(i32, i32)],
    support: &HashMap<MajorEdge, MajorEdgeSupport>,
    link_support: &HashMap<MajorJunctionLink, MajorLinkSupport>,
    link_exclusion: &HashMap<MajorJunctionLink, MajorLinkSupport>,
    span_support: &[MajorAlignedSpan],
    chain_is_phased: bool,
) -> MajorExonBuild {
    let mut exons = Vec::with_capacity(edges.len() + 1);
    let mut projected_long_block = false;
    let mut unresolved_long_block = false;
    let mut inferred_internal_block = false;
    let mut unphased_single_exon_block = false;
    let mut cursor = circ_start;
    let block_phasing = major_selected_block_phasing(edges, link_support, link_exclusion);
    let is_single_exon_path = edges.is_empty();
    for (idx, edge) in edges.iter().enumerate() {
        if cursor <= edge.donor_end {
            let status = push_major_exon_block(
                &mut exons,
                cursor,
                edge.donor_end,
                strand,
                block_phasing.get(idx).copied().unwrap_or(false),
                false,
                annotation_exons,
                support,
                span_support,
            );
            projected_long_block |= status.projected;
            unresolved_long_block |= status.unresolved;
            inferred_internal_block |= status.inferred;
            unphased_single_exon_block |= status.unphased_single_exon;
        }
        cursor = edge.acceptor_start;
    }
    if cursor <= circ_end {
        let status = push_major_exon_block(
            &mut exons,
            cursor,
            circ_end,
            strand,
            block_phasing.last().copied().unwrap_or(false),
            is_single_exon_path,
            annotation_exons,
            support,
            span_support,
        );
        projected_long_block |= status.projected;
        unresolved_long_block |= status.unresolved;
        inferred_internal_block |= status.inferred;
        unphased_single_exon_block |= status.unphased_single_exon;
    }
    if exons.is_empty() && circ_start <= circ_end {
        exons.push((circ_start, circ_end));
    }
    exons.sort_unstable();
    exons.dedup();
    MajorExonBuild {
        exons,
        projected_long_block,
        unresolved_long_block,
        unphased_junction_chain: !chain_is_phased,
        inferred_internal_block,
        unphased_single_exon_block,
    }
}

/// Returns per-block phasing status for the selected circular junction path.
fn major_selected_block_phasing(
    edges: &[MajorEdge],
    link_support: &HashMap<MajorJunctionLink, MajorLinkSupport>,
    link_exclusion: &HashMap<MajorJunctionLink, MajorLinkSupport>,
) -> Vec<bool> {
    if edges.is_empty() {
        return vec![false];
    }
    let mut junctions = Vec::with_capacity(edges.len() + 2);
    junctions.push(MajorJunction::Bsj);
    junctions.extend(edges.iter().copied().map(MajorJunction::Edge));
    junctions.push(MajorJunction::Bsj);
    junctions
        .windows(2)
        .map(|pair| {
            major_link_is_mature_supported(
                MajorJunctionLink::new(pair[0], pair[1]),
                link_support,
                link_exclusion,
            )
        })
        .collect()
}

/// Appends one transcript block after resolving unphased internal structure.
fn push_major_exon_block(
    out: &mut Vec<(i32, i32)>,
    block_start: i32,
    block_end: i32,
    strand: char,
    block_is_phased: bool,
    is_single_exon_path: bool,
    annotation_exons: &[(i32, i32)],
    support: &HashMap<MajorEdge, MajorEdgeSupport>,
    span_support: &[MajorAlignedSpan],
) -> MajorBlockBuildStatus {
    let block_len = block_end - block_start + 1;
    let annotation_edges =
        major_annotation_edges_for_block(block_start, block_end, strand, annotation_exons);
    let candidate_edges = major_infer_block_edges(
        block_start,
        block_end,
        strand,
        annotation_exons,
        support,
        span_support,
    );
    if is_single_exon_path && annotation_edges.is_empty() && candidate_edges.is_empty() {
        out.push((block_start, block_end));
        if major_single_exon_block_has_mature_support(block_start, block_end, span_support) {
            return MajorBlockBuildStatus::default();
        }
        return MajorBlockBuildStatus {
            projected: false,
            unresolved: false,
            inferred: false,
            unphased_single_exon: true,
        };
    }
    if block_is_phased
        && block_len <= MAJOR_MAX_UNSPLICED_EXON_LEN
        && candidate_edges.is_empty()
        && annotation_edges.is_empty()
    {
        out.push((block_start, block_end));
        return MajorBlockBuildStatus::default();
    }
    if !candidate_edges.is_empty() {
        return push_major_blocks_from_inferred_edges(
            out,
            block_start,
            block_end,
            &candidate_edges,
            annotation_exons,
        );
    }
    if !annotation_edges.is_empty() {
        out.push((block_start, block_end));
        return MajorBlockBuildStatus {
            projected: false,
            unresolved: block_len > MAJOR_MAX_UNSPLICED_EXON_LEN,
            inferred: false,
            unphased_single_exon: false,
        };
    }
    if annotation_exons.is_empty() {
        out.push((block_start, block_end));
        return MajorBlockBuildStatus {
            projected: false,
            unresolved: block_len > MAJOR_MAX_UNSPLICED_EXON_LEN,
            inferred: false,
            unphased_single_exon: false,
        };
    }
    let before = out.len();
    for &(exon_start, exon_end) in annotation_exons {
        if exon_end < block_start {
            continue;
        }
        if exon_start > block_end {
            break;
        }
        let start = exon_start.max(block_start);
        let end = exon_end.min(block_end);
        if start <= end {
            out.push((start, end));
        }
    }
    if out.len() == before {
        out.push((block_start, block_end));
        MajorBlockBuildStatus {
            projected: false,
            unresolved: true,
            inferred: false,
            unphased_single_exon: false,
        }
    } else {
        MajorBlockBuildStatus {
            projected: true,
            unresolved: false,
            inferred: !block_is_phased,
            unphased_single_exon: false,
        }
    }
}

/// Infers splice edges inside one unphased graph block.
///
/// Read-supported BSJ/backward/outward junctions are positive candidates.
/// Annotation-only junctions are added only when no high-confidence continuous
/// BSJ/backward alignment spans across the candidate intron, so aligned
/// junction-exclusive evidence can prevent over-splitting.
fn major_infer_block_edges(
    block_start: i32,
    block_end: i32,
    strand: char,
    annotation_exons: &[(i32, i32)],
    support: &HashMap<MajorEdge, MajorEdgeSupport>,
    span_support: &[MajorAlignedSpan],
) -> Vec<MajorEdge> {
    let mut weighted = Vec::new();
    let mut seen = HashSet::new();
    for (&edge, counts) in support {
        if edge.donor_end >= block_start
            && edge.acceptor_start <= block_end
            && edge.donor_end < edge.acceptor_start
        {
            let weight = major_edge_weight(counts);
            if weight > 0.0 && seen.insert(edge) {
                weighted.push((edge, weight));
            }
        }
    }
    for edge in major_annotation_edges_for_block(block_start, block_end, strand, annotation_exons) {
        if seen.contains(&edge) {
            continue;
        }
        if !major_edge_has_high_conf_exclusion(edge, span_support) {
            seen.insert(edge);
            weighted.push((edge, 0.01));
        }
    }
    choose_major_edge_chain(&weighted)
}

/// Returns annotation-implied splice edges inside a graph block.
fn major_annotation_edges_for_block(
    block_start: i32,
    block_end: i32,
    strand: char,
    annotation_exons: &[(i32, i32)],
) -> Vec<MajorEdge> {
    let mut clipped: Vec<(i32, i32)> = annotation_exons
        .iter()
        .filter_map(|&(start, end)| {
            let clipped_start = start.max(block_start);
            let clipped_end = end.min(block_end);
            (clipped_start <= clipped_end).then_some((clipped_start, clipped_end))
        })
        .collect();
    clipped.sort_unstable();
    clipped.dedup();
    clipped
        .windows(2)
        .filter_map(|pair| {
            let donor_end = pair[0].1;
            let acceptor_start = pair[1].0;
            (donor_end < acceptor_start).then_some(MajorEdge {
                donor_end,
                acceptor_start,
                strand,
            })
        })
        .collect()
}

/// Tests whether a single-exon circRNA block has enough evidence to be mature.
///
/// With no internal junctions, the normal neighboring-junction phasing test has
/// nothing to check. Annotation can define a plausible single-exon structure,
/// but it cannot by itself prove that read evidence excludes hidden internal
/// splicing. A single-exon call is therefore mature only when high-confidence
/// BSJ/backward continuous spans tile the complete block without gaps.
fn major_single_exon_block_has_mature_support(
    block_start: i32,
    block_end: i32,
    span_support: &[MajorAlignedSpan],
) -> bool {
    major_high_conf_spans_cover_block(block_start, block_end, span_support)
}

/// Returns whether BSJ/backward continuous spans tile the whole block.
fn major_high_conf_spans_cover_block(
    block_start: i32,
    block_end: i32,
    span_support: &[MajorAlignedSpan],
) -> bool {
    let mut spans: Vec<(i32, i32)> = span_support
        .iter()
        .filter(|span| span.bsj + span.backward >= MAJOR_MIN_MATURE_LINK_SUPPORT)
        .filter_map(|span| {
            let start = span.start.max(block_start);
            let end = span.end.min(block_end);
            (start <= end).then_some((start, end))
        })
        .collect();
    spans.sort_unstable();
    let mut cursor = block_start;
    for (start, end) in spans {
        if start > cursor {
            return false;
        }
        if end >= block_end {
            return true;
        }
        cursor = cursor.max(end + 1);
    }
    false
}

/// Returns the weighted percentage of selected exon-chain bases covered by segments.
///
/// This is a structural audit metric, not the expression-like `cov` score. BSJ
/// spans contribute full coverage, while backward/outward spans contribute their
/// circRNA assignment probability. Per-base support is capped at one so several
/// ambiguous reads cannot inflate coverage beyond the selected exon length.
fn major_segment_coverage_pct(exons: &[(i32, i32)], span_support: &[MajorAlignedSpan]) -> f64 {
    let isoform_len: i32 = exons.iter().map(|(start, end)| end - start + 1).sum();
    if isoform_len <= 0 {
        return 0.0;
    }
    let mut covered = 0.0;
    for &(exon_start, exon_end) in exons {
        let mut events: Vec<(i32, f64)> = span_support
            .iter()
            .filter_map(|span| {
                let start = span.start.max(exon_start);
                let end = span.end.min(exon_end);
                let weight = span.bsj + span.backward + span.outward;
                (start <= end && weight > 0.0).then_some([(start, weight), (end + 1, -weight)])
            })
            .flatten()
            .collect();
        events.sort_by(|a, b| a.0.cmp(&b.0));
        let mut active: f64 = 0.0;
        let mut cursor = exon_start;
        for (position, delta) in events {
            if position > cursor && active > 0.0 {
                let end = (position - 1).min(exon_end);
                if end >= cursor {
                    covered += (end - cursor + 1) as f64 * active.min(1.0);
                }
            }
            active += delta;
            if position > cursor {
                cursor = position;
            }
            if cursor > exon_end {
                break;
            }
        }
    }
    covered * 100.0 / isoform_len as f64
}

/// Appends one estimate reason while preserving stable comma-separated output.
fn major_append_estimate_reason(reasons: &mut String, reason: &str) {
    if reasons == "none" || reasons.is_empty() {
        *reasons = reason.to_string();
    } else if !reasons.split(',').any(|existing| existing == reason) {
        reasons.push(',');
        reasons.push_str(reason);
    }
}

/// Returns whether a selected isoform contains a long exon absent from GTF.
///
/// The check requires complete containment by one annotation exon. Partial
/// overlaps are not enough to trust the entire projected sequence because the
/// uncovered portion may still be an unresolved internal structure.
fn major_has_unannotated_long_exon(exons: &[(i32, i32)], annotation_exons: &[(i32, i32)]) -> bool {
    exons.iter().any(|&(start, end)| {
        end - start + 1 > MAJOR_MAX_UNSPLICED_EXON_LEN
            && !annotation_exons
                .iter()
                .any(|&(anno_start, anno_end)| anno_start <= start && anno_end >= end)
    })
}

/// Tests whether high-confidence continuous alignment excludes a splice edge.
fn major_edge_has_high_conf_exclusion(edge: MajorEdge, span_support: &[MajorAlignedSpan]) -> bool {
    span_support.iter().any(|span| {
        span.start <= edge.donor_end
            && span.end >= edge.acceptor_start
            && span.bsj + span.backward >= MAJOR_MIN_MATURE_LINK_SUPPORT
    })
}

/// Converts inferred internal edges into exon blocks.
fn push_major_blocks_from_inferred_edges(
    out: &mut Vec<(i32, i32)>,
    block_start: i32,
    block_end: i32,
    edges: &[MajorEdge],
    annotation_exons: &[(i32, i32)],
) -> MajorBlockBuildStatus {
    let before = out.len();
    let mut projected = false;
    let mut unresolved = false;
    let mut cursor = block_start;
    for edge in edges {
        if cursor <= edge.donor_end {
            let status =
                push_major_resolved_subblock(out, cursor, edge.donor_end, annotation_exons);
            projected |= status.projected;
            unresolved |= status.unresolved;
        }
        cursor = edge.acceptor_start;
    }
    if cursor <= block_end {
        let status = push_major_resolved_subblock(out, cursor, block_end, annotation_exons);
        projected |= status.projected;
        unresolved |= status.unresolved;
    }
    MajorBlockBuildStatus {
        projected: projected
            || edges.iter().any(|edge| {
                major_annotation_edges_for_block(
                    block_start,
                    block_end,
                    edge.strand,
                    annotation_exons,
                )
                .contains(edge)
            }),
        unresolved: unresolved || out.len() == before,
        inferred: true,
        unphased_single_exon: false,
    }
}

/// Appends an inferred subblock, projecting only if it is still implausibly long.
fn push_major_resolved_subblock(
    out: &mut Vec<(i32, i32)>,
    block_start: i32,
    block_end: i32,
    annotation_exons: &[(i32, i32)],
) -> MajorBlockBuildStatus {
    let block_len = block_end - block_start + 1;
    if block_len <= MAJOR_MAX_UNSPLICED_EXON_LEN {
        out.push((block_start, block_end));
        return MajorBlockBuildStatus::default();
    }
    if annotation_exons.is_empty() {
        out.push((block_start, block_end));
        return MajorBlockBuildStatus {
            projected: false,
            unresolved: true,
            inferred: false,
            unphased_single_exon: false,
        };
    }
    let before = out.len();
    for &(exon_start, exon_end) in annotation_exons {
        if exon_end < block_start {
            continue;
        }
        if exon_start > block_end {
            break;
        }
        let start = exon_start.max(block_start);
        let end = exon_end.min(block_end);
        if start <= end {
            out.push((start, end));
        }
    }
    MajorBlockBuildStatus {
        projected: out.len() > before,
        unresolved: out.len() == before,
        inferred: false,
        unphased_single_exon: false,
    }
}

/// Chooses the user-facing strand for a major isoform.
fn major_isoform_strand(circ: &CircRecord, edges: &[MajorEdge]) -> char {
    circ.strand
        .chars()
        .next()
        .filter(|strand| matches!(strand, '+' | '-'))
        .or_else(|| edges.first().map(|edge| edge.strand))
        .unwrap_or('.')
}

/// Weight used by the selected path score and coverage estimate.
fn major_edge_weight(support: &MajorEdgeSupport) -> f64 {
    support.bsj + support.backward + support.outward * 0.05
}

/// Returns a stable sample label from the output prefix.
fn major_sample_id(out_prefix: &str) -> String {
    std::path::Path::new(out_prefix)
        .file_name()
        .and_then(|name| name.to_str())
        .filter(|name| !name.is_empty())
        .unwrap_or(out_prefix)
        .to_string()
}

/// Computes a deterministic FNV-1a hash for multi-sample structure matching.
fn major_structure_hash(
    chr: &str,
    strand: char,
    start: i32,
    end: i32,
    exons: &[(i32, i32)],
) -> u64 {
    let mut hash = 0xcbf29ce484222325_u64;
    fn feed(hash: &mut u64, bytes: &[u8]) {
        for byte in bytes {
            *hash ^= u64::from(*byte);
            *hash = hash.wrapping_mul(0x100000001b3);
        }
        *hash ^= 0xff;
        *hash = hash.wrapping_mul(0x100000001b3);
    }
    feed(&mut hash, chr.as_bytes());
    feed(&mut hash, strand.to_string().as_bytes());
    feed(&mut hash, &start.to_le_bytes());
    feed(&mut hash, &end.to_le_bytes());
    for (exon_start, exon_end) in exons {
        feed(&mut hash, &exon_start.to_le_bytes());
        feed(&mut hash, &exon_end.to_le_bytes());
    }
    hash
}

/// Writes major isoforms as a GTF sidecar with coverage and structure metadata.
fn write_major_isoform_gtf(path: &str, records: &[MajorIsoformRecord]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for record in records {
        let attrs = major_gtf_attributes(record);
        writeln!(
            writer,
            "{}\tCIRI-rs\tcircRNA\t{}\t{}\t{:.3}\t{}\t.\t{}",
            record.chr, record.start, record.end, record.cov, record.strand, attrs
        )?;
        writeln!(
            writer,
            "{}\tCIRI-rs\ttranscript\t{}\t{}\t{:.3}\t{}\t.\t{}",
            record.chr, record.start, record.end, record.cov, record.strand, attrs
        )?;
        for (idx, (start, end)) in record.exons.iter().enumerate() {
            writeln!(
                writer,
                "{}\tCIRI-rs\texon\t{}\t{}\t{:.3}\t{}\t.\t{} exon_number \"{}\";",
                record.chr,
                start,
                end,
                record.cov,
                record.strand,
                attrs,
                idx + 1
            )?;
        }
    }
    writer.flush()?;
    Ok(())
}

/// Formats the stable GTF attributes shared by circRNA/transcript/exon rows.
fn major_gtf_attributes(record: &MajorIsoformRecord) -> String {
    format!(
        "gene_id \"{}\"; transcript_id \"{}\"; circRNA_id \"{}\"; isoform_id \"{}\"; sample_id \"{}\"; source_gene_id \"{}\"; isoform_rank \"1\"; isoform_class \"major\"; isoform_origin \"{}\"; estimate_reason \"{}\"; path_method \"phase_seed_union\"; cov \"{:.3}\"; segment_coverage_pct \"{:.3}\"; bsj_reads \"{}\"; path_score \"{:.3}\"; exon_count \"{}\"; isoform_len \"{}\"; structure_hash \"{:016x}\";",
        gtf_escape(&record.circ_id),
        gtf_escape(&record.isoform_id),
        gtf_escape(&record.circ_id),
        gtf_escape(&record.isoform_id),
        gtf_escape(&record.sample_id),
        gtf_escape(&record.source_gene_id),
        gtf_escape(&record.isoform_origin),
        gtf_escape(&record.estimate_reason),
        record.cov,
        record.segment_coverage_pct,
        record.bsj_reads,
        record.path_score,
        record.exons.len(),
        record.isoform_len,
        record.structure_hash
    )
}

/// Escapes double quotes in GTF attribute values.
fn gtf_escape(value: &str) -> String {
    value.replace('"', "\\\"")
}

/// Counts FASTA records emitted from the selected major isoforms.
#[derive(Debug, Clone, Copy, Default)]
struct MajorIsoformFastaSummary {
    isoforms: usize,
    circ_rnas: usize,
}

/// Writes reference-derived FASTA sequences for the selected major isoforms.
fn write_major_isoform_fasta(
    path: &str,
    records: &[MajorIsoformRecord],
    reference: &HashMap<String, String>,
) -> Result<MajorIsoformFastaSummary> {
    let mut writer = BufWriter::new(File::create(path)?);
    let mut summary = MajorIsoformFastaSummary::default();
    let mut circ_ids = HashSet::new();
    for record in records {
        if !major_isoform_should_emit_fasta(record) {
            continue;
        }
        let seq = major_isoform_sequence(record, reference)?;
        writeln!(
            writer,
            ">{} circRNA_id={} sample_id={} isoform_origin={} estimate_reason={} cov={:.3} segment_coverage_pct={:.3} bsj_reads={} path_score={:.3} exon_count={} isoform_len={} structure_hash={:016x}",
            record.isoform_id,
            record.circ_id,
            record.sample_id,
            record.isoform_origin,
            record.estimate_reason,
            record.cov,
            record.segment_coverage_pct,
            record.bsj_reads,
            record.path_score,
            record.exons.len(),
            record.isoform_len,
            record.structure_hash
        )?;
        for chunk in seq.as_bytes().chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
        summary.isoforms += 1;
        circ_ids.insert(record.circ_id.as_str());
    }
    writer.flush()?;
    summary.circ_rnas = circ_ids.len();
    Ok(summary)
}

/// Returns whether a major isoform is reliable enough for sequence FASTA output.
///
/// The GTF remains the complete audit table, including unresolved estimates.
/// FASTA is stricter because unresolved placeholders become misleading sequence
/// records. Mature isoforms are always emitted; estimates need enough direct
/// segment coverage, and unphased candidates are admitted only when weighted
/// segment coverage is high enough to make the sequence-ready subset strict.
fn major_isoform_should_emit_fasta(record: &MajorIsoformRecord) -> bool {
    if record.isoform_origin == "mature" {
        return true;
    }
    if record.isoform_origin != "estimate" {
        return false;
    }
    if record
        .estimate_reason
        .contains("low_segment_coverage_unannotated_long_exon")
    {
        return false;
    }
    if record.segment_coverage_pct >= MAJOR_MIN_CANDIDATE_SEGMENT_COVERAGE_PCT {
        return true;
    }
    if record.estimate_reason.contains("unresolved_long_block") {
        return false;
    }
    if record.estimate_reason.contains("unphased_junction_chain") {
        return major_isoform_is_trusted_unphased_candidate(record);
    }
    if record
        .estimate_reason
        .contains("unphased_single_exon_block")
    {
        return false;
    }
    if record.segment_coverage_pct < MAJOR_MIN_TRUSTED_ESTIMATE_SEGMENT_COVERAGE_PCT {
        return false;
    }
    true
}

/// Returns whether an unphased estimate is reliable enough for FASTA.
///
/// The candidate still carries `unphased_junction_chain` in the header, but the
/// coverage requirement keeps this as a high-confidence sequence set rather
/// than the full audit table.
fn major_isoform_is_trusted_unphased_candidate(record: &MajorIsoformRecord) -> bool {
    !record.estimate_reason.contains("unresolved_long_block")
        && !record
            .estimate_reason
            .contains("unphased_single_exon_block")
        && !record
            .estimate_reason
            .contains("low_segment_coverage_unannotated_long_exon")
        && record.segment_coverage_pct >= MAJOR_MIN_CANDIDATE_SEGMENT_COVERAGE_PCT
}

/// Extracts a major isoform sequence in transcript orientation.
fn major_isoform_sequence(
    record: &MajorIsoformRecord,
    reference: &HashMap<String, String>,
) -> Result<String> {
    let chr_seq = reference
        .get(&record.chr)
        .ok_or_else(|| anyhow!("missing reference sequence for {}", record.chr))?;
    let mut seq = String::with_capacity(record.isoform_len.max(0) as usize);
    for &(start, end) in &record.exons {
        if start < 1 || end < start || end as usize > chr_seq.len() {
            bail!(
                "invalid major isoform exon coordinate {}:{}-{} for {}",
                record.chr,
                start,
                end,
                record.isoform_id
            );
        }
        seq.push_str(&chr_seq[(start - 1) as usize..end as usize]);
    }
    if record.strand == '-' {
        Ok(reverse_complement(&seq))
    } else {
        Ok(seq)
    }
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

    fn test_scan_state<'a>(reference: &'a HashMap<String, String>, read_len: i32) -> ScanState<'a> {
        ScanState {
            circ_by_id: HashMap::new(),
            junction_read_to_circ: HashMap::new(),
            mate_bsj_evidence: HashMap::new(),
            reference,
            annotation: None,
            circ_spans_by_chr: HashMap::new(),
            clusters_by_chr: HashMap::new(),
            read_len,
            min_mapq: 10,
            candidates: Vec::new(),
            outward_read_ids: HashSet::new(),
            prebuilt_segment_records: Vec::new(),
            segment_groups: HashMap::new(),
        }
    }

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
        let record =
            build_backward_segment_record("read1", &records, None, &[], None, 100).unwrap();

        assert_eq!(record.type_name, "backward");
        assert_eq!(record.circ_id, "NA");
        assert_eq!(record.chrom, "chr1");
        assert_eq!(record.start, "100");
        assert_eq!(record.end, "249");
        assert_eq!(record.strand, "NA");
        assert_eq!(record.is_circular, 1);
        assert_eq!(record.r1_align_strand, "+");
        assert_eq!(record.r2_align_strand, "NA");
        assert_eq!(record.r1_segments, "200-249:?|<bsj>|100-149:?");
        assert_eq!(record.r1_cigar, "50M50B50M");
        assert_eq!(record.is_r1_bsj, 0);
    }

    #[test]
    fn outward_pair_groups_emit_linear_pair_support_without_b_marker() {
        let records = vec![
            AsAlignment {
                flag: 0x40 | 0x10,
                chr: "chr1".to_string(),
                pos: 100,
                mapq: 60,
                cigar: "1S150M21S".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x80,
                chr: "chr1".to_string(),
                pos: 200,
                mapq: 60,
                cigar: "1S50M100N50M21S".to_string(),
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
            circ_spans_by_chr: {
                let mut map = HashMap::new();
                map.insert(
                    "chr1".to_string(),
                    vec![CircSpan {
                        start: 90,
                        end: 410,
                        max_end_through: 410,
                    }],
                );
                map
            },
            clusters_by_chr: HashMap::new(),
            read_len: 100,
            min_mapq: 10,
            candidates: Vec::new(),
            outward_read_ids: HashSet::new(),
            prebuilt_segment_records: Vec::new(),
            segment_groups: HashMap::new(),
        };

        assert!(is_outward_pair_group(&records, &state));
        let record = build_outward_segment_record("read1", &records, None, 100, 10).unwrap();

        assert_eq!(record.type_name, "outward");
        assert_eq!(record.circ_id, "NA");
        assert_eq!(record.chrom, "chr1");
        assert_eq!(record.start, "100");
        assert_eq!(record.end, "399");
        assert_eq!(record.is_circular, 1);
        assert_eq!(record.is_r1_bsj, 0);
        assert_eq!(record.is_r2_bsj, 0);
        assert_eq!(record.r1_align_strand, "-");
        assert_eq!(record.r2_align_strand, "+");
        assert_eq!(record.r1_cigar, "150M1S");
        assert_eq!(record.r1_segments, "100-249:?");
        assert_eq!(record.r2_cigar, "50M100N50M1S");
        assert_eq!(record.r2_segments, "350-399:?|200-249:?");
        assert!(!record.r1_segments.contains("<bsj>"));
        assert!(!record.r2_segments.contains("<bsj>"));
        assert!(!record.r1_cigar.contains('B'));
        assert!(!record.r2_cigar.contains('B'));
    }

    #[test]
    fn outward_segments_infer_strand_from_internal_annotation() {
        let records = vec![
            AsAlignment {
                flag: 0x40 | 0x10,
                chr: "chr1".to_string(),
                pos: 100,
                mapq: 60,
                cigar: "1S150M21S".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x80,
                chr: "chr1".to_string(),
                pos: 200,
                mapq: 60,
                cigar: "1S50M100N50M21S".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
        ];
        let reference = HashMap::new();
        let mut annotation = Annotation::new();
        annotation
            .chr_exon_end_index
            .entry("chr1".to_string())
            .or_default()
            .insert(249, ("GENE1".to_string(), '+'));
        annotation
            .chr_exon_start_index
            .entry("chr1".to_string())
            .or_default()
            .insert(350, ("GENE1".to_string(), '+'));
        annotation
            .transcript_splice_index
            .entry("chr1".to_string())
            .or_default()
            .entry('+')
            .or_default()
            .insert((249, 350));
        let correction = SegmentCorrectionContext {
            reference: &reference,
            annotation: Some(&annotation),
            junction_support: None,
        };

        let record = build_outward_segment_record("read1", &records, Some(&correction), 100, 10)
            .expect("annotation-supported outward row");

        assert_eq!(record.type_name, "outward");
        assert_eq!(record.strand, "+");
        assert_eq!(record.r1_align_strand, "-");
        assert_eq!(record.r2_align_strand, "+");
        assert_eq!(record.r1_segments, "100-249:+");
        assert_eq!(record.r2_segments, "350-399:+|200-249:+");
        assert!(!record.r2_segments.contains("<bsj>"));
        assert!(!record.r2_cigar.contains('B'));
    }

    #[test]
    fn outward_pair_groups_reject_contained_primary_overlap() {
        let records = vec![
            AsAlignment {
                flag: 0x40 | 0x10,
                chr: "chr1".to_string(),
                pos: 100,
                mapq: 60,
                cigar: "150M".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x80,
                chr: "chr1".to_string(),
                pos: 100,
                mapq: 60,
                cigar: "100M".to_string(),
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
            circ_spans_by_chr: {
                let mut map = HashMap::new();
                map.insert(
                    "chr1".to_string(),
                    vec![CircSpan {
                        start: 90,
                        end: 260,
                        max_end_through: 260,
                    }],
                );
                map
            },
            clusters_by_chr: HashMap::new(),
            read_len: 100,
            min_mapq: 10,
            candidates: Vec::new(),
            outward_read_ids: HashSet::new(),
            prebuilt_segment_records: Vec::new(),
            segment_groups: HashMap::new(),
        };

        assert!(!is_outward_pair_group(&records, &state));
    }

    #[test]
    fn outward_pair_groups_reject_linear_alternative_pair() {
        let records = vec![
            AsAlignment {
                flag: 0x40 | 0x10,
                chr: "chr1".to_string(),
                pos: 100,
                mapq: 60,
                cigar: "1S100M21S".to_string(),
                seq: "A".repeat(122),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x80,
                chr: "chr1".to_string(),
                pos: 250,
                mapq: 60,
                cigar: "1S100M21S".to_string(),
                seq: "A".repeat(122),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x40 | 0x100,
                chr: "chr1".to_string(),
                pos: 1000,
                mapq: 60,
                cigar: "100M".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x80 | 0x10 | 0x100,
                chr: "chr1".to_string(),
                pos: 1120,
                mapq: 60,
                cigar: "100M".to_string(),
                seq: "A".repeat(100),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
        ];
        let reference = HashMap::new();
        let state = test_scan_state(&reference, 122);

        assert!(!is_outward_pair_group(&records, &state));
    }

    #[test]
    fn outward_pair_groups_use_main_mapq_threshold() {
        let records = vec![
            AsAlignment {
                flag: 0x40 | 0x10,
                chr: "chr1".to_string(),
                pos: 100,
                mapq: 9,
                cigar: "1S100M21S".to_string(),
                seq: "A".repeat(122),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
            AsAlignment {
                flag: 0x80,
                chr: "chr1".to_string(),
                pos: 250,
                mapq: 60,
                cigar: "1S100M21S".to_string(),
                seq: "A".repeat(122),
                from_local_clip: false,
                xa_alternatives: Vec::new(),
            },
        ];
        let reference = HashMap::new();
        let state = test_scan_state(&reference, 122);

        assert!(!is_outward_pair_group(&records, &state));
    }

    #[test]
    fn outward_geometry_accepts_unclipped_gap_without_length_cap() {
        let reverse = AsAlignment {
            flag: 0x40 | 0x10,
            chr: "chr1".to_string(),
            pos: 100,
            mapq: 60,
            cigar: "51M".to_string(),
            seq: "A".repeat(51),
            from_local_clip: false,
            xa_alternatives: Vec::new(),
        };
        let forward = AsAlignment {
            flag: 0x80,
            chr: "chr1".to_string(),
            pos: 1000,
            mapq: 60,
            cigar: "101M".to_string(),
            seq: "A".repeat(101),
            from_local_clip: false,
            xa_alternatives: Vec::new(),
        };

        assert!(has_3p_outward_pair_geometry(
            &reverse,
            (100, 150),
            &forward,
            (1000, 1100)
        ));
    }

    #[test]
    fn outward_geometry_rejects_short_unclipped_outward_length() {
        let reverse = AsAlignment {
            flag: 0x40 | 0x10,
            chr: "chr1".to_string(),
            pos: 100,
            mapq: 60,
            cigar: "51M".to_string(),
            seq: "A".repeat(51),
            from_local_clip: false,
            xa_alternatives: Vec::new(),
        };
        let forward = AsAlignment {
            flag: 0x80,
            chr: "chr1".to_string(),
            pos: 115,
            mapq: 60,
            cigar: "51M".to_string(),
            seq: "A".repeat(51),
            from_local_clip: false,
            xa_alternatives: Vec::new(),
        };

        assert!(!has_3p_outward_pair_geometry(
            &reverse,
            (100, 150),
            &forward,
            (115, 165)
        ));
    }

    #[test]
    fn outward_geometry_rejects_identical_unclipped_span() {
        let reverse = AsAlignment {
            flag: 0x40 | 0x10,
            chr: "chr1".to_string(),
            pos: 100,
            mapq: 60,
            cigar: "21M".to_string(),
            seq: "A".repeat(151),
            from_local_clip: false,
            xa_alternatives: Vec::new(),
        };
        let forward = AsAlignment {
            flag: 0x80,
            chr: "chr1".to_string(),
            pos: 100,
            mapq: 60,
            cigar: "21M".to_string(),
            seq: "A".repeat(151),
            from_local_clip: false,
            xa_alternatives: Vec::new(),
        };

        assert!(!has_3p_outward_pair_geometry(
            &reverse,
            (100, 120),
            &forward,
            (100, 120)
        ));
    }

    #[test]
    fn outward_geometry_rejects_mate_5p_clip_overlap_tail() {
        let reverse = AsAlignment {
            flag: 0x40 | 0x10,
            chr: "chr1".to_string(),
            pos: 100,
            mapq: 60,
            cigar: "15S51M30S".to_string(),
            seq: "A".repeat(96),
            from_local_clip: false,
            xa_alternatives: Vec::new(),
        };
        let forward = AsAlignment {
            flag: 0x80,
            chr: "chr1".to_string(),
            pos: 200,
            mapq: 60,
            cigar: "15S51M30S".to_string(),
            seq: "A".repeat(96),
            from_local_clip: false,
            xa_alternatives: Vec::new(),
        };

        assert!(!has_3p_outward_pair_geometry(
            &reverse,
            (100, 150),
            &forward,
            (100, 150)
        ));
    }

    #[test]
    fn outward_geometry_keeps_identical_span_with_3p_clips() {
        let reverse = AsAlignment {
            flag: 0x40 | 0x10,
            chr: "chr1".to_string(),
            pos: 100,
            mapq: 60,
            cigar: "109S21M21S".to_string(),
            seq: "A".repeat(151),
            from_local_clip: false,
            xa_alternatives: Vec::new(),
        };
        let forward = AsAlignment {
            flag: 0x80,
            chr: "chr1".to_string(),
            pos: 100,
            mapq: 60,
            cigar: "1S21M129S".to_string(),
            seq: "A".repeat(151),
            from_local_clip: false,
            xa_alternatives: Vec::new(),
        };

        assert!(has_3p_outward_pair_geometry(
            &reverse,
            (100, 120),
            &forward,
            (100, 120)
        ));
    }

    #[test]
    fn major_isoform_sequence_reverse_complements_full_chain() {
        let mut reference = HashMap::new();
        reference.insert("chrT".to_string(), "ACCGTTTAGGCC".to_string());
        let record = MajorIsoformRecord {
            circ_id: "chrT:1|12".to_string(),
            isoform_id: "chrT:1|12.major".to_string(),
            sample_id: "sample".to_string(),
            chr: "chrT".to_string(),
            start: 1,
            end: 12,
            strand: '-',
            source_gene_id: "gene".to_string(),
            exons: vec![(1, 4), (9, 12)],
            cov: 1.0,
            segment_coverage_pct: 100.0,
            path_score: 1.0,
            bsj_reads: 1,
            structure_hash: 0,
            isoform_len: 8,
            isoform_origin: "mature".to_string(),
            estimate_reason: "none".to_string(),
        };

        assert_eq!(
            major_isoform_sequence(&record, &reference).unwrap(),
            "GGCCCGGT"
        );
    }

    #[test]
    fn major_isoform_fasta_skips_obvious_unresolved_estimates() {
        let mature = MajorIsoformRecord {
            circ_id: "chrT:1|12".to_string(),
            isoform_id: "chrT:1|12.major".to_string(),
            sample_id: "sample".to_string(),
            chr: "chrT".to_string(),
            start: 1,
            end: 12,
            strand: '+',
            source_gene_id: "gene".to_string(),
            exons: vec![(1, 12)],
            cov: 1.0,
            segment_coverage_pct: 100.0,
            path_score: 1.0,
            bsj_reads: 1,
            structure_hash: 0,
            isoform_len: 12,
            isoform_origin: "mature".to_string(),
            estimate_reason: "none".to_string(),
        };
        assert!(major_isoform_should_emit_fasta(&mature));

        let mut usable_estimate = mature.clone();
        usable_estimate.isoform_origin = "estimate".to_string();
        usable_estimate.estimate_reason =
            "gtf_long_block_projection,inferred_internal_block".to_string();
        usable_estimate.exons = vec![(1, 6), (9, 12)];
        usable_estimate.isoform_len = 10;
        usable_estimate.segment_coverage_pct = 100.0;
        assert!(major_isoform_should_emit_fasta(&usable_estimate));

        let mut weakly_covered = usable_estimate.clone();
        weakly_covered.segment_coverage_pct = 49.999;
        assert!(!major_isoform_should_emit_fasta(&weakly_covered));

        let mut unphased_chain = usable_estimate.clone();
        unphased_chain.estimate_reason =
            "gtf_long_block_projection,unphased_junction_chain,inferred_internal_block".to_string();
        unphased_chain.segment_coverage_pct = 100.0;
        unphased_chain.isoform_len = 10_000;
        unphased_chain.bsj_reads = 1;
        assert!(major_isoform_should_emit_fasta(&unphased_chain));

        let mut weak_unphased_chain = unphased_chain.clone();
        weak_unphased_chain.segment_coverage_pct = 89.999;
        assert!(!major_isoform_should_emit_fasta(&weak_unphased_chain));

        let mut unresolved_unphased_chain = unphased_chain;
        unresolved_unphased_chain.estimate_reason =
            "unresolved_long_block,unphased_junction_chain".to_string();
        assert!(major_isoform_should_emit_fasta(&unresolved_unphased_chain));

        let mut unresolved = usable_estimate.clone();
        unresolved.estimate_reason = "unresolved_long_block,unphased_junction_chain".to_string();
        unresolved.segment_coverage_pct = 89.999;
        assert!(!major_isoform_should_emit_fasta(&unresolved));

        let mut unphased_single = usable_estimate.clone();
        unphased_single.exons = vec![(1, 100)];
        unphased_single.isoform_len = 100;
        unphased_single.estimate_reason = "unphased_single_exon_block".to_string();
        unphased_single.segment_coverage_pct = 89.999;
        assert!(!major_isoform_should_emit_fasta(&unphased_single));

        let mut covered_unphased_single = unphased_single.clone();
        covered_unphased_single.segment_coverage_pct = 98.0;
        assert!(major_isoform_should_emit_fasta(&covered_unphased_single));

        let mut long_annotation_guided = usable_estimate;
        long_annotation_guided.exons = vec![(1, 6000), (7000, 12050)];
        long_annotation_guided.isoform_len = 11051;
        long_annotation_guided.estimate_reason =
            "gtf_long_block_projection,inferred_internal_block".to_string();
        assert!(major_isoform_should_emit_fasta(&long_annotation_guided));

        let mut weak_unannotated_long = long_annotation_guided;
        weak_unannotated_long.estimate_reason =
            "low_segment_coverage_unannotated_long_exon".to_string();
        assert!(!major_isoform_should_emit_fasta(&weak_unannotated_long));
    }

    #[test]
    fn major_isoforms_are_rebuilt_from_segments_file() {
        let tmp_prefix = std::env::temp_dir().join(format!(
            "ciri_major_isoform_segments_{}",
            std::process::id()
        ));
        let segments_path = tmp_prefix.with_extension("segments");
        let out_prefix = tmp_prefix.to_string_lossy().to_string();
        std::fs::write(
            &segments_path,
            concat!(
                "read_id\ttype\tcirc_id\tchrom\tstart\tend\tstrand\tis_circular\tis_r1_bsj\tis_r2_bsj\tr1_align_strand\tr2_align_strand\tr1_cigar\tr1_segments\tr2_cigar\tr2_segments\n",
                "read1\tbsj\tchrT:100|300\tchrT\t100\t300\t+\t1\t1\t0\t+\tNA\t21M49N21M\t280-300:+|<bsj>|100-150:+|200-220:+\tNA\tNA\n",
                "read2\tbsj\tchrT:100|300\tchrT\t100\t300\t+\t1\t1\t0\t+\tNA\t31M49N101M\t120-150:+|200-300:+|<bsj>|100-110:+\tNA\tNA\n",
                "read3\tbsj\tchrT:100|300\tchrT\t100\t300\t+\t1\t1\t0\t+\tNA\t21M49N21M\t280-300:+|<bsj>|100-150:+|200-220:+\tNA\tNA\n",
            ),
        )
        .unwrap();
        let circ = CircRecord {
            id: "chrT:100|300".to_string(),
            chr: "chrT".to_string(),
            start: 100,
            end: 300,
            junction_read_count: "3".to_string(),
            gene_id: "geneT".to_string(),
            strand: "+".to_string(),
        };
        let mut reference = HashMap::new();
        reference.insert("chrT".to_string(), "ACGT".repeat(100));

        let summary = build_major_isoforms_from_segments_file(
            &[circ],
            segments_path.to_str().unwrap(),
            &out_prefix,
            &reference,
            None,
        )
        .unwrap();

        assert_eq!(summary.total_isoforms, 1);
        assert_eq!(summary.circ_rnas, 1);
        assert_eq!(summary.fasta_isoforms, 1);
        assert_eq!(summary.fasta_circ_rnas, 1);
        let gtf = std::fs::read_to_string(format!("{}.isoforms.gtf", out_prefix)).unwrap();
        assert!(gtf.contains("\texon\t100\t150\t"));
        assert!(gtf.contains("\texon\t200\t300\t"));
        assert!(gtf.contains("isoform_origin \"mature\";"));
        assert!(gtf.contains("estimate_reason \"none\";"));
        assert!(gtf.contains("segment_coverage_pct \"100.000\";"));
        let fasta = std::fs::read_to_string(format!("{}.isoforms.fa", out_prefix)).unwrap();
        assert!(fasta.contains("isoform_len=152"));
        assert!(fasta.contains("isoform_origin=mature"));
        assert!(fasta.contains("estimate_reason=none"));
        assert!(fasta.contains("segment_coverage_pct=100.000"));
        let _ = std::fs::remove_file(&segments_path);
        let _ = std::fs::remove_file(format!("{}.isoforms.gtf", out_prefix));
        let _ = std::fs::remove_file(format!("{}.isoforms.fa", out_prefix));
    }

    #[test]
    fn major_segment_coverage_pct_reports_selected_exon_union() {
        let exons = vec![(100, 199), (300, 399)];
        let spans = vec![
            MajorAlignedSpan {
                start: 120,
                end: 160,
                bsj: 1.0,
                backward: 0.0,
                outward: 0.0,
            },
            MajorAlignedSpan {
                start: 150,
                end: 199,
                bsj: 0.0,
                backward: 0.0,
                outward: 1.0,
            },
            MajorAlignedSpan {
                start: 350,
                end: 399,
                bsj: 0.0,
                backward: 1.0,
                outward: 0.0,
            },
        ];

        assert!((major_segment_coverage_pct(&exons, &spans) - 65.0).abs() < 1e-6);
    }

    #[test]
    fn major_segment_coverage_pct_respects_fractional_assignment() {
        let exons = vec![(100, 199)];
        let spans = vec![MajorAlignedSpan {
            start: 100,
            end: 199,
            bsj: 0.0,
            backward: 0.25,
            outward: 0.0,
        }];

        assert!((major_segment_coverage_pct(&exons, &spans) - 25.0).abs() < 1e-6);
    }

    #[test]
    fn major_non_bsj_assignment_leaves_internal_reads_mostly_unassigned() {
        let circ_records = vec![CircRecord {
            id: "chrT:100|1000".to_string(),
            chr: "chrT".to_string(),
            start: 100,
            end: 1000,
            junction_read_count: "1".to_string(),
            gene_id: "geneT".to_string(),
            strand: "+".to_string(),
        }];
        let circ_by_id: HashMap<&str, usize> = circ_records
            .iter()
            .enumerate()
            .map(|(idx, circ)| (circ.id.as_str(), idx))
            .collect();
        let circ_index_by_chr = build_major_circ_index(&circ_records);

        let boundary = assign_major_circs_from_segments(
            "backward",
            "NA",
            "chrT",
            "100",
            "1000",
            "100-150:+|950-1000:+",
            "NA",
            &circ_by_id,
            &circ_index_by_chr,
            &circ_records,
        );
        assert_eq!(boundary.len(), 1);
        assert!((boundary[0].weight - 1.0).abs() < 1e-6);

        let internal = assign_major_circs_from_segments(
            "backward",
            "NA",
            "chrT",
            "450",
            "550",
            "450-500:+|520-550:+",
            "NA",
            &circ_by_id,
            &circ_index_by_chr,
            &circ_records,
        );
        assert_eq!(internal.len(), 1);
        assert!(internal[0].weight < 0.1);
    }

    #[test]
    fn major_unannotated_long_exon_filter_requires_complete_annotation_support() {
        let exons = vec![(100, 2500), (3000, 3100)];
        assert!(major_has_unannotated_long_exon(
            &exons,
            &[(150, 2500), (3000, 3100)]
        ));
        assert!(!major_has_unannotated_long_exon(
            &exons,
            &[(100, 2500), (3000, 3100)]
        ));

        let mut reason = "gtf_long_block_projection".to_string();
        major_append_estimate_reason(&mut reason, "low_segment_coverage_unannotated_long_exon");
        assert_eq!(
            reason,
            "gtf_long_block_projection,low_segment_coverage_unannotated_long_exon"
        );
    }

    #[test]
    fn major_isoform_blocks_project_to_annotation_exons() {
        let circ = CircRecord {
            id: "chrT:100|20500".to_string(),
            chr: "chrT".to_string(),
            start: 100,
            end: 20500,
            junction_read_count: "3".to_string(),
            gene_id: "geneT".to_string(),
            strand: "+".to_string(),
        };
        let mut annotation = Annotation::new();
        annotation.gene_exon_map.insert(
            "geneT".to_string(),
            vec![(80, 150), (220, 260), (20450, 20520)],
        );

        let annotation_exons = major_annotation_exons_for_circ(&circ, Some(&annotation));
        assert_eq!(
            annotation_exons,
            vec![(100, 150), (220, 260), (20450, 20500)]
        );
        assert_eq!(
            major_exons_from_edges(
                circ.start,
                circ.end,
                '+',
                &[],
                &annotation_exons,
                &HashMap::new(),
                &HashMap::new(),
                &HashMap::new(),
                &[],
                true,
            )
            .exons,
            vec![(100, 150), (220, 260), (20450, 20500)]
        );
        let build = major_exons_from_edges(
            circ.start,
            circ.end,
            '+',
            &[],
            &annotation_exons,
            &HashMap::new(),
            &HashMap::new(),
            &HashMap::new(),
            &[],
            true,
        );
        assert_eq!(build.origin(), "estimate");
        assert_eq!(
            build.estimate_reason(),
            "gtf_long_block_projection,inferred_internal_block"
        );
    }

    #[test]
    fn major_isoform_prefers_transcript_consistent_annotation_chain() {
        let circ = CircRecord {
            id: "chrT:120|420".to_string(),
            chr: "chrT".to_string(),
            start: 120,
            end: 420,
            junction_read_count: "3".to_string(),
            gene_id: "geneT".to_string(),
            strand: "+".to_string(),
        };
        let mut annotation = Annotation::new();
        annotation.gene_exon_map.insert(
            "geneT".to_string(),
            vec![(100, 180), (120, 150), (220, 250), (300, 500), (390, 420)],
        );
        annotation.gene_transcript_exon_map.insert(
            "geneT".to_string(),
            vec![
                vec![(100, 180), (300, 500)],
                vec![(120, 150), (220, 250), (390, 420)],
            ],
        );

        assert_eq!(
            major_annotation_exons_for_circ(&circ, Some(&annotation)),
            vec![(120, 150), (220, 250), (390, 420)]
        );
    }

    #[test]
    fn major_isoform_uses_read_anchored_hybrid_annotation_estimate() {
        let circ = CircRecord {
            id: "chrT:100|1534".to_string(),
            chr: "chrT".to_string(),
            start: 100,
            end: 1534,
            junction_read_count: "2".to_string(),
            gene_id: "geneT".to_string(),
            strand: "+".to_string(),
        };
        let mut annotation = Annotation::new();
        annotation.gene_exon_map.insert(
            "geneT".to_string(),
            vec![(100, 267), (100, 1397), (1406, 1534), (1477, 1534)],
        );
        annotation.gene_transcript_exon_map.insert(
            "geneT".to_string(),
            vec![vec![(100, 267), (1406, 1534)], vec![(100, 1397)]],
        );
        let spans = vec![
            MajorAlignedSpan {
                start: 100,
                end: 199,
                bsj: 1.0,
                backward: 0.0,
                outward: 0.0,
            },
            MajorAlignedSpan {
                start: 248,
                end: 415,
                bsj: 1.0,
                backward: 0.0,
                outward: 0.0,
            },
            MajorAlignedSpan {
                start: 1465,
                end: 1534,
                bsj: 1.0,
                backward: 0.0,
                outward: 0.0,
            },
        ];

        let annotation_exons =
            major_projection_annotation_exons_for_circ(&circ, Some(&annotation), &spans);
        assert_eq!(annotation_exons, vec![(100, 1397), (1406, 1534)]);

        let build = major_exons_from_edges(
            circ.start,
            circ.end,
            '+',
            &[],
            &annotation_exons,
            &HashMap::new(),
            &HashMap::new(),
            &HashMap::new(),
            &spans,
            true,
        );
        assert_eq!(build.exons, vec![(100, 1397), (1406, 1534)]);
        assert_eq!(build.origin(), "estimate");
        assert_eq!(
            build.estimate_reason(),
            "gtf_long_block_projection,inferred_internal_block"
        );
    }

    #[test]
    fn major_isoform_infers_unphased_short_blocks_from_annotation() {
        let annotation_exons = vec![(100, 150), (220, 260), (450, 500)];
        let build = major_exons_from_edges(
            100,
            500,
            '+',
            &[],
            &annotation_exons,
            &HashMap::new(),
            &HashMap::new(),
            &HashMap::new(),
            &[],
            true,
        );
        assert_eq!(build.exons, vec![(100, 150), (220, 260), (450, 500)]);
        assert_eq!(build.origin(), "estimate");
        assert_eq!(
            build.estimate_reason(),
            "gtf_long_block_projection,inferred_internal_block"
        );
    }

    #[test]
    fn major_isoform_uses_exclusion_to_avoid_annotation_over_split() {
        let annotation_exons = vec![(100, 150), (220, 260), (450, 500)];
        let spans = vec![MajorAlignedSpan {
            start: 100,
            end: 500,
            bsj: 1.0,
            backward: 0.0,
            outward: 0.0,
        }];
        let build = major_exons_from_edges(
            100,
            500,
            '+',
            &[],
            &annotation_exons,
            &HashMap::new(),
            &HashMap::new(),
            &HashMap::new(),
            &spans,
            true,
        );
        assert_eq!(build.exons, vec![(100, 500)]);
        assert_eq!(build.origin(), "mature");
        assert_eq!(build.estimate_reason(), "none");
    }

    #[test]
    fn major_isoform_marks_uncovered_single_exon_as_estimate() {
        let spans = vec![
            MajorAlignedSpan {
                start: 100,
                end: 160,
                bsj: 1.0,
                backward: 0.0,
                outward: 0.0,
            },
            MajorAlignedSpan {
                start: 340,
                end: 500,
                bsj: 1.0,
                backward: 0.0,
                outward: 0.0,
            },
        ];
        let build = major_exons_from_edges(
            100,
            500,
            '+',
            &[],
            &[],
            &HashMap::new(),
            &HashMap::new(),
            &HashMap::new(),
            &spans,
            true,
        );
        assert_eq!(build.exons, vec![(100, 500)]);
        assert_eq!(build.origin(), "estimate");
        assert_eq!(build.estimate_reason(), "unphased_single_exon_block");
    }

    #[test]
    fn major_isoform_keeps_covered_single_exon_mature() {
        let spans = vec![
            MajorAlignedSpan {
                start: 100,
                end: 260,
                bsj: 1.0,
                backward: 0.0,
                outward: 0.0,
            },
            MajorAlignedSpan {
                start: 261,
                end: 500,
                bsj: 0.0,
                backward: 1.0,
                outward: 0.0,
            },
        ];
        let build = major_exons_from_edges(
            100,
            500,
            '+',
            &[],
            &[],
            &HashMap::new(),
            &HashMap::new(),
            &HashMap::new(),
            &spans,
            true,
        );
        assert_eq!(build.exons, vec![(100, 500)]);
        assert_eq!(build.origin(), "mature");
        assert_eq!(build.estimate_reason(), "none");
    }

    #[test]
    fn major_isoform_marks_annotation_only_single_exon_as_estimate() {
        let build = major_exons_from_edges(
            100,
            500,
            '+',
            &[],
            &[(100, 500)],
            &HashMap::new(),
            &HashMap::new(),
            &HashMap::new(),
            &[],
            true,
        );
        assert_eq!(build.exons, vec![(100, 500)]);
        assert_eq!(build.origin(), "estimate");
        assert_eq!(build.estimate_reason(), "unphased_single_exon_block");
    }

    #[test]
    fn major_isoform_rejects_clipped_annotation_single_exon_as_mature() {
        let build = major_exons_from_edges(
            100,
            400,
            '+',
            &[],
            &[(100, 400)],
            &HashMap::new(),
            &HashMap::new(),
            &HashMap::new(),
            &[],
            true,
        );
        assert_eq!(build.exons, vec![(100, 400)]);
        assert_eq!(build.origin(), "estimate");
        assert_eq!(build.estimate_reason(), "unphased_single_exon_block");
    }

    #[test]
    fn major_isoform_uses_read_supported_junction_inside_unphased_block() {
        let edge = MajorEdge {
            donor_end: 150,
            acceptor_start: 450,
            strand: '+',
        };
        let mut support = HashMap::new();
        support.insert(
            edge,
            MajorEdgeSupport {
                bsj: 1.0,
                backward: 0.0,
                outward: 0.0,
            },
        );
        let build = major_exons_from_edges(
            100,
            500,
            '+',
            &[],
            &[],
            &support,
            &HashMap::new(),
            &HashMap::new(),
            &[],
            false,
        );
        assert_eq!(build.exons, vec![(100, 150), (450, 500)]);
        assert_eq!(build.origin(), "estimate");
        assert_eq!(
            build.estimate_reason(),
            "unphased_junction_chain,inferred_internal_block"
        );
    }

    #[test]
    fn major_isoform_marks_unphased_junction_chain_as_estimate() {
        let edges = vec![
            MajorEdge {
                donor_end: 150,
                acceptor_start: 220,
                strand: '+',
            },
            MajorEdge {
                donor_end: 260,
                acceptor_start: 450,
                strand: '+',
            },
        ];
        let build = major_exons_from_edges(
            100,
            500,
            '+',
            &edges,
            &[],
            &HashMap::new(),
            &HashMap::new(),
            &HashMap::new(),
            &[],
            false,
        );
        assert_eq!(build.exons, vec![(100, 150), (220, 260), (450, 500)]);
        assert_eq!(build.origin(), "estimate");
        assert_eq!(build.estimate_reason(), "unphased_junction_chain");
    }

    #[test]
    fn major_isoform_requires_bsj_boundary_phasing_for_single_internal_edge() {
        let edge = MajorEdge {
            donor_end: 150,
            acceptor_start: 450,
            strand: '+',
        };
        let no_boundary_support = HashMap::new();
        let no_boundary_exclusion = HashMap::new();
        assert!(!major_chain_has_adjacent_link_support(
            &[edge],
            &no_boundary_support,
            &no_boundary_exclusion
        ));

        let mut boundary_support = HashMap::new();
        boundary_support.insert(
            MajorJunctionLink::new(MajorJunction::Bsj, MajorJunction::Edge(edge)),
            MajorLinkSupport {
                bsj: 0.0,
                backward: 0.0,
                outward: 1.0,
            },
        );
        assert!(!major_chain_has_adjacent_link_support(
            &[edge],
            &boundary_support,
            &no_boundary_exclusion
        ));

        boundary_support.insert(
            MajorJunctionLink::new(MajorJunction::Bsj, MajorJunction::Edge(edge)),
            MajorLinkSupport {
                bsj: 1.0,
                backward: 0.0,
                outward: 0.0,
            },
        );
        assert!(!major_chain_has_adjacent_link_support(
            &[edge],
            &boundary_support,
            &no_boundary_exclusion
        ));
        boundary_support.insert(
            MajorJunctionLink::new(MajorJunction::Edge(edge), MajorJunction::Bsj),
            MajorLinkSupport {
                bsj: 0.0,
                backward: 1.0,
                outward: 0.0,
            },
        );
        assert!(major_chain_has_adjacent_link_support(
            &[edge],
            &boundary_support,
            &no_boundary_exclusion
        ));
        let mut boundary_exclusion = HashMap::new();
        boundary_exclusion.insert(
            MajorJunctionLink::new(MajorJunction::Bsj, MajorJunction::Edge(edge)),
            MajorLinkSupport {
                bsj: 1.0,
                backward: 0.0,
                outward: 0.0,
            },
        );
        assert!(!major_chain_has_adjacent_link_support(
            &[edge],
            &boundary_support,
            &boundary_exclusion
        ));
    }

    #[test]
    fn major_segment_junction_chain_keeps_bsj_boundary_tokens() {
        let circ = CircRecord {
            id: "chrT:100|500".to_string(),
            chr: "chrT".to_string(),
            start: 100,
            end: 500,
            junction_read_count: "1".to_string(),
            gene_id: "geneT".to_string(),
            strand: "+".to_string(),
        };
        let junctions = major_segment_junction_chain("480-500:+|<bsj>|100-150:+|450-470:+", &circ);
        assert_eq!(
            junctions,
            vec![
                MajorJunction::Bsj,
                MajorJunction::Edge(MajorEdge {
                    donor_end: 150,
                    acceptor_start: 450,
                    strand: '+',
                })
            ]
        );
        assert_eq!(
            major_adjacent_junction_links(&junctions),
            vec![MajorJunctionLink::new(
                MajorJunction::Bsj,
                MajorJunction::Edge(MajorEdge {
                    donor_end: 150,
                    acceptor_start: 450,
                    strand: '+',
                })
            )]
        );
    }

    #[test]
    fn major_segment_reused_overlap_excludes_mature_phasing() {
        let circ = CircRecord {
            id: "chr18:29234877|29255341".to_string(),
            chr: "chr18".to_string(),
            start: 29234877,
            end: 29255341,
            junction_read_count: "1".to_string(),
            gene_id: "NA".to_string(),
            strand: "-".to_string(),
        };
        let reused = "29255287-29255344:-|<bsj>|29236698-29236725:-|29255273-29255337:-";
        assert!(major_segment_has_reused_overlap(reused, &circ));

        let edge = MajorEdge {
            donor_end: 29236725,
            acceptor_start: 29255273,
            strand: '-',
        };
        let mut support = HashMap::new();
        support.insert(
            MajorJunctionLink::new(MajorJunction::Bsj, MajorJunction::Edge(edge)),
            MajorLinkSupport {
                bsj: 1.0,
                backward: 0.0,
                outward: 0.0,
            },
        );
        support.insert(
            MajorJunctionLink::new(MajorJunction::Edge(edge), MajorJunction::Bsj),
            MajorLinkSupport {
                bsj: 1.0,
                backward: 0.0,
                outward: 0.0,
            },
        );
        let mut exclusion = HashMap::new();
        exclusion.insert(
            MajorJunctionLink::new(MajorJunction::Bsj, MajorJunction::Edge(edge)),
            MajorLinkSupport {
                bsj: 1.0,
                backward: 0.0,
                outward: 0.0,
            },
        );
        assert!(!major_chain_has_adjacent_link_support(
            &[edge],
            &support,
            &exclusion
        ));
    }

    #[test]
    fn major_isoform_filters_chimeric_mate_before_graph_support() {
        let tmp_prefix = std::env::temp_dir().join(format!(
            "ciri_major_isoform_chimeric_{}",
            std::process::id()
        ));
        let segments_path = tmp_prefix.with_extension("segments");
        let out_prefix = tmp_prefix.to_string_lossy().to_string();
        std::fs::write(
            &segments_path,
            concat!(
                "read_id\ttype\tcirc_id\tchrom\tstart\tend\tstrand\tis_circular\tis_r1_bsj\tis_r2_bsj\tr1_align_strand\tr2_align_strand\tr1_cigar\tr1_segments\tr2_cigar\tr2_segments\n",
                "read1\tbsj\tchrT:100|300\tchrT\t100\t300\t-\t1\t1\t1\t+\t-\t58M80B21M69N51M\t246-303:-|<bsj>|150-170:-|240-290:-\t21M69N61M109B31M\t150-170:-|240-300:-|<bsj>|100-130:-\n",
            ),
        )
        .unwrap();
        let circ = CircRecord {
            id: "chrT:100|300".to_string(),
            chr: "chrT".to_string(),
            start: 100,
            end: 300,
            junction_read_count: "1".to_string(),
            gene_id: "NA".to_string(),
            strand: "-".to_string(),
        };
        let mut reference = HashMap::new();
        reference.insert("chrT".to_string(), "ACGT".repeat(100));

        let summary = build_major_isoforms_from_segments_file(
            &[circ],
            segments_path.to_str().unwrap(),
            &out_prefix,
            &reference,
            None,
        )
        .unwrap();

        assert_eq!(summary.total_isoforms, 1);
        assert_eq!(summary.circ_rnas, 1);
        let gtf = std::fs::read_to_string(format!("{}.isoforms.gtf", out_prefix)).unwrap();
        assert!(gtf.contains("\texon\t100\t170\t"));
        assert!(gtf.contains("\texon\t240\t300\t"));
        assert!(gtf.contains("isoform_origin \"estimate\";"));
        assert!(gtf.contains("estimate_reason \"unphased_junction_chain\";"));
        let _ = std::fs::remove_file(&segments_path);
        let _ = std::fs::remove_file(format!("{}.isoforms.gtf", out_prefix));
        let _ = std::fs::remove_file(format!("{}.isoforms.fa", out_prefix));
    }

    #[test]
    fn non_bsj_local_clip_alignment_uses_cluster_window_beyond_circ_span() {
        let clip = "ACGTACGTACGTACGTACGT";
        let mut chr_seq = "T".repeat(500);
        chr_seq.replace_range(299..319, clip);
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), chr_seq);
        let records = vec![AsAlignment {
            flag: 0x40,
            chr: "chr1".to_string(),
            pos: 100,
            mapq: 60,
            cigar: "130M20S".to_string(),
            seq: format!("{}{}", "G".repeat(130), clip),
            from_local_clip: false,
            xa_alternatives: Vec::new(),
        }];
        let state = ScanState {
            circ_by_id: HashMap::new(),
            junction_read_to_circ: HashMap::new(),
            mate_bsj_evidence: HashMap::new(),
            reference: &reference,
            annotation: None,
            circ_spans_by_chr: {
                let mut map = HashMap::new();
                map.insert(
                    "chr1".to_string(),
                    vec![CircSpan {
                        start: 90,
                        end: 220,
                        max_end_through: 220,
                    }],
                );
                map
            },
            clusters_by_chr: {
                let mut map = HashMap::new();
                map.insert(
                    "chr1".to_string(),
                    vec![CircCluster {
                        chr: "chr1".to_string(),
                        start: 80,
                        end: 350,
                    }],
                );
                map
            },
            read_len: 150,
            min_mapq: 10,
            candidates: Vec::new(),
            outward_read_ids: HashSet::new(),
            prebuilt_segment_records: Vec::new(),
            segment_groups: HashMap::new(),
        };

        let enriched = add_non_bsj_local_clip_alignments(&records, &state);

        assert!(enriched.iter().any(|record| {
            record.from_local_clip && record.pos == 300 && record.cigar == "130S20M"
        }));
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
        assert_eq!(record.r1_segments, "5000-5049:+|<bsj>|4990-5039:+");
        assert_eq!(record.r1_cigar, "50M0B50M");
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
            "203745224-203745245:+|<bsj>|203741196-203741275:+|203743002-203743050:+"
        );
        assert_eq!(record.r2_cigar, "22M3948B80M1726N49M");
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
                r1_align_strand: "NA".to_string(),
                r2_align_strand: "NA".to_string(),
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
                r1_align_strand: "NA".to_string(),
                r2_align_strand: "NA".to_string(),
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
                r1_align_strand: "NA".to_string(),
                r2_align_strand: "NA".to_string(),
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
            junction_read_count: "1".to_string(),
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
            junction_read_count: "1".to_string(),
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
            junction_read_count: "1".to_string(),
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
            junction_read_count: "1".to_string(),
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
            junction_read_count: "1".to_string(),
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
                "<bsj>".to_string(),
                "93811273-93811340:-".to_string(),
            ]
        );
        assert_eq!(chain.cigar, "29M5134N55M5204B68M");
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
            junction_read_count: "1".to_string(),
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
}
