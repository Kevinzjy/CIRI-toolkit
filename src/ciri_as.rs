//! CIRI-AS sidecar reconstruction entry points.
//!
//! This module is intentionally separate from `Scan1 -> Scan2 -> Summary`.
//! The CIRI3 Rust path is already parity-sensitive, so CIRI-AS reconstruction is
//! implemented as a post-Summary sidecar that consumes the final circRNA table
//! and the original SAM/BAM alignments without feeding any evidence back into BSJ
//! detection. The first parity target is CIRI-AS v1.2 splice-junction discovery:
//! read grouping, `MSID` CIGAR classification, splice-signal adjustment, and the
//! `_splice.list` clustering contract are mirrored before cirexon/AS catalog
//! reconstruction is expanded.

use anyhow::{anyhow, bail, Context, Result};
use noodles::bam;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use crate::annotation::Annotation;
use crate::sam_bam::{detect_format, InputFormat};
use crate::utils::reverse_complement;

const MIN_INTRON: i32 = 70;
const MIN_EXON_LENGTH: i32 = 20;
const MAX_EXON_LENGTH: i32 = 2000;
const MAX_ISOFORM_PATHS_PER_CIRC: usize = 1024;
const Z_ALPHA: f64 = 1.6449;
const MAPQ_THRES: i32 = 5;
const MAPQ_UNI: i32 = 0;
const MAPQ_BOTH: i32 = 0;
const STRINGENCY: usize = 1;

/// Runtime options for the CIRI-AS sidecar.
///
/// The defaults intentionally follow `vendor/CIRI-AS/CIRI_AS_v1.2.pl`. Keeping
/// them local to this sidecar prevents CIRI-AS tuning from changing CIRI3 BSJ
/// calls or Summary filtering.
pub struct AsConfig<'a> {
    /// Original queryname-sorted SAM/BAM used by the main CIRI run.
    pub input_path: &'a str,
    /// Final CIRI3-style circRNA result table used as CIRI-AS `-C` input.
    pub circ_path: &'a str,
    /// Output prefix; files follow CIRI-AS names such as `<prefix>_splice.list`.
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
/// CIRI-AS only needs a subset of SAM fields. We keep sequence here because the
/// original Perl stores it while normalizing read orientation, even though the
/// current splice-junction stage does not consume the sequence after that point.
#[derive(Debug, Clone)]
struct AsAlignment {
    flag: i32,
    chr: String,
    pos: i32,
    mapq: i32,
    cigar: String,
    seq: String,
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
}

/// Runs the CIRI-AS sidecar splice-junction discovery and cirexon stages.
///
/// The function keeps all CIRI-AS evidence as post-Summary state. `_splice.list`
/// is still the first checkpoint; `.list` is then derived from the same splice
/// clusters plus read coverage, while `_AS.list` remains a header-only file until
/// exon path and PSI correction are implemented.
pub fn run_ciri_as(config: AsConfig<'_>) -> Result<()> {
    let (circ_records, junction_read_to_circ) = load_circ_records(config.circ_path)?;
    if circ_records.is_empty() {
        bail!("No circRNAs were loaded from {}", config.circ_path);
    }
    if junction_read_to_circ.is_empty() {
        bail!(
            "No circular junction read IDs were loaded from {}",
            config.circ_path
        );
    }

    let read_len = infer_read_length(config.input_path)?;
    if read_len < 40 {
        bail!(
            "CIRI-AS requires paired reads with inferred read length >= 40, got {}",
            read_len
        );
    }

    let clusters = build_circ_clusters(&circ_records);
    let mut state = ScanState {
        circ_by_id: circ_records
            .iter()
            .map(|r| (r.id.clone(), r.clone()))
            .collect(),
        junction_read_to_circ,
        clusters_by_chr: clusters_by_chr(clusters),
        read_len,
        candidates: Vec::new(),
        coverage: HashMap::new(),
        read_mappings: HashMap::new(),
        seen_junction_reads: HashSet::new(),
        stats: AsStats::default(),
    };
    state.stats.junction_reads_loaded = state.junction_read_to_circ.len();

    scan_alignment_groups(config.input_path, &mut state)?;
    validate_splice_motifs(&mut state.candidates, config.reference, config.annotation)?;
    state.stats.motif_validated = state.candidates.len();
    let splice_clusters = cluster_candidates(&state.candidates);
    state.stats.final_splice_clusters = splice_clusters.len();
    state.stats.junction_reads_seen = state.seen_junction_reads.len();
    let (cirexons, isoforms) = predict_cirexons(
        &state,
        &splice_clusters,
        config.reference,
        config.annotation,
    );
    state.stats.final_cirexons = cirexons.len();
    state.stats.final_isoforms = isoforms.len();

    write_splice_list(
        &format!("{}_splice.list", config.out_prefix),
        &splice_clusters,
        &state.candidates,
        &state.junction_read_to_circ,
    )?;
    write_cirexon_list(&format!("{}.list", config.out_prefix), &cirexons)?;
    write_as_header(&format!("{}_AS.list", config.out_prefix))?;
    write_isoforms(&format!("{}.isoforms", config.out_prefix), &isoforms)?;
    write_isoform_summary(
        &format!("{}.isoform_summary", config.out_prefix),
        &state,
        &cirexons,
        &isoforms,
    )?;
    write_isoform_fasta(
        &format!("{}.fa", config.out_prefix),
        &isoforms,
        config.reference,
    )?;
    write_as_log(&format!("{}.log", config.out_prefix), &state)?;
    Ok(())
}

/// Mutable state accumulated while scanning the original alignment file.
///
/// Perl stores these as package globals. Rust keeps them in one struct so later
/// BAM parallelization can shard the scan and merge candidate/state fragments
/// without changing the CIRI-AS decision routines themselves.
struct ScanState {
    circ_by_id: HashMap<String, CircRecord>,
    junction_read_to_circ: HashMap<String, String>,
    clusters_by_chr: HashMap<String, Vec<CircCluster>>,
    read_len: i32,
    candidates: Vec<PositiveCandidate>,
    coverage: HashMap<String, HashMap<i32, u32>>,
    read_mappings: HashMap<String, [Vec<ReadMapping>; 2]>,
    seen_junction_reads: HashSet<String>,
    stats: AsStats,
}

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
        group.push(AsAlignment {
            flag,
            chr,
            pos,
            mapq,
            cigar,
            seq,
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
        group.push(AsAlignment {
            flag,
            chr,
            pos,
            mapq,
            cigar: cigar_buf.clone(),
            seq: seq_buf.clone(),
        });
    }
    if !group.is_empty() {
        let id = String::from_utf8_lossy(&current_id);
        on_group(&id, &group)?;
    }
    Ok(())
}

/// Processes one queryname group through CIRI-AS known/non-known read paths.
fn process_group(read_id: &str, records: &[AsAlignment], state: &mut ScanState) -> Result<()> {
    if records.is_empty() {
        return Ok(());
    }
    record_cluster_coverage(records, state);
    if state.junction_read_to_circ.contains_key(read_id) {
        state.seen_junction_reads.insert(read_id.to_string());
        record_mapping_detail(read_id, records, state);
        mapping_check1(read_id, records, true, state)?;
    } else if records.len() > 2 && overlaps_any_circ_cluster(records, state) {
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
        let state = ScanState {
            circ_by_id: HashMap::new(),
            junction_read_to_circ: HashMap::new(),
            clusters_by_chr: HashMap::new(),
            read_len: 100,
            candidates: Vec::new(),
            coverage,
            read_mappings: HashMap::new(),
            seen_junction_reads: HashSet::new(),
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
