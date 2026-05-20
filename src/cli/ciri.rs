//! User-facing CLI entry point for the `ciri` binary.
//!
//! This module owns argument parsing and orchestration for the user pipeline.
//! Keeping it separate from `src/bin/ciri.rs` lets the binary remain a thin
//! wrapper while the behaviorally sensitive logic stays testable in the library.

use anyhow::{bail, Result};
use chrono::Local;
use clap::Parser;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::env;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::Command;
use std::time::Instant;

use crate::annotation::Annotation;
use crate::ciri_as::{
    rebuild_major_isoforms_from_segments, run_ciri_as, AsConfig, IsoformRunSummary,
};
use crate::fasta::FastaReader;
use crate::runtime::init_runtime;
use crate::sam_bam::{check_bam_sorting, detect_format, InputFormat};
use crate::scan1::Scan1;
use crate::scan2::Scan2;
use crate::summary::Summary;
use crate::utils::{
    bsj1_path_for_output, bsj2_path_for_output, bsj_path_for_output, fsj_path_for_output,
    log_path_for_output, parse_mem_str, perf_path_for_output, result_path_for_output,
    segments1_path_for_output, segments2_path_for_output, segments_non_bsj_path_for_output,
    trace_path_for_output,
};

const IGV_BSJ_COLOR: &str = "220,53,69";
const IGV_BACKWARD_COLOR: &str = "25,118,210";
const IGV_OUTWARD_COLOR: &str = "245,124,0";

/// Parsed command-line arguments for the end-to-end pipeline.
///
/// Defaults are kept aligned with CIRI3 unless there is explicit evidence that a
/// different value is required for parity.
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    disable_version_flag = true,
    about = "Run the CIRI Rust analysis pipeline",
    long_about = None
)]
struct Args {
    /// Path to the input SAM/BAM file.
    #[arg(short = 'i', long = "in")]
    in_sam: String,

    /// Output prefix; final outputs include `.out`, `.bsj`, `.segments`,
    /// `.isoforms.gtf`, `.isoforms.fa`, and IGV review sidecars.
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

    /// Stringency level (0, 1, 2). CIRI defaults to 0 to retain candidate
    /// circRNAs for downstream segment/isoform filtering; Java CIRI3 defaults
    /// to 2.
    #[arg(short = 's', long = "stringency", default_value_t = 0)]
    stringency: i32,

    /// Max spanning distance of circRNAs (Java -Max, default 200000)
    #[arg(long = "max-span", default_value_t = 200000)]
    max_span: i32,

    /// Min spanning distance of circRNAs. CIRI defaults to 50 to retain
    /// short candidate circRNAs for downstream filtering; Java CIRI3 defaults
    /// to 140.
    #[arg(long = "min-span", default_value_t = 50)]
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

    /// Comma-separated read IDs to trace through Scan1/Scan2.
    ///
    /// When set, detailed trace lines are written to `<prefix>.trace.log`.
    #[arg(long = "trace", value_name = "READS")]
    trace_reads: Option<String>,

    /// Keep internal pipeline temporary files for debugging.
    #[arg(long = "debug", default_value_t = false)]
    debug: bool,

    /// Enable profiling and write the report to `<prefix>.perf.log`.
    #[arg(long = "perf", default_value_t = false)]
    perf: bool,

    /// Resume from the completed merged segments checkpoint.
    ///
    /// The resume logic intentionally ignores shard-local `.part_*.tmp` files.
    /// If `<prefix>.segments` exists, only isoforms are rebuilt from
    /// `<prefix>.out + <prefix>.segments`. Earlier merged outputs such as
    /// `<prefix>.out + <prefix>.bsj` are intentionally not resumable because
    /// rebuilding segments still requires a full BAM/SAM rescan. All normal
    /// required inputs still must be provided so `--continue` stays an
    /// execution-mode switch, not a separate CLI shape.
    #[arg(long = "continue", default_value_t = false)]
    continue_run: bool,

    /// Print version.
    #[arg(short = 'v', long = "version", action = clap::ArgAction::SetTrue)]
    _version: bool,
}

/// Emits one aligned, timestamped progress line.
///
/// Logs are mirrored to stdout and to `<prefix>.log` so the user can inspect a
/// durable run record after the process exits.
fn log_info(log_writer: &mut BufWriter<File>, label: &str, msg: &str) -> Result<()> {
    let now = Local::now().format("%Y-%m-%d %H:%M:%S");
    let line = format!("{} [INFO] {:<18}: {}", now, label, msg);
    println!("{}", line);
    writeln!(log_writer, "{}", line)?;
    log_writer.flush()?;
    Ok(())
}

/// Sorts the mate-level BSJ rows and writes the final user-facing `.bsj`.
///
/// This runs strictly after Summary has consumed only `priority=1` rows, so the
/// user-facing sort order and `priority=0` evidence cannot perturb `.out`.
fn write_display_bsj(final_bsj: &str, bsj1_path: &str, bsj2_path: &str) -> Result<()> {
    let mut rows: Vec<(String, usize, usize, usize, String)> = Vec::new();
    for (path, stage_rank, stage_name) in
        [(bsj1_path, 0usize, "scan1"), (bsj2_path, 1usize, "scan2")]
    {
        if let Ok(file) = File::open(path) {
            let reader = BufReader::new(file);
            for line_res in reader.lines() {
                let line = line_res?;
                if line.is_empty() {
                    continue;
                }
                let mut parts = line.splitn(4, '\t');
                if let (Some(read_id), Some(mate_label), Some(priority), Some(_)) =
                    (parts.next(), parts.next(), parts.next(), parts.next())
                {
                    let mate_rank = if mate_label == "R1" { 0 } else { 1 };
                    let priority_rank = if priority == "1" { 0 } else { 1 };
                    rows.push((
                        read_id.to_string(),
                        mate_rank,
                        priority_rank,
                        stage_rank,
                        format!("{}\t{}", line, stage_name),
                    ));
                }
            }
        }
    }
    rows.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| a.1.cmp(&b.1))
            .then_with(|| a.2.cmp(&b.2))
            .then_with(|| a.3.cmp(&b.3))
    });
    let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(final_bsj)?);
    for (_, _, _, _, line) in rows {
        writeln!(writer, "{}", line)?;
    }
    writer.flush()?;
    Ok(())
}

/// Writes the `.out`-scoped BEDPE track and returns its output path and row count.
fn write_bsj_bedpe(out_prefix: &str, result_output: &str) -> Result<(String, usize)> {
    let bedpe_path = format!("{}.bedpe", out_prefix);
    let rows = write_bsj_bedpe_from_out(result_output, &bedpe_path)?;
    Ok((bedpe_path, rows))
}

/// Converts Summary-confirmed circRNA sites into a BEDPE BSJ arc track.
///
/// The BEDPE track intentionally ignores internal structure: each `.out` row
/// becomes one pair connecting the left and right BSJ anchors. BED coordinates
/// are 0-based half-open, while CIRI `.out` coordinates are 1-based closed.
fn write_bsj_bedpe_from_out(out_path: &str, bedpe_path: &str) -> Result<usize> {
    let file = File::open(out_path)?;
    let reader = BufReader::new(file);
    let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(bedpe_path)?);
    let mut rows = 0usize;
    for line_result in reader.lines() {
        let line = line_result?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 5 {
            continue;
        }
        let Ok(start) = cols[2].parse::<i32>() else {
            continue;
        };
        let Ok(end) = cols[3].parse::<i32>() else {
            continue;
        };
        if start < 1 || end < start {
            continue;
        }
        let score = cols
            .get(4)
            .copied()
            .filter(|value| value.parse::<i64>().is_ok())
            .unwrap_or("0");
        let strand = cols.get(10).copied().unwrap_or(".");
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            cols[1],
            start - 1,
            start,
            cols[1],
            end - 1,
            end,
            cols[0],
            score,
            strand,
            strand
        )?;
        rows += 1;
    }
    writer.flush()?;
    Ok(rows)
}

/// Writes the `.segments`-scoped synthetic BAM/BAI track and records the result.
fn write_and_log_segments_bam(
    log_writer: &mut BufWriter<File>,
    segments_output: &str,
    out_prefix: &str,
    reference: &HashMap<String, String>,
    threads: usize,
) -> Result<()> {
    let bam_path = format!("{}.segments.bam", out_prefix);
    write_segments_review_tracks_from_segments(
        segments_output,
        &bam_path,
        out_prefix,
        reference,
        threads,
    )?;
    let message = format!("{}", bam_path);
    log_info(log_writer, "Output segments BAM", &message)?;
    Ok(())
}

/// Formats the isoform stage summary for the CLI run log.
fn format_isoform_summary(summary: IsoformRunSummary) -> String {
    format!(
        "Assembled {} high-confidence isoforms from {} circRNAs",
        summary.fasta_isoforms, summary.fasta_circ_rnas,
    )
}

#[derive(Debug, Clone, Copy)]
struct IgvSegmentsColumns {
    read_id: usize,
    type_name: usize,
    circ_id: usize,
    chrom: usize,
    r1_cs: Option<usize>,
    r1_segments: usize,
    r2_cs: Option<usize>,
    r2_segments: usize,
}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
struct IgvChromSortKey {
    group: u8,
    rank: u32,
    suffix: String,
}

#[derive(Debug, Clone)]
struct IgvSegmentRow {
    chrom_key: IgvChromSortKey,
    chrom: String,
    chrom_start: i32,
    chrom_end: i32,
    name: String,
    sam_line: Option<String>,
}

#[derive(Debug, Clone)]
struct IgvPartPayload {
    seq: String,
    cs: String,
}

/// Converts read-level segment chains into an indexed BAM review track.
///
/// Each R1/R2 chain is split at explicit `<bsj>` markers so BSJ-crossing pieces
/// render as separate IGV alignments. Within each piece, ordinary internal
/// junctions remain CIGAR `N` gaps. Rows are coordinate-sorted before BAM
/// conversion so IGV region loading can use the generated BAI index.
fn write_segments_review_tracks_from_segments(
    segments_path: &str,
    bam_path: &str,
    out_prefix: &str,
    reference: &HashMap<String, String>,
    threads: usize,
) -> Result<()> {
    let file = File::open(segments_path)?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .transpose()?
        .ok_or_else(|| anyhow::anyhow!("empty segments file: {}", segments_path))?;
    let columns = igv_segments_columns(&header, segments_path)?;
    let mut rows = Vec::new();
    for line_result in lines {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        let Some(type_name) = cols.get(columns.type_name).copied() else {
            continue;
        };
        if !matches!(type_name, "bsj" | "backward" | "outward") {
            continue;
        }
        let Some(read_id) = cols.get(columns.read_id).copied() else {
            continue;
        };
        let Some(circ_id) = cols.get(columns.circ_id).copied() else {
            continue;
        };
        let Some(chrom) = cols.get(columns.chrom).copied() else {
            continue;
        };
        let color = igv_segment_color(type_name);
        for (mate_label, segments, cs) in [
            (
                "R1",
                cols.get(columns.r1_segments).copied().unwrap_or("NA"),
                columns
                    .r1_cs
                    .and_then(|idx| cols.get(idx).copied())
                    .unwrap_or("NA"),
            ),
            (
                "R2",
                cols.get(columns.r2_segments).copied().unwrap_or("NA"),
                columns
                    .r2_cs
                    .and_then(|idx| cols.get(idx).copied())
                    .unwrap_or("NA"),
            ),
        ] {
            push_one_segments_alignment_chain(
                &mut rows, chrom, read_id, type_name, circ_id, mate_label, segments, cs, color,
                reference,
            )?;
        }
    }
    rows.sort_by(compare_igv_segment_rows);
    let reference_lengths: HashMap<String, usize> = reference
        .iter()
        .map(|(chrom, seq)| (chrom.clone(), seq.len()))
        .collect();
    write_segments_bam_from_sorted_rows(&rows, bam_path, out_prefix, &reference_lengths, threads)?;
    Ok(())
}

/// Resolves the public `<prefix>.segments` columns needed by IGV export.
fn igv_segments_columns(header: &str, path: &str) -> Result<IgvSegmentsColumns> {
    let headers: Vec<&str> = header.split('\t').collect();
    let idx = |name: &str| -> Result<usize> {
        headers
            .iter()
            .position(|field| *field == name)
            .ok_or_else(|| anyhow::anyhow!("missing `{}` column in {}", name, path))
    };
    Ok(IgvSegmentsColumns {
        read_id: idx("read_id")?,
        type_name: idx("type")?,
        circ_id: idx("circ_id")?,
        chrom: idx("chrom")?,
        r1_cs: headers.iter().position(|field| *field == "r1_cs"),
        r1_segments: idx("r1_segments")?,
        r2_cs: headers.iter().position(|field| *field == "r2_cs"),
        r2_segments: idx("r2_segments")?,
    })
}

/// Appends one mate chain as one or more synthetic alignments split at BSJ markers.
fn push_one_segments_alignment_chain(
    rows: &mut Vec<IgvSegmentRow>,
    chrom: &str,
    read_id: &str,
    type_name: &str,
    circ_id: &str,
    mate_label: &str,
    segments: &str,
    cs: &str,
    color: &str,
    reference: &HashMap<String, String>,
) -> Result<()> {
    if segments == "NA" || segments.is_empty() {
        return Ok(());
    }
    let part_payloads = igv_part_payloads_from_cs(chrom, segments, cs, reference);
    let mut part = 1usize;
    let mut blocks: Vec<(i32, i32, char)> = Vec::new();
    for token in segments.split('|') {
        if token == "<bsj>" {
            let payload = part_payloads
                .as_ref()
                .and_then(|payloads| payloads.get(part - 1));
            if push_segments_alignment_part(
                rows,
                chrom,
                read_id,
                type_name,
                circ_id,
                mate_label,
                part,
                color,
                &mut blocks,
                payload,
            )? {
                part += 1;
            }
            continue;
        }
        if let Some((start, end, strand)) = parse_igv_segment_token(token) {
            if igv_part_coordinate_break(&blocks, start, end) {
                let payload = part_payloads
                    .as_ref()
                    .and_then(|payloads| payloads.get(part - 1));
                if push_segments_alignment_part(
                    rows,
                    chrom,
                    read_id,
                    type_name,
                    circ_id,
                    mate_label,
                    part,
                    color,
                    &mut blocks,
                    payload,
                )? {
                    part += 1;
                }
            }
            blocks.push((start, end, strand));
        }
    }
    push_segments_alignment_part(
        rows,
        chrom,
        read_id,
        type_name,
        circ_id,
        mate_label,
        part,
        color,
        &mut blocks,
        part_payloads
            .as_ref()
            .and_then(|payloads| payloads.get(part - 1)),
    )?;
    Ok(())
}

/// Returns whether a new segment must start a fresh BAM alignment part.
///
/// SAM/BAM CIGARs are reference-coordinate ordered. CIRI read chains can move
/// backward across circular or ambiguous topology even when no explicit
/// `<bsj>` token is present, so those coordinate breaks are written as separate
/// synthetic records instead of forcing the sequence-aware cs payload through a
/// reordered CIGAR.
fn igv_part_coordinate_break(blocks: &[(i32, i32, char)], start: i32, end: i32) -> bool {
    blocks.last().is_some_and(|(prev_start, prev_end, _)| {
        (*prev_start, *prev_end) > (start, end) || start <= *prev_end
    })
}

/// Appends one synthetic alignment part and clears the accumulated segment blocks.
fn push_segments_alignment_part(
    rows: &mut Vec<IgvSegmentRow>,
    chrom: &str,
    read_id: &str,
    type_name: &str,
    circ_id: &str,
    mate_label: &str,
    part: usize,
    color: &str,
    blocks: &mut Vec<(i32, i32, char)>,
    payload: Option<&IgvPartPayload>,
) -> Result<bool> {
    if blocks.is_empty() {
        return Ok(false);
    }
    blocks.retain(|(start, end, _)| *start >= 1 && *end >= *start);
    if blocks.is_empty() {
        return Ok(false);
    }
    blocks.sort_by_key(|(start, end, _)| (*start, *end));
    let chrom_start = blocks.iter().map(|(start, _, _)| *start).min().unwrap() - 1;
    let chrom_end = blocks.iter().map(|(_, end, _)| *end).max().unwrap();
    let strand = igv_bed12_strand(blocks);
    let name = format!(
        "{}|{}|{}|{}|part{}",
        read_id, type_name, circ_id, mate_label, part
    );
    rows.push(IgvSegmentRow {
        chrom_key: igv_chrom_sort_key(chrom),
        chrom: chrom.to_string(),
        chrom_start,
        chrom_end,
        name: name.clone(),
        sam_line: synthetic_segments_sam_line(
            chrom,
            chrom_start,
            name.as_str(),
            type_name,
            circ_id,
            mate_label,
            part,
            color,
            strand,
            blocks,
            payload,
        ),
    });
    blocks.clear();
    Ok(true)
}

/// Returns a human-genome-friendly chromosome key for coordinate sorting.
///
/// IGV accepts generic BED files, but review sessions here usually use hg38.
/// Putting `chr1..chr22, chrX, chrY, chrM/MT` before alternate contigs makes
/// the sidecar easier to scan while still keeping non-standard contigs sorted.
fn igv_chrom_sort_key(chrom: &str) -> IgvChromSortKey {
    let core = chrom.strip_prefix("chr").unwrap_or(chrom);
    if let Ok(rank) = core.parse::<u32>() {
        return IgvChromSortKey {
            group: 0,
            rank,
            suffix: String::new(),
        };
    }
    match core {
        "X" => IgvChromSortKey {
            group: 0,
            rank: 23,
            suffix: String::new(),
        },
        "Y" => IgvChromSortKey {
            group: 0,
            rank: 24,
            suffix: String::new(),
        },
        "M" | "MT" => IgvChromSortKey {
            group: 0,
            rank: 25,
            suffix: String::new(),
        },
        _ => IgvChromSortKey {
            group: 1,
            rank: 0,
            suffix: chrom.to_string(),
        },
    }
}

/// Orders synthetic alignment rows by coordinate, with the query name as tie-breaker.
fn compare_igv_segment_rows(left: &IgvSegmentRow, right: &IgvSegmentRow) -> Ordering {
    left.chrom_key
        .cmp(&right.chrom_key)
        .then_with(|| left.chrom_start.cmp(&right.chrom_start))
        .then_with(|| left.chrom_end.cmp(&right.chrom_end))
        .then_with(|| left.name.cmp(&right.name))
}

/// Builds a SAM alignment line from one segment row.
///
/// When the chain-level `r*_cs` payload can be checked against the part's
/// genomic span, the aligned query bases and BAM CIGAR are reconstructed from
/// cs. If the cs payload is missing or inconsistent, the writer falls back to
/// `N` bases and omits `cs:Z` rather than producing an invalid BAM.
fn synthetic_segments_sam_line(
    chrom: &str,
    chrom_start: i32,
    name: &str,
    type_name: &str,
    circ_id: &str,
    mate_label: &str,
    part: usize,
    color: &str,
    strand: char,
    blocks: &[(i32, i32, char)],
    payload: Option<&IgvPartPayload>,
) -> Option<String> {
    let (fallback_cigar, fallback_query_len) = synthetic_segments_cigar(blocks)?;
    if fallback_query_len == 0 {
        return None;
    }
    let flag = if strand == '-' { 16 } else { 0 };
    let mapq = igv_segment_mapq(type_name);
    let fallback_ref_len = cigar_query_ref_consumption(&fallback_cigar)?.1;
    let usable_payload = payload.and_then(|payload| {
        let cigar = cigar_from_part_cs(&payload.cs)?;
        let (query_len, ref_len) = cigar_query_ref_consumption(&cigar)?;
        (query_len == payload.seq.len()
            && ref_len == fallback_ref_len
            && payload
                .seq
                .bytes()
                .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N')))
        .then_some((payload, cigar))
    });
    let (cigar, seq, cs_tag) = if let Some((payload, cigar)) = usable_payload {
        (cigar, payload.seq.clone(), format!("\tcs:Z:{}", payload.cs))
    } else {
        (
            fallback_cigar,
            "N".repeat(fallback_query_len),
            String::new(),
        )
    };
    Some(format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t*\tRG:Z:{}\tYC:Z:{}\tZT:Z:{}\tCI:Z:{}\tML:Z:{}\tPT:i:{}{}\n",
        name,
        flag,
        chrom,
        chrom_start + 1,
        mapq,
        cigar,
        seq,
        type_name,
        color,
        type_name,
        circ_id,
        mate_label,
        part,
        cs_tag
    ))
}

#[derive(Debug, Clone, Copy)]
enum IgvCsOp<'a> {
    Match(usize),
    Sub { query: u8 },
    Ins(&'a str),
    Del(usize),
    Skip { raw: &'a str, len: usize, bsj: bool },
}

#[derive(Debug)]
enum IgvChainUnit {
    Segment { start: i32, end: i32 },
    Bsj,
}

/// Reconstructs per-BSJ-part sequence and cs payloads for `.segments.bam`.
///
/// The public `.segments` row stores one read-chain cs payload, while the review
/// BAM intentionally splits chains at `<bsj>`. This parser walks the segment
/// tokens and cs operations together, using reference sequence for `:` match
/// runs and starting a new part at CIRI's custom `<...>` back-splice op.
fn igv_part_payloads_from_cs(
    chrom: &str,
    segments: &str,
    cs: &str,
    reference: &HashMap<String, String>,
) -> Option<Vec<IgvPartPayload>> {
    if matches!(cs, "NA" | "*" | "") {
        return None;
    }
    let chr_seq = reference.get(chrom)?;
    let units = igv_chain_units(segments)?;
    let mut unit_idx = first_segment_unit_idx(&units)?;
    let mut ref_pos = match units.get(unit_idx)? {
        IgvChainUnit::Segment { start, .. } => *start,
        IgvChainUnit::Bsj => return None,
    };
    let mut parts = Vec::new();
    let mut current = IgvPartPayload {
        seq: String::new(),
        cs: String::new(),
    };

    let mut pos = 0usize;
    while pos < cs.len() {
        let (op, next_pos) = parse_igv_cs_op(cs, pos)?;
        pos = next_pos;
        match op {
            IgvCsOp::Match(len) => {
                let segment_end = segment_end_at(&units, unit_idx)?;
                if ref_pos + len as i32 - 1 > segment_end {
                    return None;
                }
                let start = (ref_pos - 1) as usize;
                let end = start + len;
                current.seq.push_str(chr_seq.get(start..end)?);
                current.cs.push(':');
                current.cs.push_str(&len.to_string());
                ref_pos += len as i32;
            }
            IgvCsOp::Sub { query } => {
                if ref_pos > segment_end_at(&units, unit_idx)? {
                    return None;
                }
                let raw = &cs[next_pos - 3..next_pos];
                current.seq.push(query.to_ascii_uppercase() as char);
                current.cs.push_str(raw);
                ref_pos += 1;
            }
            IgvCsOp::Ins(seq) => {
                current.seq.push_str(&seq.to_ascii_uppercase());
                current.cs.push('+');
                current.cs.push_str(seq);
            }
            IgvCsOp::Del(len) => {
                if ref_pos + len as i32 - 1 > segment_end_at(&units, unit_idx)? {
                    return None;
                }
                current.cs.push_str(&cs[next_pos - len - 1..next_pos]);
                ref_pos += len as i32;
            }
            IgvCsOp::Skip { raw, len, bsj } => {
                let next_unit_idx = next_segment_unit_idx(&units, unit_idx, bsj)?;
                let coordinate_break = match (units.get(unit_idx)?, units.get(next_unit_idx)?) {
                    (
                        IgvChainUnit::Segment {
                            start: prev_start,
                            end: prev_end,
                        },
                        IgvChainUnit::Segment { start, end },
                    ) => igv_part_coordinate_break(&[(*prev_start, *prev_end, '+')], *start, *end),
                    _ => return None,
                };
                if bsj || coordinate_break {
                    parts.push(current);
                    current = IgvPartPayload {
                        seq: String::new(),
                        cs: String::new(),
                    };
                } else {
                    current.cs.push_str(raw);
                }
                unit_idx = next_unit_idx;
                ref_pos = match units.get(unit_idx)? {
                    IgvChainUnit::Segment { start, .. } => *start,
                    IgvChainUnit::Bsj => return None,
                };
                let _ = len;
            }
        }
    }
    parts.push(current);
    Some(parts)
}

/// Parses segment tokens into the same chain units consumed by `r*_cs`.
fn igv_chain_units(segments: &str) -> Option<Vec<IgvChainUnit>> {
    let mut units = Vec::new();
    for token in segments.split('|') {
        if token == "<bsj>" {
            units.push(IgvChainUnit::Bsj);
            continue;
        }
        let (start, end, _strand) = parse_igv_segment_token(token)?;
        units.push(IgvChainUnit::Segment { start, end });
    }
    (!units.is_empty()).then_some(units)
}

fn first_segment_unit_idx(units: &[IgvChainUnit]) -> Option<usize> {
    units
        .iter()
        .position(|unit| matches!(unit, IgvChainUnit::Segment { .. }))
}

fn segment_end_at(units: &[IgvChainUnit], unit_idx: usize) -> Option<i32> {
    match units.get(unit_idx)? {
        IgvChainUnit::Segment { end, .. } => Some(*end),
        IgvChainUnit::Bsj => None,
    }
}

fn next_segment_unit_idx(
    units: &[IgvChainUnit],
    current_idx: usize,
    expect_bsj: bool,
) -> Option<usize> {
    let mut idx = current_idx + 1;
    if expect_bsj {
        if !matches!(units.get(idx)?, IgvChainUnit::Bsj) {
            return None;
        }
        idx += 1;
    } else if matches!(units.get(idx)?, IgvChainUnit::Bsj) {
        return None;
    }
    while idx < units.len() {
        if matches!(units[idx], IgvChainUnit::Segment { .. }) {
            return Some(idx);
        }
        idx += 1;
    }
    None
}

/// Parses one minimap2-style short-form cs operation, including CIRI's `<...>`.
fn parse_igv_cs_op(cs: &str, pos: usize) -> Option<(IgvCsOp<'_>, usize)> {
    let bytes = cs.as_bytes();
    let op = *bytes.get(pos)? as char;
    let mut cursor = pos + 1;
    match op {
        ':' => {
            let start = cursor;
            while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
                cursor += 1;
            }
            if cursor == start {
                return None;
            }
            Some((IgvCsOp::Match(cs[start..cursor].parse().ok()?), cursor))
        }
        '*' => {
            let _ref_base = *bytes.get(cursor)?;
            let query = *bytes.get(cursor + 1)?;
            Some((IgvCsOp::Sub { query }, cursor + 2))
        }
        '+' => {
            let start = cursor;
            while bytes
                .get(cursor)
                .is_some_and(|base| base.is_ascii_alphabetic())
            {
                cursor += 1;
            }
            (cursor > start).then_some((IgvCsOp::Ins(&cs[start..cursor]), cursor))
        }
        '-' => {
            let start = cursor;
            while bytes
                .get(cursor)
                .is_some_and(|base| base.is_ascii_alphabetic())
            {
                cursor += 1;
            }
            (cursor > start).then_some((IgvCsOp::Del(cursor - start), cursor))
        }
        '~' | '<' => {
            cursor += 2;
            let len_start = cursor;
            while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
                cursor += 1;
            }
            if cursor == len_start {
                return None;
            }
            cursor += 2;
            if cursor > cs.len() {
                return None;
            }
            Some((
                IgvCsOp::Skip {
                    raw: &cs[pos..cursor],
                    len: cs[len_start..cursor - 2].parse().ok()?,
                    bsj: op == '<',
                },
                cursor,
            ))
        }
        _ => None,
    }
}

/// Checks that a part-level cs payload consumes the same query and reference as
/// the CIGAR that will be written to BAM.
#[cfg(test)]
fn cs_consumption_matches_cigar(cs: &str, cigar: &str) -> bool {
    let Some((cs_query, cs_ref)) = cs_query_ref_consumption(cs) else {
        return false;
    };
    let Some((cigar_query, cigar_ref)) = cigar_query_ref_consumption(cigar) else {
        return false;
    };
    cs_query == cigar_query && cs_ref == cigar_ref
}

/// Converts a part-level cs payload into a BAM CIGAR.
///
/// `:` runs and `*` substitutions both become `M`, while `+`, `-`, and `~`
/// become `I`, `D`, and `N`. CIRI's custom `<...>` should have been consumed
/// when the chain was split into BAM parts, so it is rejected here.
fn cigar_from_part_cs(cs: &str) -> Option<String> {
    let mut out = String::new();
    let mut pending: Option<(char, usize)> = None;
    let mut pos = 0usize;
    while pos < cs.len() {
        let (op, next) = parse_igv_cs_op(cs, pos)?;
        match op {
            IgvCsOp::Match(len) => push_cigar_run(&mut out, &mut pending, 'M', len),
            IgvCsOp::Sub { .. } => push_cigar_run(&mut out, &mut pending, 'M', 1),
            IgvCsOp::Ins(seq) => push_cigar_run(&mut out, &mut pending, 'I', seq.len()),
            IgvCsOp::Del(len) => push_cigar_run(&mut out, &mut pending, 'D', len),
            IgvCsOp::Skip { len, bsj, .. } => {
                if bsj {
                    return None;
                }
                push_cigar_run(&mut out, &mut pending, 'N', len);
            }
        }
        pos = next;
    }
    flush_cigar_run(&mut out, &mut pending);
    (!out.is_empty()).then_some(out)
}

fn push_cigar_run(out: &mut String, pending: &mut Option<(char, usize)>, op: char, len: usize) {
    if len == 0 {
        return;
    }
    match pending.as_mut() {
        Some((pending_op, pending_len)) if *pending_op == op => *pending_len += len,
        _ => {
            flush_cigar_run(out, pending);
            *pending = Some((op, len));
        }
    }
}

fn flush_cigar_run(out: &mut String, pending: &mut Option<(char, usize)>) {
    if let Some((op, len)) = pending.take() {
        out.push_str(&len.to_string());
        out.push(op);
    }
}

#[cfg(test)]
fn cs_query_ref_consumption(cs: &str) -> Option<(usize, usize)> {
    let mut query = 0usize;
    let mut reference = 0usize;
    let mut pos = 0usize;
    while pos < cs.len() {
        let (op, next) = parse_igv_cs_op(cs, pos)?;
        match op {
            IgvCsOp::Match(len) => {
                query += len;
                reference += len;
            }
            IgvCsOp::Sub { .. } => {
                query += 1;
                reference += 1;
            }
            IgvCsOp::Ins(seq) => query += seq.len(),
            IgvCsOp::Del(len) => reference += len,
            IgvCsOp::Skip { len, .. } => reference += len,
        }
        pos = next;
    }
    Some((query, reference))
}

fn cigar_query_ref_consumption(cigar: &str) -> Option<(usize, usize)> {
    let mut query = 0usize;
    let mut reference = 0usize;
    let mut number = String::new();
    for op in cigar.chars() {
        if op.is_ascii_digit() {
            number.push(op);
            continue;
        }
        if number.is_empty() {
            return None;
        }
        let len = number.parse::<usize>().ok()?;
        match op {
            'M' | '=' | 'X' => {
                query += len;
                reference += len;
            }
            'I' | 'S' => query += len,
            'D' | 'N' => reference += len,
            'H' | 'P' => {}
            _ => return None,
        }
        number.clear();
    }
    number.is_empty().then_some((query, reference))
}

/// Converts sorted segment blocks into a splice-aware CIGAR string.
fn synthetic_segments_cigar(blocks: &[(i32, i32, char)]) -> Option<(String, usize)> {
    let mut cigar = String::new();
    let mut query_len = 0usize;
    let mut prev_end0: Option<i32> = None;
    for &(start, end, _) in blocks {
        if start < 1 || end < start {
            return None;
        }
        let start0 = start - 1;
        if let Some(prev) = prev_end0 {
            let gap = start0 - prev;
            if gap < 0 {
                return None;
            }
            if gap > 0 {
                cigar.push_str(&format!("{}N", gap));
            }
        }
        let len = (end - start + 1) as usize;
        cigar.push_str(&format!("{}M", len));
        query_len += len;
        prev_end0 = Some(end);
    }
    Some((cigar, query_len))
}

/// Returns a synthetic MAPQ that keeps segment classes visually separable.
fn igv_segment_mapq(type_name: &str) -> u8 {
    match type_name {
        "bsj" => 60,
        "backward" => 45,
        "outward" => 30,
        _ => 0,
    }
}

/// Writes a coordinate-sorted synthetic BAM and BAI from already sorted rows.
///
/// This intentionally requires `samtools` because the repo does not yet own a
/// native BAI writer and `.segments.bam/.bai` is part of the normal review
/// output contract. The SAM temp is deleted after successful conversion so the
/// indexed BAM is the durable high-performance IGV artifact.
fn write_segments_bam_from_sorted_rows(
    rows: &[IgvSegmentRow],
    bam_path: &str,
    out_prefix: &str,
    reference_lengths: &HashMap<String, usize>,
    threads: usize,
) -> Result<()> {
    let samtools = require_samtools()?;
    let sam_path = format!("{}.segments.sam.tmp", out_prefix);
    write_segments_sam_from_sorted_rows(rows, &sam_path, reference_lengths)?;
    let view_threads = threads.max(1).to_string();
    let status = Command::new(&samtools)
        .arg("view")
        .arg("-@")
        .arg(&view_threads)
        .arg("-b")
        .arg("-o")
        .arg(bam_path)
        .arg(&sam_path)
        .status()?;
    if !status.success() {
        bail!("samtools view failed while creating {}", bam_path);
    }
    let status = Command::new(&samtools)
        .arg("index")
        .arg("-@")
        .arg(&view_threads)
        .arg(bam_path)
        .status()?;
    if !status.success() {
        bail!("samtools index failed while creating {}.bai", bam_path);
    }
    let _ = fs::remove_file(&sam_path);
    Ok(())
}

/// Writes sorted SAM records used as the conversion source for the review BAM.
fn write_segments_sam_from_sorted_rows(
    rows: &[IgvSegmentRow],
    sam_path: &str,
    reference_lengths: &HashMap<String, usize>,
) -> Result<()> {
    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, File::create(sam_path)?);
    writeln!(writer, "@HD\tVN:1.6\tSO:coordinate")?;
    for (chrom, len) in igv_bam_reference_lengths(rows, reference_lengths) {
        writeln!(writer, "@SQ\tSN:{}\tLN:{}", chrom, len)?;
    }
    writeln!(
        writer,
        "@RG\tID:bsj\tSM:CIRI_segments\tDS:back-spliced junction read segments"
    )?;
    writeln!(
        writer,
        "@RG\tID:backward\tSM:CIRI_segments\tDS:backward read segments"
    )?;
    writeln!(
        writer,
        "@RG\tID:outward\tSM:CIRI_segments\tDS:outward-facing read segments"
    )?;
    writeln!(
        writer,
        "@PG\tID:CIRI\tPN:CIRI\tVN:{}",
        env!("CARGO_PKG_VERSION")
    )?;
    for row in rows {
        if let Some(line) = &row.sam_line {
            writer.write_all(line.as_bytes())?;
        }
    }
    writer.flush()?;
    Ok(())
}

/// Returns reference lengths for all contigs touched by the synthetic rows.
fn igv_bam_reference_lengths(
    rows: &[IgvSegmentRow],
    reference_lengths: &HashMap<String, usize>,
) -> Vec<(String, usize)> {
    let mut lengths: HashMap<String, usize> = HashMap::new();
    for row in rows {
        let observed_len = usize::try_from(row.chrom_end.max(1)).unwrap_or(1);
        let len = reference_lengths
            .get(row.chrom.as_str())
            .copied()
            .unwrap_or(observed_len)
            .max(observed_len);
        lengths
            .entry(row.chrom.clone())
            .and_modify(|current| *current = (*current).max(len))
            .or_insert(len);
    }
    let mut entries: Vec<_> = lengths.into_iter().collect();
    entries
        .sort_by(|(left, _), (right, _)| igv_chrom_sort_key(left).cmp(&igv_chrom_sort_key(right)));
    entries
}

/// Finds the `samtools` executable used for BAM/BAI review-track generation.
fn find_samtools() -> Option<String> {
    if let Ok(path) = env::var("SAMTOOLS") {
        if samtools_is_usable(&path) {
            return Some(path);
        }
    }
    if samtools_is_usable("samtools") {
        return Some("samtools".to_string());
    }
    None
}

/// Requires `samtools` before the pipeline starts touching large input files.
///
/// CIRI-toolkit writes `<prefix>.segments.bam/.bai` as a standard IGV review
/// sidecar. Failing early avoids spending minutes or hours on Scan1/Scan2 only
/// to discover that the final indexed review artifact cannot be produced.
fn require_samtools() -> Result<String> {
    find_samtools().ok_or_else(|| {
        anyhow::anyhow!(
            "`samtools` is not installed. Please install samtools in PATH or set SAMTOOLS=/path/to/samtools"
        )
    })
}

/// Checks whether a candidate `samtools` command can be executed.
fn samtools_is_usable(command: &str) -> bool {
    Command::new(command)
        .arg("--version")
        .output()
        .map(|output| output.status.success())
        .unwrap_or(false)
}

/// Parses a `start-end:strand` segment token from `<prefix>.segments`.
fn parse_igv_segment_token(token: &str) -> Option<(i32, i32, char)> {
    let (range, strand_text) = token.rsplit_once(':')?;
    let (start, end) = range.split_once('-')?;
    let strand = strand_text.chars().next().unwrap_or('.');
    Some((start.parse().ok()?, end.parse().ok()?, strand))
}

/// Returns the synthetic alignment strand when all blocks agree, otherwise `.`.
fn igv_bed12_strand(blocks: &[(i32, i32, char)]) -> char {
    let mut strand = None;
    for &(_, _, current) in blocks {
        if !matches!(current, '+' | '-') {
            return '.';
        }
        match strand {
            Some(previous) if previous != current => return '.',
            Some(_) => {}
            None => strand = Some(current),
        }
    }
    strand.unwrap_or('.')
}

/// Returns a stable RGB color for one segment type.
fn igv_segment_color(type_name: &str) -> &'static str {
    match type_name {
        "bsj" => IGV_BSJ_COLOR,
        "backward" => IGV_BACKWARD_COLOR,
        "outward" => IGV_OUTWARD_COLOR,
        _ => "128,128,128",
    }
}

/// Loads inputs, runs Scan1 -> Scan2 -> Summary, and writes outputs.
///
/// The current `ciri` entry keeps the historical direct CIRI3-style arguments
/// instead of forcing a `detect` subcommand. After Summary finishes, the CLI now
/// always performs the read-level circRNA segments pass that feeds future
/// full-length reconstruction.
pub fn main() -> Result<()> {
    let run_started = Instant::now();
    if env::args()
        .skip(1)
        .any(|arg| arg == "-v" || arg == "--version")
    {
        println!("{} {}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    let args = Args::parse();
    let _ = require_samtools()?;
    let mem_limit = parse_mem_str(&args.mem_per_thread);
    let result_output = result_path_for_output(&args.out_prefix);
    let log_output = log_path_for_output(&args.out_prefix);
    let trace_output = trace_path_for_output(&args.out_prefix);
    let perf_output = perf_path_for_output(&args.out_prefix);
    let bsj1_output = bsj1_path_for_output(&args.out_prefix);
    let bsj_output = bsj_path_for_output(&args.out_prefix);
    let bsj2_output = bsj2_path_for_output(&args.out_prefix);
    let segments1_output = segments1_path_for_output(&args.out_prefix);
    let segments2_output = segments2_path_for_output(&args.out_prefix);
    let segments_non_bsj_output = segments_non_bsj_path_for_output(&args.out_prefix);
    let fsj_output = fsj_path_for_output(&args.out_prefix);
    let mut log_writer = BufWriter::new(File::create(&log_output)?);

    init_runtime(
        args.trace_reads.as_deref(),
        args.trace_reads.as_ref().map(|_| trace_output.as_str()),
        args.perf.then_some(perf_output.as_str()),
    )?;

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()?;
    }

    log_info(&mut log_writer, "Reference FASTA", &args.ref_fasta)?;
    let mut fasta = FastaReader::new();
    fasta.read_fasta(&args.ref_fasta)?;

    log_info(
        &mut log_writer,
        "Annotation GTF",
        args.gtf.as_deref().unwrap_or("None provided"),
    )?;
    let mut annotation = Annotation::new();
    if let Some(gtf_path) = &args.gtf {
        annotation.read_gtf(gtf_path)?;
    }

    if args.trace_reads.is_some() {
        log_info(&mut log_writer, "Read trace", &trace_output)?;
    }
    if args.debug {
        log_info(
            &mut log_writer,
            "Debug temp",
            "Keeping internal temporary files after successful completion",
        )?;
    }
    if args.perf {
        log_info(&mut log_writer, "Perf report", &perf_output)?;
    }

    let input_path = args.in_sam.as_str();
    let format = detect_format(input_path)?;
    let format_str = match format {
        InputFormat::Bam => {
            check_bam_sorting(input_path)?;
            "BAM (queryname-sorted)"
        }
        InputFormat::Sam => "SAM (text-based)",
    };
    log_info(&mut log_writer, "Input format", format_str)?;

    let segments_output = format!("{}.segments", args.out_prefix);
    if args.continue_run {
        if Path::new(&segments_output).is_file() {
            if !Path::new(&result_output).is_file() {
                bail!(
                    "--continue found {} but missing required circRNA table {}",
                    segments_output,
                    result_output
                );
            }
            log_info(
                &mut log_writer,
                "Resume from checkpoint",
                "Start building circRNA isoforms...",
            )?;
            let isoform_summary = rebuild_major_isoforms_from_segments(
                &result_output,
                &segments_output,
                &args.out_prefix,
                &fasta.chr_tcga_map,
                args.gtf.as_ref().map(|_| &annotation),
            )?;
            log_info(
                &mut log_writer,
                "Isoforms output",
                &format!(
                    "{}.isoforms.gtf, {}.isoforms.fa",
                    args.out_prefix, args.out_prefix
                ),
            )?;
            log_info(
                &mut log_writer,
                "Isoforms summary",
                &format_isoform_summary(isoform_summary),
            )?;
            log_info(
                &mut log_writer,
                "Total runtime",
                &format!("{:.2} seconds", run_started.elapsed().as_secs_f64()),
            )?;
            return Ok(());
        }

        bail!(
            "--continue did not find a resumable checkpoint: expected completed {}",
            segments_output
        );
    }

    // Stage boundaries are logged explicitly because most benchmarking and parity
    // work is reasoned about in terms of Scan1 / Scan2 / Summary timings.
    // 3. Scan 1
    log_info(
        &mut log_writer,
        "=== STAGE 1/3 ===",
        "Back-splicing junction identification...",
    )?;
    log_info(
        &mut log_writer,
        "Running scan 1",
        "Identifying BSJ sites...",
    )?;
    let mut scan1 = Scan1::new(
        args.min_mapq,
        args.min_span,
        args.max_span,
        args.linear_range_size_min,
    );
    scan1.set_mem_limit(mem_limit);
    scan1.run_with_priority_and_segments(
        input_path,
        &bsj1_output,
        Some(&segments1_output),
        &fasta.chr_tcga_map,
        &annotation,
    )?;
    log_info(
        &mut log_writer,
        "Scan 1 summary",
        &format!(
            "{} mapped reads, {} BSJ1 reads",
            scan1.mapped_reads, scan1.bsj1_reads
        ),
    )?;

    // 4. Indexing
    log_info(
        &mut log_writer,
        "Loading BSJ sites",
        "Generating candidate BSJ index...",
    )?;
    let scan2_seq_len = (scan1.read_len - 12).max(1);
    // Scan1 is no longer needed once its output file and derived read length have
    // been materialized, so drop it before Scan2 to keep whole-genome RSS lower.
    drop(scan1);
    let mut scan2 = Scan2::new(args.min_mapq, args.linear_range_size_min, scan2_seq_len);
    scan2.set_mem_limit(mem_limit);
    scan2.build_index(&bsj1_output)?;

    // 5. Scan 2
    log_info(
        &mut log_writer,
        "Running scan 2",
        "Curating splicing signals & counting FSJs...",
    )?;
    let scan2_segment_artifacts = scan2.run_with_display_and_segments(
        input_path,
        &bsj2_output,
        &fsj_output,
        Some(&bsj1_output),
        None,
        Some(&segments2_output),
        Some(&segments_non_bsj_output),
        &fasta.chr_tcga_map,
    )?;
    log_info(
        &mut log_writer,
        "Scan 2 summary",
        &format!("{} BSJ2 reads rescued", scan2.rescued_reads),
    )?;
    scan2.release_working_set();

    // 6. Finalization
    log_info(
        &mut log_writer,
        "Post-processing",
        "Clustering sites and filtering results...",
    )?;
    let mut summary = Summary::new(args.stringency);
    summary.run_from_bsj_files(
        &[&bsj1_output, &bsj2_output],
        &result_output,
        &scan2.fsj_map,
        &fasta.chr_tcga_map,
        &annotation,
    )?;
    scan2.release_working_set();
    write_display_bsj(&bsj_output, &bsj1_output, &bsj2_output)?;

    log_info(
        &mut log_writer,
        "BSJ summary",
        &format!(
            // This count comes from Summary's retained circ/read assignments,
            // not from `Scan1 + Scan2` raw BSJ accumulation. Keeping the final
            // user-facing metric here prevents stage-local bookkeeping from
            // being mistaken for the clustered output size.
            "{} circRNAs, {} BSJ reads detected",
            summary.circ_count, summary.final_bsj_reads
        ),
    )?;
    log_info(&mut log_writer, "Output BSJ file", &result_output)?;
    let (bedpe_path, _bedpe_rows) = write_bsj_bedpe(&args.out_prefix, &result_output)?;
    log_info(
        &mut log_writer,
        "Output BEDPE file",
        &format!("{}", bedpe_path),
    )?;

    log_info(
        &mut log_writer,
        "=== STAGE 2/3 ===",
        "Internal splice junction identification...",
    )?;
    log_info(
        &mut log_writer,
        "Processing segments",
        "Generating read-level junction paths...",
    )?;
    let mut segment_progress_log =
        |label: &str, message: &str| log_info(&mut log_writer, label, message);
    let non_bsj_segment_evidence_paths: Vec<&str> = scan2_segment_artifacts
        .non_bsj_segment_evidence_paths
        .iter()
        .map(String::as_str)
        .collect();
    let segment_summary_result = run_ciri_as(AsConfig {
        input_path,
        circ_path: &result_output,
        bsj_path: Some(&bsj_output),
        segment_evidence_paths: vec![&segments1_output, &segments2_output],
        non_bsj_segment_evidence_paths,
        out_prefix: &args.out_prefix,
        keep_temp_files: args.debug,
        reference: &fasta.chr_tcga_map,
        annotation: args.gtf.as_ref().map(|_| &annotation),
        min_mapq: args.min_mapq,
        progress_log: Some(&mut segment_progress_log),
    });
    let segment_summary = segment_summary_result?;
    log_info(
        &mut log_writer,
        "Output segments file",
        &format!("{}.segments", args.out_prefix),
    )?;
    write_and_log_segments_bam(
        &mut log_writer,
        &segments_output,
        &args.out_prefix,
        &fasta.chr_tcga_map,
        args.threads,
    )?;
    log_info(
        &mut log_writer,
        "Segments summary",
        &format!(
            "{} reads ({} BSJ, {} backward, {} outward)",
            segment_summary.total_segments,
            segment_summary.bsj_segments,
            segment_summary.backward_segments,
            segment_summary.outward_segments
        ),
    )?;

    log_info(
        &mut log_writer,
        "=== STAGE 3/3 ===",
        "Full-length isoform reconstruction...",
    )?;
    log_info(
        &mut log_writer,
        "Reconstructing isoforms",
        "Calling circRNA isoforms from splice graphs...",
    )?;
    let isoform_summary = rebuild_major_isoforms_from_segments(
        &result_output,
        &segments_output,
        &args.out_prefix,
        &fasta.chr_tcga_map,
        args.gtf.as_ref().map(|_| &annotation),
    )?;
    log_info(
        &mut log_writer,
        "Output isoform files",
        &format!(
            "{}.isoforms.gtf, {}.isoforms.fa",
            args.out_prefix, args.out_prefix
        ),
    )?;
    log_info(
        &mut log_writer,
        "Isoform summary",
        &format_isoform_summary(isoform_summary),
    )?;

    log_info(
        &mut log_writer,
        "Total runtime",
        &format!("{:.2} seconds", run_started.elapsed().as_secs_f64()),
    )?;

    if !args.debug {
        cleanup_pipeline_temp_files(
            &[
                &bsj1_output,
                &bsj2_output,
                &segments1_output,
                &segments2_output,
                &segments_non_bsj_output,
                &fsj_output,
            ],
            &scan2_segment_artifacts.non_bsj_segment_evidence_paths,
        );
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn igv_part_payloads_reconstruct_sequence_and_split_bsj_cs() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), "NNNNACGTCCCCGGGG".to_string());

        let payloads =
            igv_part_payloads_from_cs("chr1", "5-8:+|<bsj>|13-16:+", ":4<cc4gg:2*ga:1", &reference)
                .unwrap();

        assert_eq!(payloads.len(), 2);
        assert_eq!(payloads[0].seq, "ACGT");
        assert_eq!(payloads[0].cs, ":4");
        assert_eq!(payloads[1].seq, "GGAG");
        assert_eq!(payloads[1].cs, ":2*ga:1");
    }

    #[test]
    fn igv_part_payloads_split_coordinate_break_without_bsj_marker() {
        let mut reference = HashMap::new();
        reference.insert(
            "chr1".to_string(),
            "AAAACCCCGGGGTTTTAAAACCCCGGGG".to_string(),
        );

        let payloads =
            igv_part_payloads_from_cs("chr1", "21-24:+|10-13:+", ":4~gg7tt:2*gc:1", &reference)
                .unwrap();

        assert_eq!(payloads.len(), 2);
        assert_eq!(payloads[0].seq, "CCCC");
        assert_eq!(payloads[0].cs, ":4");
        assert_eq!(payloads[1].seq, "GGCT");
        assert_eq!(payloads[1].cs, ":2*gc:1");
    }

    #[test]
    fn cs_consumption_must_match_synthetic_cigar() {
        assert!(cs_consumption_matches_cigar(":4~gt3ag:2", "4M3N2M"));
        assert!(!cs_consumption_matches_cigar(":2+a:2", "4M"));
        assert!(!cs_consumption_matches_cigar(":2-a:2", "4M"));
    }

    #[test]
    fn cigar_from_part_cs_preserves_indels() {
        assert_eq!(cigar_from_part_cs(":2+a:1*ag-g:3").unwrap(), "2M1I2M1D3M");
        assert_eq!(cigar_from_part_cs(":4~gt3ag:2").unwrap(), "4M3N2M");
        assert!(cigar_from_part_cs(":4<gt3ag:2").is_none());
    }

    #[test]
    fn synthetic_segments_sam_line_writes_real_seq_and_cs_when_valid() {
        let payload = IgvPartPayload {
            seq: "ACGT".to_string(),
            cs: ":4".to_string(),
        };
        let line = synthetic_segments_sam_line(
            "chr1",
            4,
            "read1|bsj|circ|R1|part1",
            "bsj",
            "circ",
            "R1",
            1,
            IGV_BSJ_COLOR,
            '+',
            &[(5, 8, '+')],
            Some(&payload),
        )
        .unwrap();

        assert!(line.contains("\t4M\t"));
        assert!(line.contains("\tACGT\t"));
        assert!(line.contains("\tcs:Z::4\n"));
    }

    #[test]
    fn synthetic_segments_sam_line_falls_back_when_cs_and_cigar_disagree() {
        let payload = IgvPartPayload {
            seq: "AACGT".to_string(),
            cs: "+a:4".to_string(),
        };
        let line = synthetic_segments_sam_line(
            "chr1",
            4,
            "read1|bsj|circ|R1|part1",
            "bsj",
            "circ",
            "R1",
            1,
            IGV_BSJ_COLOR,
            '+',
            &[(5, 8, '+')],
            Some(&payload),
        )
        .unwrap();

        assert!(line.contains("\t1I4M\t"));
        assert!(line.contains("\tAACGT\t"));
        assert!(line.contains("\tcs:Z:+a:4\n"));
    }

    #[test]
    fn synthetic_segments_sam_line_falls_back_when_cs_ref_span_disagrees() {
        let payload = IgvPartPayload {
            seq: "ACGT".to_string(),
            cs: ":4-a".to_string(),
        };
        let line = synthetic_segments_sam_line(
            "chr1",
            4,
            "read1|bsj|circ|R1|part1",
            "bsj",
            "circ",
            "R1",
            1,
            IGV_BSJ_COLOR,
            '+',
            &[(5, 8, '+')],
            Some(&payload),
        )
        .unwrap();

        assert!(line.contains("\t4M\t"));
        assert!(line.contains("\tNNNN\t"));
        assert!(!line.contains("\tcs:Z:"));
    }
}

/// Removes internal stage sidecars after the final user-facing outputs exist.
///
/// Official CLI outputs should remain focused on the final result files. The
/// intermediate `.bsj1/.bsj2/.segments1/.segments2` stems are implementation
/// details kept only when `--debug` is set, while shard paths are removed here
/// because Scan1/Scan2/segments may leave either merged stems or `.part_*.tmp`
/// files depending on which fast path was active.
fn cleanup_pipeline_temp_files(stems: &[&str], extra_paths: &[String]) {
    for path in extra_paths {
        let _ = std::fs::remove_file(path);
    }
    for stem in stems {
        let _ = std::fs::remove_file(stem);
        for idx in 0..rayon::current_num_threads().max(1) {
            let _ = std::fs::remove_file(crate::utils::part_path(stem, idx));
        }
    }
}
