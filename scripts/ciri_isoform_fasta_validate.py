#!/usr/bin/env python3
"""Validate CIRI isoform FASTA sequences against GTF exon chains.

The primary check is deterministic and does not require a mapper: every FASTA
record is reconstructed from `<prefix>.isoforms.gtf` and the reference FASTA,
including reverse-complementing negative-strand transcripts. An optional BLAT
sample check can be enabled to audit genomic alignment strand without making
full BLAT a routine dependency.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    """Return the reverse-complemented DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1].upper()


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    """Read FASTA records as `{id: (full_header, sequence)}`."""
    records: dict[str, tuple[str, str]] = {}
    name: str | None = None
    header = ""
    seq: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = (header, "".join(seq).upper())
                header = line[1:]
                name = header.split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        records[name] = (header, "".join(seq).upper())
    return records


def parse_fasta_header_attrs(header: str) -> dict[str, str]:
    """Parse whitespace-separated `key=value` fields from one FASTA header."""
    attrs: dict[str, str] = {}
    for token in header.split()[1:]:
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        attrs[key] = value
    return attrs


def parse_gtf_attrs(text: str) -> dict[str, str]:
    """Parse one GTF attribute column into a dictionary."""
    attrs: dict[str, str] = {}
    for item in text.strip().split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def read_isoform_gtf(
    path: Path,
) -> tuple[dict[str, dict[str, str]], dict[str, list[tuple[int, int]]]]:
    """Read transcript metadata and exon chains from CIRI isoform GTF."""
    transcripts: dict[str, dict[str, str]] = {}
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = parse_gtf_attrs(fields[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue
            if fields[2] == "transcript":
                transcripts[transcript_id] = {
                    "chrom": fields[0],
                    "strand": fields[6],
                    "type": attrs.get("type", "NA"),
                    "evidence": attrs.get("evidence", "NA"),
                    "isoform_len": attrs.get("isoform_len", "NA"),
                    "score": attrs.get("score", "NA"),
                }
            elif fields[2] == "exon":
                exons[transcript_id].append((int(fields[3]), int(fields[4])))
    return transcripts, exons


def expected_sequence(
    transcript_id: str,
    transcript: dict[str, str],
    exons: dict[str, list[tuple[int, int]]],
    reference: dict[str, tuple[str, str]],
) -> str:
    """Reconstruct one transcript sequence from sorted GTF exons."""
    chrom = transcript["chrom"]
    if chrom not in reference:
        raise KeyError(f"{transcript_id}: chromosome {chrom!r} not found in reference")
    pieces = []
    ref_seq = reference[chrom][1]
    for start, end in sorted(exons.get(transcript_id, [])):
        pieces.append(ref_seq[start - 1 : end])
    seq = "".join(pieces).upper()
    if transcript["strand"] == "-":
        seq = reverse_complement(seq)
    return seq


def expected_cirexon(transcript: dict[str, str], exon_chain: list[tuple[int, int]]) -> str:
    """Return the compact `cirexon` header value expected from GTF exons."""
    strand = transcript["strand"]
    return ",".join(f"{start}-{end}:{strand}" for start, end in sorted(exon_chain))


def stratified_sample(
    fasta_records: dict[str, tuple[str, str]],
    transcripts: dict[str, dict[str, str]],
    sample_size: int,
) -> list[str]:
    """Select a deterministic sample balanced by strand and mature/estimate type."""
    names = sorted(fasta_records)
    if sample_size <= 0 or sample_size >= len(names):
        return names
    buckets: dict[tuple[str, str], list[str]] = defaultdict(list)
    for name in names:
        meta = transcripts.get(name)
        if meta is None:
            continue
        buckets[(meta["strand"], meta["type"])].append(name)
    selected: list[str] = []
    keys = sorted(buckets)
    base = sample_size // max(1, len(keys))
    remainder = sample_size % max(1, len(keys))
    for idx, key in enumerate(keys):
        limit = base + (1 if idx < remainder else 0)
        values = buckets[key]
        if limit >= len(values):
            selected.extend(values)
            continue
        if limit == 1:
            selected.append(values[len(values) // 2])
            continue
        step = (len(values) - 1) / (limit - 1)
        selected.extend(values[round(i * step)] for i in range(limit))
    return sorted(set(selected), key=names.index)


def write_fasta_subset(
    path: Path,
    names: list[str],
    fasta_records: dict[str, tuple[str, str]],
) -> None:
    """Write selected FASTA records for BLAT sampling."""
    with path.open("w") as handle:
        for name in names:
            header, seq = fasta_records[name]
            handle.write(f">{header}\n")
            for idx in range(0, len(seq), 80):
                handle.write(seq[idx : idx + 80] + "\n")


def run_blat_sample(
    blat: Path,
    reference: Path,
    fasta_records: dict[str, tuple[str, str]],
    transcripts: dict[str, dict[str, str]],
    sample_size: int,
    chunk_size: int,
    min_identity: int,
    max_examples: int,
) -> tuple[Counter[str], list[str]]:
    """Run BLAT on a deterministic FASTA sample and summarize best-hit strand."""
    selected = stratified_sample(fasta_records, transcripts, sample_size)
    stats: Counter[str] = Counter()
    examples: list[str] = []
    with tempfile.TemporaryDirectory(prefix="ciri_isoform_blat_") as tmp:
        tmp_dir = Path(tmp)
        psl_path = tmp_dir / "sample.psl"
        for chunk_idx in range(0, len(selected), chunk_size):
            chunk_names = selected[chunk_idx : chunk_idx + chunk_size]
            query_path = tmp_dir / f"chunk_{chunk_idx // chunk_size:04d}.fa"
            chunk_psl = tmp_dir / f"chunk_{chunk_idx // chunk_size:04d}.psl"
            write_fasta_subset(query_path, chunk_names, fasta_records)
            cmd = [
                str(blat),
                "-t=dna",
                "-q=rna",
                f"-minIdentity={min_identity}",
                "-minScore=30",
                "-noHead",
                str(reference),
                str(query_path),
                str(chunk_psl),
            ]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
            with chunk_psl.open() as src, psl_path.open("a") as dst:
                shutil.copyfileobj(src, dst)

        best: dict[str, tuple[int, float, float, str, int, int, int, int]] = {}
        malformed = 0
        if psl_path.exists():
            with psl_path.open() as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 21:
                        malformed += 1
                        continue
                    matches = int(fields[0])
                    mismatches = int(fields[1])
                    rep_matches = int(fields[2])
                    query_name = fields[9]
                    query_size = int(fields[10])
                    query_start = int(fields[11])
                    query_end = int(fields[12])
                    strand = fields[8][0]
                    target_start = int(fields[15])
                    target_end = int(fields[16])
                    score = matches + rep_matches - mismatches
                    coverage = (query_end - query_start) / query_size if query_size else 0.0
                    identity = (
                        matches / (matches + mismatches)
                        if matches + mismatches
                        else 0.0
                    )
                    record = (
                        score,
                        coverage,
                        identity,
                        strand,
                        matches,
                        mismatches,
                        target_start,
                        target_end,
                    )
                    if query_name not in best or record[:3] > best[query_name][:3]:
                        best[query_name] = record
        stats["sampled"] = len(selected)
        stats["best_hits"] = len(best)
        stats["malformed_psl_rows"] = malformed
        for name in selected:
            meta = transcripts.get(name)
            record = best.get(name)
            if meta is None:
                continue
            if record is None:
                stats["no_hit"] += 1
                if len(examples) < max_examples:
                    examples.append(f"{name}: no BLAT hit")
                continue
            _score, coverage, identity, strand, matches, mismatches, t_start, t_end = record
            if strand == meta["strand"]:
                stats["strand_match"] += 1
            else:
                stats["strand_mismatch"] += 1
                if len(examples) < max_examples:
                    examples.append(
                        f"{name}: GTF strand {meta['strand']} but BLAT best hit {strand}"
                    )
            if coverage >= 0.99:
                stats["coverage_ge_99"] += 1
            else:
                stats["coverage_lt_99"] += 1
                if len(examples) < max_examples:
                    examples.append(
                        f"{name}: BLAT coverage {coverage:.3f}, identity {identity:.4f}, "
                        f"matches {matches}, mismatches {mismatches}, target {t_start}-{t_end}"
                    )
            if identity >= 0.999:
                stats["identity_ge_999"] += 1
            else:
                stats["identity_lt_999"] += 1
                if len(examples) < max_examples:
                    examples.append(
                        f"{name}: BLAT identity {identity:.4f}, coverage {coverage:.3f}"
                    )
    return stats, examples


def validate(args: argparse.Namespace) -> int:
    """Run deterministic FASTA validation and optional BLAT strand auditing."""
    fasta_records = read_fasta(args.fasta)
    reference = read_fasta(args.reference)
    transcripts, exons = read_isoform_gtf(args.gtf)

    stats: Counter[str] = Counter()
    by_type: Counter[str] = Counter()
    by_strand: Counter[str] = Counter()
    by_evidence: Counter[str] = Counter()
    examples: list[str] = []

    for name, (header, seq) in fasta_records.items():
        transcript = transcripts.get(name)
        if transcript is None:
            stats["missing_gtf_transcript"] += 1
            if len(examples) < args.max_examples:
                examples.append(f"{name}: missing transcript row in GTF")
            continue
        by_type[transcript["type"]] += 1
        by_strand[transcript["strand"]] += 1
        for evidence in transcript["evidence"].split(","):
            by_evidence[evidence] += 1
        try:
            expected = expected_sequence(name, transcript, exons, reference)
        except KeyError as exc:
            stats["missing_reference_chrom"] += 1
            if len(examples) < args.max_examples:
                examples.append(str(exc))
            continue
        if seq == expected:
            stats["sequence_exact"] += 1
        else:
            exon_chain = sorted(exons.get(name, []))
            chrom = transcript["chrom"]
            plus_seq = "".join(reference[chrom][1][start - 1 : end] for start, end in exon_chain)
            opposite = plus_seq.upper()
            if transcript["strand"] == "+":
                opposite = reverse_complement(opposite)
            if seq == opposite:
                stats["opposite_strand_sequence"] += 1
            else:
                stats["sequence_mismatch"] += 1
            if len(examples) < args.max_examples:
                examples.append(
                    f"{name}: sequence mismatch, observed len {len(seq)}, expected len {len(expected)}"
                )
        isoform_len = transcript["isoform_len"]
        if isoform_len != "NA" and len(seq) == int(isoform_len):
            stats["gtf_len_match"] += 1
        else:
            stats["gtf_len_mismatch"] += 1
            if len(examples) < args.max_examples:
                examples.append(
                    f"{name}: FASTA length {len(seq)} != GTF isoform_len {isoform_len}"
                )
        header_attrs = parse_fasta_header_attrs(header)
        if "len" in header_attrs:
            if len(seq) == int(header_attrs["len"]):
                stats["header_len_match"] += 1
            else:
                stats["header_len_mismatch"] += 1
        if "type" in header_attrs:
            if header_attrs["type"] == transcript["type"]:
                stats["header_type_match"] += 1
            else:
                stats["header_type_mismatch"] += 1
        if "evidence" in header_attrs:
            if header_attrs["evidence"] == transcript["evidence"]:
                stats["header_evidence_match"] += 1
            else:
                stats["header_evidence_mismatch"] += 1
        if "cirexon" in header_attrs:
            if header_attrs["cirexon"] == expected_cirexon(transcript, exons.get(name, [])):
                stats["header_cirexon_match"] += 1
            else:
                stats["header_cirexon_mismatch"] += 1

    print(f"FASTA records: {len(fasta_records)}")
    print(f"GTF transcripts: {len(transcripts)}")
    print(f"Strand counts in FASTA: {dict(sorted(by_strand.items()))}")
    print(f"Type counts in FASTA: {dict(sorted(by_type.items()))}")
    print(f"Evidence counts in FASTA: {dict(sorted(by_evidence.items()))}")
    print(f"Sequence validation: {dict(sorted(stats.items()))}")

    exit_code = 0
    failure_keys = [
        "missing_gtf_transcript",
        "missing_reference_chrom",
        "opposite_strand_sequence",
        "sequence_mismatch",
        "gtf_len_mismatch",
        "header_evidence_mismatch",
        "header_len_mismatch",
        "header_type_mismatch",
        "header_cirexon_mismatch",
    ]
    if any(stats[key] for key in failure_keys):
        exit_code = 1

    if args.blat is not None:
        blat_stats, blat_examples = run_blat_sample(
            args.blat,
            args.reference,
            fasta_records,
            transcripts,
            args.blat_sample_size,
            args.blat_chunk_size,
            args.blat_min_identity,
            args.max_examples,
        )
        print(f"BLAT sample validation: {dict(sorted(blat_stats.items()))}")
        examples.extend(blat_examples)
        if blat_stats["strand_mismatch"] or blat_stats["malformed_psl_rows"]:
            exit_code = 1

    if examples:
        print("Examples:")
        for example in examples[: args.max_examples]:
            print(f"  - {example}")
    return exit_code


def build_arg_parser() -> argparse.ArgumentParser:
    """Build CLI parser for the FASTA validation script."""
    parser = argparse.ArgumentParser(
        description="Validate CIRI isoform FASTA against GTF exon chains and reference."
    )
    parser.add_argument("--fasta", required=True, type=Path, help="CIRI isoforms FASTA")
    parser.add_argument("--gtf", required=True, type=Path, help="CIRI isoforms GTF")
    parser.add_argument("--reference", required=True, type=Path, help="Reference FASTA")
    parser.add_argument(
        "--blat",
        type=Path,
        help="Optional BLAT executable for sampled genomic strand validation",
    )
    parser.add_argument(
        "--blat-sample-size",
        type=int,
        default=100,
        help="Number of FASTA records to sample when --blat is set",
    )
    parser.add_argument(
        "--blat-chunk-size",
        type=int,
        default=100,
        help="Records per BLAT invocation when --blat is set",
    )
    parser.add_argument(
        "--blat-min-identity",
        type=int,
        default=95,
        help="BLAT -minIdentity value for sampled checks",
    )
    parser.add_argument(
        "--max-examples",
        type=int,
        default=10,
        help="Maximum diagnostic examples to print",
    )
    return parser


def main() -> int:
    """CLI entry point."""
    parser = build_arg_parser()
    args = parser.parse_args()
    return validate(args)


if __name__ == "__main__":
    sys.exit(main())
