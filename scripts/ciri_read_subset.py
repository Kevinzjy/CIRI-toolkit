#!/usr/bin/env python3
"""Build read subsets from paired CIRI `.out`/`.result` files.

Use this when parity debugging needs a compact read subset:

1. `all-bsj` mode keeps every BSJ-supporting read reported by either result.
2. `diff` mode keeps only reads implicated in current circ/read/assignment diffs,
   with optional full circRNA supporting-read context.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path


def load_result(
    path: Path,
) -> tuple[set[str], set[str], set[tuple[str, str]], dict[str, set[str]], dict[str, set[str]]]:
    """Load circ/read/assignment sets plus lookup tables from one CIRI result."""
    circ_ids: set[str] = set()
    read_ids: set[str] = set()
    assignments: set[tuple[str, str]] = set()
    circ_to_reads: dict[str, set[str]] = defaultdict(set)
    read_to_reasons: dict[str, set[str]] = defaultdict(set)

    with path.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle):
            if line_no == 0:
                continue
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            circ_id = fields[0]
            circ_ids.add(circ_id)
            if len(fields) <= 11 or not fields[11]:
                continue
            for read_id in fields[11].split(","):
                if not read_id:
                    continue
                read_ids.add(read_id)
                assignments.add((circ_id, read_id))
                circ_to_reads[circ_id].add(read_id)
                read_to_reasons[read_id].add(circ_id)

    return circ_ids, read_ids, assignments, circ_to_reads, read_to_reasons


def add_reads(reason_map: dict[str, set[str]], reads: set[str], reason: str) -> None:
    """Add a batch of reads under one inclusion reason."""
    for read_id in reads:
        reason_map[read_id].add(reason)


def collect_all_bsj_reads(
    java_read_ids: set[str],
    rust_read_ids: set[str],
) -> dict[str, set[str]]:
    """Collect every BSJ-supporting read reported by Java or Rust."""
    reason_map: dict[str, set[str]] = defaultdict(set)
    add_reads(reason_map, java_read_ids, "bsj_java")
    add_reads(reason_map, rust_read_ids, "bsj_rust")
    return reason_map


def collect_diff_reads(
    java_circ_ids: set[str],
    rust_circ_ids: set[str],
    java_read_ids: set[str],
    rust_read_ids: set[str],
    java_assignments: set[tuple[str, str]],
    rust_assignments: set[tuple[str, str]],
    java_circ_to_reads: dict[str, set[str]],
    rust_circ_to_reads: dict[str, set[str]],
    include_circ_context: bool,
) -> tuple[
    dict[str, set[str]],
    set[str],
    set[str],
    set[str],
    set[str],
    set[tuple[str, str]],
    set[tuple[str, str]],
]:
    """Collect only the reads directly involved in current Java/Rust differences."""
    circ_only_java = java_circ_ids - rust_circ_ids
    circ_only_rust = rust_circ_ids - java_circ_ids
    read_only_java = java_read_ids - rust_read_ids
    read_only_rust = rust_read_ids - java_read_ids
    assignment_only_java = java_assignments - rust_assignments
    assignment_only_rust = rust_assignments - java_assignments

    reason_map: dict[str, set[str]] = defaultdict(set)
    add_reads(reason_map, read_only_java, "read_only_java")
    add_reads(reason_map, read_only_rust, "read_only_rust")
    add_reads(reason_map, {read_id for _, read_id in assignment_only_java}, "assignment_only_java")
    add_reads(reason_map, {read_id for _, read_id in assignment_only_rust}, "assignment_only_rust")

    if include_circ_context:
        touched_java_circs = set(circ_only_java)
        touched_java_circs.update(circ_id for circ_id, _ in assignment_only_java)
        touched_rust_circs = set(circ_only_rust)
        touched_rust_circs.update(circ_id for circ_id, _ in assignment_only_rust)

        for circ_id in touched_java_circs:
            add_reads(reason_map, java_circ_to_reads.get(circ_id, set()), f"circ_context_java:{circ_id}")
        for circ_id in touched_rust_circs:
            add_reads(reason_map, rust_circ_to_reads.get(circ_id, set()), f"circ_context_rust:{circ_id}")

    return (
        reason_map,
        circ_only_java,
        circ_only_rust,
        read_only_java,
        read_only_rust,
        assignment_only_java,
        assignment_only_rust,
    )


def write_read_list(path: Path, reads: list[str]) -> None:
    """Write one read ID per line for `samtools view -N`."""
    with path.open("w", encoding="utf-8") as handle:
        for read_id in reads:
            handle.write(f"{read_id}\n")


def write_manifest(path: Path, reason_map: dict[str, set[str]]) -> None:
    """Write a TSV manifest for the selected read subset."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["read_id", "reasons"])
        for read_id in sorted(reason_map):
            writer.writerow([read_id, ",".join(sorted(reason_map[read_id]))])


def build_subset(input_alignment: Path, read_list: Path, output_alignment: Path) -> None:
    """Materialize a SAM/BAM subset with only the selected read IDs."""
    output_format = output_alignment.suffix.lower()
    cmd = ["samtools", "view", "-h", "-N", str(read_list), "-@", "16"]
    if output_format == ".bam":
        cmd.extend(["-b", "-o", str(output_alignment)])
    elif output_format == ".sam":
        cmd.extend(["-o", str(output_alignment)])
    else:
        raise ValueError("subset output must end with .bam or .sam")
    cmd.append(str(input_alignment))
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build BSJ-wide or diff-focused read subsets from two CIRI results."
    )
    parser.add_argument("java_result", type=Path, help="Reference result file (Java/CIRI3)")
    parser.add_argument("rust_result", type=Path, help="Current result file (Rust/CIRI-toolkit)")
    parser.add_argument(
        "-o",
        "--output-prefix",
        type=Path,
        required=True,
        help="Prefix for generated files; emits <prefix>.reads.txt and <prefix>.manifest.tsv",
    )
    parser.add_argument(
        "-i",
        "--input-alignment",
        type=Path,
        help="Optional source SAM/BAM. When set, also emit a subset alignment.",
    )
    parser.add_argument(
        "--subset-output",
        type=Path,
        help="Optional subset SAM/BAM path. Defaults to <prefix>.subset.bam when --input-alignment is set.",
    )
    parser.add_argument(
        "--mode",
        choices=["all-bsj", "diff"],
        default="all-bsj",
        help=(
            "Subset selection mode. 'all-bsj' keeps every BSJ-supporting read "
            "reported by either result. 'diff' keeps only directly differing "
            "reads plus optional circ context."
        ),
    )
    parser.add_argument(
        "--no-circ-context",
        action="store_true",
        help="In diff mode, skip full supporting-read context for touched circRNAs.",
    )
    args = parser.parse_args()

    (
        java_circ_ids,
        java_read_ids,
        java_assignments,
        java_circ_to_reads,
        _java_read_to_circs,
    ) = load_result(args.java_result)
    (
        rust_circ_ids,
        rust_read_ids,
        rust_assignments,
        rust_circ_to_reads,
        _rust_read_to_circs,
    ) = load_result(args.rust_result)

    circ_only_java: set[str] = set()
    circ_only_rust: set[str] = set()
    read_only_java: set[str] = set()
    read_only_rust: set[str] = set()
    assignment_only_java: set[tuple[str, str]] = set()
    assignment_only_rust: set[tuple[str, str]] = set()

    if args.mode == "all-bsj":
        reason_map = collect_all_bsj_reads(java_read_ids, rust_read_ids)
    else:
        (
            reason_map,
            circ_only_java,
            circ_only_rust,
            read_only_java,
            read_only_rust,
            assignment_only_java,
            assignment_only_rust,
        ) = collect_diff_reads(
            java_circ_ids,
            rust_circ_ids,
            java_read_ids,
            rust_read_ids,
            java_assignments,
            rust_assignments,
            java_circ_to_reads,
            rust_circ_to_reads,
            include_circ_context=not args.no_circ_context,
        )

    selected_reads = sorted(reason_map)
    read_list_path = args.output_prefix.with_suffix(".reads.txt")
    manifest_path = args.output_prefix.with_suffix(".manifest.tsv")
    write_read_list(read_list_path, selected_reads)
    write_manifest(manifest_path, reason_map)

    print(f"mode={args.mode}")
    print(f"circ_only_java={len(circ_only_java)}")
    print(f"circ_only_rust={len(circ_only_rust)}")
    print(f"read_only_java={len(read_only_java)}")
    print(f"read_only_rust={len(read_only_rust)}")
    print(f"assignment_only_java={len(assignment_only_java)}")
    print(f"assignment_only_rust={len(assignment_only_rust)}")
    print(f"selected_reads={len(selected_reads)}")
    print(f"read_list={read_list_path}")
    print(f"manifest={manifest_path}")

    if args.input_alignment is not None:
        subset_output = (
            args.subset_output
            if args.subset_output is not None
            else args.output_prefix.with_suffix(".subset.bam")
        )
        build_subset(args.input_alignment, read_list_path, subset_output)
        print(f"subset_alignment={subset_output}")


if __name__ == "__main__":
    main()
