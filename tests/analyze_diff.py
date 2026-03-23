#!/usr/bin/env python3
"""Quickly compare CIRI results at circRNA and read levels."""

from __future__ import annotations

import argparse
from pathlib import Path


def load_result_sets(path: Path) -> tuple[set[str], set[str], set[tuple[str, str]]]:
    """Load circRNA IDs, read IDs, and (circRNA_ID, read_ID) assignments."""
    circ_ids: set[str] = set()
    read_ids: set[str] = set()
    read_assignments: set[tuple[str, str]] = set()

    with path.open("r", encoding="utf-8") as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            circ_id = fields[0]
            circ_ids.add(circ_id)
            if len(fields) > 11 and fields[11]:
                for read_id in fields[11].split(","):
                    if not read_id:
                        continue
                    read_ids.add(read_id)
                    read_assignments.add((circ_id, read_id))

    return circ_ids, read_ids, read_assignments


def print_pr(prefix: str, common_n: int, pred_n: int, truth_n: int) -> None:
    """Print precision/recall for a label prefix."""
    print(f"{prefix}_java={truth_n}")
    print(f"{prefix}_rust={pred_n}")
    print(f"{prefix}_common={common_n}")
    print(f"{prefix}_only_java={truth_n - common_n}")
    print(f"{prefix}_only_rust={pred_n - common_n}")
    if pred_n > 0:
        print(f"{prefix}_precision={common_n / pred_n:.4f}")
    if truth_n > 0:
        print(f"{prefix}_recall={common_n / truth_n:.4f}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare two CIRI .result files.")
    parser.add_argument("java_result", type=Path, help="Reference result file (e.g. Java)")
    parser.add_argument("rust_result", type=Path, help="Current result file (e.g. Rust)")
    parser.add_argument(
        "--show-ids",
        action="store_true",
        help="Print full ID lists for only_java / only_rust.",
    )
    parser.add_argument(
        "--show-read-ids",
        action="store_true",
        help="Print full read ID lists for only_java / only_rust.",
    )
    parser.add_argument(
        "--show-read-assignments",
        action="store_true",
        help="Print full (circRNA_ID, read_ID) assignment differences.",
    )
    args = parser.parse_args()

    java_circ_ids, java_read_ids, java_assignments = load_result_sets(args.java_result)
    rust_circ_ids, rust_read_ids, rust_assignments = load_result_sets(args.rust_result)

    circ_common = java_circ_ids & rust_circ_ids
    circ_only_java = java_circ_ids - rust_circ_ids
    circ_only_rust = rust_circ_ids - java_circ_ids

    read_common = java_read_ids & rust_read_ids
    read_only_java = java_read_ids - rust_read_ids
    read_only_rust = rust_read_ids - java_read_ids

    assignment_common = java_assignments & rust_assignments
    assignment_only_java = java_assignments - rust_assignments
    assignment_only_rust = rust_assignments - java_assignments

    print_pr(
        "circ",
        common_n=len(circ_common),
        pred_n=len(rust_circ_ids),
        truth_n=len(java_circ_ids),
    )
    print_pr(
        "read",
        common_n=len(read_common),
        pred_n=len(rust_read_ids),
        truth_n=len(java_read_ids),
    )
    print_pr(
        "read_assignment",
        common_n=len(assignment_common),
        pred_n=len(rust_assignments),
        truth_n=len(java_assignments),
    )

    if args.show_ids:
        print("ONLY_JAVA_IDS")
        for item in sorted(circ_only_java):
            print(item)
        print("ONLY_RUST_IDS")
        for item in sorted(circ_only_rust):
            print(item)

    if args.show_read_ids:
        print("ONLY_JAVA_READ_IDS")
        for item in sorted(read_only_java):
            print(item)
        print("ONLY_RUST_READ_IDS")
        for item in sorted(read_only_rust):
            print(item)

    if args.show_read_assignments:
        print("ONLY_JAVA_READ_ASSIGNMENTS")
        for circ_id, read_id in sorted(assignment_only_java):
            print(f"{circ_id}\t{read_id}")
        print("ONLY_RUST_READ_ASSIGNMENTS")
        for circ_id, read_id in sorted(assignment_only_rust):
            print(f"{circ_id}\t{read_id}")


if __name__ == "__main__":
    main()
