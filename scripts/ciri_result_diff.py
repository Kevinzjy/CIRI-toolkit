#!/usr/bin/env python3
"""Compare two CIRI `.out`/`.result` files at circ/read/assignment/FSJ levels.

Use this when validating baseline CIRI3 parity between a reference result
(`Java`, old Rust baseline, etc.) and the current Rust output. The script is not
for segments truth evaluation; see `ciri_segments_eval.py` for that.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def load_result(
    path: Path,
) -> tuple[set[str], set[str], set[tuple[str, str]], dict[str, int]]:
    """Load circRNA IDs, read IDs, assignments, and per-circ FSJ counts."""
    circ_ids: set[str] = set()
    read_ids: set[str] = set()
    read_assignments: set[tuple[str, str]] = set()
    fsj_counts: dict[str, int] = {}

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
            if len(fields) > 6:
                fsj_counts[circ_id] = int(fields[6])
            if len(fields) > 11 and fields[11]:
                for read_id in fields[11].split(","):
                    if not read_id:
                        continue
                    read_ids.add(read_id)
                    read_assignments.add((circ_id, read_id))

    return circ_ids, read_ids, read_assignments, fsj_counts


def print_pr(prefix: str, common_n: int, pred_n: int, truth_n: int) -> None:
    """Print precision/recall for one comparison layer."""
    print(f"{prefix}_java={truth_n}")
    print(f"{prefix}_rust={pred_n}")
    print(f"{prefix}_common={common_n}")
    print(f"{prefix}_only_java={truth_n - common_n}")
    print(f"{prefix}_only_rust={pred_n - common_n}")
    if pred_n > 0:
        print(f"{prefix}_precision={common_n / pred_n:.4f}")
    if truth_n > 0:
        print(f"{prefix}_recall={common_n / truth_n:.4f}")


def print_fsj_stats(
    java_fsj: dict[str, int],
    rust_fsj: dict[str, int],
    shared_circs: set[str],
) -> list[tuple[str, int, int, int]]:
    """Print FSJ agreement statistics for shared circRNAs."""
    diffs = []
    for circ_id in sorted(shared_circs):
        java_count = java_fsj.get(circ_id, 0)
        rust_count = rust_fsj.get(circ_id, 0)
        if java_count != rust_count:
            diffs.append((circ_id, java_count, rust_count, rust_count - java_count))

    print(f"fsj_shared_circ={len(shared_circs)}")
    print(f"fsj_same={len(shared_circs) - len(diffs)}")
    print(f"fsj_diff={len(diffs)}")
    print(f"fsj_java_total_shared={sum(java_fsj.get(circ_id, 0) for circ_id in shared_circs)}")
    print(f"fsj_rust_total_shared={sum(rust_fsj.get(circ_id, 0) for circ_id in shared_circs)}")
    return diffs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare two CIRI .out/.result files for circ/read/assignment/FSJ parity."
    )
    parser.add_argument("java_result", type=Path, help="Reference result file (e.g. Java)")
    parser.add_argument("rust_result", type=Path, help="Current result file (e.g. Rust)")
    parser.add_argument(
        "--show-ids",
        action="store_true",
        help="Print full circRNA ID lists for only_java / only_rust.",
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
    parser.add_argument(
        "--show-fsj-diff",
        action="store_true",
        help="Print per-circ FSJ differences for shared circRNAs.",
    )
    args = parser.parse_args()

    java_circ_ids, java_read_ids, java_assignments, java_fsj = load_result(args.java_result)
    rust_circ_ids, rust_read_ids, rust_assignments, rust_fsj = load_result(args.rust_result)

    circ_common = java_circ_ids & rust_circ_ids
    circ_only_java = java_circ_ids - rust_circ_ids
    circ_only_rust = rust_circ_ids - java_circ_ids

    read_common = java_read_ids & rust_read_ids
    read_only_java = java_read_ids - rust_read_ids
    read_only_rust = rust_read_ids - java_read_ids

    assignment_common = java_assignments & rust_assignments
    assignment_only_java = java_assignments - rust_assignments
    assignment_only_rust = rust_assignments - java_assignments

    print_pr("circ", len(circ_common), len(rust_circ_ids), len(java_circ_ids))
    print_pr("read", len(read_common), len(rust_read_ids), len(java_read_ids))
    print_pr(
        "read_assignment",
        len(assignment_common),
        len(rust_assignments),
        len(java_assignments),
    )
    fsj_diffs = print_fsj_stats(java_fsj, rust_fsj, circ_common)

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

    if args.show_fsj_diff:
        print("FSJ_DIFF")
        for circ_id, java_count, rust_count, delta in fsj_diffs:
            print(f"{circ_id}\t{java_count}\t{rust_count}\t{delta}")


if __name__ == "__main__":
    main()
