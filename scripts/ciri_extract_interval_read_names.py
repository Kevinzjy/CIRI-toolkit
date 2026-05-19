#!/usr/bin/env python3
"""Collect read names whose alignments overlap genomic intervals.

This helper scans `samtools view -h` output sequentially, records every read
name whose mapped alignment overlaps one of the requested intervals, and writes
a plain read-name list that can be fed back into `samtools view -N`.
"""

from __future__ import annotations

import argparse
import subprocess
from collections import defaultdict
from pathlib import Path


def load_intervals(path: Path, limit: int) -> dict[str, list[tuple[int, int]]]:
    """Load BED-like TSV intervals grouped by chromosome."""
    by_chr: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as handle:
        next(handle)
        for idx, line in enumerate(handle):
            if limit and idx >= limit:
                break
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            by_chr[fields[0]].append((int(fields[1]), int(fields[2])))
    for chr_name in by_chr:
        by_chr[chr_name].sort()
    return by_chr


def overlaps(intervals: list[tuple[int, int]], start: int, end: int) -> bool:
    """Return whether a 1-based closed reference span overlaps any interval."""
    for iv_start, iv_end in intervals:
        if end < iv_start:
            return False
        if start <= iv_end and end >= iv_start:
            return True
    return False


def reference_length_from_cigar(cigar: str) -> int:
    """Return the reference-consuming CIGAR length."""
    ref_len = 0
    num = 0
    for ch in cigar:
        if ch.isdigit():
            num = num * 10 + int(ch)
        else:
            if ch in {"M", "D", "N", "=", "X"}:
                ref_len += num
            num = 0
    return ref_len


def main() -> int:
    """Parse CLI arguments and write the sorted overlapping read-name list."""
    parser = argparse.ArgumentParser(description="Extract read names overlapping interval regions.")
    parser.add_argument("bam", type=Path, help="Input BAM/SAM readable by samtools view")
    parser.add_argument("intervals", type=Path, help="TSV with chr/start/end columns and a header")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output read-name list")
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Optional number of interval rows to use from the input TSV (0 = all).",
    )
    args = parser.parse_args()

    interval_map = load_intervals(args.intervals, args.limit)
    wanted: set[str] = set()

    proc = subprocess.Popen(["samtools", "view", "-h", str(args.bam)], stdout=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        if not line or line[0] == "@":
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue
        flag = int(fields[1])
        if flag & 0x4:
            continue
        chr_name = fields[2]
        if chr_name not in interval_map:
            continue
        pos = int(fields[3])
        end = pos + max(reference_length_from_cigar(fields[5]) - 1, 0)
        if overlaps(interval_map[chr_name], pos, end):
            wanted.add(fields[0])
    ret = proc.wait()
    if ret != 0:
        return ret

    with args.output.open("w", encoding="utf-8") as handle:
        for read_id in sorted(wanted):
            handle.write(f"{read_id}\n")
    print(args.output)
    print(f"read_count={len(wanted)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
