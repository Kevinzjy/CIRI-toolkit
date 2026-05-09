#!/usr/bin/env python3
"""Evaluate `<prefix>.segments` against simulator `*.reads.tsv` truth.

Use this after running the default post-Summary segments phase on simulator
fixtures. The script focuses on the current read-level contract and prints two
extra views that are useful during the current development stage:

- predicted `type` (`bsj` / `backward` / `outward`)
- `circ_id`
- `is_circular`
- mate-level genomic segment sets from `r1_segments` / `r2_segments`
- mate-level read-chain segment strings from `r1_segments` / `r2_segments`
- mate-level junction chains derived from retained read-chain segments
- mate-level `is_r1_bsj` / `is_r2_bsj` for confirmed `type=bsj` rows only

Additional summaries:

- `BSJ_ONLY`: read-pair and mate-level metrics restricted to predicted `type=bsj`
- `BACKWARD_TRUTH_COMPOSITION`: how predicted `type=backward` rows overlap
  simulator truth. In backward rows, `<bsj>` / `B` is treated only as a topology
  marker and chain break, not as confirmed mate-level BSJ evidence.
- `OUTWARD_TRUTH_COMPOSITION`: how predicted `type=outward` rows overlap
  simulator truth. Outward rows carry pair-level circular topology support but
  do not encode a mate-level BSJ junction.

The default output format is a compact human-readable report. Use `--format tsv`
to emit tidy rows that are easier to redirect into downstream shell tooling.
By default, segment-set comparison ignores fragments shorter than 10 bp because
these short clipped pieces are not treated as high-confidence CIRI-AS splice
templates.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


PAIR_KEYS_BY_TYPE = {
    "bsj": [
        "circ_id",
        "is_circular",
    ],
    "backward": [
        "is_circular",
    ],
    "outward": [
        "is_circular",
    ],
}
PRED_TYPES = ("bsj", "backward", "outward")

MATE_NAMES = ("r1", "r2")
SEGMENT_ERROR_CLASSES = (
    "exact",
    "count_diff",
    "shift_1bp",
    "shift_2bp",
    "shift_gt2",
    "strand_diff",
    "other",
)


def pred_mate_bsj_key(mate: str) -> str:
    """Return the current `<prefix>.segments` mate-level BSJ column name."""
    return f"is_{mate}_bsj"


def truth_mate_bsj_key(mate: str) -> str:
    """Return the simulator truth mate-level BSJ column name."""
    return f"{mate}_is_bsj"


def example_keys_for_type(pred_type: str) -> list[str]:
    """Return example metric keys that are meaningful for one predicted type.

    `type=backward` rows can contain `<bsj>` / `B` markers and `type=outward`
    rows carry pair-level circular topology, but neither marker is
    Summary-confirmed mate BSJ evidence. For that reason, non-BSJ examples
    deliberately exclude `is_r*_bsj`.
    """
    keys = (
        [*PAIR_KEYS_BY_TYPE[pred_type]]
        + [f"{mate}_segment_set" for mate in MATE_NAMES]
        + [f"{mate}_segment_chain" for mate in MATE_NAMES]
        + [f"{mate}_junction_chain" for mate in MATE_NAMES]
    )
    if pred_type == "bsj":
        keys += [pred_mate_bsj_key(mate) for mate in MATE_NAMES]
    return keys


def parse_segment_token(token: str) -> tuple[int, int, str] | None:
    """Parse one `start-end:strand` segment token."""
    try:
        span, strand = token.rsplit(":", 1)
        start_text, end_text = span.split("-", 1)
        start = int(start_text)
        end = int(end_text)
    except ValueError:
        return None
    return start, end, strand


def segment_len(segment: tuple[int, int, str]) -> int:
    """Return the inclusive genomic length of one parsed segment."""
    return segment[1] - segment[0] + 1


def infer_truth_type(row: dict[str, str]) -> str:
    """Map simulator truth labels to the current segments `type` contract."""
    if row["is_bsj"] == "1":
        return "bsj"
    if row["is_circular"] == "1":
        return "backward"
    return "forward"


def parse_segments(text: str, min_seg_len: int) -> tuple[tuple[int, int, str], ...]:
    """Parse and length-filter one mate's segment tokens."""
    if text == "NA":
        return tuple()
    segments = []
    for token in text.split("|"):
        if not token or token == "<bsj>":
            continue
        segment = parse_segment_token(token)
        if segment is None:
            continue
        if segment_len(segment) >= min_seg_len:
            segments.append(segment)
    return tuple(sorted(segments))


def normalize_segment_chain(text: str, min_seg_len: int) -> tuple[str, ...]:
    """Normalize one mate's segment tokens while preserving read-chain order.

    Short fragments are filtered with the same threshold used by segment-set
    metrics. `<bsj>` is retained only when it separates two retained segment
    tokens, so a filtered short fragment next to the BSJ does not leave dangling
    or duplicate marker tokens.
    """
    if text == "NA":
        return tuple()
    chain: list[str] = []
    pending_bsj = False
    for token in text.split("|"):
        if not token:
            continue
        if token == "<bsj>":
            pending_bsj = bool(chain)
            continue
        segment = parse_segment_token(token)
        if segment is None or segment_len(segment) < min_seg_len:
            continue
        if pending_bsj and chain and chain[-1] != "<bsj>":
            chain.append("<bsj>")
        start, end, strand = segment
        chain.append(f"{start}-{end}:{strand}")
        pending_bsj = False
    if chain and chain[-1] == "<bsj>":
        chain.pop()
    return tuple(chain)


def normalize_segment_set(text: str, min_seg_len: int) -> tuple[str, ...]:
    """Normalize one mate's segment tokens into an order-insensitive tuple.

    `<bsj>` is intentionally excluded here because BSJ presence is tracked by the
    dedicated `is_r1_bsj` / `is_r2_bsj` fields. The current evaluation goal is
    to ask whether the same genomic segments were recovered, regardless of read
    order or mapper-specific token order.
    """
    return tuple(
        f"{start}-{end}:{strand}"
        for start, end, strand in parse_segments(text, min_seg_len)
    )


def normalize_junction_chain(text: str, min_seg_len: int) -> tuple[str, ...]:
    """Normalize one mate into a read-chain junction path.

    Full-length reconstruction should eventually use splice/BSJ connections as
    graph edges, not require every terminal exon boundary to be exact. This
    metric therefore keeps only connections between retained segments: linear
    junctions are represented by their donor/acceptor coordinates, while a BSJ
    marker is represented as `B` without re-checking the adjacent segment ends.
    """
    chain = normalize_segment_chain(text, min_seg_len)
    if not chain:
        return tuple()

    junctions: list[str] = []
    previous_segment: tuple[int, int, str] | None = None
    pending_bsj = False
    for token in chain:
        if token == "<bsj>":
            pending_bsj = previous_segment is not None
            continue

        segment = parse_segment_token(token)
        if segment is None:
            previous_segment = None
            pending_bsj = False
            continue

        if previous_segment is not None:
            if pending_bsj:
                junctions.append("B")
            else:
                _prev_start, prev_end, prev_strand = previous_segment
                start, _end, strand = segment
                if prev_strand == strand:
                    junctions.append(f"N:{prev_end}>{start}:{strand}")
                else:
                    junctions.append(f"N:{prev_end}>{start}:{prev_strand}/{strand}")
        previous_segment = segment
        pending_bsj = False
    return tuple(junctions)


def classify_segment_error(
    pred_segments: tuple[tuple[int, int, str], ...],
    truth_segments: tuple[tuple[int, int, str], ...],
) -> str:
    """Classify one segment-set comparison after length filtering."""
    if pred_segments == truth_segments:
        return "exact"
    if len(pred_segments) != len(truth_segments):
        return "count_diff"
    if not pred_segments:
        return "other"

    pred_no_strand = tuple((start, end) for start, end, _strand in pred_segments)
    truth_no_strand = tuple((start, end) for start, end, _strand in truth_segments)
    if pred_no_strand == truth_no_strand:
        return "strand_diff"

    max_shift = max(
        max(abs(pred_start - truth_start), abs(pred_end - truth_end))
        for (pred_start, pred_end), (truth_start, truth_end) in zip(
            pred_no_strand, truth_no_strand
        )
    )
    if max_shift <= 1:
        return "shift_1bp"
    if max_shift <= 2:
        return "shift_2bp"
    return "shift_gt2"


def load_truth(path: Path) -> tuple[dict[str, dict[str, str]], Counter[str]]:
    """Load simulator read-level truth keyed by read ID."""
    truth: dict[str, dict[str, str]] = {}
    counts: Counter[str] = Counter()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row["type"] = infer_truth_type(row)
            truth[row["read_id"]] = row
            counts[row["type"]] += 1
    return truth, counts


def emit_summary_row_tsv(
    section: str,
    metric: str,
    matched: int,
    total: int | None = None,
) -> None:
    """Emit one tidy TSV summary row."""
    rate = "NA" if total in (None, 0) else f"{matched / total:.4f}"
    total_text = "NA" if total is None else str(total)
    print(f"summary\t{section}\t{metric}\t{matched}\t{total_text}\t{rate}")


def emit_example_row_tsv(
    section: str,
    metric: str,
    read_id: str,
    pred_value: str,
    truth_value: str,
) -> None:
    """Emit one tidy TSV example row."""
    print(
        f"example\t{section}\t{metric}\t{read_id}\t{pred_value}\t{truth_value}"
    )


def format_rate(matched: int, total: int) -> str:
    """Return a guarded floating-point rate string."""
    if total == 0:
        return "NA"
    return f"{matched / total:.4f}"


def print_table(title: str, header: list[str], rows: list[list[str]]) -> None:
    """Print a simple aligned text table."""
    print(title)
    widths = [len(col) for col in header]
    for row in rows:
        for idx, cell in enumerate(row):
            widths[idx] = max(widths[idx], len(cell))
    print("  ".join(cell.ljust(widths[idx]) for idx, cell in enumerate(header)))
    for row in rows:
        print("  ".join(cell.ljust(widths[idx]) for idx, cell in enumerate(row)))
    print()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate CIRI `<prefix>.segments` against simulator `*.reads.tsv` truth."
    )
    parser.add_argument("truth_reads_tsv", type=Path, help="Simulator truth `.reads.tsv` file")
    parser.add_argument("pred_segments_tsv", type=Path, help="CIRI `<prefix>.segments` file")
    parser.add_argument(
        "--examples",
        type=int,
        default=5,
        help="Maximum number of mismatch examples to print per category (default: 5).",
    )
    parser.add_argument(
        "--format",
        choices=["text", "tsv"],
        default="text",
        help="Output format. 'text' prints report-style tables; 'tsv' is tidy and shell-friendly.",
    )
    parser.add_argument(
        "--min-seg-len",
        type=int,
        default=10,
        help="Minimum inclusive segment length used in segment-set metrics (default: 10).",
    )
    args = parser.parse_args()

    if args.format == "tsv":
        emit_summary = emit_summary_row_tsv
        emit_example = emit_example_row_tsv
        print("kind\tsection\tmetric\tvalue_or_read_id\ttotal_or_pred\ttruth_or_rate")

    truth, truth_counts = load_truth(args.truth_reads_tsv)
    pred_counts: Counter[str] = Counter()
    confusion: Counter[tuple[str, str]] = Counter()
    exact: Counter[tuple[str, str]] = Counter()
    segment_errors: Counter[tuple[str, str, str]] = Counter()
    examples: dict[tuple[str, str], list[tuple[str, str, str]]] = defaultdict(list)
    seen_pred: set[str] = set()

    with args.pred_segments_tsv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            read_id = row["read_id"]
            pred_type = row["type"]
            pred_counts[pred_type] += 1
            seen_pred.add(read_id)

            truth_row = truth.get(read_id)
            if truth_row is None:
                confusion[(pred_type, "missing_truth")] += 1
                continue

            truth_type = truth_row["type"]
            confusion[(pred_type, truth_type)] += 1

            for key in PAIR_KEYS_BY_TYPE.get(pred_type, []):
                if row[key] == truth_row[key]:
                    exact[(pred_type, key)] += 1
                elif len(examples[(pred_type, key)]) < args.examples:
                    examples[(pred_type, key)].append((read_id, row[key], truth_row[key]))
            for mate in MATE_NAMES:
                set_key = f"{mate}_segment_set"
                chain_key = f"{mate}_segment_chain"
                junction_key = f"{mate}_junction_chain"
                pred_segments = parse_segments(row[f"{mate}_segments"], args.min_seg_len)
                truth_segments = parse_segments(
                    truth_row[f"{mate}_segments"], args.min_seg_len
                )
                pred_chain = normalize_segment_chain(
                    row[f"{mate}_segments"], args.min_seg_len
                )
                truth_chain = normalize_segment_chain(
                    truth_row[f"{mate}_segments"], args.min_seg_len
                )
                pred_set = normalize_segment_set(
                    row[f"{mate}_segments"], args.min_seg_len
                )
                truth_set = normalize_segment_set(
                    truth_row[f"{mate}_segments"], args.min_seg_len
                )
                pred_junction_chain = normalize_junction_chain(
                    row[f"{mate}_segments"], args.min_seg_len
                )
                truth_junction_chain = normalize_junction_chain(
                    truth_row[f"{mate}_segments"], args.min_seg_len
                )
                error_class = classify_segment_error(pred_segments, truth_segments)
                segment_errors[(pred_type, mate, error_class)] += 1
                if pred_set == truth_set:
                    exact[(pred_type, set_key)] += 1
                elif len(examples[(pred_type, set_key)]) < args.examples:
                    examples[(pred_type, set_key)].append(
                        (read_id, ",".join(pred_set), ",".join(truth_set))
                    )
                if pred_chain == truth_chain:
                    exact[(pred_type, chain_key)] += 1
                elif len(examples[(pred_type, chain_key)]) < args.examples:
                    examples[(pred_type, chain_key)].append(
                        (read_id, "|".join(pred_chain), "|".join(truth_chain))
                    )
                if pred_junction_chain == truth_junction_chain:
                    exact[(pred_type, junction_key)] += 1
                elif len(examples[(pred_type, junction_key)]) < args.examples:
                    examples[(pred_type, junction_key)].append(
                        (
                            read_id,
                            "|".join(pred_junction_chain),
                            "|".join(truth_junction_chain),
                        )
                    )
                if pred_type == "bsj":
                    bsj_key = pred_mate_bsj_key(mate)
                    truth_bsj_key = truth_mate_bsj_key(mate)
                    if row[bsj_key] == truth_row[truth_bsj_key]:
                        exact[(pred_type, bsj_key)] += 1
                    elif len(examples[(pred_type, bsj_key)]) < args.examples:
                        examples[(pred_type, bsj_key)].append(
                            (read_id, row[bsj_key], truth_row[truth_bsj_key])
                        )

    recall_counts: Counter[str] = Counter()
    for read_id, truth_row in truth.items():
        if truth_row["type"] in ("bsj", "backward") and read_id in seen_pred:
            recall_counts[truth_row["type"]] += 1

    backward_total = pred_counts["backward"]
    backward_truth_bsj = confusion[("backward", "bsj")]
    backward_truth_backward = confusion[("backward", "backward")]
    backward_truth_forward = confusion[("backward", "forward")]
    outward_total = pred_counts["outward"]
    outward_truth_bsj = confusion[("outward", "bsj")]
    outward_truth_backward = confusion[("outward", "backward")]
    outward_truth_forward = confusion[("outward", "forward")]

    if args.format == "tsv":
        emit_summary("GLOBAL", "truth_bsj", truth_counts["bsj"])
        emit_summary("GLOBAL", "truth_backward", truth_counts["backward"])
        emit_summary("GLOBAL", "truth_forward", truth_counts["forward"])
        emit_summary("GLOBAL", "pred_bsj", pred_counts["bsj"])
        emit_summary("GLOBAL", "pred_backward", pred_counts["backward"])
        emit_summary("GLOBAL", "pred_outward", pred_counts["outward"])
        for pred_type in PRED_TYPES:
            for truth_type in ("bsj", "backward", "forward"):
                emit_summary(
                    "CONFUSION",
                    f"{pred_type}_vs_{truth_type}",
                    confusion[(pred_type, truth_type)],
                )
        emit_summary("GLOBAL", "bsj_recall", recall_counts["bsj"], truth_counts["bsj"])
        emit_summary("GLOBAL", "backward_recall", recall_counts["backward"], truth_counts["backward"])
        emit_summary("GLOBAL", "bsj_precision", confusion[("bsj", "bsj")], pred_counts["bsj"])
        emit_summary(
            "GLOBAL",
            "backward_precision",
            confusion[("backward", "backward")],
            pred_counts["backward"],
        )
        for key in PAIR_KEYS_BY_TYPE["bsj"]:
            emit_summary("BSJ_ONLY", f"bsj_{key}_exact", exact[("bsj", key)], pred_counts["bsj"])
        for mate in MATE_NAMES:
            emit_summary(
                "BSJ_ONLY",
                f"bsj_{mate}_segment_set_exact",
                exact[("bsj", f"{mate}_segment_set")],
                pred_counts["bsj"],
            )
            for error_class in SEGMENT_ERROR_CLASSES:
                emit_summary(
                    "BSJ_SEGMENT_ERRORS",
                    f"{mate}_{error_class}",
                    segment_errors[("bsj", mate, error_class)],
                    pred_counts["bsj"],
                )
            emit_summary(
                "BSJ_ONLY",
                f"bsj_{mate}_segment_chain_exact",
                exact[("bsj", f"{mate}_segment_chain")],
                pred_counts["bsj"],
            )
            emit_summary(
                "BSJ_ONLY",
                f"bsj_{mate}_junction_chain_exact",
                exact[("bsj", f"{mate}_junction_chain")],
                pred_counts["bsj"],
            )
            emit_summary(
                "BSJ_ONLY",
                f"bsj_is_{mate}_bsj_exact",
                exact[("bsj", pred_mate_bsj_key(mate))],
                pred_counts["bsj"],
            )
        for key in PAIR_KEYS_BY_TYPE["backward"]:
            emit_summary(
                "BACKWARD_ONLY",
                f"backward_{key}_exact",
                exact[("backward", key)],
                pred_counts["backward"],
            )
        for mate in MATE_NAMES:
            emit_summary(
                "BACKWARD_ONLY",
                f"backward_{mate}_segment_set_exact",
                exact[("backward", f"{mate}_segment_set")],
                pred_counts["backward"],
            )
            for error_class in SEGMENT_ERROR_CLASSES:
                emit_summary(
                    "BACKWARD_SEGMENT_ERRORS",
                    f"{mate}_{error_class}",
                    segment_errors[("backward", mate, error_class)],
                    pred_counts["backward"],
                )
            emit_summary(
                "BACKWARD_ONLY",
                f"backward_{mate}_segment_chain_exact",
                exact[("backward", f"{mate}_segment_chain")],
                pred_counts["backward"],
            )
            emit_summary(
                "BACKWARD_ONLY",
                f"backward_{mate}_junction_chain_exact",
                exact[("backward", f"{mate}_junction_chain")],
                pred_counts["backward"],
            )
        for key in PAIR_KEYS_BY_TYPE["outward"]:
            emit_summary(
                "OUTWARD_ONLY",
                f"outward_{key}_exact",
                exact[("outward", key)],
                pred_counts["outward"],
            )
        for mate in MATE_NAMES:
            emit_summary(
                "OUTWARD_ONLY",
                f"outward_{mate}_segment_set_exact",
                exact[("outward", f"{mate}_segment_set")],
                pred_counts["outward"],
            )
            for error_class in SEGMENT_ERROR_CLASSES:
                emit_summary(
                    "OUTWARD_SEGMENT_ERRORS",
                    f"{mate}_{error_class}",
                    segment_errors[("outward", mate, error_class)],
                    pred_counts["outward"],
                )
            emit_summary(
                "OUTWARD_ONLY",
                f"outward_{mate}_segment_chain_exact",
                exact[("outward", f"{mate}_segment_chain")],
                pred_counts["outward"],
            )
            emit_summary(
                "OUTWARD_ONLY",
                f"outward_{mate}_junction_chain_exact",
                exact[("outward", f"{mate}_junction_chain")],
                pred_counts["outward"],
            )
        emit_summary("BACKWARD_TRUTH_COMPOSITION", "backward_truth_is_bsj", backward_truth_bsj, backward_total)
        emit_summary(
            "BACKWARD_TRUTH_COMPOSITION",
            "backward_truth_is_backward",
            backward_truth_backward,
            backward_total,
        )
        emit_summary(
            "BACKWARD_TRUTH_COMPOSITION",
            "backward_truth_is_forward",
            backward_truth_forward,
            backward_total,
        )
        emit_summary(
            "BACKWARD_TRUTH_COMPOSITION",
            "backward_truth_is_circular",
            backward_truth_bsj + backward_truth_backward,
            backward_total,
        )
        emit_summary("OUTWARD_TRUTH_COMPOSITION", "outward_truth_is_bsj", outward_truth_bsj, outward_total)
        emit_summary(
            "OUTWARD_TRUTH_COMPOSITION",
            "outward_truth_is_backward",
            outward_truth_backward,
            outward_total,
        )
        emit_summary(
            "OUTWARD_TRUTH_COMPOSITION",
            "outward_truth_is_forward",
            outward_truth_forward,
            outward_total,
        )
        emit_summary(
            "OUTWARD_TRUTH_COMPOSITION",
            "outward_truth_is_circular",
            outward_truth_bsj + outward_truth_backward,
            outward_total,
        )
        for pred_type in PRED_TYPES:
            for key in example_keys_for_type(pred_type):
                example_rows = examples.get((pred_type, key))
                if not example_rows:
                    continue
                for read_id, pred_value, truth_value in example_rows:
                    emit_example("EXAMPLES", f"{pred_type}_{key}", read_id, pred_value, truth_value)
        return

    print_table(
        "GLOBAL",
        ["source", "bsj", "backward", "outward", "forward"],
        [
            ["truth", str(truth_counts["bsj"]), str(truth_counts["backward"]), "NA", str(truth_counts["forward"])],
            ["pred", str(pred_counts["bsj"]), str(pred_counts["backward"]), str(pred_counts["outward"]), "NA"],
            [
                "recall/precision",
                f"{format_rate(recall_counts['bsj'], truth_counts['bsj'])} / {format_rate(confusion[('bsj', 'bsj')], pred_counts['bsj'])}",
                f"{format_rate(recall_counts['backward'], truth_counts['backward'])} / {format_rate(confusion[('backward', 'backward')], pred_counts['backward'])}",
                f"NA / {format_rate(confusion[('outward', 'forward')], pred_counts['outward'])}",
                "NA",
            ],
        ],
    )

    print_table(
        "CONFUSION",
        ["pred\\truth", "bsj", "backward", "forward"],
        [
            ["bsj", str(confusion[("bsj", "bsj")]), str(confusion[("bsj", "backward")]), str(confusion[("bsj", "forward")])],
            [
                "backward",
                str(confusion[("backward", "bsj")]),
                str(confusion[("backward", "backward")]),
                str(confusion[("backward", "forward")]),
            ],
            [
                "outward",
                str(confusion[("outward", "bsj")]),
                str(confusion[("outward", "backward")]),
                str(confusion[("outward", "forward")]),
            ],
        ],
    )

    print_table(
        "BSJ",
        ["group", "metric", "matched", "total", "rate"],
        [
            ["pair", "circ_id", str(exact[("bsj", "circ_id")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "circ_id")], pred_counts["bsj"])],
            ["pair", "is_circular", str(exact[("bsj", "is_circular")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "is_circular")], pred_counts["bsj"])],
            ["r1", "is_bsj", str(exact[("bsj", "is_r1_bsj")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "is_r1_bsj")], pred_counts["bsj"])],
            ["r2", "is_bsj", str(exact[("bsj", "is_r2_bsj")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "is_r2_bsj")], pred_counts["bsj"])],
            ["r1", "segment_set", str(exact[("bsj", "r1_segment_set")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "r1_segment_set")], pred_counts["bsj"])],
            ["r2", "segment_set", str(exact[("bsj", "r2_segment_set")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "r2_segment_set")], pred_counts["bsj"])],
            ["r1", "segment_chain", str(exact[("bsj", "r1_segment_chain")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "r1_segment_chain")], pred_counts["bsj"])],
            ["r2", "segment_chain", str(exact[("bsj", "r2_segment_chain")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "r2_segment_chain")], pred_counts["bsj"])],
            ["r1", "junction_chain", str(exact[("bsj", "r1_junction_chain")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "r1_junction_chain")], pred_counts["bsj"])],
            ["r2", "junction_chain", str(exact[("bsj", "r2_junction_chain")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "r2_junction_chain")], pred_counts["bsj"])],
        ],
    )

    print_table(
        f"BSJ_SEGMENT_ERRORS min_seg_len={args.min_seg_len}",
        ["mate", "error_class", "count", "total", "rate"],
        [
            [
                mate,
                error_class,
                str(segment_errors[("bsj", mate, error_class)]),
                str(pred_counts["bsj"]),
                format_rate(segment_errors[("bsj", mate, error_class)], pred_counts["bsj"]),
            ]
            for mate in MATE_NAMES
            for error_class in SEGMENT_ERROR_CLASSES
        ],
    )

    print_table(
        "BACKWARD",
        ["group", "metric", "matched", "total", "rate"],
        [
            ["pair", "is_circular", str(exact[("backward", "is_circular")]), str(pred_counts["backward"]), format_rate(exact[("backward", "is_circular")], pred_counts["backward"])],
            ["r1", "segment_set", str(exact[("backward", "r1_segment_set")]), str(pred_counts["backward"]), format_rate(exact[("backward", "r1_segment_set")], pred_counts["backward"])],
            ["r2", "segment_set", str(exact[("backward", "r2_segment_set")]), str(pred_counts["backward"]), format_rate(exact[("backward", "r2_segment_set")], pred_counts["backward"])],
            ["r1", "segment_chain", str(exact[("backward", "r1_segment_chain")]), str(pred_counts["backward"]), format_rate(exact[("backward", "r1_segment_chain")], pred_counts["backward"])],
            ["r2", "segment_chain", str(exact[("backward", "r2_segment_chain")]), str(pred_counts["backward"]), format_rate(exact[("backward", "r2_segment_chain")], pred_counts["backward"])],
            ["r1", "junction_chain", str(exact[("backward", "r1_junction_chain")]), str(pred_counts["backward"]), format_rate(exact[("backward", "r1_junction_chain")], pred_counts["backward"])],
            ["r2", "junction_chain", str(exact[("backward", "r2_junction_chain")]), str(pred_counts["backward"]), format_rate(exact[("backward", "r2_junction_chain")], pred_counts["backward"])],
        ],
    )

    print_table(
        f"BACKWARD_SEGMENT_ERRORS min_seg_len={args.min_seg_len}",
        ["mate", "error_class", "count", "total", "rate"],
        [
            [
                mate,
                error_class,
                str(segment_errors[("backward", mate, error_class)]),
                str(pred_counts["backward"]),
                format_rate(
                    segment_errors[("backward", mate, error_class)],
                    pred_counts["backward"],
                ),
            ]
            for mate in MATE_NAMES
            for error_class in SEGMENT_ERROR_CLASSES
        ],
    )

    print_table(
        "OUTWARD",
        ["group", "metric", "matched", "total", "rate"],
        [
            ["pair", "is_circular", str(exact[("outward", "is_circular")]), str(pred_counts["outward"]), format_rate(exact[("outward", "is_circular")], pred_counts["outward"])],
            ["r1", "segment_set", str(exact[("outward", "r1_segment_set")]), str(pred_counts["outward"]), format_rate(exact[("outward", "r1_segment_set")], pred_counts["outward"])],
            ["r2", "segment_set", str(exact[("outward", "r2_segment_set")]), str(pred_counts["outward"]), format_rate(exact[("outward", "r2_segment_set")], pred_counts["outward"])],
            ["r1", "segment_chain", str(exact[("outward", "r1_segment_chain")]), str(pred_counts["outward"]), format_rate(exact[("outward", "r1_segment_chain")], pred_counts["outward"])],
            ["r2", "segment_chain", str(exact[("outward", "r2_segment_chain")]), str(pred_counts["outward"]), format_rate(exact[("outward", "r2_segment_chain")], pred_counts["outward"])],
            ["r1", "junction_chain", str(exact[("outward", "r1_junction_chain")]), str(pred_counts["outward"]), format_rate(exact[("outward", "r1_junction_chain")], pred_counts["outward"])],
            ["r2", "junction_chain", str(exact[("outward", "r2_junction_chain")]), str(pred_counts["outward"]), format_rate(exact[("outward", "r2_junction_chain")], pred_counts["outward"])],
        ],
    )

    print_table(
        f"OUTWARD_SEGMENT_ERRORS min_seg_len={args.min_seg_len}",
        ["mate", "error_class", "count", "total", "rate"],
        [
            [
                mate,
                error_class,
                str(segment_errors[("outward", mate, error_class)]),
                str(pred_counts["outward"]),
                format_rate(
                    segment_errors[("outward", mate, error_class)],
                    pred_counts["outward"],
                ),
            ]
            for mate in MATE_NAMES
            for error_class in SEGMENT_ERROR_CLASSES
        ],
    )

    print_table(
        "BACKWARD_TRUTH_COMPOSITION",
        ["metric", "matched", "total", "rate"],
        [
            ["truth_is_bsj", str(backward_truth_bsj), str(backward_total), format_rate(backward_truth_bsj, backward_total)],
            ["truth_is_backward", str(backward_truth_backward), str(backward_total), format_rate(backward_truth_backward, backward_total)],
            ["truth_is_forward", str(backward_truth_forward), str(backward_total), format_rate(backward_truth_forward, backward_total)],
            [
                "truth_is_circular",
                str(backward_truth_bsj + backward_truth_backward),
                str(backward_total),
                format_rate(backward_truth_bsj + backward_truth_backward, backward_total),
            ],
        ],
    )

    print_table(
        "OUTWARD_TRUTH_COMPOSITION",
        ["metric", "matched", "total", "rate"],
        [
            ["truth_is_bsj", str(outward_truth_bsj), str(outward_total), format_rate(outward_truth_bsj, outward_total)],
            ["truth_is_backward", str(outward_truth_backward), str(outward_total), format_rate(outward_truth_backward, outward_total)],
            ["truth_is_forward", str(outward_truth_forward), str(outward_total), format_rate(outward_truth_forward, outward_total)],
            [
                "truth_is_circular",
                str(outward_truth_bsj + outward_truth_backward),
                str(outward_total),
                format_rate(outward_truth_bsj + outward_truth_backward, outward_total),
            ],
        ],
    )

    for pred_type in PRED_TYPES:
        for key in example_keys_for_type(pred_type):
            example_rows = examples.get((pred_type, key))
            if not example_rows:
                continue
            print_table(
                f"EXAMPLES {pred_type}_{key}",
                ["read_id", "pred", "truth"],
                [[read_id, pred_value, truth_value] for read_id, pred_value, truth_value in example_rows],
            )


if __name__ == "__main__":
    main()
