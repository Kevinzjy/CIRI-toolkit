#!/usr/bin/env python3
"""Evaluate `<prefix>.segments` against simulator `*.reads.tsv` truth.

Use this after running the default post-Summary segments phase on simulator
fixtures. The script focuses on the current read-level contract and prints two
extra views that are useful during the current development stage:

- predicted `type` (`bsj` / `backward`)
- `circ_id`
- `is_circular`
- mate-level genomic segment sets from `r1_segments` / `r2_segments`
- mate-level `r1_is_bsj` / `r2_is_bsj`

Additional summaries:

- `BSJ_ONLY`: read-pair and mate-level metrics restricted to predicted `type=bsj`
- `BACKWARD_AS_MISSED_BSJ`: how many predicted `type=backward` rows are actually
  BSJ reads in simulator truth, which is useful while backward output is still
  acting as a pool of missed-BSJ candidates

The default output format is a compact human-readable report. Use `--format tsv`
to emit tidy rows that are easier to redirect into downstream shell tooling.
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
        "is_bsj",
    ],
    "backward": [
        "is_circular",
        "is_bsj",
    ],
}

MATE_NAMES = ("r1", "r2")


def infer_truth_type(row: dict[str, str]) -> str:
    """Map simulator truth labels to the current segments `type` contract."""
    if row["is_bsj"] == "1":
        return "bsj"
    if row["is_circular"] == "1":
        return "backward"
    return "forward"


def normalize_segment_set(text: str) -> tuple[str, ...]:
    """Normalize one mate's segment tokens into an order-insensitive tuple.

    `<bsj>` is intentionally excluded here because BSJ presence is tracked by the
    dedicated `r1_is_bsj` / `r2_is_bsj` fields. The current evaluation goal is
    to ask whether the same genomic segments were recovered, regardless of read
    order or mapper-specific token order.
    """
    if text == "NA":
        return tuple()
    tokens = [token for token in text.split("|") if token and token != "<bsj>"]
    return tuple(sorted(tokens))


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
    args = parser.parse_args()

    if args.format == "tsv":
        emit_summary = emit_summary_row_tsv
        emit_example = emit_example_row_tsv
        print("kind\tsection\tmetric\tvalue_or_read_id\ttotal_or_pred\ttruth_or_rate")

    truth, truth_counts = load_truth(args.truth_reads_tsv)
    pred_counts: Counter[str] = Counter()
    confusion: Counter[tuple[str, str]] = Counter()
    exact: Counter[tuple[str, str]] = Counter()
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
                pred_set = normalize_segment_set(row[f"{mate}_segments"])
                truth_set = normalize_segment_set(truth_row[f"{mate}_segments"])
                if pred_set == truth_set:
                    exact[(pred_type, set_key)] += 1
                elif len(examples[(pred_type, set_key)]) < args.examples:
                    examples[(pred_type, set_key)].append(
                        (read_id, ",".join(pred_set), ",".join(truth_set))
                    )
                bsj_key = f"{mate}_is_bsj"
                if row[bsj_key] == truth_row[bsj_key]:
                    exact[(pred_type, bsj_key)] += 1
                elif len(examples[(pred_type, bsj_key)]) < args.examples:
                    examples[(pred_type, bsj_key)].append(
                        (read_id, row[bsj_key], truth_row[bsj_key])
                    )

    recall_counts: Counter[str] = Counter()
    for read_id, truth_row in truth.items():
        if truth_row["type"] in ("bsj", "backward") and read_id in seen_pred:
            recall_counts[truth_row["type"]] += 1

    backward_total = pred_counts["backward"]
    backward_truth_bsj = confusion[("backward", "bsj")]
    backward_truth_backward = confusion[("backward", "backward")]
    backward_truth_forward = confusion[("backward", "forward")]

    if args.format == "tsv":
        emit_summary("GLOBAL", "truth_bsj", truth_counts["bsj"])
        emit_summary("GLOBAL", "truth_backward", truth_counts["backward"])
        emit_summary("GLOBAL", "truth_forward", truth_counts["forward"])
        emit_summary("GLOBAL", "pred_bsj", pred_counts["bsj"])
        emit_summary("GLOBAL", "pred_backward", pred_counts["backward"])
        for pred_type in ("bsj", "backward"):
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
            emit_summary(
                "BSJ_ONLY",
                f"bsj_{mate}_is_bsj_exact",
                exact[("bsj", f"{mate}_is_bsj")],
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
            emit_summary(
                "BACKWARD_ONLY",
                f"backward_{mate}_is_bsj_exact",
                exact[("backward", f"{mate}_is_bsj")],
                pred_counts["backward"],
            )
        emit_summary("BACKWARD_AS_MISSED_BSJ", "backward_truth_is_bsj", backward_truth_bsj, backward_total)
        emit_summary(
            "BACKWARD_AS_MISSED_BSJ",
            "backward_truth_is_backward",
            backward_truth_backward,
            backward_total,
        )
        emit_summary(
            "BACKWARD_AS_MISSED_BSJ",
            "backward_truth_is_forward",
            backward_truth_forward,
            backward_total,
        )
        emit_summary(
            "BACKWARD_AS_MISSED_BSJ",
            "backward_truth_is_circular",
            backward_truth_bsj + backward_truth_backward,
            backward_total,
        )
        for pred_type in ("bsj", "backward"):
            for key in (
                [*PAIR_KEYS_BY_TYPE[pred_type]]
                + [f"{mate}_segment_set" for mate in MATE_NAMES]
                + [f"{mate}_is_bsj" for mate in MATE_NAMES]
            ):
                example_rows = examples.get((pred_type, key))
                if not example_rows:
                    continue
                for read_id, pred_value, truth_value in example_rows:
                    emit_example("EXAMPLES", f"{pred_type}_{key}", read_id, pred_value, truth_value)
        return

    print_table(
        "GLOBAL",
        ["source", "bsj", "backward", "forward"],
        [
            ["truth", str(truth_counts["bsj"]), str(truth_counts["backward"]), str(truth_counts["forward"])],
            ["pred", str(pred_counts["bsj"]), str(pred_counts["backward"]), "NA"],
            [
                "recall/precision",
                f"{format_rate(recall_counts['bsj'], truth_counts['bsj'])} / {format_rate(confusion[('bsj', 'bsj')], pred_counts['bsj'])}",
                f"{format_rate(recall_counts['backward'], truth_counts['backward'])} / {format_rate(confusion[('backward', 'backward')], pred_counts['backward'])}",
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
        ],
    )

    print_table(
        "BSJ",
        ["group", "metric", "matched", "total", "rate"],
        [
            ["pair", "circ_id", str(exact[("bsj", "circ_id")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "circ_id")], pred_counts["bsj"])],
            ["pair", "is_circular", str(exact[("bsj", "is_circular")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "is_circular")], pred_counts["bsj"])],
            ["pair", "is_bsj", str(exact[("bsj", "is_bsj")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "is_bsj")], pred_counts["bsj"])],
            ["r1", "segment_set", str(exact[("bsj", "r1_segment_set")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "r1_segment_set")], pred_counts["bsj"])],
            ["r1", "is_bsj", str(exact[("bsj", "r1_is_bsj")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "r1_is_bsj")], pred_counts["bsj"])],
            ["r2", "segment_set", str(exact[("bsj", "r2_segment_set")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "r2_segment_set")], pred_counts["bsj"])],
            ["r2", "is_bsj", str(exact[("bsj", "r2_is_bsj")]), str(pred_counts["bsj"]), format_rate(exact[("bsj", "r2_is_bsj")], pred_counts["bsj"])],
        ],
    )

    print_table(
        "BACKWARD",
        ["group", "metric", "matched", "total", "rate"],
        [
            ["pair", "is_circular", str(exact[("backward", "is_circular")]), str(pred_counts["backward"]), format_rate(exact[("backward", "is_circular")], pred_counts["backward"])],
            ["pair", "is_bsj", str(exact[("backward", "is_bsj")]), str(pred_counts["backward"]), format_rate(exact[("backward", "is_bsj")], pred_counts["backward"])],
            ["r1", "segment_set", str(exact[("backward", "r1_segment_set")]), str(pred_counts["backward"]), format_rate(exact[("backward", "r1_segment_set")], pred_counts["backward"])],
            ["r1", "is_bsj", str(exact[("backward", "r1_is_bsj")]), str(pred_counts["backward"]), format_rate(exact[("backward", "r1_is_bsj")], pred_counts["backward"])],
            ["r2", "segment_set", str(exact[("backward", "r2_segment_set")]), str(pred_counts["backward"]), format_rate(exact[("backward", "r2_segment_set")], pred_counts["backward"])],
            ["r2", "is_bsj", str(exact[("backward", "r2_is_bsj")]), str(pred_counts["backward"]), format_rate(exact[("backward", "r2_is_bsj")], pred_counts["backward"])],
        ],
    )

    print_table(
        "BACKWARD_AS_MISSED_BSJ",
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

    for pred_type in ("bsj", "backward"):
        for key in (
            [*PAIR_KEYS_BY_TYPE[pred_type]]
            + [f"{mate}_segment_set" for mate in MATE_NAMES]
            + [f"{mate}_is_bsj" for mate in MATE_NAMES]
        ):
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
