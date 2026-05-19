---
name: ciri-segments-isoform-eval
description: Use for CIRI-toolkit project-local evaluation of read-level segments, BSJ/backward/outward semantics, outward accuracy, segment graph behavior, major isoform reconstruction, mature/estimate classification, FASTA confidence, and chr1 simulator or real-data isoform accuracy.
---

# CIRI Segments And Isoform Evaluation

Use this skill inside `CIRI-toolkit` when the task asks whether `bsj`,
`backward`, `outward`, `mature`, `estimate`, or FASTA isoforms are correct.

## Boundaries

- Treat `.out/.bsj` CIRI3 parity as fixed unless the user explicitly reopens
  BSJ detection.
- Treat `<prefix>.segments` as the isoform-stage input boundary. Isoform
  evaluation should parse the written `.segments` file, matching the production
  `--continue` behavior.
- Do not use CIRI-AS/CIRI-full/RO remap as a target contract. They are only
  historical references for risk and algorithm ideas.
- Keep truth labels and prediction labels separate:
  - simulator truth is generated-model truth, normally `bsj/outward/forward`
  - prediction may contain `bsj/backward/outward`
  - `backward` prediction is BSJ-like evidence without a precise final BSJ site

## Standard Workflow

1. Identify the artifact set:
   - simulated chr1: BAM, annotation, reference, simulator truth, `.segments`,
     `.isoforms.gtf`, `.isoforms.fa`
   - real data: `.out`, `.segments`, `.isoforms.gtf`, `.isoforms.fa`,
     optionally `.segments.bam/.bai` for IGV review
2. Use existing scripts before writing ad hoc parsers:
   - `scripts/ciri_segments_eval.py` for read-level truth/prediction evaluation
   - `scripts/ciri_result_diff.py` when checking `.out/.bsj` parity
   - `scripts/ciri_read_subset.py` and
     `scripts/ciri_extract_interval_read_names.py` for focused examples
3. Split metrics instead of blending them:
   - `bsj`, `backward`, `outward` read classes
   - `mature` vs `estimate`
   - FASTA-retained vs GTF-only isoforms
   - BSJ-only, BSJ+backward, and BSJ+backward+outward evidence tiers
4. Always inspect representative examples after a metric shift:
   - read ID
   - circRNA coordinate
   - R1/R2 segment strings and align strands
   - junction-chain support
   - exclusive evidence showing an internal junction was not used
   - annotation support and segment coverage percentage

## Segment Semantics

- `bsj` means a read/mate explicitly identifies a BSJ site.
- `backward` means BSJ-like circular evidence exists, but the exact BSJ site is
  not fully resolved.
- `outward` means an outward-facing pair supports a circular origin but cannot
  define the BSJ range by itself.
- `R1/R2 segments` should be interpreted as reconstructed aligned pieces; use
  separate mate align-strand fields when auditing mapper orientation.
- Do not count a linear-compatible read as outward unless linear mapping has
  been conservatively rejected by the implemented rules.

## Isoform Evaluation Rules

- A mature isoform requires adjacent junction-chain support. A long block is
  not mature just because a read covers one side of it.
- Phasing support and exclusive evidence have equal weight when deciding
  whether a block is mature or must be downgraded to estimate.
- Outward and backward reads may fill missing path evidence, but BSJ/backward
  phasing has higher priority than unphased outward support.
- For estimate isoforms, check the estimate reason before trusting length:
  annotated internal exons, aligned junction evidence, and segment coverage are
  more informative than a hard length cutoff.
- FASTA should contain only high-confidence isoforms. When reviewing exclusions,
  first check whether segment coverage is high enough or whether annotation and
  read evidence make the estimate biologically plausible.

## Useful Commands

Prefer release mode for real fixture or performance-sensitive checks:

```bash
cargo run --release -- -i tests/chr1/test.bam -o tmp/eval.ciri \
  -r tests/chr1/chr1.fa -a tests/chr1/test.annotation.gtf -t 16 -s 0
```

For fast isoform-only iteration from completed segments:

```bash
cargo run --release -- -i tests/chr1/test.bam -o tmp/eval.ciri \
  -r tests/chr1/chr1.fa -a tests/chr1/test.annotation.gtf -t 16 -s 0 --continue
```

## Reporting

Report counts and rates with the denominator stated. When judging a failure,
include one or two concrete read/circRNA examples with enough fields for IGV or
text inspection. Do not present one aggregate accuracy as the full answer when
error modes differ by evidence tier.
