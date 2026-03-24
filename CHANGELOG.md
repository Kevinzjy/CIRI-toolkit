# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, adapted to the current release flow of
`CIRI-toolkit`.

## [0.1.1] - 2026-03-24

### Added
- Added `tests/analyze_diff.py` FSJ parity reporting for shared circRNAs.
- Added `<prefix>.log` as a first-class pipeline output so stage summaries and
  runtime information are preserved outside the terminal session.
- Added CLI `--debug <read_id_list>` to trace selected reads through the
  pipeline and write the detailed trace to `<prefix>.debug.log`.
- Added CLI `--perf` to emit Scan1/Scan2 profiling summaries to
  `<prefix>.perf.log`.
- Added source labels to the final `<prefix>.bsj` output so each BSJ-supporting
  read records whether it was emitted by `scan1` or `scan2`.

### Changed
- Normalized output naming around one prefix:
  - final result: `<prefix>.out`
  - final BSJ reads: `<prefix>.bsj`
  - intermediate Scan1 output: `<prefix>.bsj1`
  - shard-local temp files: `<prefix>.bsj1.part_0001.tmp`,
    `<prefix>.bsj2.part_0001.tmp`, `<prefix>.fsj.part_0001.tmp`
- Removed the final `.fsj` output file. FSJ counts remain available in the
  final `.out`, while shard-local FSJ spill files are still used internally and
  cleaned up automatically.
- Expanded runtime logging with stage summaries:
  - Scan1 mapped reads and BSJ1 reads
  - Scan2 rescued BSJ2 reads
  - final circRNA count, final BSJ read count, and total runtime
- Increased comment coverage around the parity-sensitive Scan2/FSJ path and the
  output/logging conventions that were established through repeated validation.
- Updated final circRNA sorting so canonical `chr*` chromosomes are emitted
  before scaffold/contig names such as `GL*` and `KI*`.
- Refined Scan2 BAM parity with Java CIRI3:
  - aligned FSJ bucket-gating with Java `GetFSJClass.getFSJ(...)`
  - aligned BAM mate-run overwrite semantics for per-mate alignment groups
  - aligned Scan2 representative-sequence overwrite behavior with Java
- Updated project documentation to use the verified hg38 whole-genome parity
  baseline:
  - FASTA: `/data/public/database/gencode/hg38/_BWAindex/hg38.fa`
  - GTF: `/data/public/database/gencode/hg38/gencode.v44.annotation.gtf`

### Fixed
- Fixed Scan2 FSJ counting drift on whole-genome BAM inputs. With the verified
  hg38 FASTA/GTF pairing, shared-circ `#non_junction_reads` now match Java
  exactly (`fsj_diff=0` in `tests/analyze_diff.py`).
- Fixed Scan2 whole-genome parity regressions caused by BAM mate-switch handling
  within supplementary-heavy read groups.
- Fixed documentation drift so debugging, parity checking, and release-state
  notes reflect the current code and validation results.

### Notes
- `tests/chr1` remains fully aligned at circ/read/read-assignment/FSJ levels.
- Whole-genome hg38 validation now reaches full FSJ parity and near-complete
  circ/read parity; the remaining tiny set of family-level residual cases has
  been manually reviewed and judged to be more reasonable on the Rust side than
  on the Java side.
