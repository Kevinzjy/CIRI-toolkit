---
name: ciri-release-and-doc-sync
description: Use for CIRI-toolkit project-local release preparation, README/AGENTS/docs synchronization, output-contract documentation, stale CIRI-rs wording checks, release-note drafting, and pre-commit validation after CIRI CLI, segments, isoform, BEDPE, or segments BAM behavior changes.
---

# CIRI Release And Documentation Sync

Use this skill inside `CIRI-toolkit` when preparing a release or synchronizing
documentation after output, CLI, logging, segments, isoform, or workflow changes.

## Scope

- Keep public user docs focused on what users need to run and interpret CIRI.
- Keep internal design, validation, and future development details in `AGENTS.md`
  or `docs/`.
- Keep CIRI3 parity boundaries explicit: `vendor/CIRI3` defines core BSJ
  behavior; segments, isoforms, BEDPE, and segments BAM are CIRI extensions.
- Do not revive CIRI-AS/CIRI-full/RO remap as current roadmap language. Mention
  them only as historical references when needed.

## Standard Workflow

1. Inspect changed behavior:
   - CLI defaults and flags
   - logging order and messages
   - output files and temp files
   - segments/isoform evidence contracts
   - release-visible performance or dependency requirements
2. Sync the project-facing docs:
   - `README.md` for user-facing install, run, outputs, options, and examples
   - `AGENTS.md` for agent workflow, constraints, and current development stage
   - `docs/01-development-status.md` for status and roadmap
   - `docs/07-full-length-reconstruction.md` for segments/isoform contracts
   - `docs/00-index.md` if new or renamed docs must be discoverable
3. Search for stale wording:
   - `CIRI-rs`, unless intentionally describing old files or historical data
   - abandoned `CIRI-AS`, `CIRI-full`, or `RO remap` roadmap language
   - obsolete outputs such as `.segments.bed`
   - old defaults such as Java-only `-S 2` or `--min-span 140` when describing
     the Rust CLI defaults
4. Validate formatting and build state:
   - run `cargo fmt` after code or Rust-doc edits
   - run `cargo build` or `cargo check` when the docs sync accompanies code
   - run targeted tests when behavior changed

## Output Contract Checklist

Current public outputs are:

- `<prefix>.out`
- `<prefix>.bsj`
- `<prefix>.bedpe`
- `<prefix>.segments`
- `<prefix>.segments.bam`
- `<prefix>.segments.bam.bai`
- `<prefix>.isoforms.gtf`
- `<prefix>.isoforms.fa`

Internal temp outputs include `.bsj1`, `.bsj2`, `.segments1`, `.segments2`,
`.segments.non_bsj`, and `.part_XXXX.tmp` shards. They should remain internal
unless `--debug` is being documented.

## Release Checklist

- Confirm `samtools` is documented as required because segments BAM indexing
  depends on it.
- Confirm `--continue` is documented as isoform-only rebuild from completed
  `<prefix>.segments`.
- Confirm GTF `source` and FASTA naming use `CIRI` consistently.
- Confirm README examples use current defaults: `-s 0` and `--min-span 50`.
- Confirm release notes separate core BSJ parity from extension features.
- Before commit, inspect `git diff --stat` and avoid mixing unrelated cleanup
  with release-contract changes.

## Reporting

Summaries should list what user-facing contract changed, which docs were
updated, and what validation was run. Keep implementation history short unless
the user asks for a full development record.
