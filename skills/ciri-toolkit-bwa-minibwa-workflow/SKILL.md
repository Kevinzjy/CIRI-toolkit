---
name: ciri-toolkit-bwa-minibwa-workflow
description: Use when running CIRI-toolkit on paired-end short-read data aligned with BWA-MEM or minibwa, including command templates, input requirements, expected outputs, and basic result review.
---

# CIRI-toolkit BWA/minibwa User Workflow

Use this skill when helping a user run CIRI-toolkit from paired-end FASTQ data
aligned by `bwa mem` or `minibwa map`.

## Before Running

- Prepare paired-end FASTQ files, a reference FASTA, and optionally a GTF
  annotation.
- Install `samtools` because the workflow writes BAM and CIRI-toolkit writes an
  indexed segments BAM for IGV review.
- Keep the BAM unsorted by coordinate. Do not run `samtools sort` before
  CIRI-toolkit.
- Use the same reference FASTA for mapping and for `ciri -r`.

## Choose A Mapper

Use BWA-MEM when you want the conservative, widely used alignment route:

```bash
bwa mem -t <threads> -T 19 <ref.fa> <R1.fq.gz> <R2.fq.gz> \
  | samtools view -bS -@ <threads> -o <sample>.bam -
```

Use minibwa when you want a faster BWA-compatible route:

```bash
minibwa map -t <threads> -x sr --adap=no -m19 -s19 <ref.fa> <R1.fq.gz> <R2.fq.gz> \
  | samtools view -bS -@ <threads> -o <sample>.bam -
```

Use one mapper consistently across samples unless the goal is to compare mappers.

## Run CIRI-toolkit

With annotation:

```bash
ciri \
  -i <sample>.bam \
  -o <prefix> \
  -r <ref.fa> \
  -a <annotation.gtf> \
  -t <threads>
```

Without annotation:

```bash
ciri \
  -i <sample>.bam \
  -o <prefix> \
  -r <ref.fa> \
  -t <threads>
```

CIRI-toolkit defaults are intended to keep candidates for downstream segments
and isoform reconstruction. If the user asks for stricter CIRI3-style filtering,
run:

```bash
ciri -i <sample>.bam -o <prefix> -r <ref.fa> -a <annotation.gtf> \
  -t <threads> -s 2 --min-span 140
```

## Expected Outputs

A normal run writes:

- `<prefix>.out`: circRNA result table
- `<prefix>.bsj`: BSJ evidence table
- `<prefix>.bedpe`: BSJ anchor track for IGV
- `<prefix>.segments`: read-level segment evidence
- `<prefix>.segments.bam` and `<prefix>.segments.bam.bai`: IGV review BAM
- `<prefix>.isoforms.gtf`: reconstructed major isoform structures
- `<prefix>.isoforms.fa`: high-confidence isoform sequences
- `<prefix>.log`: run log

Use `--continue` only when `<prefix>.out` and `<prefix>.segments` already exist
and the user wants to rebuild isoform outputs without rerunning the full BAM
scan.

## Basic Checks

- Confirm `<prefix>.log` ends without an error.
- Confirm `<prefix>.out` is non-empty when circRNAs are expected.
- Load `<prefix>.bedpe` and `<prefix>.segments.bam` in IGV for manual review of
  BSJ anchors and read-level segment structure.
- If no circRNA is detected, CIRI-toolkit may still finish successfully with
  header-only or empty outputs.

## Multi-Sample Use

For a shared circRNA catalog, run first-pass detection per sample, merge the
catalog, then run second-pass analysis:

```bash
ciri -i <sample1>.bam -o <sample1>.first -r <ref.fa> -a <annotation.gtf> \
  -t <threads> --1st-pass
ciri -i <sample2>.bam -o <sample2>.first -r <ref.fa> -a <annotation.gtf> \
  -t <threads> --1st-pass

ciri-merge -o <catalog_prefix> <sample1>.first.out <sample2>.first.out

ciri -i <sample1>.bam -o <sample1>.second -r <ref.fa> -a <annotation.gtf> \
  -t <threads> --2nd-pass --circ <catalog_prefix>.circ.bed
ciri -i <sample2>.bam -o <sample2>.second -r <ref.fa> -a <annotation.gtf> \
  -t <threads> --2nd-pass --circ <catalog_prefix>.circ.bed

ciri-assemble -o <assemble_prefix> \
  <sample1>.second.out <sample1>.second.segments \
  <sample2>.second.out <sample2>.second.segments
```

