# CIRI-toolkit

`CIRI-toolkit` is a high-performance circRNA detection and isoform reconstruction toolkit for large-scale transcriptome data. It provides CIRI3-compatible BSJ detection, read-level circRNA segment reconstruction, high-confidence isoform sequence output, and IGV-ready review tracks.

## Performance Benchmarking

Benchmark dataset: ~80 GB BAM file (~300 GB SAM).

| Tool / workflow | Output scope | Runtime |
|---|---|---:|
| CIRI3 | BSJ detection | ~80 min |
| CIRI-toolkit | BSJ detection | ~6 min |
| CIRI-toolkit | BSJ detection + isoform reconstruction | <20 min |

## Prerequisites

- Rust >= 1.85.0, if building from source
- GCC >= 5
- `pigz` or `gzip`
- `samtools` in `PATH`, or set `SAMTOOLS=/path/to/samtools`

## Installation

Download a prebuilt binary from the [release page](https://bioinfo.ioz.ac.cn/git/zhangjy/CIRI-toolkit/releases), or build from source:

```bash
git clone https://bioinfo.ioz.ac.cn/git/zhangjy/CIRI-toolkit.git
cd CIRI-toolkit
cargo build --release

# Run the compiled binary
./target/release/ciri --help
```

## Quick Start

```bash
# Step 1. Align reads with BWA-MEM and write a BAM file.
bwa mem -t <threads> -T 19 <bwa_index> <R1.fastq.gz> <R2.fastq.gz> \
  | samtools view -bS -@ <threads> -o <bam_file> -

# Step 2. Run CIRI-toolkit.
ciri \
  -i <bam_file> \
  -o <prefix> \
  -r <reference_fasta> \
  -a <annotation_gtf> \
  -t <threads>
```

## Usage

```text
CIRI-toolkit: fast circRNA detection and isoform reconstruction

Usage: ciri [OPTIONS] --in <IN_SAM> --out <OUT_PREFIX> --ref <REF_FASTA>

Required arguments:
  -i, --in <SAM/BAM>              Input SAM/BAM file; BAM must not be coordinate-sorted
  -o, --out <PREFIX>              Output prefix
  -r, --ref <FASTA>               Reference genome FASTA

Optional arguments:
  -a, --anno <GTF>                GTF annotation file
  -m, --mapq <MIN_MAPQ>           Minimum MAPQ for candidate BSJ reads (default: 10)
  -s, --stringency <LEVEL>        Summary filter level [0/1/2] (default: 0)
  --min-span <SIZE>               Minimum circRNA span (default: 50)
  --max-span <SIZE>               Maximum circRNA span (default: 200000)
  --linear-range-size-min <SIZE>  Linear competition search range size (default: 50000)
  -t, --threads <THREADS>         Number of worker threads (default: auto)
  -M, --mem-per-thread <MEM>      Maximum memory per thread (default: 512M; e.g., 2G)

Review and debugging arguments:
  --trace <READS>                 Comma-separated read IDs to trace in Scan1/Scan2
  --debug                         Keep internal pipeline temporary files
  --perf                          Write profiling report to <prefix>.perf.log
  --continue                      Rebuild isoforms only from existing <prefix>.out + <prefix>.segments
  -h, --help                      Print help
  -v, --version                   Print version
```

**Note:** CIRI-toolkit keeps the CIRI3 BSJ detection, Scan2 rescue, and Summary stringency formulas, but its default command-line policy is tuned for downstream segment and isoform reconstruction. CIRI-toolkit defaults to `-s 0 --min-span 50`, while CIRI3 defaults to `-S 2 -Min 140`. To use the same filtering policy as CIRI3, run CIRI-toolkit with `-s 2 --min-span 140`.

## Outputs

The main output files are:

- `<prefix>.out`: CIRI3-compatible circRNA result table
- `<prefix>.bsj`: mate-level BSJ evidence display used for review/debugging
- `<prefix>.segments`: read-level BSJ/backward/outward segments
- `<prefix>.isoforms.gtf`: major isoform structure audit table for all reported circRNAs; attributes use `gene_id`, `transcript_id`, `type`, `evidence`, `weakness`, `score`, `bsj_reads`, `weight`, `exon_count`, and `isoform_len`
- `<prefix>.isoforms.fa`: high-confidence circRNA isoform sequences; FASTA IDs use `<circRNA_id>.iso1` and headers keep only compact sequence fields such as `type`, `len`, and `cirexon`
- `<prefix>.bedpe`: IGV-compatible BSJ anchor track
- `<prefix>.segments.bam` and `<prefix>.segments.bam.bai`: IGV-compatible segment alignments
- `<prefix>.log`: run log

## Validate isoform FASTA

Use `scripts/ciri_isoform_fasta_validate.py` to verify that every FASTA record
matches the exon chain in `<prefix>.isoforms.gtf` and the supplied reference,
including negative-strand reverse-complement handling:

```bash
python scripts/ciri_isoform_fasta_validate.py \
  --fasta <prefix>.isoforms.fa \
  --gtf <prefix>.isoforms.gtf \
  --reference <reference_fasta>
```

An optional sampled BLAT strand check can be added when the UCSC BLAT binary is
available:

```bash
python scripts/ciri_isoform_fasta_validate.py \
  --fasta <prefix>.isoforms.fa \
  --gtf <prefix>.isoforms.gtf \
  --reference <reference_fasta> \
  --blat /data/public/software/UCSC_utility/blat/blat \
  --blat-sample-size 100
```

## `.out` Format

`<prefix>.out` uses the CIRI3-compatible 13-column format:

```text
1.  circRNA_ID
2.  chr
3.  circRNA_start
4.  circRNA_end
5.  #junction_reads
6.  SM_MS_SMS
7.  #non_junction_reads
8.  junction_reads_ratio
9.  circRNA_type
10. gene_id
11. strand
12. junction_reads_ID
13. Score
```

## Generate a simulation dataset with `ciri-simulator`

CIRI-toolkit also includes `ciri-simulator`, which generates paired-end circRNA reads and structured truth tables from a reference genome and annotation.

```bash
ciri-simulator \
  -r <reference_fasta> \
  -a <source_gtf> \
  -o <prefix> \
  --circ-count 100 \
  --circ-coverage 10 \
  --linear-coverage 0.1 \
  --exon-exclusive-rate 0.25 \
  --seed 5
```

Simulator outputs:

- `<prefix>_1.fq.gz` and `<prefix>_2.fq.gz`: simulated paired-end reads
- `<prefix>.annotation.gtf`: masked linear annotation for de novo evaluation
- `<prefix>.isoforms.tsv`: isoform-level truth table
- `<prefix>.reads.tsv`: read-pair-level truth table
