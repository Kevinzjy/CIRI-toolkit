# CIRI-toolkit (Rust)

A high-performance Rust reimplementation of the CIRI3 circular RNA identification algorithm.

## Features
- **Pixel-level Parity**: Achieves >97% read-level alignment accuracy with original CIRI3 Java implementation.
- **Split-mapping Identification (Scan 1)**: Efficient detection of Back-Spliced Junction (BSJ) signals.
- **PEM/SMS Rescue (Scan 2)**: Advanced rescue mechanism for junction reads with ambiguous alignments.
- **Global Unique Read Assignment**: Resolves alignment ambiguity by assigning each read to the single most likely circRNA site.
- **Automated Annotation**: GTF-based circular RNA type (exon/intron/intergenic) and gene ID annotation.
- **高性能**: Significant speedup and reduced memory footprint compared to the original Java implementation.

## Usage

### Build
```bash
cargo build --release
```

### Execution
```bash
./target/release/ciri-toolkit \
    -i input.sam \
    -o output_prefix \
    -r reference.fa \
    -a annotation.gtf \
    -m 10 \
    -s 1
```

#### Parameters:
- `-i, --in`: Path to the input SAM file (mapped with BWA-MEM).
- `-o, --out`: Prefix for output files (`.BSJ1`, `.result`).
- `-r, --ref`: Path to the reference genome FASTA file.
- `-a, --anno`: (Optional) Path to the GTF annotation file.
- `-m, --mapq`: Minimum Mapping Quality (default: 10).
- `-s, --stringency`: Stringency level (0, 1, or 2; default: 1).

## Output Format
The tool generates a 13-column `.result` file compatible with CIRI3 downstream analysis:
1. `circRNA_ID`: Chromosome and 1-based coordinates.
2. `chr`: Chromosome name.
3. `circRNA_start`: 1-based start position.
4. `circRNA_end`: 1-based end position.
5. `#junction_reads`: Number of BSJ reads assigned to this site.
6. `SM_MS_SMS`: Distribution of junction read CIGAR types.
7. `#non_junction_reads`: Number of linear reads (FSJ) overlapping the junction.
8. `junction_reads_ratio`: Ratio of BSJ reads to total junction-crossing reads.
9. `circRNA_type`: Genomic feature type (e.g., exon, intergenic).
10. `gene_id`: Associated gene ID from GTF.
11. `strand`: Strand information.
12. `junction_reads_ID`: List of unique Read IDs providing evidence for the BSJ.
13. `Score`: Confidence score based on signal strength.

## Implementation Notes
- **Coordinate Systems**: Handles conversion between Java's 0-indexed substring logic and genomic 1-based coordinates.
- **Smith-Waterman Parameters**: Uses `(1, -1, -3)` for site clustering and `(1, -1, -1)` for sequence validation to match CIRI3 behavior.
- **Ambiguity Resolution**: Unlike the non-deterministic `HashMap` iteration in Java, the Rust implementation uses a deterministic priority ranking system.
