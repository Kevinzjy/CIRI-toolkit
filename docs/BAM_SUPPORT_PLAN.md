# BAM Support Implementation Plan

This document serves as a handover for the next development phase: adding native BAM support to `CIRI-toolkit`.

## 1. Architectural Strategy: The "Abstract Record" Pattern
To support both SAM and BAM without duplicating the core BSJ identification logic, we should transition to a trait-based abstraction.

### Proposed Trait: `Alignment`
```rust
pub trait Alignment {
    fn read_id(&self) -> &str;
    fn flag(&self) -> i32;
    fn chrom(&self) -> &str;
    fn pos(&self) -> i32;
    fn mapq(&self) -> i32;
    fn cigar(&self) -> &str;
    fn seq(&self) -> &str;
}
```

- **SAM Implementation**: Currently handled by `AlignmentView` (zero-copy from Mmap).
- **BAM Implementation**: Wrapper around `noodles_bam::Record`.

## 2. Key Technical Challenges for BAM

### A. Parallel I/O (Sharding)
- **Current SAM Logic**: Uses physical file offsets (`file_size / threads`).
- **BAM Constraint**: BAM is BGZF-compressed. Physical offsets might break in the middle of a compressed block.
- **Solution**: 
    1. **Option 1 (Indexed)**: Use `.bai` index to get virtual offsets (`VirtualPosition`) for genomic regions.
    2. **Option 2 (Streaming)**: Use a single producer thread to decompress BGZF blocks and multiple worker threads to process decoded records.

### B. CIGAR Parsing
- BAM stores CIGAR as a series of 32-bit integers (`op_len << 4 | op_kind`). 
- **Optimization**: Our current `misd` function parses strings. For BAM, we can implement a `misd_raw(cigar: &[u32])` to avoid string conversion entirely, further boosting performance.

## 3. Preservation of Parity
Any implementation of BAM support **MUST** maintain the following logic found in `is_bsj_hg2.rs` and `scan1.rs`:
1. **Coordinate Systems**: Ensure 1-based genomic coordinates are correctly mapped to 0-based sequence extraction.
2. **Mate-Check (`str4_ok`)**: The logic for verifying if the mate read falls within the circular range must be preserved.
3. **Priority Ranking**: The `Summary` module's ranking system (ScanType > Signal > MQ > Length) is file-format agnostic and should remain the final arbiter.

## 4. Suggested Libraries
- **`noodles-bam`**: Preferred for its pure-Rust implementation and alignment with the modern Rust ecosystem.
- **`noodles-bgzf`**: For handling the underlying compressed stream if custom sharding is needed.

---
*Prepared for the next Agent on March 10, 2026.*
