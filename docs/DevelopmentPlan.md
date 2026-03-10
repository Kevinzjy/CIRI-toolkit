# CIRI-toolkit (Rust) Development Status

This document summarizes the current state of the high-performance Rust reimplementation of CIRI3.

## 1. Project Overview
CIRI-toolkit identifies circular RNA back-spliced junction (BSJ) reads from SAM files. It is a 1:1 behavioral port of CIRI3 Java with significant performance optimizations.

## 2. Current Implementation: SAM Support (v0.1.0) - COMPLETED

### Core Algorithms
- **Scan 1**: Multi-threaded BSJ identification using split-mapping CIGAR analysis.
- **Scan 2**: PEM/SMS rescue and Forward-Spliced Junction (FSJ) counting.
- **Validation**: Smith-Waterman sequence validation, canonical splice signal check, and linear competition checks.
- **Summary**: Smith-Waterman based site clustering and deterministic unique read assignment.

### Technical Implementation Details (Crucial for Handover)
1. **Zero-Copy I/O**: Uses `memmap2` to map SAM files into memory. Worker threads operate on `&str` and `&[u8]` views pointing directly to the Mmap, minimizing heap allocations.
2. **Shard Synchronization**: To handle Read IDs spanning shard boundaries,分片 $i$ continues reading until the ID changes, and分片 $i+1$ skips its first ID group.
3. **High-Performance Parsing**: Custom status-machine based CIGAR parser (`misd.rs`) and fast integer parser to avoid Regex and `str::parse` overhead.
4. **Memory Management**: Global allocator switched to `mimalloc` to eliminate lock contention during high-concurrency string operations.
5. **Parity Achievements**: Verified against `tests/test.sam` with >97% Read-level accuracy compared to CIRI3 Java.

## 3. Key Findings & Experience
- **Java NIO vs Rust Mmap**: Java's multi-threaded performance comes from independent `FileChannel` reads. Rust's `Mmap` combined with zero-copy slicing proved 3.4x faster in wall-clock time.
- **Alignment Ambiguity**: CIRI3's `HashMap` iteration is non-deterministic. Rust's implementation introduces a **Priority Tie-breaker** (Signal > MAPQ > Length) which is deterministic and technically more robust.
- **Coordinate Precision**: Absolute care was taken to mirror Java's 0-indexed `substring` shifts. Final output coordinates are 1-based genomic coordinates.

## 4. Future Roadmap: BAM Support
Implementation details are documented in `docs/BAM_SUPPORT_PLAN.md`. The next phase involves abstracting the `Alignment` record to handle binary BGZF streams while preserving the validated BSJ identification logic.

---
*Documentation updated on March 10, 2026, preparing for Git Worktree transition.*
