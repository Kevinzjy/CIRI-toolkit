//! CIRI-toolkit: High-performance circular RNA identification.
//!
//! This crate provides a high-performance Rust implementation of the CIRI3 algorithm
//! for identifying circular RNA back-spliced junctions (BSJ) from SAM or BAM alignment files.
//! It supports multi-threaded processing, zero-copy parsing, and deterministic read assignment.

pub mod annotation;
pub mod fasta;
pub mod index_compare;
pub mod is_bsj_hg2;
pub mod misd;
pub mod scan1;
pub mod scan2;
pub mod summary;
pub mod sam_bam;
pub mod utils;
pub mod tests;
