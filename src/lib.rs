//! CIRI-toolkit: High-performance circular RNA identification.
//!
//! This crate provides a high-performance Rust implementation of the CIRI3 algorithm
//! for identifying circular RNA back-spliced junctions (BSJ) from SAM or BAM alignment files.
//! It supports multi-threaded processing, deterministic read assignment, and
//! parity-oriented optimization against the Java reference implementation.

pub mod annotation;
pub mod circ_catalog;
pub mod ciri_as;
pub mod cli;
pub mod fasta;
pub mod index_compare;
pub mod is_bsj_hg2;
pub mod misd;
pub mod runtime;
pub mod sam_bam;
pub mod scan1;
pub mod scan2;
pub mod simulator;
pub mod summary;
pub mod tests;
pub mod utils;
