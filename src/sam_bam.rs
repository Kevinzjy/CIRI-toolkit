//! Multi-format I/O module: SAM and BAM support.
//!
//! This module provides utilities for detecting the format of genomic alignment files
//! and validating their sort order. CIRI-toolkit requires Name-sorted or Unsorted
//! inputs for correct BSJ identification.

use anyhow::{anyhow, Result};
use noodles::bam;
use std::fs::File;
use std::io::Read;

/// Supported input alignment formats.
pub enum InputFormat {
    /// Sequence Alignment Map (Text-based)
    Sam,
    /// Binary Alignment Map (BGZF-compressed)
    Bam,
}

/// Detects the alignment-file format from the leading magic bytes.
///
/// This check is intentionally lightweight because it runs before any expensive
/// validation or parser setup.
///
/// # Arguments
/// * `path` - Path to the input file.
///
/// # Returns
/// `Ok(InputFormat::Bam)` for BAM/BGZF, `Ok(InputFormat::Sam)` otherwise.
pub fn detect_format(path: &str) -> Result<InputFormat> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 4];
    // Attempt to read the first 4 bytes.
    if file.read_exact(&mut magic).is_err() {
        return Ok(InputFormat::Sam); // Small file or header-only
    }

    // BAM files start with "BAM\x01". BGZF files start with GZIP magic + extra field.
    if &magic == b"BAM\x01" || &magic == b"\x1f\x8b\x08\x04" {
        Ok(InputFormat::Bam)
    } else {
        Ok(InputFormat::Sam)
    }
}

/// Validates that a BAM file is not coordinate-sorted.
///
/// CIRI-toolkit's logic depends on Read 1 and Read 2 being close to each other in the file.
/// Coordinate-sorted files scatter read pairs across the genome, breaking the streaming logic.
/// The check is conservative: false negatives are preferable to silently running
/// the pipeline on a definitely unsupported ordering.
///
/// # Arguments
/// * `path` - Path to the BAM file.
pub fn check_bam_sorting(path: &str) -> Result<()> {
    let file = File::open(path)?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header()?;

    // Convert header to string for robust sort-order check.
    let header_str = format!("{:?}", header);
    if header_str.contains("SO:coordinate") || header_str.contains("SortOrder(coordinate)") {
        return Err(anyhow!("Coordinate-sorted BAM detected. CIRI-toolkit requires Name-sorted or Unsorted BAM. Please use 'samtools sort -n' first."));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_detect_sam() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "@HD\tVN:1.6\tSO:queryname").unwrap();
        assert!(matches!(
            detect_format(file.path().to_str().unwrap()).unwrap(),
            InputFormat::Sam
        ));
    }

    #[test]
    fn test_detect_bam_bgzf() {
        let mut file = NamedTempFile::new().unwrap();
        // Standard BGZF magic bytes.
        file.write_all(&[0x1f, 0x8b, 0x08, 0x04]).unwrap();
        assert!(matches!(
            detect_format(file.path().to_str().unwrap()).unwrap(),
            InputFormat::Bam
        ));
    }

    #[test]
    fn test_detect_bam_magic() {
        let mut file = NamedTempFile::new().unwrap();
        // Standard BAM magic bytes.
        file.write_all(b"BAM\x01").unwrap();
        assert!(matches!(
            detect_format(file.path().to_str().unwrap()).unwrap(),
            InputFormat::Bam
        ));
    }
}
