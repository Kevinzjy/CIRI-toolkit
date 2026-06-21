//! FASTA module: High-speed sequence loading with Needletail.
//!
//! This module provides the `FastaReader` struct for loading reference genome sequences
//! into memory for rapid retrieval during junction validation.

use anyhow::Result;
use needletail::parse_fastx_file;
use std::collections::HashMap;

/// In-memory reference FASTA store.
///
/// The full chromosome sequence is kept as uppercase text because both Scan1 and
/// Scan2 perform many short substring comparisons, and the reference size used by
/// CIRI is small enough that this tradeoff is simpler and faster than on-demand I/O.
pub struct FastaReader {
    /// Maps from chromosome name to its sequence length.
    pub chr_len_map: HashMap<String, usize>,
    /// Maps from chromosome name (e.g., "chr1") to its full uppercase sequence string.
    pub chr_tcga_map: HashMap<String, String>,
}

impl FastaReader {
    /// Initializes an empty FASTA store.
    pub fn new() -> Self {
        Self {
            chr_len_map: HashMap::new(),
            chr_tcga_map: HashMap::new(),
        }
    }

    /// Reads reference sequences from a FASTA file.
    ///
    /// `needletail` handles parsing, while this wrapper normalizes names and case
    /// so later junction-validation code can use cheap byte/substring comparisons.
    ///
    /// # Arguments
    /// * `fasta_file` - Path to the reference genome FASTA file.
    pub fn read_fasta(&mut self, fasta_file: &str) -> Result<()> {
        let mut reader = parse_fastx_file(fasta_file)?;
        while let Some(record) = reader.next() {
            let seq_record = record?;
            // Use the first whitespace-separated part of the ID as the chromosome name.
            let chr_name = String::from_utf8_lossy(
                seq_record.id().split(|&b| b == b' ').next().unwrap_or(b"*"),
            )
            .to_string();
            let seq = String::from_utf8_lossy(seq_record.seq().as_ref())
                .to_string()
                .to_uppercase();

            self.chr_len_map.insert(chr_name.clone(), seq.len());
            self.chr_tcga_map.insert(chr_name, seq);
        }
        Ok(())
    }
}

impl Default for FastaReader {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_fasta_reader() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1 test\nATGCATGC\n>chr2\nGGGGCCCC").unwrap();

        let mut reader = FastaReader::new();
        reader.read_fasta(file.path().to_str().unwrap()).unwrap();

        assert_eq!(reader.chr_tcga_map.get("chr1").unwrap(), "ATGCATGC");
        assert_eq!(reader.chr_tcga_map.get("chr2").unwrap(), "GGGGCCCC");
        assert_eq!(reader.chr_tcga_map.len(), 2);
    }
}
