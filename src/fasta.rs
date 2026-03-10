use std::collections::HashMap;
use anyhow::Result;
use needletail::parse_fastx_file;

/// FASTA Reader Module.
/// Responsible for parsing reference genome sequences and storing them in memory 
/// for rapid retrieval during junction validation and sequence competition checks.
pub struct FastaReader {
    /// Maps from chromosome name (e.g., "chr1") to its full sequence string.
    pub chr_tcga_map: HashMap<String, String>,
    /// Maps from chromosome name to its sequence length.
    pub chr_len_map: HashMap<String, usize>,
}

impl FastaReader {
    /// Initializes an empty FastaReader.
    pub fn new() -> Self {
        Self {
            chr_tcga_map: HashMap::new(),
            chr_len_map: HashMap::new(),
        }
    }

    /// Reads reference sequences from a FASTA file.
    /// Uses the `needletail` library for high-performance parsing of potentially multi-gigabyte reference files.
    /// 
    /// - `fa_file`: Path to the reference genome FASTA file.
    pub fn read_fa(&mut self, fa_file: &str) -> Result<()> {
        let mut reader = parse_fastx_file(fa_file)?;
        
        while let Some(record_res) = reader.next() {
            let record = record_res?;
            
            // Convert byte array record IDs to Strings.
            let id_str = String::from_utf8(record.id().to_vec())?;
            // Use only the first whitespace-separated part of the ID as the chromosome name (standard practice).
            let chr_name = id_str.split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            
            // Normalize sequence: convert to uppercase to ensure uniform comparisons across the tool.
            let seq = String::from_utf8(record.seq().to_vec())?.to_uppercase();
            
            // Store the sequence and its length.
            self.chr_len_map.insert(chr_name.clone(), seq.len());
            self.chr_tcga_map.insert(chr_name, seq);
        }
        Ok(())
    }
}
