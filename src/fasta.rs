/// FASTA module: High-speed sequence loading with Needletail.
/// Handles multi-line FASTA and produces a global reference map.

use std::collections::HashMap;
use anyhow::Result;
use needletail::parse_fastx_file;

pub struct FastaReader {
    pub chr_len_map: HashMap<String, usize>,
    pub chr_tcga_map: HashMap<String, String>,
}

impl FastaReader {
    pub fn new() -> Self {
        Self {
            chr_len_map: HashMap::new(),
            chr_tcga_map: HashMap::new(),
        }
    }

    pub fn read_fasta(&mut self, fasta_file: &str) -> Result<()> {
        let mut reader = parse_fastx_file(fasta_file)?;
        while let Some(record) = reader.next() {
            let seq_record = record?;
            let chr_name = String::from_utf8_lossy(seq_record.id().split(|&b| b == b' ').next().unwrap_or(b"*")).to_string();
            let seq = String::from_utf8_lossy(seq_record.seq().as_ref()).to_string().to_uppercase();
            
            // Store the sequence and its length.
            self.chr_len_map.insert(chr_name.clone(), seq.len());
            self.chr_tcga_map.insert(chr_name, seq);
        }
        Ok(())
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
