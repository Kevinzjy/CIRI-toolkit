use std::fs::File;
use std::io::{Read};
use anyhow::{Result, anyhow};
use noodles::bam;
use noodles::sam;

pub enum InputFormat {
    Sam,
    Bam,
}

pub fn detect_format(path: &str) -> Result<InputFormat> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 4];
    file.read_exact(&mut magic)?;
    
    if &magic == b"BAM\x01" || &magic == b"\x1f\x8b\x08\x04" {
        Ok(InputFormat::Bam)
    } else {
        Ok(InputFormat::Sam)
    }
}

pub fn check_bam_sorting(path: &str) -> Result<()> {
    let file = File::open(path)?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header()?;
    
    // Check Sort Order in header using string search if structured access is failing
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
        assert!(matches!(detect_format(file.path().to_str().unwrap()).unwrap(), InputFormat::Sam));
    }

    #[test]
    fn test_detect_bam_bgzf() {
        let mut file = NamedTempFile::new().unwrap();
        // Standard BGZF magic bytes
        file.write_all(&[0x1f, 0x8b, 0x08, 0x04]).unwrap();
        assert!(matches!(detect_format(file.path().to_str().unwrap()).unwrap(), InputFormat::Bam));
    }

    #[test]
    fn test_detect_bam_magic() {
        let mut file = NamedTempFile::new().unwrap();
        // Standard BAM magic bytes
        file.write_all(b"BAM\x01").unwrap();
        assert!(matches!(detect_format(file.path().to_str().unwrap()).unwrap(), InputFormat::Bam));
    }
}
