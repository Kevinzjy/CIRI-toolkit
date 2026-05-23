//! External circRNA catalog parsing for multi-sample second-pass runs.
//!
//! The catalog is intentionally narrower than CIRI `.out`: `ciri-merge` writes
//! BED6 sites that define cohort-level BSJ candidate boundaries, while each
//! sample still rescans its own BAM/SAM to decide whether those candidates have
//! read support. Keeping this as a typed module avoids leaking BED coordinate
//! conventions into Scan2's Java-compatible candidate payload.

use anyhow::{bail, Context, Result};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// One external BSJ candidate from a BED6 circRNA catalog.
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub struct CircCatalogRecord {
    /// Reference sequence name.
    pub chrom: String,
    /// CIRI-style 1-based BSJ start coordinate.
    pub start: i32,
    /// CIRI-style 1-based BSJ end coordinate.
    pub end: i32,
    /// Strand label carried by BED6 (`+`, `-`, or `NA` when unknown).
    pub strand: String,
}

/// Reads a BED6 circRNA catalog produced by `ciri-merge`.
///
/// BED uses 0-based half-open coordinates, while Scan2 candidate records use the
/// same 1-based coordinates as CIRI `.out`. The conversion is performed exactly
/// once here so downstream code can work in CIRI coordinates.
pub fn read_circ_bed6(path: &str) -> Result<Vec<CircCatalogRecord>> {
    let file = File::open(path).with_context(|| format!("open circ catalog {path}"))?;
    let reader = BufReader::with_capacity(1024 * 1024, file);
    let mut records = Vec::new();
    let mut seen = HashSet::new();
    for (line_idx, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track")
            || trimmed.starts_with("browser")
        {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 6 {
            bail!(
                "{}:{} has {} columns, expected BED6",
                path,
                line_idx + 1,
                fields.len()
            );
        }
        let chrom_start = fields[1]
            .parse::<i64>()
            .with_context(|| format!("parse chromStart at {}:{}", path, line_idx + 1))?;
        let chrom_end = fields[2]
            .parse::<i64>()
            .with_context(|| format!("parse chromEnd at {}:{}", path, line_idx + 1))?;
        if chrom_start < 0 || chrom_end <= chrom_start {
            bail!(
                "{}:{} has invalid BED interval {}-{}",
                path,
                line_idx + 1,
                chrom_start,
                chrom_end
            );
        }
        let start = chrom_start + 1;
        if start > i32::MAX as i64 || chrom_end > i32::MAX as i64 {
            bail!(
                "{}:{} exceeds supported i32 coordinate range",
                path,
                line_idx + 1
            );
        }
        let record = CircCatalogRecord {
            chrom: fields[0].to_string(),
            start: start as i32,
            end: chrom_end as i32,
            strand: normalize_bed_strand(fields[5]).to_string(),
        };
        if seen.insert((
            record.chrom.clone(),
            record.start,
            record.end,
            record.strand.clone(),
        )) {
            records.push(record);
        }
    }
    Ok(records)
}

/// Normalizes BED strand labels without inventing biological direction.
fn normalize_bed_strand(raw: &str) -> &str {
    match raw {
        "+" | "-" => raw,
        _ => "NA",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn bed6_catalog_converts_to_ciri_coordinates() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "chr1\t9\t20\tchr1:10|20:+\t7\t+").unwrap();
        writeln!(file, "chr1\t9\t20\tduplicate\t7\t+").unwrap();
        writeln!(file, "chr2\t29\t40\tunknown\t3\t.").unwrap();

        let records = read_circ_bed6(file.path().to_str().unwrap()).unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(
            records[0],
            CircCatalogRecord {
                chrom: "chr1".to_string(),
                start: 10,
                end: 20,
                strand: "+".to_string(),
            }
        );
        assert_eq!(records[1].strand, "NA");
    }
}
