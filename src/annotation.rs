use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use anyhow::Result;

/// GTF Annotation Module.
/// Parses reference GTF (General Transfer Format) files and provides efficient 
/// lookup for exons and gene coordinates.
/// 
/// This module is crucial for rescuing BSJ candidates using known transcript boundaries, 
/// especially when splicing motifs are not clearly identifiable.
pub struct Annotation {
    /// Maps from genomic coordinate string ("Chr\tStart") to its exon metadata ("GeneID\tStrand").
    pub chr_exon_start_map: HashMap<String, String>,
    /// Maps from genomic coordinate string ("Chr\tEnd") to its exon metadata ("GeneID\tStrand").
    pub chr_exon_end_map: HashMap<String, String>,
    /// Maps from Gene ID to a vector of its exons (as start and end pairs).
    pub gene_exon_map: HashMap<String, Vec<(i32, i32)>>,
}

impl Annotation {
    /// Initializes an empty Annotation instance.
    pub fn new() -> Self {
        Self {
            chr_exon_start_map: HashMap::new(),
            chr_exon_end_map: HashMap::new(),
            gene_exon_map: HashMap::new(),
        }
    }

    /// Reads a GTF file and populates the internal maps.
    /// It filters for lines with the 'exon' feature and extracts relevant gene information.
    /// 
    /// - `annotation_file`: Path to the input GTF file.
    pub fn read_gtf(&mut self, annotation_file: &str) -> Result<()> {
        let file = File::open(annotation_file)?;
        let reader = BufReader::new(file);
        for line_res in reader.lines() {
            let line = line_res?;
            // Skip comments and header lines
            if line.starts_with('#') { continue; } 
            
            let parts: Vec<&str> = line.split('\t').collect();
            // Expected GTF format has 9 columns; column 3 (index 2) is the feature type.
            if parts.len() < 9 || parts[2] != "exon" { continue; }
            
            let chr = parts[0];
            let start = parts[3].parse::<i32>()?;
            let end = parts[4].parse::<i32>()?;
            let strand = parts[6];
            let attributes = parts[8];
            
            // Extract 'gene_id' from the semicolon-separated attributes list (column 9).
            let gene_id = attributes.split(';')
                .filter(|s| s.trim().starts_with("gene_id"))
                .next()
                .and_then(|s| s.split('"').nth(1))
                .unwrap_or("NA")
                .to_string();
            
            let start_key = format!("{}\t{}", chr, start);
            let end_key = format!("{}\t{}", chr, end);
            let value = format!("{}\t{}", gene_id, strand);
            
            // Store coordinates for rapid boundary lookup.
            self.chr_exon_start_map.insert(start_key, value.clone());
            self.chr_exon_end_map.insert(end_key, value);
            // Group exons by gene ID.
            self.gene_exon_map.entry(gene_id).or_insert_with(Vec::new).push((start, end));
        }
        Ok(())
    }

    /// Verifies if a given junction site pair matches known exon start and end boundaries 
    /// within a small shift offset. This allows "rescuing" BSJs that align perfectly to known exons.
    /// 
    /// - `chr`: Chromosome name.
    /// - `tmp_site1`: Predicted upstream junction coordinate.
    /// - `tmp_site2`: Predicted downstream junction coordinate.
    /// - `adjt_bp`: Range of shift offsets to check (boundary flexibility).
    /// 
    /// Returns: `Some(shift, strand, donor_placeholder, acceptor_placeholder)` if a match is found.
    #[allow(dead_code)]
    pub fn find_exon_match(&self, chr: &str, tmp_site1: i32, tmp_site2: i32, adjt_bp: i32) -> Option<(i32, String, String, String)> {
        for i in 0..=adjt_bp {
            let s1 = tmp_site1 + i;
            let s2 = tmp_site2 + i;
            let start_key = format!("{}\t{}", chr, s1);
            let end_key = format!("{}\t{}", chr, s2);
            
            // Check if both coordinates are existing exon boundaries in the same gene.
            if let (Some(v1), Some(v2)) = (self.chr_exon_start_map.get(&start_key), self.chr_exon_end_map.get(&end_key)) {
                let p1: Vec<&str> = v1.split('\t').collect();
                let p2: Vec<&str> = v2.split('\t').collect();
                
                // Confirm both boundaries belong to the same Gene ID.
                if p1[0] == p2[0] {
                    // Boundary match confirmed. Use standard GT/AG signals as placeholders.
                    return Some((i, p1[1].to_string(), "GT".to_string(), "AG".to_string()));
                }
            }
        }
        None
    }
}
