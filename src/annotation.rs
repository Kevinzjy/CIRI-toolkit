//! Annotation module: GTF parsing and exon boundary mapping.
//!
//! This module provides the `Annotation` struct for parsing reference GTF (General Transfer Format)
//! files and providing efficient lookup for exons and gene coordinates.
//!
//! It is used to validate BSJ candidates against known transcript boundaries.

use anyhow::Result;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};

/// One gene span on a chromosome, reconstructed from all exon intervals seen in
/// the GTF.
///
/// Java `Summary.java` uses a chromosome-grouped gene interval list as a
/// fallback when a circRNA does not land exactly on exon boundaries. Rust keeps
/// the same information so final annotation labels can mirror that fallback.
pub struct GeneSpan {
    pub gene_id: String,
    pub start: i32,
    pub end: i32,
}

/// Parsed GTF-derived lookup tables used during BSJ validation and final labeling.
///
/// The current representation is intentionally string-keyed because Scan1 hot
/// paths already build Java-shaped `"chr\tpos"` keys; keeping that contract avoids
/// repeated format conversions during parity-sensitive checks.
pub struct Annotation {
    /// Maps from genomic coordinate string ("Chr\tStart") to its exon metadata ("GeneID\tStrand").
    pub chr_exon_start_map: HashMap<String, String>,
    /// Maps from genomic coordinate string ("Chr\tEnd") to its exon metadata ("GeneID\tStrand").
    pub chr_exon_end_map: HashMap<String, String>,
    /// Numeric exon-start lookup for post-Summary CIRI-AS boundary correction.
    ///
    /// Scan1/Scan2 keep using the Java-shaped string maps above for parity.
    /// Segments reconstruction tests thousands of nearby splice-site candidates,
    /// so it uses this side index to avoid formatting `"chr\tpos"` strings in
    /// the hot loop.
    pub chr_exon_start_index: HashMap<String, HashMap<i32, (String, char)>>,
    /// Numeric exon-end lookup for post-Summary CIRI-AS boundary correction.
    pub chr_exon_end_index: HashMap<String, HashMap<i32, (String, char)>>,
    /// Maps from Gene ID to a vector of its exons (as start and end pairs).
    pub gene_exon_map: HashMap<String, Vec<(i32, i32)>>,
    /// Transcript-consistent introns keyed as `chr\tleft_exon_end\tright_exon_start\tstrand`.
    ///
    /// This is intentionally a side lookup rather than a replacement for the
    /// Java-shaped exon-boundary maps above: CIRI3 parity paths still need the
    /// boundary-only behavior, while CIRI-AS read-level segments can use adjacent
    /// exon pairs from the same transcript to resolve ambiguous splice offsets.
    pub transcript_splice_map: HashSet<String>,
    /// Numeric transcript-consistent splice-pair index for CIRI-AS hot loops.
    pub transcript_splice_index: HashMap<String, HashMap<char, HashSet<(i32, i32)>>>,
    /// Groups genes by chromosome with their outer exon span for Summary's
    /// exon/intron/intergenic fallback labeling.
    pub chr_gene_map: HashMap<String, Vec<GeneSpan>>,
}

impl Annotation {
    /// Initializes an empty annotation store.
    pub fn new() -> Self {
        Self {
            chr_exon_start_map: HashMap::new(),
            chr_exon_end_map: HashMap::new(),
            chr_exon_start_index: HashMap::new(),
            chr_exon_end_index: HashMap::new(),
            gene_exon_map: HashMap::new(),
            transcript_splice_map: HashSet::new(),
            transcript_splice_index: HashMap::new(),
            chr_gene_map: HashMap::new(),
        }
    }

    /// Reads a GTF file and populates the exon-boundary lookup tables.
    ///
    /// Only `exon` features are retained because the rest of the pipeline only
    /// needs exon boundary membership and gene IDs for parity with CIRI3.
    ///
    /// # Arguments
    /// * `annotation_file` - Path to the input GTF file.
    pub fn read_gtf(&mut self, annotation_file: &str) -> Result<()> {
        let file = File::open(annotation_file)?;
        let reader = BufReader::new(file);
        let mut gene_bounds: HashMap<(String, String), (i32, i32)> = HashMap::new();
        let mut transcript_exons: HashMap<(String, String), Vec<(i32, i32, String)>> =
            HashMap::new();
        for line_res in reader.lines() {
            let line = line_res?;
            if line.starts_with('#') {
                continue;
            }

            let parts: Vec<&str> = line.split('\t').collect();
            // Feature type is in the 3rd column (index 2).
            if parts.len() < 9 || parts[2] != "exon" {
                continue;
            }

            let chr = parts[0];
            let start = parts[3].parse::<i32>()?;
            let end = parts[4].parse::<i32>()?;
            let strand = parts[6];
            let attributes = parts[8];

            // Extract 'gene_id' from the semicolon-separated attributes list.
            let gene_id = attributes
                .split(';')
                .filter(|s| s.trim().starts_with("gene_id"))
                .next()
                .and_then(|s| s.split('"').nth(1))
                .unwrap_or("NA")
                .to_string();
            let transcript_id = attributes
                .split(';')
                .find(|s| s.trim().starts_with("transcript_id"))
                .and_then(|s| s.split('"').nth(1))
                .map(str::to_string);

            let gene_id_key = gene_id.clone();
            let start_key = format!("{}\t{}", chr, start);
            let end_key = format!("{}\t{}", chr, end);
            let value = format!("{}\t{}", gene_id, strand);

            self.chr_exon_start_map.insert(start_key, value.clone());
            self.chr_exon_end_map.insert(end_key, value);
            let strand_char = strand.chars().next().unwrap_or('?');
            self.chr_exon_start_index
                .entry(chr.to_string())
                .or_default()
                .insert(start, (gene_id_key.clone(), strand_char));
            self.chr_exon_end_index
                .entry(chr.to_string())
                .or_default()
                .insert(end, (gene_id_key.clone(), strand_char));
            self.gene_exon_map
                .entry(gene_id)
                .or_insert_with(Vec::new)
                .push((start, end));
            if let Some(transcript_id) = transcript_id {
                transcript_exons
                    .entry((chr.to_string(), transcript_id))
                    .or_default()
                    .push((start, end, strand.to_string()));
            }
            let gene_key = (chr.to_string(), gene_id_key);
            gene_bounds
                .entry(gene_key)
                .and_modify(|bounds| {
                    bounds.0 = bounds.0.min(start);
                    bounds.1 = bounds.1.max(end);
                })
                .or_insert((start, end));
        }
        self.transcript_splice_map.clear();
        self.transcript_splice_index.clear();
        for ((chr, _transcript_id), mut exons) in transcript_exons {
            exons.sort_by_key(|(start, end, _strand)| (*start, *end));
            for pair in exons.windows(2) {
                let left = &pair[0];
                let right = &pair[1];
                if left.2 != right.2 {
                    continue;
                }
                self.transcript_splice_map
                    .insert(format!("{}\t{}\t{}\t{}", chr, left.1, right.0, left.2));
                let strand = left.2.chars().next().unwrap_or('?');
                self.transcript_splice_index
                    .entry(chr.clone())
                    .or_default()
                    .entry(strand)
                    .or_default()
                    .insert((left.1, right.0));
            }
        }
        self.chr_gene_map.clear();
        for ((chr, gene_id), (start, end)) in gene_bounds {
            self.chr_gene_map.entry(chr).or_default().push(GeneSpan {
                gene_id,
                start,
                end,
            });
        }
        for spans in self.chr_gene_map.values_mut() {
            spans.sort_by_key(|span| span.start);
        }
        Ok(())
    }

    /// Checks whether a shifted junction pair matches known exon boundaries.
    ///
    /// This helper is mainly useful for tests and experiments; Scan1's hot path
    /// now performs the same lookup inline to stay closer to Java branch order.
    ///
    /// # Arguments
    /// * `chr` - Chromosome name.
    /// * `tmp_site1` - Predicted upstream junction coordinate.
    /// * `tmp_site2` - Predicted downstream junction coordinate.
    /// * `adjt_bp` - Range of shift offsets to check.
    ///
    /// # Returns
    /// `Some(shift, strand, donor_placeholder, acceptor_placeholder)` if a match is found.
    #[allow(dead_code)]
    pub fn find_exon_match(
        &self,
        chr: &str,
        tmp_site1: i32,
        tmp_site2: i32,
        adjt_bp: i32,
    ) -> Option<(i32, String, String, String)> {
        for i in 0..=adjt_bp {
            let s1 = tmp_site1 + i;
            let s2 = tmp_site2 + i;
            let start_key = format!("{}\t{}", chr, s1);
            let end_key = format!("{}\t{}", chr, s2);

            if let (Some(v1), Some(v2)) = (
                self.chr_exon_start_map.get(&start_key),
                self.chr_exon_end_map.get(&end_key),
            ) {
                let p1: Vec<&str> = v1.split('\t').collect();
                let p2: Vec<&str> = v2.split('\t').collect();

                if p1[0] == p2[0] {
                    return Some((i, p1[1].to_string(), "GT".to_string(), "AG".to_string()));
                }
            }
        }
        None
    }
}
