use std::collections::BTreeMap;

/// Splicing Signal Comparator.
/// This module is responsible for identifying biological splicing signals (motifs) 
/// in the genomic sequences flanking a predicted Back-Spliced Junction (BSJ).
pub struct IndexCompare;

impl IndexCompare {
    /// Compares sequences for standard (canonical and semi-canonical) splicing signals.
    /// This is used for nuclear chromosomes.
    /// 
    /// - `end_string1`: Sequence flanking the predicted upstream site (acceptor side).
    /// - `end_string2`: Sequence flanking the predicted downstream site (donor side).
    /// 
    /// Returns a map of shift offsets to a formatted string containing:
    /// `offset, strand, donor_signal, acceptor_signal`.
    pub fn index_compare(end_string1: &str, end_string2: &str) -> BTreeMap<i32, String> {
        let mut index_strand_map = BTreeMap::new();
        let up_s1 = end_string1.to_uppercase();
        let up_s2 = end_string2.to_uppercase();
        
        // Canonical splicing motifs:
        // For '+' strand: Donor is GT, Acceptor is AG.
        // For '-' strand: Donor (complement) is CT, Acceptor (complement) is AC.
        // bibases[0] contains upstream motifs, bibases[1] contains downstream motifs.
        let bibases = [
            ["AC", "AG"], // motifs for '-' and '+' strands respectively
            ["CT", "GT"]  
        ];
        let strand_index = ["-", "+"];

        for i in 0..=1 {
            let mut pre_index: i32 = -1;
            loop {
                let target = bibases[0][i];
                let search_start = (pre_index + 1) as usize;
                if search_start >= up_s1.len() { break; }
                
                // Search for the first motif in the upstream flanking sequence.
                if let Some(index) = up_s1[search_start..].find(target) {
                    let actual_index = (search_start + index) as i32;
                    // Verify if the matching second motif exists in the downstream sequence at the same relative position.
                    if up_s2.len() >= (actual_index + 2) as usize {
                        let s2_part = &up_s2[actual_index as usize..(actual_index + 2) as usize];
                        if s2_part.eq_ignore_ascii_case(bibases[1][i]) {
                            // Motif pair found! Store it with the relative shift (offset).
                            index_strand_map.insert(actual_index, format!("{}\t{}\t{}\t{}", actual_index, strand_index[i], target, bibases[1][i]));
                        }
                    } else {
                        break;
                    }
                    pre_index = actual_index;
                } else {
                    break;
                }
            }
        }
        index_strand_map
    }

    /// Compares sequences for mitochondrial or specialized splicing signals.
    /// Mitochondrial DNA (MT) has less constrained splicing, so a wider variety of motifs is checked.
    /// 
    /// Returns the same map format as `index_compare`.
    pub fn index_compare_chrm(end_string1: &str, end_string2: &str) -> BTreeMap<i32, String> {
        let mut index_strand_map = BTreeMap::new();
        let up_s1 = end_string1.to_uppercase();
        let up_s2 = end_string2.to_uppercase();
        
        // Expanded set of splicing motifs for mitochondrial DNA or specialized cases.
        // Includes non-canonical pairs like GC/AG, AT/AC, etc.
        let bibases_mut = [
            ["AC", "AG", "GC", "AG", "AT", "AC", "AT", "AG"],
            ["CT", "GT", "CT", "GC", "GT", "AT", "CT", "AT"]
        ];
        let strand_index = ["-", "+"];

        for i in 0..=7 {
            let mut pre_index: i32 = -1;
            let strand = i % 2; // Alternates between '-' (even i) and '+' (odd i) strands in this array.
            loop {
                let target = bibases_mut[0][i];
                let search_start = (pre_index + 1) as usize;
                if search_start >= up_s1.len() { break; }
                
                if let Some(index) = up_s1[search_start..].find(target) {
                    let actual_index = (search_start + index) as i32;
                    if up_s2.len() >= (actual_index + 2) as usize {
                        let s2_part = &up_s2[actual_index as usize..(actual_index + 2) as usize];
                        if s2_part.eq_ignore_ascii_case(bibases_mut[1][i]) {
                            index_strand_map.insert(actual_index, format!("{}\t{}\t{}\t{}", actual_index, strand_index[strand], target, bibases_mut[1][i]));
                        }
                    } else {
                        break;
                    }
                    pre_index = actual_index;
                } else {
                    break;
                }
            }
        }
        index_strand_map
    }
}
