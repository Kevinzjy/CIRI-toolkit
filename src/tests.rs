#[cfg(test)]
mod tests {
    use crate::scan2::Scan2;
    use crate::is_bsj_hg2::IsBSJHg2;
    use std::collections::HashMap;

    fn mock_aln(id: &str, flag: &str, chr: &str, pos: &str, mq: &str, cigar: &str, seq: &str) -> Vec<String> {
        vec![id.to_string(), flag.to_string(), chr.to_string(), pos.to_string(), mq.to_string(), 
             cigar.to_string(), "*".to_string(), "*".to_string(), "0".to_string(), seq.to_string(), "*".to_string()]
    }

    #[test]
    fn test_scan2_rescue_simulate_10025() {
        let mut scan2 = Scan2::new(10, 5, 100);
        
        // 构造索引： site2 = 155252631
        let bsj1_line = "fake\t88M12S\t1\tchr1\t155252202\t155252631\t+\tAG\tGT\t1\n";
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp, bsj1_line.as_bytes()).unwrap();
        scan2.build_index(tmp.path().to_str().unwrap()).unwrap();

        let mut chr_map = HashMap::new();
        // simulate:10025 chr1 155252202 155252631 (+) AG-GT
        // site1=155252202, site2=155252631
        let mut seq = "N".repeat(155252190);
        seq.push_str("NNNNNNNNNNAG"); // site1 (155252202) -> index 155252201
        seq.push_str(&"N".repeat(420));
        seq.push_str("NNNNNNNNNNGT"); // site2 (155252631) -> index 155252630
        seq.push_str(&"N".repeat(1000));
        chr_map.insert("chr1".to_string(), seq);

        // Read 序列: GCAT...GAGG TTTTACTGGGACC
        // 88M 结束于 index 87 (POS 155252633). 
        // 信号 AG 在 POS 155252631. 
        // 所以我们要在 index 85 截取。 bias = -2.
        let read_seq = "GCATCGTTCGCTTCACCAAGATCCTAAGCCTGATGAGGCTGCTCCGCCTCTCCCGCCTCATCCGCTACATACACCAGTGGGAGGAGGTTTTACTGGGACC";
        let alignments = vec![
            mock_aln("simulate:10025", "147", "chr1", "155252546", "60", "88M12S", read_seq),
        ];
        let mut stand_map = HashMap::new();
        stand_map.insert(0, "0".to_string() + read_seq); // Pos strand flag 0

        let mut output = Vec::new();
        let mut is_bsj_hg2 = IsBSJHg2::new(5, 10);
        scan2.process_read("simulate:10025", &alignments, &stand_map, &mut output, &chr_map, &mut is_bsj_hg2).unwrap();
        let result = String::from_utf8(output).unwrap();
        
        assert!(result.contains("simulate:10025"), "Rescue failed for 10025. Result was: {}", result);
    }
}
