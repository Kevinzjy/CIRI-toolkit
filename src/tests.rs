#[cfg(test)]
mod tests {
    use crate::scan2::Scan2;
    use crate::is_bsj_hg2::IsBSJHg2;
    use crate::utils::AlignmentRecord;
    use std::collections::HashMap;
    use std::borrow::Cow;
    use std::fs::File;
    use tempfile::tempdir;

    fn mock_record<'a>(flag: i32, chr: &'a str, pos: i32, mq: i32, cigar: &'a str, seq: &'a str) -> AlignmentRecord<'a> {
        AlignmentRecord {
            flag,
            chrom: Cow::Borrowed(chr),
            pos,
            mapq: mq,
            cigar: Cow::Borrowed(cigar),
            seq: Cow::Borrowed(seq),
        }
    }

    #[test]
    fn test_scan2_rescue_simulate_10025() {
        let mut scan2 = Scan2::new(10, 5, 100);
        let bsj1_line = "fake\t88M12S\t1\tchr1\t155252202\t155252631\t+\tAG\tGT\t1\n";
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp, bsj1_line.as_bytes()).unwrap();
        scan2.build_index(tmp.path().to_str().unwrap()).unwrap();

        let mut chr_map = HashMap::new();
        let mut seq = "N".repeat(155252190);
        seq.push_str("NNNNNNNNNNAG");
        seq.push_str(&"N".repeat(420));
        seq.push_str("NNNNNNNNNNGT");
        seq.push_str(&"N".repeat(1000));
        chr_map.insert("chr1".to_string(), seq);

        let read_seq = "N".repeat(100);
        let alignments = vec![
            mock_record(99, "chr1", 155252546, 60, "88M12S", &read_seq),
        ];
        let mut stand_map = HashMap::new();
        stand_map.insert(0, ('0', Cow::Borrowed(read_seq.as_str())));

        let mut results = Vec::new();
        let mut local_fsj = HashMap::new();
        let mut is_bsj_hg2 = IsBSJHg2::new(5, 10);
        
        scan2.process_group_view("simulate:10025", &alignments, &stand_map, &mut results, &mut local_fsj, &chr_map, &mut is_bsj_hg2).unwrap();
        
        assert!(!results.is_empty(), "Rescue failed for 10025.");
    }

    #[test]
    fn test_full_pipeline_mock() {
        use crate::fasta::FastaReader;
        use crate::annotation::Annotation;
        use crate::scan1::Scan1;
        use crate::summary::Summary;
        use std::io::Write;

        let dir = tempdir().unwrap();
        let fa_path = dir.path().join("test.fa");
        let gtf_path = dir.path().join("test.gtf");
        let sam_path = dir.path().join("test.sam");
        let out_prefix = dir.path().join("out").to_str().unwrap().to_string();

        {
            let mut fa_file = File::create(&fa_path).unwrap();
            writeln!(fa_file, ">chr1\n{}", "N".repeat(1000)).unwrap();
            
            let mut gtf_file = File::create(&gtf_path).unwrap();
            writeln!(gtf_file, "chr1\tTEST\texon\t100\t200\t.\t+\t.\tgene_id \"G1\";").unwrap();

            let mut sam_file = File::create(&sam_path).unwrap();
            writeln!(sam_file, "@HD\tVN:1.6\tSO:queryname").unwrap();
        }
        
        let mut fasta = FastaReader::new();
        fasta.read_fasta(fa_path.to_str().unwrap()).unwrap();
        
        let mut annotation = Annotation::new();
        annotation.read_gtf(gtf_path.to_str().unwrap()).unwrap();

        let mut scan1 = Scan1::new(10, 100, 100000, 5);
        let _ = scan1.run(sam_path.to_str().unwrap(), &out_prefix, &fasta.chr_tcga_map, &annotation).unwrap();
        
        let mut scan2 = Scan2::new(10, 5, 100);
        let bsj1_path = format!("{}.BSJ1", out_prefix);
        // Ensure the file exists if scan1 didn't find anything
        if !std::path::Path::new(&bsj1_path).exists() {
            File::create(&bsj1_path).unwrap();
        }
        scan2.build_index(&bsj1_path).unwrap();
        scan2.run(sam_path.to_str().unwrap(), &format!("{}.BSJ2", out_prefix), &fasta.chr_tcga_map).unwrap();

        let mut summary = Summary::new(0);
        summary.run(&bsj1_path, &out_prefix, &scan2.fsj_map, &fasta.chr_tcga_map, &annotation).unwrap();
        
        assert!(std::path::Path::new(&format!("{}.result", out_prefix)).exists());
    }

    #[test]
    #[ignore]
    fn test_integration_real_bam() {
        use crate::fasta::FastaReader;
        use crate::annotation::Annotation;
        use crate::scan1::Scan1;
        use crate::summary::Summary;

        let bam_path = "tests/test.bam";
        let fa_path = "tests/chr1.fa";
        let gtf_path = "tests/chr1.gtf";
        
        let dir = tempdir().unwrap();
        let out_prefix = dir.path().join("integration").to_str().unwrap().to_string();

        if !std::path::Path::new(bam_path).exists() { return; }

        let mut fasta = FastaReader::new();
        fasta.read_fasta(fa_path).unwrap();
        let mut annotation = Annotation::new();
        annotation.read_gtf(gtf_path).unwrap();

        let mut scan1 = Scan1::new(10, 140, 200000, 5);
        scan1.run(bam_path, &out_prefix, &fasta.chr_tcga_map, &annotation).unwrap();
        
        let mut scan2 = Scan2::new(10, 5, 100);
        let bsj1_path = format!("{}.BSJ1", out_prefix);
        scan2.build_index(&bsj1_path).unwrap();
        scan2.run(bam_path, &format!("{}.BSJ2", out_prefix), &fasta.chr_tcga_map).unwrap();

        let mut summary = Summary::new(0);
        summary.run(&bsj1_path, &out_prefix, &scan2.fsj_map, &fasta.chr_tcga_map, &annotation).unwrap();

        assert!(std::path::Path::new(&format!("{}.result", out_prefix)).exists());
    }
}
