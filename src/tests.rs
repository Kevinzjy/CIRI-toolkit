//! High-level integration-oriented tests for the public pipeline stages.

#[cfg(test)]
mod tests {
    use crate::scan2::Scan2;
    use crate::utils::{
        bsj1_path_for_output, bsj2_path_for_output, bsj_path_for_output, fsj_path_for_output,
        result_path_for_output,
    };
    use std::fs::File;
    use tempfile::tempdir;

    #[test]
    /// Verifies that the three pipeline stages can run end-to-end on a tiny mock input.
    fn test_full_pipeline_mock() {
        use crate::annotation::Annotation;
        use crate::fasta::FastaReader;
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
            writeln!(
                gtf_file,
                "chr1\tTEST\texon\t100\t200\t.\t+\t.\tgene_id \"G1\";"
            )
            .unwrap();

            let mut sam_file = File::create(&sam_path).unwrap();
            writeln!(sam_file, "@HD\tVN:1.6\tSO:queryname").unwrap();
        }

        let mut fasta = FastaReader::new();
        fasta.read_fasta(fa_path.to_str().unwrap()).unwrap();

        let mut annotation = Annotation::new();
        annotation.read_gtf(gtf_path.to_str().unwrap()).unwrap();

        let mut scan1 = Scan1::new(10, 100, 100000, 5);
        let _ = scan1
            .run(
                sam_path.to_str().unwrap(),
                &out_prefix,
                &fasta.chr_tcga_map,
                &annotation,
            )
            .unwrap();

        let mut scan2 = Scan2::new(10, 5, 100);
        let bsj1_path = bsj1_path_for_output(&out_prefix);
        let bsj2_path = bsj2_path_for_output(&out_prefix);
        let bsj_path = bsj_path_for_output(&out_prefix);
        let fsj_path = fsj_path_for_output(&out_prefix);
        let result_path = result_path_for_output(&out_prefix);
        // Ensure the file exists if scan1 didn't find anything
        if !std::path::Path::new(&bsj1_path).exists() {
            File::create(&bsj1_path).unwrap();
        }
        scan2.build_index(&bsj1_path).unwrap();
        scan2
            .run(
                sam_path.to_str().unwrap(),
                &bsj1_path,
                &bsj_path,
                &bsj2_path,
                &fsj_path,
                &fasta.chr_tcga_map,
            )
            .unwrap();

        let mut summary = Summary::new(0);
        summary
            .run(
                &bsj_path,
                &result_path,
                &scan2.fsj_map,
                &fasta.chr_tcga_map,
                &annotation,
            )
            .unwrap();

        assert!(std::path::Path::new(&result_path).exists());
    }

    #[test]
    #[ignore]
    /// Smoke test for a real BAM input when local test fixtures are available.
    fn test_integration_real_bam() {
        use crate::annotation::Annotation;
        use crate::fasta::FastaReader;
        use crate::scan1::Scan1;
        use crate::summary::Summary;

        let bam_path = "tests/test.bam";
        let fa_path = "tests/chr1.fa";
        let gtf_path = "tests/chr1.gtf";

        let dir = tempdir().unwrap();
        let out_prefix = dir.path().join("integration").to_str().unwrap().to_string();

        if !std::path::Path::new(bam_path).exists() {
            return;
        }

        let mut fasta = FastaReader::new();
        fasta.read_fasta(fa_path).unwrap();
        let mut annotation = Annotation::new();
        annotation.read_gtf(gtf_path).unwrap();

        let mut scan1 = Scan1::new(10, 140, 200000, 5);
        scan1
            .run(bam_path, &out_prefix, &fasta.chr_tcga_map, &annotation)
            .unwrap();

        let mut scan2 = Scan2::new(10, 5, 100);
        let bsj1_path = bsj1_path_for_output(&out_prefix);
        let bsj2_path = bsj2_path_for_output(&out_prefix);
        let bsj_path = bsj_path_for_output(&out_prefix);
        let fsj_path = fsj_path_for_output(&out_prefix);
        let result_path = result_path_for_output(&out_prefix);
        scan2.build_index(&bsj1_path).unwrap();
        scan2
            .run(
                bam_path,
                &bsj1_path,
                &bsj_path,
                &bsj2_path,
                &fsj_path,
                &fasta.chr_tcga_map,
            )
            .unwrap();

        let mut summary = Summary::new(0);
        summary
            .run(
                &bsj_path,
                &result_path,
                &scan2.fsj_map,
                &fasta.chr_tcga_map,
                &annotation,
            )
            .unwrap();

        assert!(std::path::Path::new(&result_path).exists());
    }
}
