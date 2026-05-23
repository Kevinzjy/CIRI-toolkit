use pretty_assertions::assert_eq;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use super::{
    format_segments, is_outward_facing_truth, run, sample_circular_insert_len, sample_insert_len,
    select_fastq_compressor_from_availability, FastqCompressor, Lcg64, SimulateArgs, SourceBase,
};

const ISOFORM_HEADER: &str = "circ_id\tchrom\tstart\tend\tstrand\tgene_id\ttranscript_id\tcoverage\tread_cnt\tbsj_read_cnt\tisoform_cnt\tisoform_exons\tisoform_len\tisoform_read_cnt\tisoform_bsj_read_cnt";
const READS_HEADER: &str = "read_id\tcirc_id\tchrom\tstart\tend\tstrand\tisoform_id\tis_circular\tis_bsj\tr1_segments\tr1_is_bsj\tr2_segments\tr2_is_bsj\ttype";
const USAGE_HEADER: &str =
    "circ_id\tchrom\tstart\tend\tstrand\tisoform_id\tusage\tread_cnt\tbsj_read_cnt";

fn base_args() -> SimulateArgs {
    SimulateArgs {
        ref_fasta: "unused.fa".to_string(),
        gtf: "unused.gtf".to_string(),
        out_prefixes: vec!["unused".to_string()],
        chrom: None,
        circ_count: 1,
        switching_event: 0,
        circ_coverage: 10.0,
        linear_coverage: 0.1,
        scale: 0.5,
        read_len: 150,
        insert_len: 260,
        insert_sd: 40.0,
        insert_len_minor: 420,
        insert_sd_minor: 60.0,
        minor_insert_fraction: 0.10,
        error_rate: 0.002,
        exon_exclusive_rate: 0.25,
        seed: 5,
    }
}

#[test]
fn selects_gzip_when_pigz_is_unavailable() {
    assert_eq!(
        select_fastq_compressor_from_availability(false, true).unwrap(),
        FastqCompressor::Gzip
    );
}

#[test]
fn fails_when_no_fastq_compressor_exists() {
    let err = select_fastq_compressor_from_availability(false, false).unwrap_err();
    assert!(err.to_string().contains("neither pigz nor gzip"));
}

#[test]
fn short_circ_sampling_biases_insert_toward_circ_length() {
    let args = base_args();
    let mut circ_rng = Lcg64::new(7);
    let mut generic_rng = Lcg64::new(7);
    let short_circ_len = 220usize;
    let samples = 2048usize;

    let circ_mean = (0..samples)
        .map(|_| sample_circular_insert_len(short_circ_len, &args, &mut circ_rng) as f64)
        .sum::<f64>()
        / samples as f64;
    let generic_mean = (0..samples)
        .map(|_| sample_insert_len(&args, &mut generic_rng) as f64)
        .sum::<f64>()
        / samples as f64;

    assert!(
        circ_mean < generic_mean - 20.0,
        "short circ insert mean should shift below generic PE mixture: circ_mean={circ_mean}, generic_mean={generic_mean}"
    );
    assert!(
        (circ_mean - short_circ_len as f64).abs() < 30.0,
        "short circ insert mean should stay close to circ length: circ_mean={circ_mean}, circ_len={short_circ_len}"
    );
}

#[test]
fn simulator_segments_keep_read_chain_bsj_order() {
    let source_map = vec![
        SourceBase {
            coord: 100,
            exon_idx: 0,
        },
        SourceBase {
            coord: 101,
            exon_idx: 0,
        },
        SourceBase {
            coord: 200,
            exon_idx: 1,
        },
        SourceBase {
            coord: 201,
            exon_idx: 1,
        },
    ];

    let plus = format_segments(&source_map, '+', 2, 3, true);
    assert_eq!(plus.text, "200-201:+|<bsj>|100-100:+");
    assert!(plus.is_bsj);

    let minus_source_map = vec![
        SourceBase {
            coord: 201,
            exon_idx: 1,
        },
        SourceBase {
            coord: 200,
            exon_idx: 1,
        },
        SourceBase {
            coord: 101,
            exon_idx: 0,
        },
        SourceBase {
            coord: 100,
            exon_idx: 0,
        },
    ];
    let minus = format_segments(&minus_source_map, '-', 2, 3, true);
    assert_eq!(minus.text, "100-101:-|<bsj>|201-201:-");
    assert!(minus.is_bsj);
}

#[test]
fn simulator_outward_truth_requires_clear_pair_offset() {
    let source_map: Vec<SourceBase> = (100..300)
        .map(|coord| SourceBase { coord, exon_idx: 0 })
        .collect();

    assert!(is_outward_facing_truth(&source_map, '-', 40, 70, 100));
    assert!(is_outward_facing_truth(&source_map, '-', 0, 100, 50));
    assert!(!is_outward_facing_truth(&source_map, '-', 40, 55, 100));
    assert!(!is_outward_facing_truth(&source_map, '-', 40, 40, 100));
    assert!(!is_outward_facing_truth(&source_map, '+', 0, 140, 100));
}

#[test]
fn long_circ_sampling_matches_generic_insert_distribution() {
    let args = base_args();
    let mut circ_rng = Lcg64::new(11);
    let mut generic_rng = Lcg64::new(11);

    for _ in 0..256 {
        assert_eq!(
            sample_circular_insert_len(480, &args, &mut circ_rng),
            sample_insert_len(&args, &mut generic_rng)
        );
    }
}

#[test]
fn simulator_output_contract_is_stable() {
    let temp_dir = tempfile::tempdir().expect("create temporary simulator output directory");
    let output_dir = temp_dir.path().join("out");
    let prefix = output_dir.join("sim_contract");
    let (ref_fasta, gtf) = write_simulator_contract_fixture(temp_dir.path());

    let summary = run(SimulateArgs {
        ref_fasta: ref_fasta.to_string_lossy().into_owned(),
        gtf: gtf.to_string_lossy().into_owned(),
        out_prefixes: vec![prefix.to_string_lossy().into_owned()],
        chrom: Some("chrTest".to_string()),
        circ_count: 10,
        switching_event: 0,
        circ_coverage: 8.0,
        linear_coverage: 0.01,
        scale: 0.5,
        read_len: 80,
        insert_len: 320,
        insert_sd: 60.0,
        insert_len_minor: 550,
        insert_sd_minor: 80.0,
        minor_insert_fraction: 0.10,
        error_rate: 0.002,
        exon_exclusive_rate: 0.25,
        seed: 5,
    })
    .expect("run simulator");

    assert_output_file_set(&output_dir);

    let isoforms = read_tsv(&prefix.with_extension("isoforms.tsv"), ISOFORM_HEADER);
    let reads = read_tsv(&prefix.with_extension("reads.tsv"), READS_HEADER);
    let annotation = fs::read_to_string(prefix.with_extension("annotation.gtf"))
        .expect("read filtered annotation");

    assert_eq!(isoforms.len(), 10);
    assert!(annotation.lines().all(|line| {
        let fields: Vec<&str> = line.split('\t').collect();
        fields.len() >= 3 && matches!(fields[2], "gene" | "transcript" | "exon")
    }));

    let circ_read_pairs_from_isoforms: usize = isoforms.iter().map(|row| parse_usize(row, 8)).sum();
    let bsj_read_pairs_from_isoforms: usize = isoforms.iter().map(|row| parse_usize(row, 9)).sum();
    let isoform_count: usize = isoforms.iter().map(|row| parse_usize(row, 10)).sum();
    let circ_read_pairs = reads.iter().filter(|row| row[1] != "NA").count();
    let linear_read_pairs = reads.iter().filter(|row| row[1] == "NA").count();
    let bsj_read_pairs = reads.iter().filter(|row| row[8] == "1").count();
    let outward_read_pairs = reads.iter().filter(|row| row[13] == "outward").count();
    let backward_read_pairs = reads.iter().filter(|row| row[13] == "backward").count();
    let bsj_reads = reads
        .iter()
        .map(|row| parse_usize(row, 10) + parse_usize(row, 12))
        .sum::<usize>();

    assert_eq!(circ_read_pairs_from_isoforms, circ_read_pairs);
    assert_eq!(bsj_read_pairs_from_isoforms, bsj_read_pairs);
    assert_eq!(linear_read_pairs, reads.len() - circ_read_pairs);
    for row in &reads {
        assert_eq!(
            parse_usize(row, 8),
            parse_usize(row, 10) | parse_usize(row, 12)
        );
        if row[13] == "bsj" {
            assert_eq!(row[8], "1");
        }
        if row[13] == "outward" {
            assert_eq!(row[7], "1");
            assert_eq!(row[8], "0");
            assert_eq!(row[10], "0");
            assert_eq!(row[12], "0");
        }
        if row[13] == "forward" {
            assert_eq!(row[8], "0");
        }
        if row[1] == "NA" {
            assert_eq!(row[3], "NA");
            assert_eq!(row[4], "NA");
            assert_eq!(row[6], "NA");
            assert_eq!(row[7], "0");
            assert_eq!(row[8], "0");
            assert_eq!(row[13], "forward");
        }
        if row[10] == "1" {
            assert!(row[9].contains("<bsj>"));
        }
        if row[12] == "1" {
            assert!(row[11].contains("<bsj>"));
        }
        assert!(!row[9].contains("chr"));
        assert!(!row[11].contains("chr"));
    }

    assert_eq!(
        fastq_record_count(&prefix.with_file_name("sim_contract_1.fq.gz")),
        reads.len()
    );
    assert_eq!(
        fastq_record_count(&prefix.with_file_name("sim_contract_2.fq.gz")),
        reads.len()
    );

    assert_eq!(summary.circ_count, isoforms.len());
    assert_eq!(summary.isoform_count, isoform_count);
    assert_eq!(summary.total_read_pairs, reads.len());
    assert_eq!(summary.circ_read_pairs, circ_read_pairs);
    assert_eq!(summary.linear_read_pairs, linear_read_pairs);
    assert_eq!(summary.bsj_reads, bsj_reads);
    assert_eq!(summary.bsj_read_pairs, bsj_read_pairs);
    assert!(
        outward_read_pairs > 0,
        "fixed simulator contract should contain pair-level outward truth"
    );
    assert_eq!(
        backward_read_pairs, 0,
        "simulator truth should not split BSJ-crossing reads into a separate backward class"
    );

    let expected_summary = format!(
        "Simulated {} circRNAs, {isoform_count} circular isoforms\n\
         Total {} read pairs, {circ_read_pairs} circRNA read pairs, {linear_read_pairs} linear read pairs\n\
         BSJ feature: {bsj_reads} reads / {bsj_read_pairs} read pairs",
        isoforms.len(),
        reads.len(),
    );
    assert_eq!(summary.to_string(), expected_summary);
}

#[test]
fn multi_sample_simulator_emits_isoform_switching_truth() {
    let temp_dir = tempfile::tempdir().expect("create temporary simulator output directory");
    let output_dir = temp_dir.path().join("out");
    let sample_a = output_dir.join("sampleA");
    let sample_b = output_dir.join("sampleB");
    let (ref_fasta, gtf) = write_simulator_contract_fixture(temp_dir.path());

    let summary = run(SimulateArgs {
        ref_fasta: ref_fasta.to_string_lossy().into_owned(),
        gtf: gtf.to_string_lossy().into_owned(),
        out_prefixes: vec![
            sample_a.to_string_lossy().into_owned(),
            sample_b.to_string_lossy().into_owned(),
        ],
        chrom: Some("chrTest".to_string()),
        circ_count: 8,
        switching_event: 3,
        circ_coverage: 80.0,
        linear_coverage: 0.0,
        scale: 0.0,
        read_len: 80,
        insert_len: 320,
        insert_sd: 60.0,
        insert_len_minor: 550,
        insert_sd_minor: 80.0,
        minor_insert_fraction: 0.10,
        error_rate: 0.002,
        exon_exclusive_rate: 0.25,
        seed: 7,
    })
    .expect("run multi-sample simulator");

    assert_eq!(summary.sample_count, 2);
    assert_eq!(summary.circ_count, 8);
    assert_eq!(
        read_tsv(&sample_a.with_extension("isoforms.tsv"), ISOFORM_HEADER).len(),
        read_tsv(&sample_b.with_extension("isoforms.tsv"), ISOFORM_HEADER).len()
    );
    assert!(sample_a.with_extension("usage.tsv").exists());
    assert!(sample_b.with_extension("usage.tsv").exists());

    let usage_a = read_tsv(&sample_a.with_extension("usage.tsv"), USAGE_HEADER);
    let usage_b = read_tsv(&sample_b.with_extension("usage.tsv"), USAGE_HEADER);
    assert_eq!(usage_a.len(), usage_b.len());
    let switching = inferred_switching_circs(&usage_a, &usage_b);
    assert_eq!(switching.len(), 3);
}

#[test]
fn multi_sample_switching_event_clamps_to_available_two_isoform_circs() {
    let temp_dir = tempfile::tempdir().expect("create temporary simulator output directory");
    let output_dir = temp_dir.path().join("out");
    let sample_a = output_dir.join("sampleA");
    let sample_b = output_dir.join("sampleB");
    let (ref_fasta, gtf) = write_simulator_contract_fixture(temp_dir.path());

    run(SimulateArgs {
        ref_fasta: ref_fasta.to_string_lossy().into_owned(),
        gtf: gtf.to_string_lossy().into_owned(),
        out_prefixes: vec![
            sample_a.to_string_lossy().into_owned(),
            sample_b.to_string_lossy().into_owned(),
        ],
        chrom: Some("chrTest".to_string()),
        circ_count: 3,
        switching_event: 100,
        circ_coverage: 40.0,
        linear_coverage: 0.0,
        scale: 0.0,
        read_len: 80,
        insert_len: 320,
        insert_sd: 60.0,
        insert_len_minor: 550,
        insert_sd_minor: 80.0,
        minor_insert_fraction: 0.10,
        error_rate: 0.002,
        exon_exclusive_rate: 0.25,
        seed: 17,
    })
    .expect("run multi-sample simulator with oversized switching target");

    let usage_a = read_tsv(&sample_a.with_extension("usage.tsv"), USAGE_HEADER);
    let usage_b = read_tsv(&sample_b.with_extension("usage.tsv"), USAGE_HEADER);
    let isoforms = read_tsv(&sample_a.with_extension("isoforms.tsv"), ISOFORM_HEADER);
    let available_two_isoform = isoforms
        .iter()
        .filter(|row| parse_usize(row, 10) >= 2)
        .count();
    let switching = inferred_switching_circs(&usage_a, &usage_b);
    assert_eq!(switching.len(), available_two_isoform);
    assert!(switching.len() <= 3);
}

fn inferred_switching_circs(sample_a: &[Vec<String>], sample_b: &[Vec<String>]) -> Vec<String> {
    use std::collections::BTreeMap;

    let mut by_sample: [BTreeMap<String, BTreeMap<String, f64>>; 2] =
        [BTreeMap::new(), BTreeMap::new()];
    for (target, rows) in by_sample.iter_mut().zip([sample_a, sample_b]) {
        for row in rows {
            target
                .entry(row[0].clone())
                .or_default()
                .insert(row[5].clone(), row[6].parse::<f64>().unwrap());
        }
    }
    by_sample[0]
        .iter()
        .filter_map(|(circ_id, usage_a)| {
            let usage_b = by_sample[1].get(circ_id)?;
            let major_a = usage_a
                .iter()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())?;
            let major_b = usage_b
                .iter()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())?;
            (major_a.0 != major_b.0).then(|| circ_id.clone())
        })
        .collect()
}

fn write_simulator_contract_fixture(dir: &Path) -> (PathBuf, PathBuf) {
    let ref_fasta = dir.join("contract.fa");
    let gtf = dir.join("contract.gtf");
    let reference = (0..20_000)
        .map(|idx| match idx % 4 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        })
        .collect::<String>();
    fs::write(&ref_fasta, format!(">chrTest\n{reference}\n")).expect("write test FASTA");

    let mut gtf_text = String::new();
    for tx_idx in 0..12 {
        let gene_id = format!("gene{tx_idx}");
        let transcript_id = format!("tx{tx_idx}");
        let strand = if tx_idx % 2 == 0 { '+' } else { '-' };
        let tx_start = 100 + tx_idx * 1_200;
        let exon_count = 5usize;
        let exon_len = 120usize;
        let intron_len = 30usize;
        let tx_end = tx_start + (exon_count - 1) * (exon_len + intron_len) + exon_len - 1;
        let attrs = format!("gene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";");
        gtf_text.push_str(&format!(
            "chrTest\tfixture\tgene\t{tx_start}\t{tx_end}\t.\t{strand}\t.\tgene_id \"{gene_id}\";\n"
        ));
        gtf_text.push_str(&format!(
            "chrTest\tfixture\ttranscript\t{tx_start}\t{tx_end}\t.\t{strand}\t.\t{attrs}\n"
        ));
        for exon_idx in 0..exon_count {
            let start = tx_start + exon_idx * (exon_len + intron_len);
            let end = start + exon_len - 1;
            gtf_text.push_str(&format!(
                "chrTest\tfixture\texon\t{start}\t{end}\t.\t{strand}\t.\t{attrs} exon_number \"{}\";\n",
                exon_idx + 1
            ));
        }
    }
    fs::write(&gtf, gtf_text).expect("write test GTF");
    (ref_fasta, gtf)
}

fn assert_output_file_set(output_dir: &Path) {
    let mut names: Vec<String> = fs::read_dir(output_dir)
        .expect("read simulator output directory")
        .map(|entry| {
            entry
                .expect("read simulator output entry")
                .file_name()
                .to_string_lossy()
                .into_owned()
        })
        .collect();
    names.sort();

    let expected = vec![
        "sim_contract.annotation.gtf",
        "sim_contract.isoforms.tsv",
        "sim_contract.reads.tsv",
        "sim_contract_1.fq.gz",
        "sim_contract_2.fq.gz",
    ];
    assert_eq!(names, expected);
}

fn read_tsv(path: &Path, expected_header: &str) -> Vec<Vec<String>> {
    let content =
        fs::read_to_string(path).unwrap_or_else(|err| panic!("read TSV {}: {err}", path.display()));
    let mut lines = content.lines();
    assert_eq!(
        lines.next().unwrap_or_default(),
        expected_header,
        "unexpected header in {}",
        path.display()
    );
    lines
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect()
}

fn parse_usize(row: &[String], idx: usize) -> usize {
    row[idx]
        .parse()
        .unwrap_or_else(|err| panic!("parse usize from column {idx} value {:?}: {err}", row[idx]))
}

fn fastq_record_count(path: &Path) -> usize {
    let content = read_fastq_text(path);
    let line_count = content.lines().count();
    assert_eq!(line_count % 4, 0, "FASTQ line count must be divisible by 4");
    line_count / 4
}

fn read_fastq_text(path: &Path) -> String {
    if path.extension().and_then(|ext| ext.to_str()) == Some("gz") {
        let output = Command::new("gzip")
            .arg("-dc")
            .arg(path)
            .output()
            .unwrap_or_else(|err| panic!("run gzip -dc {}: {err}", path.display()));
        assert!(
            output.status.success(),
            "gzip -dc {} exited with {}",
            path.display(),
            output.status
        );
        String::from_utf8(output.stdout)
            .unwrap_or_else(|err| panic!("decode FASTQ {} as UTF-8: {err}", path.display()))
    } else {
        fs::read_to_string(path)
            .unwrap_or_else(|err| panic!("read FASTQ {}: {err}", path.display()))
    }
}
