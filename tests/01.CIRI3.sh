# Run CIRI3
# bwa mem -t 16 -T 19 ./chr1/chr1.fa ./chr1/test_1.fq.gz ./chr1/test_2.fq.gz > ./chr1/test.sam
# time java -jar ../vendor/CIRI3/CIRI3_Java_1.8.0.jar -T 16 -I ./chr1/test.sam -O ./chr1/CIRI3_result.txt -F ./chr1/chr1.fa -A ./chr1/chr1.gtf

# Run CIRI-toolkit
# samtools view -bS -@ 16 -o ./chr1/test.bam ./chr1/test.sam

# Baseline sampling (optional):
# CIRI_PROFILE_SCAN1=1 /usr/bin/time -f "elapsed=%E cpu=%P mem=%MKB" \
#   ../target/release/ciri-toolkit -i ./chr1/test.bam -o ./chr1/CIRI3_result.txt -r ./chr1/chr1.fa -a ./chr1/chr1.gtf -t 16

time cargo run --release -- -i ./chr1/test.sam -o ./chr1/CIRI-rs.ciri -r ./chr1/chr1.fa -a ./chr1/chr1.gtf -t 16

# Keep for validation
python analyze_diff.py ./chr1/CIRI3_result.txt ./chr1/CIRI-rs.ciri.result
