# 1. Run ciri-simulator
#cargo run --bin ciri-simulator --release -- -r ./chr1/chr1.fa -a ./chr1/chr1.gtf -o ./chr1/test --circ-coverage 100 --linear-coverage 10 --circ-count 20000

# 2. Run CIRI3
#bwa mem -t 16 -T 19 ./chr1/chr1.fa ./chr1/test_1.fq.gz ./chr1/test_2.fq.gz > ./chr1/test.sam
#time java -jar ../vendor/CIRI3/CIRI3_Java_1.8.0.jar -T 16 -I ./chr1/test.sam -O ./chr1/CIRI3_result.txt -F ./chr1/chr1.fa -A ./chr1/test.annotation.gtf -S 0

# Run CIRI-toolkit
#samtools view -bS -@ 16 -o ./chr1/test.bam ./chr1/test.sam

# Baseline sampling (optional):
# CIRI_PROFILE_SCAN1=1 /usr/bin/time -f "elapsed=%E cpu=%P mem=%MKB" \
#   ../target/release/ciri -i ./chr1/test.bam -o ./chr1/CIRI3_result.txt -r ./chr1/chr1.fa -a ./chr1/test.annotation.gtf -t 16

time cargo run --release -- -i ./chr1/test.bam -o ./chr1/simulate.ciri \
    -r ./chr1/chr1.fa -a ./chr1/test.annotation.gtf -t 16 -s 0

# CIRI3 BSJ-level validation
python ../scripts/ciri_result_diff.py ./chr1/CIRI3_result.txt ./chr1/simulate.ciri.out

# Read-level evaluation
#python ../scripts/ciri_segments_eval.py chr1/test.reads.tsv chr1/simulate.ciri.segments
