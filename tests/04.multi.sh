# Run CIRI-simulator
cargo run --bin ciri-simulator --release -- -r ./chr1/chr1.fa -a ./chr1/chr1.gtf -o ./multi/sample1 ./multi/sample2 --circ-coverage 100 --linear-coverage 10 --circ-count 10000 --switching-event 2000

# Run bwa-mem
bwa mem -t 16 -T 19 ./chr1/chr1.fa ./multi/sample1_1.fq.gz ./multi/sample1_2.fq.gz | samtools view -bS -@ 16 -o ./multi/sample1.bam -
bwa mem -t 16 -T 19 ./chr1/chr1.fa ./multi/sample2_1.fq.gz ./multi/sample2_2.fq.gz | samtools view -bS -@ 16 -o ./multi/sample2.bam -

# 1st pass BSJ detection
mkdir -p ./multi/1st_pass
cargo run --release -- -i ./multi/sample1.bam -o ./multi/1st_pass/sample1 -r ./chr1/chr1.fa -a ./chr1/chr1.gtf -s 0 --1st-pass
cargo run --release -- -i ./multi/sample2.bam -o ./multi/1st_pass/sample2 -r ./chr1/chr1.fa -a ./chr1/chr1.gtf -s 0 --1st-pass

# Merge BSJ detection results
cargo run --release --bin ciri-merge -- -i ./multi/1st_pass/sample1.out ./multi/1st_pass/sample2.out -o ./multi/1st_pass.bed

# 2nd pass segment detection
mkdir -p ./multi/2nd_pass
cargo run --release -- -i ./multi/sample1.bam -o ./multi/2nd_pass/sample1 -r ./chr1/chr1.fa -a ./chr1/chr1.gtf -s 0 --2nd-pass --circ ./multi/1st_pass.bed
cargo run --release -- -i ./multi/sample2.bam -o ./multi/2nd_pass/sample2 -r ./chr1/chr1.fa -a ./chr1/chr1.gtf -s 0 --2nd-pass --circ ./multi/1st_pass.bed

# Integrative sequence assembly
printf "sample1\t./multi/2nd_pass/sample1\nsample2\t./multi/2nd_pass/sample2\n" > ./multi/sample_list.tsv
cargo run --release --bin ciri-assemble -- -i ./multi/sample_list.tsv -o ./multi/merged -r ./chr1/chr1.fa -a ./chr1/chr1.gtf 
