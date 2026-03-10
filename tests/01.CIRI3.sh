bwa mem -t 16 -T 19 chr1.fa test_1.fq.gz test_2.fq.gz > test.sam
java -jar ../vendor/CIRI3/CIRI3_Java_1.8.0.jar -T 16 -I test.sam -O CIRI3_result.txt -F chr1.fa -A chr1.gtf -S 0
../target/release/ciri-toolkit -i ./test.bam -o ./test.ciri -r ./chr1.fa -a ./chr1.gtf -s 0 -t 16

