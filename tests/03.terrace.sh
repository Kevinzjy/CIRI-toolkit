STAR --genomeDir chr1/STAR --runThreadN 16 --readFilesIn chr1/test_1.fq.gz chr1/test_2.fq.gz --readFilesCommand zcat --outFileNamePrefix chr1/test_star_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within

time terrace -i chr1/test_star_Aligned.sortedByCoord.out.bam -o chr1/test_star_terrace.gtf -fa chr1/chr1.fa --read_length 150 -r chr1/chr1.gtf -fe chr1/test_star_terrace.csv
