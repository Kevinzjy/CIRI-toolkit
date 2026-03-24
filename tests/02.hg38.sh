#python extract_diff_reads.py /home/zhangjy/data/01.ZhaoLab/02.CIRI-rs/CIRI3_result.txt /home/zhangjy/data/01.ZhaoLab/02.CIRI-rs/CIRI-rs.ciri.result -o ./hg38/diff -i /home/zhangjy/data/01.ZhaoLab/02.CIRI-rs/RNA015434_S1.bam
#python extract_diff_reads.py /home/zhangjy/data/01.ZhaoLab/02.CIRI-rs/CIRI3_result.txt /home/zhangjy/data/01.ZhaoLab/02.CIRI-rs/CIRI-rs.ciri.result --mode diff --no-circ-context -o hg38/diff_only

java -jar ../vendor/CIRI3/CIRI3_Java_1.8.0.jar -T 16 -I ./hg38/diff.subset.bam -O ./hg38/diff.subset.CIRI3_result.txt -F /data/public/database/gencode/hg38/_BWAindex/hg38.fa -A /data/public/database/gencode/hg38/gencode.v44.annotation.gtf -S 0

../target/release/ciri-toolkit -i hg38/diff.subset.bam -o ./hg38/diff.subset.CIRI-rs.ciri -r /data/public/database/gencode/hg38/_BWAindex/hg38.fa -a /data/public/database/gencode/hg38/gencode.v44.annotation.gtf -s 0 -t 16

python analyze_diff.py ./hg38/diff.subset.CIRI3_result.txt ./hg38/diff.subset.CIRI-rs.ciri.result 
