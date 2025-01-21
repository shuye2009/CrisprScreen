#!/bin/bash
wd=$1
cd $wd
if [ -f fastq_raw_reads_count.tab ]; then
	rm fastq_raw_reads_count.tab
fi

for f in *R1_001.fastq.gz; do 
	count=$(zcat $f |grep "@" | wc -l); echo -e "$f\t$count" >> fastq_raw_reads_count.tab; 
done
