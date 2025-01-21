#!/bin/bash

for f in *_R1_001.fastq.gz; do
	basef=$(basename "$f" "_R1_001.fastq.gz")
	echo -e "\n$f" >> preprocess_and_alignment_stats.txt
	zcat ${basef}_unmatched.fastq.gz | tail -n 20 >> preprocess_and_alignment_stats.txt
	cat ${basef}_bowtie.log >> preprocess_and_alignment_stats.txt
done
