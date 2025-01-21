#!/bin/bash

wd1=$HOME/Edyta/dropoutScreen/Nov18_2020/Greenblatt_CRISPR_12
wd2=$HOME/Edyta/dropoutScreen/Nov27_2020/Greenblatt
wd3=$HOME/Edyta/dropoutScreen/Mar18_2021/GREENBLATT
wd4=$HOME/Edyta/dropoutScreen/May28_2021/GREENBLATT
wd5=$HOME/Edyta/dropoutScreen/Jun24_2021/GREENBLATT
wd6=$HOME/Edyta/dropoutScreen/Jul05_2021/GREENBLATT


for wd in $wd6; do
	cd $wd
	rm *_unmatched.fastq.gz
	rm *_preProcessed.fastq.gz
	rm *_unaligned.fastq
	rm *_aligned.tab
done
