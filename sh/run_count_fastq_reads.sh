#!/bin/bash

cmd=$HOME/Edyta/dropoutScreen/scripts/count_fastq_reads.sh
wd1=$HOME/Edyta/dropoutScreen/Nov18_2020/Greenblatt_CRISPR_12
wd2=$HOME/Edyta/dropoutScreen/Nov27_2020/Greenblatt
wd3=$HOME/Edyta/dropoutScreen/Mar18_2021/GREENBLATT
wd4=$HOME/Edyta/dropoutScreen/May28_2021/GREENBLATT
wd5=$HOME/Edyta/dropoutScreen/Jun24_2021/GREENBLATT


for wd in $wd1 $wd2 $wd3 $wd4 $wd5; do
	submitjob -m 10 $cmd $wd
done
