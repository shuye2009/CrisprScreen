#!/bin/bash

cmd=$HOME/Zuyao/crispr/scripts/pipeline_for_Mar02_2020.sh
#wd=$HOME/Zuyao/crispr/Mar02_2020/Greenblatt
#wd=$HOME/Zuyao/crispr/Nov18_2020/Greenblatt_CRISPR_24
wd=$HOME/Zuyao/crispr/Mar25_2021/GREENBLATT
cd $wd
for f in *_R1_001.fastq.gz; do
	submitjob -m 10 $cmd $f $wd
done
