#!/bin/bash

cmd=$HOME/Edyta/dropoutScreen/scripts/collect_preprocess_and_alignment_stats.sh
wd1=$HOME/Edyta/dropoutScreen/Nov18_2020/Greenblatt_CRISPR_12
wd2=$HOME/Edyta/dropoutScreen/Nov27_2020/Greenblatt
wd3=$HOME/Edyta/dropoutScreen/Mar18_2021/GREENBLATT
wd4=$HOME/Edyta/dropoutScreen/May28_2021/GREENBLATT
wd5=$HOME/Edyta/dropoutScreen/Jun24_2021/GREENBLATT
wd6=$HOME/Edyta/dropoutScreen/Jul05_2021/GREENBLATT


for wd in $wd6; do
	cd $wd
	if [ -f preprocess_and_alignment_stats.txt ]; then
		rm preprocess_and_alignment_stats.txt
	fi
	$cmd
done
