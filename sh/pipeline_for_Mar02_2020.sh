#!/bin/bash
## used for preprocessing reads with stagger sequence at both ends

f=$1  # input is in the format of *_001.fastq.gz
wd=$2
cd $wd
TKOdir=$HOME/Zuyao/crispr/TKO
bf=$(basename $f _R1_001.fastq.gz)

if [ 1 -eq 1 ]; then
echo "1.0 preprocessing $bf"
zcat $f | perl $HOME/perl_script/preprocessReads_trimBothends.pl - 2> >(gzip > ${bf}_unmatched.fastq.gz) | gzip - > ${bf}_preProcessed.fastq.gz
fi

if [ 1 -eq 1 ]; then
echo "2.0 aligning to TKOv2"
trim5=0
trim3=0
logp=${bf}_bowtie.log
zcat ${bf}_preProcessed.fastq.gz | bowtie -v2 -m1 -p4 -5 "$trim5" -3 "$trim3" --un ${bf}_unaligned.fastq --sam-nohead ${TKOdir}/TKOv2.1 - 2> "$logp" > ${bf}_aligned.tab
#zcat ${f} | bowtie -v2 -m1 -p4 -5 "$trim5" -3 "$trim3" --un ${bf}_unaligned.fastq --sam-nohead ${TKOdir}/TKOv2.1 - 2> "$logp" > ${bf}_aligned.tab

cut -f3 ${bf}_aligned.tab | sort | uniq -c > ${bf}_readCount.tab

aligned_read_count=$(wc -l ${bf}_aligned.tab | cut -d " " -f1)
echo "aligned read count ${aligned_read_count}"
fi

if [ 1 -eq 1 ]; then
echo "3.0 counting reads per guide"
#perl $HOME/perl_script/countReads_perGuide.pl ${bf}_readCount.tab ${bf}_readCount_bf.tab T > ${bf}_geneCount.tab  #normalized to per million
perl $HOME/perl_script/countReads_perGuide.pl ${bf}_readCount.tab ${bf}_guideRawCount.tab F > ${bf}_geneRawCount.tab
fi

echo "Finish all"
