#!/bin/bash

cpus=$1
reftargz=$2
t_bamfile=$3
i_bamfile=$4
t_ppqt=$5

genome=`echo $reftargz|awk -F "/" '{print $NF}'|sed "s/.tar.gz//g"`
t_bam=`echo $t_bamfile|awk -F "/" '{print $NF}'`
i_bam=`echo $i_bamfile|awk -F "/" '{print $NF}'`

if [ "$genome" == "mm10" ]; then
genomesize="2652783500"
elif [ "$genome" == "hg19" ]; then
genomesize="2864785220"
elif [ "$genome" == "hg38" ]; then
genomesize="2913022398"
elif [ "$genome" == "mm9" ]; then
genomesize="2620345972"
fi

treatmentbw=`echo $t_bam|sed "s/.bam/.bw/g"`
inputbw=`echo $i_bam|sed "s/.bam/.bw/g"`
treatmentextsize=`cat $t_ppqt|awk -F"\t" '{print $3}'|awk -F"," '{print $1}'`

treatmentsortedbam=`echo $t_bam|sed "s/.bam/.sorted.bam/g"`
inputsortedbam=`echo $i_bam|sed "s/.bam/.sorted.bam/g"`

samtools sort -@ $cpus -o $treatmentsortedbam $t_bamfile
samtools sort -@ $cpus -o $inputsortedbam $i_bamfile

samtools index $treatmentsortedbam
samtools index $inputsortedbam

bamCoverage --bam $treatmentsortedbam -o $treatmentbw --binSize 25 --smoothLength 75 --ignoreForNormalization chrX chrY chrM --numberOfProcessors $cpus --normalizeUsing RPGC --effectiveGenomeSize $genomesize --extendReads $treatmentextsize
bamCoverage --bam $inputsortedbam -o $inputbw --binSize 25 --smoothLength 75 --ignoreForNormalization chrX chrY chrM --numberOfProcessors $cpus --normalizeUsing RPGC --effectiveGenomeSize $genomesize --extendReads 200
