#!/bin/bash

cpus=$1
treatmentbamfile=$2
inputbamfile=$3

treatmentbam=`echo $treatmentbamfile|awk -F "/" '{print $NF}'`
inputbam=`echo $inputbamfile|awk -F "/" '{print $NF}'`

t_label=${treatmentbam%%.*}
i_label=${inputbam%%.*}
outpdf=${treatmentbam%.*}_vs_${inputbam%.*}.fingerprint.pdf

treatmentsortedbam=`echo $treatmentbam|sed "s/.bam/.sorted.bam/g"`
inputsortedbam=`echo $inputbam|sed "s/.bam/.sorted.bam/g"`

samtools sort -@ $cpus -o $treatmentsortedbam $treatmentbamfile
samtools sort -@ $cpus -o $inputsortedbam $inputbamfile

samtools index $treatmentsortedbam
samtools index $inputsortedbam

plotFingerprint -b $treatmentsortedbam $inputsortedbam --labels $t_label $i_label -p $cpus --skipZeros --plotFile $outpdf
