#!/bin/bash
bamfilename=$1
bam=`echo $bamfilename|awk -F "/" '{print $NF}'`

tagalignfile=`echo $bam|sed "s/.bam/.DD.tagAlign/g"`
outbamfile=`echo $bam|sed "s/.bam/.DD.bam/g"`
outbamflagstat=${outbamfile}.flagstat

macs2 filterdup -i $bam -o TmpTagAlign
awk -F"\t" -v OFS="\t" '{if ($2>0 && $3>0) {print}}' TmpTagAlign > TmpTagAlign2
samtools view -H $bam|grep "^@SQ"|cut -f2,3|sed "s/SN://"|sed "s/LN://g" > GenomeFile
awk -F"\t" -v OFS="\t" '{print $1,1,$2}' GenomeFile |sort -k1,1 -k2,2n > GenomeFileBed
bedtools intersect -wa -f 1.0 -a TmpTagAlign2 -b GenomeFileBed > $tagalignfile
bedtools bedtobam -i $tagalignfile -g GenomeFile > $outbamfile
samtools flagstat $outbamfile > $outbamflagstat
