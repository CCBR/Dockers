#!/bin/bash
cpus=$1
infq=$2
reftargz=$3

tar xzvf $reftargz
ref=`echo $reftargz|awk -F "/" '{print $NF}'|sed "s/.tar.gz//g"`

basefilename=`echo $infq|awk -F "/" '{print $NF}'|sed "s/.fastq.gz//g"`
bamfile=${basefilename}.bam
sortedbamfile=${basefilename}.sorted.bam
sortedbambaifile=${sortedbamfile}.bai
sortedq5bamfile=${basefilename}.sorted.Q5.bam
sortedq5bambaifile=${basefilename}.sorted.Q5.bam.bai
sortedbamflagstatfile=${sortedbamfile}.flagstat
sortedq5bamflagstatfile=${sortedq5bamfile}.flagstat

bwa mem -t $cpus $ref $infq > $bamfile
samtools sort -@ $cpus -o $sortedbamfile $bamfile
samtools index $sortedbamfile
samtools view -@ $cpus -b -q 6 $sortedbamfile -o $sortedq5bamfile
samtools index $sortedq5bamfile
samtools flagstat $sortedbamfile > $sortedbamflagstatfile
samtools flagstat $sortedq5bamfile > $sortedq5bamflagstatfile
