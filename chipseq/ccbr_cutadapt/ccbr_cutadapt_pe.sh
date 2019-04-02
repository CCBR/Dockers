#!/bin/bash
set -e -x -o pipefail
ncpus=`nproc`

infastq1=$1
infastq2=$2
outfastq1=`echo $infastq1|sed "s/.fastq.gz/.trim.fastq.gz/g"`
outfastq2=`echo $infastq2|sed "s/.fastq.gz/.trim.fastq.gz/g"`


cutadapt \
--pair-filter=any \
--nextseq-trim=2 \
--trim-n \
-n 5 -O 5 \
-q 10,10 \
-m 35 \
-b file:/opt/TruSeq_and_nextera_adapters.consolidated.fa \
-B file:/opt/TruSeq_and_nextera_adapters.consolidated.fa \
-j $ncpus \
-o $outfastq1 \
-p $outfastq2 \
$infastq1 \
$infastq2
