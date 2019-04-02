#!/bin/bash
set -e -x -o pipefail
ncpus=`nproc`

infastq=$1
outfastq=`echo $infastq|sed "s/.fastq.gz/.trim.fastq.gz/g"`

cutadapt \
--nextseq-trim=2 \
--trim-n \
-n 5 -O 5 \
-q 10,10 \
-m 35 \
-b file:/opt/TruSeq_and_nextera_adapters.consolidated.fa \
-j $ncpus \
-o $outfastq \
$infastq
