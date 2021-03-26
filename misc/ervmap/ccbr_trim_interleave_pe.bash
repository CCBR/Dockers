#!/bin/bash

. /opt2/conda/etc/profile.d/conda.sh
conda activate python3

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="Trim PE reads using cutadapt"      # this is optional
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--infastq1',required=True, help='input R1 fastq.gz file')
parser.add_argument('--infastq2',required=True, help='input R2 fastq.gz file')
parser.add_argument('--threads',required=True, help='number of threads')
parser.add_argument('--outinterleavedfastq',required=True, help='output fastq.gz file')
EOF

infastq1=$INFASTQ1
infastq2=$INFASTQ2
outfastq=$OUTINTERLEAVEDFASTQ
ncpus=$THREADS

cutadapt \
--pair-filter=any \
--nextseq-trim=2 \
--trim-n \
-n 5 -O 5 \
-q 10,10 \
-m 35:35 \
-b file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa \
-B file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa \
-j $ncpus \
-o ${infastq1%.*}.trimmed.fastq.gz \
-p ${infastq2%.*}.trimmed.fastq.gz \
$infastq1 $infastq2

reformat.sh -Xmx20g in1=${infastq1%.*}.trimmed.fastq.gz in2=${infastq2%.*}.trimmed.fastq.gz out=$outfastq

rm -f ${infastq1%.*}.trimmed.fastq.gz ${infastq2%.*}.trimmed.fastq.gz

conda deactivate
