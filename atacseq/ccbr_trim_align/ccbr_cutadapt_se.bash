#!/bin/bash

. /opt2/conda/etc/profile.d/conda.sh
conda activate python3

set -e -x -o pipefail
ncpus=`nproc`
ARGPARSE_DESCRIPTION="Trim SE reads using cutadapt"      # this is optional
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--infastq',required=True, help='input fastq.gz file')
EOF

infastq=$INFASTQ
outfastq=`echo $infastq|sed "s/.fastq.gz/.trim.fastq.gz/g"`
samplename=`echo $infastq|sed "s/.fastq.gz//g"|sed "s/.R1//g"`

cutadapt \
--nextseq-trim=2 \
--trim-n \
-n 5 -O 5 \
-q 10,10 \
-m 35 \
-b file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa \
-j $ncpus \
-o $outfastq \
$infastq

nreads=`zcat $infastq|wc -l`
nreadstrim=`zcat $outfastq|wc -l`
echo "$nreads $nreadstrim"|awk '{printf("%d\tInput Nreads\n%d\tAfter trimming\n",$1/4,$2/4)}' > ${samplename}.nreads.txt

conda deactivate
