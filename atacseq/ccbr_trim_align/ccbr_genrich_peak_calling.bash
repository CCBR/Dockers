#!/bin/bash
set -e -x -o pipefail
ARGPARSE_DESCRIPTION="call atac-seq peaks using genrich v0.6"      # this is optional
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--bam',required=True, help='input name sorted bam file')
parser.add_argument('--peakfile',required=True, help='output narrowPeak file')
EOF

excludelist=$(samtools view -H $BAM|grep ^@SQ|cut -f2|sed "s/SN://g"|awk "\$1 !~ /[0-9]$/"|grep -v "chrX"|tr "\n" ","|awk "{print substr(\$_,1,length(\$_)-1)}")
echo $excludelist
Genrich  -t $BAM -o $PEAKFILE  -j -r -e $excludelist -v -s 20 -m 6 -x -y
