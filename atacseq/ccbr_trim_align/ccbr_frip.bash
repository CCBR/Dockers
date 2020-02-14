#!/bin/bash

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="fraction of reads in peaks/DHS/Enhancers/Promoters" 
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--narrowpeak',required=True, help='input narrowPeak file')
parser.add_argument('--tagalign',required=True, help='input tagAlign.gz file')
parser.add_argument('--samplename',required=True, help='sample name')
parser.add_argument('--genome',required=True, help='hg19/hg38/mm9/mm10')
parser.add_argument('--out',required=True, help='output txt file')
parser.add_argument('--scriptsfolder',required=False, default='/opt2', help='folder where the scripts are... used for debuging without rebuilding the docker')
EOF

set -e -x -o pipefail
nreads=$(zcat $TAGALIGN|wc -l)

frip_reads=$(bedtools intersect -wa -a $TAGALIGN -b $NARROWPEAK|wc -l)
frip=$(echo "scale=3;$frip_reads/$nreads"|bc)

fridhs_reads=$(bedtools intersect -wa -a $TAGALIGN -b /opt2/db/${GENOME}.DHS.bed.gz|wc -l)
fridhs=$(echo "scale=3;$fridhs_reads/$nreads"|bc)

fripro_reads=$(bedtools intersect -wa -a $TAGALIGN -b /opt2/db/${GENOME}.promoters.bed.gz|wc -l)
fripro=$(echo "scale=3;$fripro_reads/$nreads"|bc)

frienh_reads=$(bedtools intersect -wa -a $TAGALIGN -b /opt2/db/${GENOME}.enhancers.bed.gz|wc -l)
frienh=$(echo "scale=3;$frienh_reads/$nreads"|bc)

echo -ne "# $SAMPLENAME\tFRiP\t$frip\n" > $OUT
echo -ne "# $SAMPLENAME\tFRiDHS\t$fridhs\n" >> $OUT
echo -ne "# $SAMPLENAME\tFRiPro\t$fripro\n" >> $OUT
echo -ne "# $SAMPLENAME\tFRiEnh\t$frienh\n" >> $OUT
