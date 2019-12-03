#!/bin/bash
set -e -x -o pipefail
ARGPARSE_DESCRIPTION="call atac-seq peaks using genrich v0.6"      # this is optional
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--bamrep1',required=True, help='input name sorted bam file for replicate 1')
parser.add_argument('--bamrep2',required=True, help='input name sorted bam file for replicate 2')
parser.add_argument('--peakfile1',required=True, help='output narrowPeak file for replicate 1')
parser.add_argument('--peakfile2',required=True, help='output narrowPeak file for replicate 1')
parser.add_argument('--mergedpeakfile',required=True, help='output narrowPeak file for both replicates combined')
# only required if filtering
parser.add_argument('--filterpeaks',required=False, default="False", help='filterpeaks by qvalue: True or False')
parser.add_argument('--qfilter',required=False, default=2, help='default qfiltering value is 2 for q=0.01')
EOF

excludelist=$(samtools view -H $BAMREP1|grep ^@SQ|cut -f2|sed "s/SN://g"|awk "\$1 !~ /[0-9]$/"|grep -v "chrX"|grep -v "chrY"|tr "\n" ","|awk "{print substr(\$_,1,length(\$_)-1)}")
# echo $excludelist
Genrich  -t $BAMREP1 -o $PEAKFILE1  -j -r -e $excludelist -v -s 20 -m 6 -x -y -p 0.01 -q 1
Genrich  -t $BAMREP2 -o $PEAKFILE2  -j -r -e $excludelist -v -s 20 -m 6 -x -y -p 0.01 -q 1
Genrich  -t $BAMREP1,$BAMREP2 -o $MERGEDPEAKFILE  -j -r -e $excludelist -v -s 20 -m 6 -x -y -p 0.01 -q 1

if [ $FILTERPEAKS == "True" ];then
  qvalue=$QFILTER
  for f in $PEAKFILE1 $PEAKFILE2 $MERGEDPEAKFILE;do
    filteredpeakfile=$(echo $f|sed "s/.narrowPeak/.qfilter.narrowPeak/g")
    awk -F"\t" -v q=$qvalue "{if (\$9>q){print}}" $f > $filteredpeakfile
  done
fi
