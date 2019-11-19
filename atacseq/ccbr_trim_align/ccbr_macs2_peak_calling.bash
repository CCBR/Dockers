#!/bin/bash
. /opt/conda/etc/profile.d/conda.sh
conda activate python3

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="call atac-seq peaks using macs2"      # this is optional
source /opt/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--tagalign',required=True, help='input tagAlign.gz file')
parser.add_argument('--extsize',required=False, default=73, help='extsize')
parser.add_argument('--shiftsize',required=False, default=37, help='shiftsize')
parser.add_argument('--samplename',required=True, help='samplename')
parser.add_argument('--genomename',required=True, help='hg19/hg38/mm9/mm10')
# only required if filtering
parser.add_argument('--filterpeaks',required=False, default="False", help='filterpeaks by qvalue: True or False')
parser.add_argument('--qfilter',required=False, default=2, help='default qfiltering value is 2 for q=0.01')
# only required if saving bigwigs
parser.add_argument('--dedupbam',required=False, help='dedupbam')
parser.add_argument('--savebigwig',required=False, default="False", help='save bigwig file: True or False')

EOF

extsize=$EXTSIZE
shiftsize=$SHIFTSIZE
tagAlign=$TAGALIGN
samplename=$SAMPLENAME
prefix=${samplename}.macs2
genome=$GENOMENAME
dedupbam=$DEDUPBAM

if [ $genome == "hg19" ]; then g="hs";fi
if [ $genome == "hg38" ]; then g="hs";fi
if [ $genome == "mm10" ]; then g="mm";fi
if [ $genome == "mm9" ]; then g="mm";fi

if [ $SAVEBIGWIG == "True" ]; then
  macs2 callpeak -t $tagAlign -f BED -n $prefix -g $g -p 0.01 --shift -$shiftsize --extsize $extsize --keep-dup all --call-summits --nomodel --SPMR -B
else
  macs2 callpeak -t $tagAlign -f BED -n $prefix -g $g -p 0.01 --shift -$shiftsize --extsize $extsize --keep-dup all --call-summits --nomodel
fi

# remove duplicate peak calls
mv ${prefix}_peaks.narrowPeak ${prefix}_peaks.narrowPeak.tmp
sort -k9,9gr ${prefix}_peaks.narrowPeak.tmp|awk -F"\t" '!NF || !seen[$1":"$2"-"$3]++'|sort -k1,1 -k2,2n > ${prefix}_peaks.narrowPeak

# qvalue filter of 0.05
# awk -F"\t" '{if ($9>1.302){print}}' ${prefix}_peaks.narrowPeak > ${prefix}.qfilter.narrowPeak
# qvalue filter of 0.01
if [ $FILTERPEAKS == "True" ];then
  qvalue=$QFILTER
  awk -F"\t" -v q=$qvalue '{if ($9>q){print}}' ${prefix}_peaks.narrowPeak > ${prefix}.qfilter.narrowPeak
fi

if [ $SAVEBIGWIG == "True" ];then
  if [ ! -f $dedupbam ];then
    exit "Dedupbam file:\"$dedupbam\" not accessible"
  fi
  samtools view -H $dedupbam|grep "^@SQ"|cut -f2,3|sed "s/SN://g"|sed "s/LN://g" > ${samplename}.genome
  bedSort ${prefix}_treat_pileup.bdg ${prefix}_treat_pileup.bdg
  bedGraphToBigWig ${prefix}_treat_pileup.bdg ${samplename}.genome ${prefix}.bw
  rm -f ${samplename}.genome
fi

rm -f ${prefix}_treat_pileup.bdg
rm -f ${prefix}_control_lambda.bdg

conda deactivate
