#!/bin/bash
. /opt2/conda/etc/profile.d/conda.sh
conda activate python3

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="call atac-seq peaks using macs2"      # this is optional
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--tagalign1',required=True, help='input tagAlign.gz file for replicate 1')
parser.add_argument('--tagalign2',required=False, help='input tagAlign.gz file for replicate 2')
parser.add_argument('--tagalign3',required=False, help='input tagAlign.gz file for replicate 3')
parser.add_argument('--tagalign4',required=False, help='input tagAlign.gz file for replicate 4')

parser.add_argument('--rep1name',required=True, help='Name for replicate 1')
parser.add_argument('--rep2name',required=False, help='Name for replicate 2')
parser.add_argument('--rep3name',required=False, help='Name for replicate 3')
parser.add_argument('--rep4name',required=False, help='Name for replicate 4')

parser.add_argument('--genome',required=True,help="hg19/38 or mm9/10")
parser.add_argument('--pooledpeakfile',required=False, help='output narrowPeak file for both replicates combined')
parser.add_argument('--concensusbedfile',required=False, help='consensus overlapping MAX majority peaks in bed format')

parser.add_argument('--extsize',required=False, default=200, help='extsize')
parser.add_argument('--shiftsize',required=False, default=100, help='shiftsize')
parser.add_argument('--samplename',required=True, help='samplename,i.e., all replicates belong to this sample')

# only required if filtering
parser.add_argument('--filterpeaks',required=False, default="True", help='filterpeaks by qvalue: True or False')
parser.add_argument('--qfilter',required=False, default=0.693147, help='default qfiltering value is 0.693147 (-log10 of 0.5) for q=0.5')

# only required if saving bigwigs/tn5knicks.bed/bam etc.
parser.add_argument('--genomefile',required=True, help='dedupbam based genome file ... required by bedGraphToBigWig')
# parser.add_argument('--savebigwig',required=False, default="False", help='save bigwig file: True or False')
# parser.add_argument('--savetn5knicksbed',required=False, default="False", help='save tn5knicks bed and bam file: True or False')

parser.add_argument('--scriptsfolder',required=False, default='/opt2', help='folder where the scripts are... used for debuging without rebuilding the docker')

EOF


nreplicates=1
if [ $TAGALIGN2 ];then
	nreplicates=2
	if [ $TAGALIGN3 ];then
		nreplicates=3
		if [ $TAGALIGN4 ];then
			nreplicates=4
		fi
	fi
fi

if [ $REP2NAME ];then
	nreplicates=2
	if [ $REP3NAME ];then
		nreplicates=3
		if [ $REP4NAME ];then
			nreplicates=4
		fi
	fi
fi

if [ "$nreplicates" -eq "2" ];then
	if [ ! $REP2NAME ]; then
		echo "Name for replicate 2 is required!"
		exit
	fi
fi
if [ "$nreplicates" -eq "3" ];then
	if [ ! $REP3NAME ]; then
		echo "Name for replicate 3 is required!"
		exit
	fi
fi
if [ "$nreplicates" -eq "4" ];then
	if [ ! $REP4NAME ]; then
		echo "Name for replicate 4 is required!"
		exit
	fi
fi

if [ "$nreplicates" -eq "2" ];then
	if [ ! $TAGALIGN2 ]; then
		echo "sorted bam file for replicate 2 is required!"
		exit
	fi
fi
if [ "$nreplicates" -eq "3" ];then
	if [ ! $TAGALIGN3 ]; then
		echo "sorted bam file for replicate 3 is required!"
		exit
	fi
fi
if [ "$nreplicates" -eq "4" ];then
	if [ ! $TAGALIGN4 ]; then
		echo "sorted bam file for replicate 4 is required!"
		exit
	fi
fi

if [ "$nreplicates" -ge "2" ];then
	if [ ! $POOLEDPEAKFILE ]; then
		echo "Pooled peak file is required if replicates are present!"
		exit
	fi
	if [ ! $CONCENSUSBEDFILE ]; then
		echo "Concensus bed file is required if replicates are present!"
		exit
	fi
fi

extsize=$EXTSIZE
shiftsize=$SHIFTSIZE
tagAlign=$TAGALIGN
samplename=$SAMPLENAME
prefix=${samplename}.macs2
genome=$GENOME
genomefile=$GENOMEFILE

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
rm -f ${prefix}_peaks.narrowPeak.tmp
mv ${prefix}_peaks.narrowPeak ${prefix}.narrowPeak

# qvalue filter of 0.05
# awk -F"\t" '{if ($9>1.302){print}}' ${prefix}_peaks.narrowPeak > ${prefix}.qfilter.narrowPeak
# qvalue filter of 0.01
if [ $FILTERPEAKS == "True" ];then
  qvalue=$QFILTER
  awk -F"\t" -v q=$qvalue '{if ($9>q){print}}' ${prefix}.narrowPeak > ${prefix}.qfilter.narrowPeak
fi

if [ $SAVEBIGWIG == "True" ];then
  bedSort ${prefix}_treat_pileup.bdg ${prefix}_treat_pileup.bdg
#   sort-bed --max-mem 2G ${prefix}_treat_pileup.bdg > ${prefix}_treat_pileup.bdg.tmp
#   mv ${prefix}_treat_pileup.bdg.tmp ${prefix}_treat_pileup.bdg
  bedGraphToBigWig ${prefix}_treat_pileup.bdg $genomefile ${prefix}.bw
fi
if [ $SAVETN5KNICKSBED == "True" ]; then
  KNICKSBED=${samplename}.tn5knicks.bed
  KNICKSBAM=${samplename}.tn5knicks.bam
  zcat $tagAlign|awk -F"\t" -v OFS="\t" '{if ($6=="+") {print $1,$2,$2+1} else {print $1,$3,$3+1}}'> $KNICKSBED
  bedSort $KNICKSBED $KNICKSBED
#   sort-bed --max-mem 2G $KNICKBED > ${KNICKBED}.tmp
#   mv ${KNICKBED}.tmp $KNICKBED
  awk -F"\t" -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' $KNICKSBED > ${KNICKSBED%.*}.tmp.bed
  bedToBam -i ${KNICKSBED%.*}.tmp.bed -g $genomefile > $KNICKSBAM
  samtools sort -@4 -o ${KNICKSBAM%.*}.sorted.bam $KNICKSBAM
  mv ${KNICKSBAM%.*}.sorted.bam $KNICKSBAM
  rm -f ${KNICKSBED%.*}.tmp.bed 
  pigz -f -p4 $KNICKSBED
  samtools index $KNICKSBAM
fi
rm -f ${prefix}_treat_pileup.bdg
rm -f ${prefix}_control_lambda.bdg

conda deactivate
