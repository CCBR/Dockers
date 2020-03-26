#!/bin/bash
. /opt2/conda/etc/profile.d/conda.sh
conda activate python3

callPeaks(){
# inputs
	TAGALIGN=$1
	PREFIX=$2
	GENOME=$3
	SHIFTSIZE=$4
	EXTSIZE=$5
	FILTERPEAKS=$6
	QFILTER=$7
	GENOMEFILE=$8
	SCRIPTSFOLDER=$9

# call peaks
	if [ $GENOME == "hg19" ]; then g="hs";fi
	if [ $GENOME == "hg38" ]; then g="hs";fi
	if [ $GENOME == "mm10" ]; then g="mm";fi
	if [ $GENOME == "mm9" ]; then g="mm";fi
	macs2 callpeak -t $TAGALIGN -f BED -n $PREFIX -g $g -p 0.01 --shift -$SHIFTSIZE --extsize $EXTSIZE --keep-dup all --call-summits --nomodel --SPMR -B

# remove duplicate peak calls
	mv ${PREFIX}_peaks.narrowPeak ${PREFIX}_peaks.narrowPeak.tmp
	sort -k9,9gr ${PREFIX}_peaks.narrowPeak.tmp|awk -F"\t" '!NF || !seen[$1":"$2"-"$3]++'|sort -k1,1 -k2,2n > ${PREFIX}_peaks.narrowPeak
	rm -f ${PREFIX}_peaks.narrowPeak.tmp
	mv ${PREFIX}_peaks.narrowPeak ${PREFIX}.narrowPeak

# qvalue filter of 0.5 ... very lenient
	if [ $FILTERPEAKS == "True" ];then
	  qvalue=$QFILTER
	  awk -F"\t" -v q=$qvalue '{if ($9>q){print}}' ${PREFIX}.narrowPeak > ${PREFIX}.qfilter.narrowPeak
	fi

# annotate
	for f in ${PREFIX}.narrowPeak ${PREFIX}.qfilter.narrowPeak;do
		Rscript ${SCRIPTSFOLDER}/ccbr_annotate_peaks.R -n $f  -a ${f}.annotated -g $GENOME -l ${f}.genelist -f ${f}.annotation_summary 
		cut -f1,2 ${f}.annotation_summary > ${f}.annotation_distribution
	done

# save bigwig
# MAY NEED TO BE NORMALIZED FOR GENOME SIZE ... DONT KNOW FOR SURE
# MACS says
#
# MACS will save fragment pileup signal per million reads
#
# Genrich normalization is as per 1x genome coverage... and macs2 normalization is per million reads 
	bedSort ${PREFIX}_treat_pileup.bdg ${PREFIX}_treat_pileup.bdg
	bedGraphToBigWig ${PREFIX}_treat_pileup.bdg $GENOMEFILE ${PREFIX}.bw
	rm -f ${PREFIX}_treat_pileup.bdg
	rm -f ${PREFIX}_control_lambda.bdg	

# save knicks bam and bed
  KNICKSBED=${PREFIX}.tn5knicks.bed
  KNICKSBAM=${PREFIX}.tn5knicks.bam
  zcat $TAGALIGN|awk -F"\t" -v OFS="\t" '{if ($6=="+") {print $1,$2,$2+1} else {print $1,$3,$3+1}}'> $KNICKSBED
  bedSort $KNICKSBED $KNICKSBED
  awk -F"\t" -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' $KNICKSBED > ${KNICKSBED%.*}.tmp.bed
  bedToBam -i ${KNICKSBED%.*}.tmp.bed -g $GENOMEFILE > $KNICKSBAM
  samtools sort -@4 -o ${KNICKSBAM%.*}.sorted.bam $KNICKSBAM
  mv ${KNICKSBAM%.*}.sorted.bam $KNICKSBAM
  rm -f ${KNICKSBED%.*}.tmp.bed 
  pigz -f -p4 $KNICKSBED
  samtools index $KNICKSBAM
}

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
# parser.add_argument('--pooledpeakfile',required=False, help='output narrowPeak file for all replicates combined') # pooled peakfile will be named <samplename>.macs.narrowPeak
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
	if [ ! $CONCENSUSBEDFILE ]; then
		echo "Concensus bed file is required if replicates are present!"
		exit
	fi
fi


# replicate 1 peak calling
callPeaks $TAGALIGN1 ${REP1NAME}.macs2 $GENOME $SHIFTSIZE $EXTSIZE $FILTERPEAKS $QFILTER $GENOMEFILE $SCRIPTSFOLDER
# 
# replicate 2 peak calling
if [ "$nreplicates" -ge 2 ]; then
callPeaks $TAGALIGN2 ${REP2NAME}.macs2 $GENOME $SHIFTSIZE $EXTSIZE $FILTERPEAKS $QFILTER $GENOMEFILE $SCRIPTSFOLDER
fi
# 
# replicate 3 peak calling
if [ "$nreplicates" -ge 3 ]; then
callPeaks $TAGALIGN3 ${REP3NAME}.macs2 $GENOME $SHIFTSIZE $EXTSIZE $FILTERPEAKS $QFILTER $GENOMEFILE $SCRIPTSFOLDER
fi
# 
# replicate 4 peak calling
if [ "$nreplicates" -ge 4 ]; then
callPeaks $TAGALIGN4 ${REP4NAME}.macs2 $GENOME $SHIFTSIZE $EXTSIZE $FILTERPEAKS $QFILTER $GENOMEFILE $SCRIPTSFOLDER
fi

if [ "$nreplicates" -eq 2 ];then
pooled="${REP1NAME}_${REP2NAME}.pooled.tagAlign"
zcat $TAGALIGN1 $TAGALIGN2 > ${pooled}
fi

if [ "$nreplicates" -eq 3 ];then
pooled="${REP1NAME}_${REP2NAME}_${REP3NAME}.pooled.tagAlign"
zcat $TAGALIGN1 $TAGALIGN2 $TAGALIGN3 > ${pooled}
fi

if [ "$nreplicates" -eq 4 ];then
pooled="${REP1NAME}_${REP2NAME}_${REP3NAME}_${REP4NAME}.pooled.tagAlign"
zcat $TAGALIGN1 $TAGALIGN2 $TAGALIGN3 $TAGALIGN4 > ${pooled}
fi

# pooled peak calling
bedSort ${pooled} ${pooled}
pigz -f -p4 ${pooled}
callPeaks ${pooled}.gz ${SAMPLENAME}.macs2 $GENOME $SHIFTSIZE $EXTSIZE $FILTERPEAKS $QFILTER $GENOMEFILE $SCRIPTSFOLDER

# concensus peak calling
if [ "$nreplicates" -eq 2 ];then
python ${SCRIPTSFOLDER}/ccbr_get_concensus_peaks.py --peakfiles ${REP1NAME}.macs2.narrowPeak ${REP2NAME}.macs2.narrowPeak --outbed $CONCENSUSBEDFILE
fi

if [ "$nreplicates" -eq 3 ];then
python ${SCRIPTSFOLDER}/ccbr_get_concensus_peaks.py --peakfiles ${REP1NAME}.macs2.narrowPeak ${REP2NAME}.macs2.narrowPeak ${REP3NAME}.macs2.narrowPeak --outbed $CONCENSUSBEDFILE
fi

if [ "$nreplicates" -eq 4 ];then
python ${SCRIPTSFOLDER}/ccbr_get_concensus_peaks.py --peakfiles ${REP1NAME}.macs2.narrowPeak ${REP2NAME}.macs2.narrowPeak ${REP3NAME}.macs2.narrowPeak ${REP4NAME}.macs2.narrowPeak --outbed $CONCENSUSBEDFILE
fi

Rscript ${SCRIPTSFOLDER}/ccbr_annotate_bed.R -b $CONCENSUSBEDFILE -a ${f}.annotated -g $GENOME -l ${f}.genelist -f ${f}.annotation_summary
cut -f1,2 ${CONCENSUSBEDFILE}.annotation_summary > ${CONCENSUSBEDFILE}.annotation_distribution

conda deactivate
