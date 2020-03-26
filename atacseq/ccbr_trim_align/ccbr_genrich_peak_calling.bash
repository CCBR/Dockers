#!/bin/bash
set -e -x -o pipefail
. /opt2/conda/etc/profile.d/conda.sh
conda activate python3

callGenrichPeaks(){
BAM=$1 
PEAKFILE=$2 
EXCLUDELIST=$3
READSBEDFILE=${PEAKFILE%.*}.reads.bed
READSBWFILE=${PEAKFILE%.*}.reads.bw
KNICKSBEDFILE=${PEAKFILE%.*}.tn5knicks.bed
KNICKSBAMFILE=${PEAKFILE%.*}.tn5knicks.bam 
Genrich -t $BAM -o $PEAKFILE -j -r -e $EXCLUDELIST -v -s 5 -m 6 -b $READSBEDFILE -q 1 -l 200 -g 200
processReadsBed $BAM $READSBEDFILE $READSBWFILE $GENOME $KNICKSBEDFILE $KNICKSBAMFILE
}

processReadsBed() {
BAM=$1
BED=$2
BW=$3
GENOME=$4
KNICKSBED=$5
KNICKSBAM=$6
if [ "$GENOME" == "hg19" ]; then
effectiveSize=2700000000
elif [ "$GENOME" == "hg38" ]; then
effectiveSize=2700000000
elif [ "$GENOME" == "mm9" ]; then
effectiveSize=2400000000
elif [ "$GENOME" == "mm10" ]; then
effectiveSize=2400000000
fi
samtools view -H $BAM|grep ^@SQ|cut -f2,3|sed "s/SN://g"|sed "s/LN://g" > ${BAM}.genome
bedSort $BED $BED
sf=$(awk -F"\t" -v size=$effectiveSize '{sum=sum+$3-$2}END{print sum/size}' $BED)
genomeCoverageBed -bg -i $BED -g ${BAM}.genome |awk -F"\t" -v OFS="\t" -v sf=$sf '{print $1,$2,$3,$4/sf}' > ${BED}.bg
bedGraphToBigWig ${BED}.bg ${BAM}.genome $BW

awk -F"\t" -v OFS="\t" '{if ($3-$2==100) {print $1,$2+50,$2+51} }' $BED > $KNICKSBED
awk -F"\t" -v OFS="\t" '{if ($3-$2!=100) {print $1,$2+50,$2+51} }' $BED >> $KNICKSBED
awk -F"\t" -v OFS="\t" '{if ($3-$2!=100) {print $1,$3-51,$3-50} }' $BED >> $KNICKSBED
bedSort $KNICKSBED $KNICKSBED
awk -F"\t" -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' $KNICKSBED > ${KNICKSBED%.*}.tmp.bed
bedToBam -i ${KNICKSBED%.*}.tmp.bed -g ${BAM}.genome > $KNICKSBAM
samtools sort -@4 -o ${KNICKSBAM%.*}.sorted.bam $KNICKSBAM
mv ${KNICKSBAM%.*}.sorted.bam $KNICKSBAM
rm -f ${KNICKSBED%.*}.tmp.bed ${BED}.bg
pigz -f -p4 $BED
pigz -f -p4 $KNICKSBED
samtools index $KNICKSBAM
# exit
}

ARGPARSE_DESCRIPTION="call atac-seq peaks using genrich v0.6 with 2 replicates" 
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--bamrep1',required=True, help='input name sorted bam file for replicate 1')
parser.add_argument('--bamrep2',required=False, help='input name sorted bam file for replicate 2')
parser.add_argument('--bamrep3',required=False, help='input name sorted bam file for replicate 3')
parser.add_argument('--bamrep4',required=False, help='input name sorted bam file for replicate 4')

parser.add_argument('--peakfile1',required=True, help='output narrowPeak file for replicate 1')
parser.add_argument('--peakfile2',required=False, help='output narrowPeak file for replicate 2')
parser.add_argument('--peakfile3',required=False, help='output narrowPeak file for replicate 3')
parser.add_argument('--peakfile4',required=False, help='output narrowPeak file for replicate 4')

parser.add_argument('--genome',required=True,help="hg19/38 or mm9/10")
parser.add_argument('--pooledpeakfile',required=False, help='output narrowPeak file for both replicates combined')
parser.add_argument('--concensusbedfile',required=False, help='consensus overlapping MAX majority peaks in bed format')

# only required if filtering
parser.add_argument('--filterpeaks',required=False, default="True", help='filterpeaks by qvalue: True or False')
parser.add_argument('--qfilter',required=False, default=0.693147, help='default qfiltering value is 0.693147 (-log10 of 0.5) for q=0.5')

parser.add_argument('--scriptsfolder',required=False, default='/opt2', help='folder where the scripts are... used for debuging without rebuilding the docker')
EOF

nreplicates=1
if [ $BAMREP2 ];then
	nreplicates=2
	if [ $BAMREP3 ];then
		nreplicates=3
		if [ $BAMREP4 ];then
			nreplicates=4
		fi
	fi
fi

if [ $PEAKFILE2 ];then
	nreplicates=2
	if [ $PEAKFILE3 ];then
		nreplicates=3
		if [ $PEAKFILE4 ];then
			nreplicates=4
		fi
	fi
fi

if [ "$nreplicates" -eq "2" ];then
	if [ ! $PEAKFILE2 ]; then
		echo "output narrowPeak file for replicate 2 is required!"
		exit
	fi
fi
if [ "$nreplicates" -eq "3" ];then
	if [ ! $PEAKFILE3 ]; then
		echo "output narrowPeak file for replicate 3 is required!"
		exit
	fi
fi
if [ "$nreplicates" -eq "4" ];then
	if [ ! $PEAKFILE4 ]; then
		echo "output narrowPeak file for replicate 4 is required!"
		exit
	fi
fi

if [ "$nreplicates" -eq "2" ];then
	if [ ! $BAMREP2 ]; then
		echo "sorted bam file for replicate 2 is required!"
		exit
	fi
fi
if [ "$nreplicates" -eq "3" ];then
	if [ ! $BAMREP3 ]; then
		echo "sorted bam file for replicate 3 is required!"
		exit
	fi
fi
if [ "$nreplicates" -eq "4" ];then
	if [ ! $BAMREP4 ]; then
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

excludelist=$(samtools view -H $BAMREP1|grep ^@SQ|cut -f2|sed "s/SN://g"|awk "\$1 !~ /[0-9]$/"|grep -v "chrX"|grep -v "chrY"|tr "\n" ","|awk "{print substr(\$_,1,length(\$_)-1)}")
# echo $excludelist

# replicate 1 peak calling
callGenrichPeaks $BAMREP1 $PEAKFILE1 "$excludelist"
# 
# replicate 2 peak calling
if [ "$nreplicates" -ge 2 ]; then
callGenrichPeaks $BAMREP2 $PEAKFILE2 "$excludelist"
fi
# 
# replicate 3 peak calling
if [ "$nreplicates" -ge 3 ]; then
callGenrichPeaks $BAMREP3 $PEAKFILE3 "$excludelist"
fi
# 
# replicate 4 peak calling
if [ "$nreplicates" -ge 4 ]; then
callGenrichPeaks $BAMREP4 $PEAKFILE4 "$excludelist"
fi


if [ "$nreplicates" -ge 2 ]; then

	# merged peak calling and concensus peak calling
	READSBEDFILE=${POOLEDPEAKFILE%.*}.reads.bed
	READSBWFILE=${POOLEDPEAKFILE%.*}.reads.bw
	KNICKSBEDFILE=${POOLEDPEAKFILE%.*}.tn5knicks.bed
	KNICKSBAMFILE=${POOLEDPEAKFILE%.*}.tn5knicks.bam

	if [ "$nreplicates" -eq "2" ];then
	Genrich  -t $BAMREP1,$BAMREP2 -o $POOLEDPEAKFILE  -j -r -e $excludelist -v -s 5 -m 6 -b $READSBEDFILE -q 1 -l 200 -g 200
	python ${SCRIPTSFOLDER}/ccbr_get_concensus_peaks.py --peakfiles $PEAKFILE1 $PEAKFILE2 --outbed $CONCENSUSBEDFILE
	elif [ "$nreplicates" -eq "3" ];then
	Genrich  -t $BAMREP1,$BAMREP2,$BAMREP3 -o $POOLEDPEAKFILE  -j -r -e $excludelist -v -s 5 -m 6 -b $READSBEDFILE -q 1 -l 200 -g 200
	python ${SCRIPTSFOLDER}/ccbr_get_concensus_peaks.py --peakfiles $PEAKFILE1 $PEAKFILE2 $PEAKFILE3 --outbed $CONCENSUSBEDFILE
	elif [ "$nreplicates" -eq "4" ];then
	Genrich  -t $BAMREP1,$BAMREP2,$BAMREP3,$BAMREP4 -o $POOLEDPEAKFILE  -j -r -e $excludelist -v -s 5 -m 6 -b $READSBEDFILE -q 1 -l 200 -g 200
	python ${SCRIPTSFOLDER}/ccbr_get_concensus_peaks.py --peakfiles $PEAKFILE1 $PEAKFILE2 $PEAKFILE3 $PEAKFILE4 --outbed $CONCENSUSBEDFILE
	fi

	processReadsBed $BAMREP2 $READSBEDFILE $READSBWFILE $GENOME $KNICKSBEDFILE $KNICKSBAMFILE
fi

files="$PEAKFILE1"
if [ "$nreplicates" -eq "2" ];then files="$PEAKFILE1 $PEAKFILE2 $POOLEDPEAKFILE"; fi
if [ "$nreplicates" -eq "3" ];then files="$PEAKFILE1 $PEAKFILE2 $PEAKFILE3 $POOLEDPEAKFILE"; fi
if [ "$nreplicates" -eq "4" ];then files="$PEAKFILE1 $PEAKFILE2 $PEAKFILE3 $PEAKFILE4 $POOLEDPEAKFILE"; fi

Rscript ${SCRIPTSFOLDER}/ccbr_annotate_bed.R -b $CONCENSUSBEDFILE -a ${f}.annotated -g $GENOME -l ${f}.genelist -f ${f}.annotation_summary
cut -f1,2 ${CONCENSUSBEDFILE}.annotation_summary > ${CONCENSUSBEDFILE}.annotation_distribution

if [ $FILTERPEAKS == "True" ];then
  qvalue=$QFILTER
  for f in $files;do
	filteredpeakfile=$(echo $f|sed "s/.narrowPeak/.qfilter.narrowPeak/g")
	awk -F"\t" -v q=$qvalue "{if (\$9>q){print}}" $f > $filteredpeakfile
	Rscript ${SCRIPTSFOLDER}/ccbr_annotate_peaks.R -n $f -a ${f}.annotated -g $GENOME -l ${f}.genelist -f ${f}.annotation_summary 
	Rscript ${SCRIPTSFOLDER}/ccbr_annotate_peaks.R -n $filteredpeakfile -a ${filteredpeakfile}.annotated -g $GENOME -l ${filteredpeakfile}.genelist -f ${filteredpeakfile}.annotation_summary
	cut -f1,2 ${f}.annotation_summary > ${f}.annotation_distribution
	cut -f1,2 ${filteredpeakfile}.annotation_summary > ${filteredpeakfile}.annotation_distribution
  done
fi

conda deactivate