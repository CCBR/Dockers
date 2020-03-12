#!/bin/bash
set -e -x -o pipefail

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
pigz -p4 $BED
pigz -p4 $KNICKSBED
samtools index $KNICKSBAM
# exit
}

ARGPARSE_DESCRIPTION="call atac-seq peaks using genrich v0.6 with 2 replicates" 
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--bamrep1',required=True, help='input name sorted bam file for replicate 1')
parser.add_argument('--bamrep2',required=True, help='input name sorted bam file for replicate 2')
parser.add_argument('--peakfile1',required=True, help='output narrowPeak file for replicate 1')
parser.add_argument('--peakfile2',required=True, help='output narrowPeak file for replicate 1')
parser.add_argument('--genome',required=True,help="hg19/38 or mm9/10")
parser.add_argument('--mergedpeakfile',required=True, help='output narrowPeak file for both replicates combined')
# only required if filtering
parser.add_argument('--filterpeaks',required=False, default="True", help='filterpeaks by qvalue: True or False')
parser.add_argument('--qfilter',required=False, default=0.693147, help='default qfiltering value is 0.693147 (-log10 of 0.5) for q=0.5')
parser.add_argument('--scriptsfolder',required=False, default='/opt2', help='folder where the scripts are... used for debuging without rebuilding the docker')
EOF



excludelist=$(samtools view -H $BAMREP1|grep ^@SQ|cut -f2|sed "s/SN://g"|awk "\$1 !~ /[0-9]$/"|grep -v "chrX"|grep -v "chrY"|tr "\n" ","|awk "{print substr(\$_,1,length(\$_)-1)}")
# echo $excludelist

# replicate 1 peak calling
READSBEDFILE1=${PEAKFILE1%.*}.reads.bed
READSBWFILE1=${PEAKFILE1%.*}.reads.bw
KNICKSBEDFILE1=${PEAKFILE1%.*}.tn5knicks.bed
KNICKSBAMFILE1=${PEAKFILE1%.*}.tn5knicks.bam

Genrich  -t $BAMREP1 -o $PEAKFILE1  -j -r -e $excludelist -v -s 5 -m 6 -b $READSBEDFILE1 -q 1 -l 200 -g 200

processReadsBed $BAMREP1 $READSBEDFILE1 $READSBWFILE1 $GENOME $KNICKSBEDFILE1 $KNICKSBAMFILE1

# replicate 2 peak calling
READSBEDFILE2=${PEAKFILE2%.*}.reads.bed
READSBWFILE2=${PEAKFILE2%.*}.reads.bw
KNICKSBEDFILE2=${PEAKFILE2%.*}.tn5knicks.bed
KNICKSBAMFILE2=${PEAKFILE2%.*}.tn5knicks.bam

Genrich -t $BAMREP2 -o $PEAKFILE2  -j -r -e $excludelist -v -s 5 -m 6 -b $READSBEDFILE2 -q 1 -l 200 -g 200

processReadsBed $BAMREP2 $READSBEDFILE2 $READSBWFILE2 $GENOME $KNICKSBEDFILE2 $KNICKSBAMFILE2


# merged peak calling
READSBEDFILE3=${MERGEDPEAKFILE%.*}.reads.bed
READSBWFILE3=${MERGEDPEAKFILE%.*}.reads.bw
KNICKSBEDFILE3=${MERGEDPEAKFILE%.*}.tn5knicks.bed
KNICKSBAMFILE3=${MERGEDPEAKFILE%.*}.tn5knicks.bam

Genrich  -t $BAMREP1,$BAMREP2 -o $MERGEDPEAKFILE  -j -r -e $excludelist -v -s 5 -m 6 -b $READSBEDFILE3 -q 1 -l 200 -g 200

processReadsBed $BAMREP2 $READSBEDFILE3 $READSBWFILE3 $GENOME $KNICKSBEDFILE3 $KNICKSBAMFILE3

if [ $FILTERPEAKS == "True" ];then
  qvalue=$QFILTER
  for f in $PEAKFILE1 $PEAKFILE2 $MERGEDPEAKFILE;do
    filteredpeakfile=$(echo $f|sed "s/.narrowPeak/.qfilter.narrowPeak/g")
    awk -F"\t" -v q=$qvalue "{if (\$9>q){print}}" $f > $filteredpeakfile
    Rscript ${SCRIPTSFOLDER}/ccbr_annotate_peaks.R -n $f -a ${f}.annotated -g $GENOME -l ${f}.genelist -f ${f}.annotation_summary 
    Rscript ${SCRIPTSFOLDER}/ccbr_annotate_peaks.R -n $filteredpeakfile -a ${filteredpeakfile}.annotated -g $GENOME -l ${filteredpeakfile}.genelist -f ${filteredpeakfile}.annotation_summary
    cut -f1,2 ${f}.annotation_summary > ${f}.annotation_distribution
    cut -f1,2 ${filteredpeakfile}.annotation_summary > ${filteredpeakfile}.annotation_distribution
  done
fi
