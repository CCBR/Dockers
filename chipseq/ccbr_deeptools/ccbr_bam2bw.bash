#!/bin/bash
set -e -x -o pipefail
ARGPARSE_DESCRIPTION="convert treatment & input bams to bigwigs"      # this is optional
source /opt/argparse.bash || exit 1
# source /home/kopardev/Desktop/Dockers/chipseq/ccbr_macs_peakcalling/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--treatmentbam',required=True, help='treatment bam file')
parser.add_argument('--treatmentppqt',required=True, help='treatment bam file')
parser.add_argument('--inputbam',required=True, help='IgG or Input bam file')
parser.add_argument('--effectivegenomesize',required=True, type=int, help='Effective genome size')
parser.add_argument('--ncpus',required=False, type=int, default=1, help='nproc')

EOF


treatmentbw=`echo $TREATMENTBAM|sed "s/.bam/.bw/g"`
inputbw=`echo $INPUTBAM|sed "s/.bam/.bw/g"`
inputnormbw=`echo $TREATMENTBAM|sed "s/.bam/.inputnorm.bw/g"`

treatmentextsize=`cat $TREATMENTPPQT|awk -F"\t" '{print $3}'|awk -F"," '{print $1}'`
if [ $treatmentextsize -lt 0 ];then
    treatmentextsize=`cat $TREATMENTPPQT|awk -F"\t" '{print $3}'|awk -F"," '{print $2}'`
fi  
if [ $treatmentextsize -lt 0 ];then
    treatmentextsize="200"
fi  

treatmentsortedbam=`echo $TREATMENTBAM|sed "s/.bam/.sorted.bam/g"`
inputsortedbam=`echo $INPUTBAM|sed "s/.bam/.sorted.bam/g"`

samtools sort -@ $NCPUS -o $treatmentsortedbam $TREATMENTBAM
samtools sort -@ $NCPUS -o $inputsortedbam $INPUTBAM

samtools index $treatmentsortedbam
samtools index $inputsortedbam

bamCoverage --bam $treatmentsortedbam -o $treatmentbw --binSize 25 --smoothLength 75 --ignoreForNormalization chrX chrY chrM --numberOfProcessors $NCPUS --normalizeUsing RPGC --effectiveGenomeSize $EFFECTIVEGENOMESIZE --extendReads $treatmentextsize

bamCoverage --bam $inputsortedbam -o $inputbw --binSize 25 --smoothLength 75 --ignoreForNormalization chrX chrY chrM --numberOfProcessors $NCPUS --normalizeUsing RPGC --effectiveGenomeSize $EFFECTIVEGENOMESIZE --extendReads 200

bigwigCompare --binSize 25 --outFileName $inputnormbw --outFileFormat 'bigwig' --bigwig1 $treatmentbw --bigwig2 $inputbw --operation 'subtract' --skipNonCoveredRegions --numberOfProcessors $NCPUS
