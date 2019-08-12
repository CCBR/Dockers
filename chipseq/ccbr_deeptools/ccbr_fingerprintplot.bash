#!/bin/bash
set -e -x -o pipefail
ARGPARSE_DESCRIPTION="Produce fingerprintplot and other metagene plots for ChiPSeq data"      # this is optional
source /opt/argparse.bash || exit 1
# source /home/kopardev/Desktop/Dockers/chipseq/ccbr_macs_peakcalling/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--treatmentbam',required=True, help='treatment bam file')
parser.add_argument('--treatmentbigwig',required=True, help='treatment bigwig file')
parser.add_argument('--inputbam',required=True, help='IgG or Input bam file')
parser.add_argument('--inputbigwig',required=True, help='IgG or Input bigwig file')
parser.add_argument('--genesbed',required=True, help='metagene bed file')
parser.add_argument('--ncpus',required=False, type=int, default=1, help='nproc')
EOF

t_label=${TREATMENTBAM%%.*}
i_label=${INPUTBAM%%.*}
fingerprintplot=${TREATMENTBAM%.*}_vs_${INPUTBAM%.*}.fingerprint.pdf
metageneheatmap=${TREATMENTBAM%.*}_vs_${INPUTBAM%.*}.metageneheatmap.pdf
metageneprofileplot=${TREATMENTBAM%.*}_vs_${INPUTBAM%.*}.metageneprofile.pdf

cpus=`nproc`

treatmentsortedbam=`echo $TREATMENTBAM|sed "s/.bam/.sorted.bam/g"`
inputsortedbam=`echo $INPUTBAM|sed "s/.bam/.sorted.bam/g"`

samtools sort -@ $cpus -o $treatmentsortedbam $TREATMENTBAM
samtools sort -@ $cpus -o $inputsortedbam $INPUTBAM

samtools index $treatmentsortedbam
samtools index $inputsortedbam

plotFingerprint -b $treatmentsortedbam $inputsortedbam --labels $t_label $i_label -p $cpus --skipZeros --plotFile $fingerprintplot

computeMatrix scale-regions -S $TREATMENTBIGWIG $INPUTBIGWIG -R $GENESBED -p $cpus --upstream 3000 --regionBodyLength 2000 --downstream 1000 --skipZeros -o matrix.txt --smartLabels

plotHeatmap -m matrix.txt -out $metageneheatmap --colorMap 'PuOr_r' --yAxisLabel 'average RPGC' --regionsLabel 'genes' --legendLocation 'none'

plotProfile -m matrix.txt -out $metageneprofileplot --plotHeight 15 --plotWidth 15 --perGroup --yAxisLabel 'average RPGC' --plotType 'se' --legendLocation upper-right