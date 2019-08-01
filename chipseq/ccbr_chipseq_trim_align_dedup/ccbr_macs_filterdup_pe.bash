#!/bin/bash

. /opt/conda/etc/profile.d/conda.sh
conda activate macs

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="Filter duplicates with macs2"      # this is optional
source /opt/argparse.bash || exit 1
# source /home/kopardev/Desktop/Dockers/chipseq/ccbr_macs_peakcalling/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--samplename',required=True, help='Sample name')
parser.add_argument('--bam',required=True, help='input bam file')
parser.add_argument('--genomesize',required=True, help='effective genome size',type=int)
parser.add_argument('--outTagAlign',required=True, help='gzipped Q5DD output tagAlign')
parser.add_argument('--outBam',required=True, help='output bam file')

EOF

cpus=`nproc`

macs2 filterdup -f BAMPE -i $BAM -g $GENOMESIZE -o TmpTagAlign
awk -F"\t" -v OFS="\t" '{if ($2>0 && $3>0) {print}}' TmpTagAlign > TmpTagAlign2

samtools view -@ $cpus -H $BAM|grep "^@SQ"|cut -f2,3|sed "s/SN://"|sed "s/LN://g" > GenomeFile
awk -F"\t" -v OFS="\t" '{print $1,1,$2}' GenomeFile |sort -k1,1 -k2,2n > GenomeFileBed
bedtools intersect -wa -f 1.0 -a TmpTagAlign2 -b GenomeFileBed > TagAlign

awk -F"\t" '{x++;printf("%s\tr%d\n",$0,x)}' TagAlign |bedtools bedtobam -g GenomeFile -i - > OutBam

pigz -p $cpus -f TagAlign
mv TagAlign.gz $OUTTAGALIGN

samtools sort -@ $cpus -o $OUTBAM OutBam
samtools flagstat $OUTBAM > ${OUTBAM}.flagstat

rm -f TmpTagAlign TmpTagAlign2 GenomeFile GenomeFileBed OutBam

nreads_q5dd=`zcat $OUTTAGALIGN|wc -l`
echo "$nreads_q5dd"|awk '{printf("%d\tQ5DD reads\n",2 * $1)}' >> ${SAMPLENAME}.nreads.txt

conda deactivate
