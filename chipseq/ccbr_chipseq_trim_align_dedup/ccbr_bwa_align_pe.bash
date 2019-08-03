#!/bin/bash
set -e -x -o pipefail
ncpus=`nproc`
ARGPARSE_DESCRIPTION="Align PE reads using BWA and Q5 filter the bam file"      # this is optional
source /opt/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--samplename',required=True, help='Sample name')
parser.add_argument('--infastq1',required=True, help='input R1 fastq.gz file')
parser.add_argument('--infastq2',required=True, help='input R2 fastq.gz file')
parser.add_argument('--genome',required=True, help='hg19/hg38/mm9/mm10')
parser.add_argument('--outq5bam',required=True, help='output q5 bam file')

EOF

infastq1=$INFASTQ1
infastq2=$INFASTQ2
samplename=$SAMPLENAME
genome=$GENOME
log="${samplename}.bwa.log"
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

bwa mem -t $ncpus /index/$genome $infastq1 $infastq2 > ${samplename}.bwa.sam 2> $log

samtools view -@ $ncpus -bS -o ${samplename}.bwa.bam ${samplename}.bwa.sam
samtools sort -@ $ncpus -T ${samplename}.tmp -o ${samplename}.bwa.sorted.bam ${samplename}.bwa.bam
mv ${samplename}.bwa.sorted.bam ${samplename}.bwa.bam
samtools index ${samplename}.bwa.bam

samtools view -b ${samplename}.bwa.bam $CHROMOSOMES > ${samplename}.bwa.chrs.bam
samtools index ${samplename}.bwa.chrs.bam

python /opt/bam_filter_by_mapq.py -i ${samplename}.bwa.chrs.bam -o $OUTQ5BAM -q 6
samtools index $OUTQ5BAM

samtools flagstat ${samplename}.bwa.bam > ${samplename}.bwa.bam.flagstat
samtools flagstat $OUTQ5BAM > ${OUTQ5BAM}.flagstat

nreads_mapped=`grep -m1 mapped ${samplename}.bwa.bam.flagstat|awk '{print $1}'`
echo "$nreads_mapped"|awk '{printf("%d\tMapped reads\n",$1)}' >> ${samplename}.nreads.txt

nreads_mapped_q5=`grep -m1 mapped ${OUTQ5BAM}.flagstat|awk '{print $1}'`
echo "$nreads_mapped_q5"|awk '{printf("%d\tMapped Q5 reads\n",$1)}' >> ${samplename}.nreads.txt

rm -f ${samplename}.bwa.sam ${samplename}.bwa.chrs.bam*