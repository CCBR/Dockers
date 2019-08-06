#!/bin/bash
set -e -x -o pipefail
ncpus=`nproc`
ARGPARSE_DESCRIPTION="Trim PE reads using cutadapt"      # this is optional
source /opt/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--infastq1',required=True, help='input R1 fastq.gz file')
parser.add_argument('--infastq2',required=True, help='input R2 fastq.gz file')
parser.add_argument('--samplename',required=True, help='samplename')
parser.add_argument('--genomename',required=True, help='hg19/hg38/mm9/mm10')
parser.add_argument('--multimapping',required=False, default=4, help='hg19/hg38/mm9/mm10')
EOF

infastq1=$INFASTQ1
infastq2=$INFASTQ2
samplename=$SAMPLENAME
genome=$GENOMENAME
multimapping=$MULTIMAPPING
log="${samplename}.bowtie2.log"
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"


bowtie2 -X2000 -k $multimapping --local --mm --threads $ncpus -x /index/$genome \
 -1 $infastq1 -2 $infastq2 > ${samplename}.bowtie2.sam \
 2> $log


samtools view -@ $ncpus -bS -o ${samplename}.bowtie2.bam ${samplename}.bowtie2.sam

samtools sort -@ $ncpus ${samplename}.bowtie2.bam ${samplename}.bowtie2.sorted

samtools flagstat ${samplename}.bowtie2.bam > ${samplename}.bowtie2.bam.flagstat

samtools index ${samplename}.bowtie2.sorted.bam
samtools view -@ $ncpus -F 516 -u ${samplename}.bowtie2.sorted.bam $CHROMOSOMES > ${samplename}.tmp1.bam
samtools sort -@ $ncpus -n ${samplename}.tmp1.bam ${samplename}.tmp1.sorted
samtools view -@ $ncpus -h ${samplename}.tmp1.sorted.bam > ${samplename}.tmp1.sorted.sam
cat ${samplename}.tmp1.sorted.sam | \
atac_assign_multimappers.py -k 4 > ${samplename}.tmp2.sorted.sam
samtools view -@ $ncpus -bS -o ${samplename}.tmp3.bam ${samplename}.tmp2.sorted.sam

samtools view -@ $ncpus -F 256 -u ${samplename}.tmp3.bam > ${samplename}.tmp4.bam 
samtools sort -@ $ncpus ${samplename}.tmp4.bam ${samplename}.dup
samtools view -@ $ncpus -F 1796 -u ${samplename}.tmp3.bam > ${samplename}.tmp5.bam 
samtools sort -@ $ncpus ${samplename}.tmp5.bam ${samplename}.filt
samtools index ${samplename}.filt.bam
samtools flagstat ${samplename}.filt.bam > ${samplename}.filt.bam.flagstat

java -Xmx4G -jar /opt/picardcloud.jar MarkDuplicates \
INPUT=${samplename}.filt.bam \
OUTPUT=${samplename}.dupmark.bam \
METRICS_FILE=${samplename}.dupmetric \
VALIDATION_STRINGENCY=LENIENT \
ASSUME_SORTED=true \
REMOVE_DUPLICATES=false

samtools view -F 1796 -b -@ $ncpus -o ${samplename}.dedup.bam ${samplename}.dupmark.bam
samtools index ${samplename}.dedup.bam
samtools flagstat ${samplename}.dedup.bam > ${samplename}.dedup.bam.flagstat
samtools sort -@ $ncpus -n ${samplename}.dedup.bam ${samplename}.dedup.qsorted

bedtools bamtobed -i ${samplename}.dedup.bam | \
awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}'| \
awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | \
gzip -c > ${samplename}.tagAlign.gz

nreads=`grep -m1 total ${samplename}.filt.bam.flagstat|awk '{print $1}'`
echo "$nreads"|awk '{printf("%d\tMapped reads\n",$1)}' >> ${samplename}.nreads.txt
nreads=`grep -m1 total ${samplename}.dedup.bam.flagstat|awk '{print $1}'`
echo "$nreads"|awk '{printf("%d\tAfter deduplication\n",$1)}' >> ${samplename}.nreads.txt

rm -f ${samplename}.tmp1.bam \
${samplename}.tmp1.sorted.sam \
${samplename}.tmp1.sorted.bam \
${samplename}.tmp2.sorted.sam \
${samplename}.tmp3.bam \
${samplename}.tmp4.bam \
${samplename}.tmp5.bam \
${samplename}.dup.bam \
${samplename}.dupmark.bam \
${samplename}.dedup.qsorted.bam \
${samplename}.bowtie2.bam \
${samplename}.bowtie2.sam \
${samplename}.bowtie2.sorted.bam* 

