#!/bin/bash


set -e -x -o pipefail
ARGPARSE_DESCRIPTION="Adapter trimmed (and blacklist filtered) fastqs are aligned to genome using bowtie2, multimappers are properly assigned, deduplicated using picard, filtered based on mapq, bams converted to tagAlign files."      # this is optional
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--infastq1',required=True, help='input R1 fastq.gz file')
parser.add_argument('--infastq2',required=True, help='input R2 fastq.gz file')
parser.add_argument('--samplename',required=True, help='samplename')
parser.add_argument('--threads',required=True, help='number of threads')
parser.add_argument('--genomename',required=True, help='hg19/hg38/mm9/mm10')
parser.add_argument('--multimapping',required=False, default=4, help='hg19/hg38/mm9/mm10')
parser.add_argument('--scriptsfolder',required=False, default='/opt2', help='folder where the scripts are... used for debuging without rebuilding the docker')
parser.add_argument('--keepfiles',required=False,default='False',help='to keep intermediate files, set this to True')
EOF

infastq1=$INFASTQ1
infastq2=$INFASTQ2
samplename=$SAMPLENAME
genome=$GENOMENAME
multimapping=$MULTIMAPPING
log="${samplename}.bowtie2.log"
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
ncpus=$THREADS

bowtie2 -X2000 -k $multimapping --very-sensitive --threads $ncpus -x /index/$genome \
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
${SCRIPTSFOLDER}/atac_assign_multimappers.py -k $multimapping --paired-end > ${samplename}.tmp2.sorted.sam
samtools view -@ $ncpus -bS -o ${samplename}.tmp3.bam ${samplename}.tmp2.sorted.sam
samtools sort -@ $ncpus ${samplename}.tmp3.bam ${samplename}.tmp3.sorted


bash ${SCRIPTSFOLDER}/ccbr_bam2nrf.bash \
--bam ${samplename}.tmp3.sorted.bam \
--preseq ${samplename}.preseq \
--preseqlog ${samplename}.preseq.log \
--nrf {samplename}.nrf \
--scriptsfolder $SCRIPTSFOLDER

samtools view -@ $ncpus -F 256 -u ${samplename}.tmp3.bam > ${samplename}.tmp4.bam
samtools sort -@ $ncpus ${samplename}.tmp4.bam ${samplename}.dup
samtools view -@ $ncpus -F 1796 -u ${samplename}.dup.bam > ${samplename}.tmp5.bam
samtools sort -@ $ncpus ${samplename}.tmp5.bam ${samplename}.filt
samtools index ${samplename}.filt.bam
samtools flagstat ${samplename}.filt.bam > ${samplename}.filt.bam.flagstat

java -Xmx4G -jar ${SCRIPTSFOLDER}/picardcloud.jar MarkDuplicates \
INPUT=${samplename}.filt.bam \
OUTPUT=${samplename}.dupmark.bam \
METRICS_FILE=${samplename}.dupmetric \
VALIDATION_STRINGENCY=LENIENT \
ASSUME_SORTED=true \
REMOVE_DUPLICATES=false

samtools view -F 1796 -b -@ $ncpus -o ${samplename}.dedup.tmp.bam ${samplename}.dupmark.bam

samtools index ${samplename}.dedup.tmp.bam

. /opt2/conda/etc/profile.d/conda.sh
conda activate python3
python ${SCRIPTSFOLDER}/ccbr_bam_filter_by_mapq.py -i ${samplename}.dedup.tmp.bam -o ${samplename}.dedup.bam -q 6
conda deactivate


samtools index ${samplename}.dedup.bam
samtools flagstat ${samplename}.dedup.bam > ${samplename}.dedup.bam.flagstat
samtools view -H ${samplename}.dedup.bam|grep "^@SQ"|cut -f2,3|sed "s/SN://g"|sed "s/LN://g" > ${samplename}.genome
# samtools sort -@ $ncpus -n ${samplename}.dedup.bam ${samplename}.dedup.qsorted

bedtools bamtobed -i ${samplename}.dedup.bam | \
awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}'| \
awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | \
gzip -c > ${samplename}.tagAlign.gz

nreads=`grep -m1 total ${samplename}.filt.bam.flagstat|awk '{print $1}'`
echo "$nreads"|awk '{printf("%d\tMapped reads\n",$1)}' >> ${samplename}.nreads.txt
nreads=`grep -m1 total ${samplename}.dedup.bam.flagstat|awk '{print $1}'`
echo "$nreads"|awk '{printf("%d\tAfter deduplication\n",$1)}' >> ${samplename}.nreads.txt

mv ${samplename}.tmp1.sorted.bam ${samplename}.qsorted.bam

if [ $KEEPFILES == "False" ];then
rm -f ${samplename}.tmp1.bam \
${samplename}.tmp1.sorted.sam \
${samplename}.tmp2.sorted.sam \
${samplename}.tmp3.sorted.bam \
${samplename}.tmp3.bam \
${samplename}.tmp4.bam \
${samplename}.tmp5.bam \
${samplename}.dup.bam \
${samplename}.dedup.tmp.bam* \
${samplename}.dupmark.bam \
${samplename}.bowtie2.bam \
${samplename}.bowtie2.sam \
${samplename}.filt.bam \
${samplename}.bowtie2.sorted.bam*
fi

# rm -f ${samplename}.dedup.qsorted.bam
