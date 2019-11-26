#!/bin/bash
set -e -x -o pipefail
ARGPARSE_DESCRIPTION="""Trim PE reads using cutadapt and align for ATACSeq peakcalling.
This script runs 3 scripts in series:
1. ccbr_cutadapt_pe
2. ccbr_remove_blacklisted_reads_pe
3. ccbr_bowtie2_align_pe
It expects:
1. access to input fastq files (make sure the correct binding of folders;symlinks will not work
2. bowtie2 index files need to be bound/mounted at /index eg. for mm10 on biowulf /data/CCBR_Pipeliner/db/PipeDB/Indices/mm10_basic/indexes/ needs to be bound to /index"""
source /opt/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--infastq1',required=True, help='input R1 fastq.gz file')
parser.add_argument('--infastq2',required=True, help='input R2 fastq.gz file')
parser.add_argument('--threads',required=True, help='number of threads')
parser.add_argument('--genome',required=True, help='hg19/hg38/mm9/mm10')
parser.add_argument('--scriptsfolder',required=False, default='/opt', help='folder where the scripts are... used for debuging without rebuilding the docker')
parser.add_argument('--keepfiles',required=False,default='False',help='to keep intermediate files, set this to True')
EOF

infastq1=$INFASTQ1
infastq2=$INFASTQ2
samplename=`echo $infastq1|sed "s/.fastq.gz//g"|sed "s/.R1//g"|awk -F"/" '{print $NF}'`
trimmedfastq1="${samplename}.R1.trim.fastq.gz"
trimmedfastq2="${samplename}.R2.trim.fastq.gz"
noBLfastq1="${samplename}.R1.noBL.fastq.gz"
noBLfastq2="${samplename}.R2.noBL.fastq.gz"
ncpus=$THREADS

bash ${SCRIPTSFOLDER}/ccbr_cutadapt_pe.bash \
--infastq1 $infastq1 \
--infastq2 $infastq2 \
--samplename $samplename \
--threads $ncpus \
--outfastq1 $trimmedfastq1 \
--outfastq2 $trimmedfastq2

bash ${SCRIPTSFOLDER}/ccbr_remove_blacklisted_reads_pe.bash \
--infastq1 $trimmedfastq1 \
--infastq2 $trimmedfastq2 \
--samplename $samplename \
--threads $ncpus \
--outfastq1 $noBLfastq1 \
--outfastq2 $noBLfastq2

bash ${SCRIPTSFOLDER}/ccbr_bowtie2_align_pe.bash \
--samplename $samplename \
--threads $ncpus \
--infastq1 $noBLfastq1 \
--infastq2 $noBLfastq2 \
--genome $GENOME \
--scriptsfolder $SCRIPTSFOLDER \
--keepfiles $KEEPFILES

fastqc --threads $ncpus --format fastq $infastq1 $infastq2 $noBLfastq1 $noBLfastq2

if [ $KEEPFILES == "False" ];then
rm -f $trimmedfastq1 \
$trimmedfastq2 \
$noBLfastq1 \
$noBLfastq2
fi
