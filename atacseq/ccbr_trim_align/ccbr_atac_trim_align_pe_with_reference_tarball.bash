#!/bin/bash
set -e -x -o pipefail
ARGPARSE_DESCRIPTION="""Trim PE reads using cutadapt and align for ATACSeq peakcalling.
This script runs 3 scripts in series:
1. ccbr_cutadapt_pe
2. ccbr_remove_blacklisted_reads_pe
3. ccbr_bowtie2_align_pe
It expects:
1. access to input fastq files (make sure the correct binding of folders;symlinks will not work
2. bowtie2 index files will be mounted at /index from a tarball provided"""
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--infastq1',required=True, help='input R1 fastq.gz file')
parser.add_argument('--infastq2',required=True, help='input R2 fastq.gz file')
parser.add_argument('--threads',required=True, help='number of threads')
parser.add_argument('--genome',required=True, help='hg19/hg38/mm9/mm10')
parser.add_argument('--samplename',required=True,help='sample name')
parser.add_argument('--indextarball',required=True, help='eg. mm10_bowtie2_index.tar.gz or a URL to the tarball')
parser.add_argument('--scriptsfolder',required=False, default='/opt2', help='folder where the scripts are... used for debuging without rebuilding the docker')
parser.add_argument('--keepfiles',required=False,default='False',help='to keep intermediate files, set this to True')
EOF

infastq1=$INFASTQ1
infastq2=$INFASTQ2
# samplename=`echo $infastq1|sed "s/.fastq.gz//g"|sed "s/.R1//g"|awk -F"/" '{print $NF}'`
trimmedfastq1="${SAMPLENAME}.R1.trim.fastq.gz"
trimmedfastq2="${SAMPLENAME}.R2.trim.fastq.gz"
noBLfastq1="${SAMPLENAME}.R1.noBL.fastq.gz"
noBLfastq2="${SAMPLENAME}.R2.noBL.fastq.gz"
ncpus=$THREADS

working_dir=$(pwd)
mkdir -p /index
y=$(echo $INDEXTARBALL|grep "://"|wc -l)
if [ $y == 1 ]; then 
	echo "Index tarball is a URL"
	cd /index
	wget $INDEXTARBALL
	filename=$(echo $INDEXTARBALL|awk -F"/" '{print $NF}')
	tar xzvf $filename
else 
	echo "Index tarball is Not a URL"
	cp $INDEXTARBALL /index
	cd /index
	tar xzvf $INDEXTARBALL
fi
cd $working_dir

bash ${SCRIPTSFOLDER}/ccbr_cutadapt_pe.bash \
--infastq1 $infastq1 \
--infastq2 $infastq2 \
--samplename $SAMPLENAME \
--threads $ncpus \
--outfastq1 $trimmedfastq1 \
--outfastq2 $trimmedfastq2

bash ${SCRIPTSFOLDER}/ccbr_remove_blacklisted_reads_pe.bash \
--infastq1 $trimmedfastq1 \
--infastq2 $trimmedfastq2 \
--samplename $SAMPLENAME \
--threads $ncpus \
--outfastq1 $noBLfastq1 \
--outfastq2 $noBLfastq2

bash ${SCRIPTSFOLDER}/ccbr_bowtie2_align_pe.bash \
--samplename $SAMPLENAME \
--threads $ncpus \
--infastq1 $noBLfastq1 \
--infastq2 $noBLfastq2 \
--genome $GENOME \
--scriptsfolder $SCRIPTSFOLDER \
--keepfiles $KEEPFILES

fastqc --threads $ncpus --format fastq $infastq1 $infastq2 $noBLfastq1 $noBLfastq2 --outdir $working_dir

if [ $KEEPFILES == "False" ];then
rm -f $trimmedfastq1 \
$trimmedfastq2 \
$noBLfastq1 \
$noBLfastq2
fi
