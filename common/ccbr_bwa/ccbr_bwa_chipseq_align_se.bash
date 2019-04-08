#!/bin/bash
set -e -x -o pipefail
cpus=`nproc`
ARGPARSE_DESCRIPTION="Remove blacklist reads and align SE ChIPSeq data to genome"      # this is optional
source /opt/argparse.bash || exit 1
# source /home/kopardev/Desktop/Dockers/chipseq/ccbr_macs_peakcalling/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--genome',required=True, help='mm9/10 or hg19/38')
parser.add_argument('--infastq',required=True, help='Input fastq file')
EOF

basefilename=`echo $INFASTQ|sed "s/.fastq.gz//g"`
bamfile=${basefilename}.bam
sortedbamfile=${basefilename}.sorted.bam
sortedbambaifile=${sortedbamfile}.bai
sortedq5bamfile=${basefilename}.sorted.Q5.bam
sortedq5bambaifile=${basefilename}.sorted.Q5.bam.bai
sortedbamflagstatfile=${sortedbamfile}.flagstat
sortedq5bamflagstatfile=${sortedq5bamfile}.flagstat

bwa mem -t $cpus ${GENOME}_blacklist $INFASTQ | samtools view -@ $cpus -f4 -b -o notAlignedToBlacklist.bam

bedtools bamtofastq -i notAlignedToBlacklist.bam -fq notAlignedToBlacklist.fastq

bwa mem -t $cpus $GENOME notAlignedToBlacklist.fastq > $bamfile

samtools sort -@ $cpus -o $sortedbamfile $bamfile

samtools index $sortedbamfile

samtools view -@ $cpus -b -q 6 $sortedbamfile -o $sortedq5bamfile

samtools index $sortedq5bamfile

samtools flagstat $sortedbamfile > $sortedbamflagstatfile

samtools flagstat $sortedq5bamfile > $sortedq5bamflagstatfile