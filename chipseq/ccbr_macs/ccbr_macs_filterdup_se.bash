#!/bin/bash
set -e -x -o pipefail
ARGPARSE_DESCRIPTION="Filter duplicates with macs2"      # this is optional
source /opt/argparse.bash || exit 1
# source /home/kopardev/Desktop/Dockers/chipseq/ccbr_macs_peakcalling/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--bam',required=True, help='input bam file')
parser.add_argument('--genomesize',required=True, help='effective genome size',type=int)
EOF

tagalignfile=`echo $BAM|sed "s/.bam/.DD.tagAlign/g"`
outbamfile=`echo $BAM|sed "s/.bam/.DD.bam/g"`
cpus=`nproc`

macs2 filterdup -i $BAM -g $GENOMESIZE -o TmpTagAlign
awk -F"\t" -v OFS="\t" '{if ($2>0 && $3>0) {print}}' TmpTagAlign > TmpTagAlign2

samtools view -@ $cpus -H $BAM|grep "^@SQ"|cut -f2,3|sed "s/SN://"|sed "s/LN://g" > GenomeFile
awk -F"\t" -v OFS="\t" '{print $1,1,$2}' GenomeFile |sort -k1,1 -k2,2n > GenomeFileBed
bedtools intersect -wa -f 1.0 -a TmpTagAlign2 -b GenomeFileBed > TagAlign
cat TagAlign|wc -l |awk '{printf("%d\tAfter Deduplication\n",$1)}' > nreads.txt

bedtools bedtobam -i TagAlign -g GenomeFile > OutBam

mv TagAlign $tagalignfile
mv OutBam $outbamfile

pigz -p $cpus -f $tagalignfile