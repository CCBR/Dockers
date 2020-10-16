#!/bin/bash

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="Filter duplicates with picard"      # this is optional
source /opt/argparse.bash || exit 1
# source /home/kopardev/Desktop/Dockers/chipseq/ccbr_macs_peakcalling/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--samplename',required=True, help='Sample name')
parser.add_argument('--bam',required=True, help='input bam file')
parser.add_argument('--outBam',required=True, help='output bam file')

EOF

java -Xmx4G -jar /opt/picardcloud.jar MarkDuplicates \
INPUT=$BAM \
OUTPUT=$OUTBAM \
METRICS_FILE=${SAMPLENAME}.dupmetric \
VALIDATION_STRINGENCY=LENIENT \
ASSUME_SORTED=true \
REMOVE_DUPLICATES=true

samtools flagstat $OUTBAM > ${OUTBAM}.flagstat

nreads_mapped=`grep -m1 mapped ${OUTBAM}.flagstat|awk '{print $1}'`
echo "$nreads_mapped"|awk '{printf("%d\tQ5DD reads\n",$1)}' >> ${SAMPLENAME}.nreads.txt

python /opt/bam_pe_2_bedgraph.py -i $OUTBAM 