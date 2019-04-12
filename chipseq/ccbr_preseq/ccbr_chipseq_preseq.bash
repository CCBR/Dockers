#!/bin/bash
set -e -x -o pipefail
cpus=`nproc`
ARGPARSE_DESCRIPTION="ccurve and nrf calculation"      # this is optional
source /opt/argparse.bash || exit 1
# source /home/kopardev/Desktop/Dockers/chipseq/ccbr_macs_peakcalling/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--bam',required=True, help='bam file')
parser.add_argument('--ncpus',required=False, default=1, help='number of processors')
parser.add_argument('--paired', action='store_true',
                    default=False, help='Are reads paired?')
EOF

BAMBASE=`echo ${BAM%.*}`

if [ $PAIRED ]; then
	paired="-P"
else 
	paired=""
fi

preseq c_curve $paired -B -o ${BAMBASE}.c_curve $BAM
preseq lc_extrap $paired -B -D -o ${BAMBASE}.preseq $BAM -seed 12345 -v -l 100000000000 2> ${BAMBASE}.preseqlog
# cat ${BAMBASE}.preseqlog
# python /opt/ccbr_nrf.py ${BAMBASE}.preseqlog
python /opt/ccbr_nrf.py ${BAMBASE}.preseqlog > ${BAMBASE}.nrf.txt
samtools sort -@ $NCPUS -o sorted.bam $BAM
samtools index sorted.bam
samtools idxstats sorted.bam > ${BAM}.idxstats
rm -f sorted.bam*