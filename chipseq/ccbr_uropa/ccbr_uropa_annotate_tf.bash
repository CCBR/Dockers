#!/bin/bash
set -e -x -o pipefail
cpus=`nproc`
ARGPARSE_DESCRIPTION="annotate bed file using uropa"      # this is optional
source /opt/argparse.bash || exit 1
# source /home/kopardev/Desktop/Dockers/chipseq/ccbr_macs_peakcalling/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--bed',required=True, help='called peaks in bed format')
parser.add_argument('--gtf',required=True, help='annotation file')
EOF

cat tf.json |sed "s/BEDFILEGOESHERE/${BED}/g" |sed "s/GTFFILEGOESHERE/${GTF}/g" > run_uropa.json
uropa -i run_uropa.json -s --threads `nproc`

