#!/bin/bash

set -e -x -o pipefail
. /opt2/conda/etc/profile.d/conda.sh
conda activate python3
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

ARGPARSE_DESCRIPTION="use preseq to find duplication levels and calculate NRF/PBC1/PBC2"      # this is optional
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--bam',required=True, help='sored bam file')
parser.add_argument('--preseq',required=True, help='preseq output file')
parser.add_argument('--preseqlog',required=True, help='preseqlog file')
parser.add_argument('--nrf',required=True, help='NRF output file')
parser.add_argument('--scriptsfolder',required=False, default='/opt2', help='folder where the scripts are... used for debuging without rebuilding the docker')
EOF

preseq lc_extrap -B -D -o $PRESEQ $BAM -seed 12345 -v -l 100000000000 2> $PRESEQLOG
python ${SCRIPTSFOLDER}/nrf.py $PRESEQLOG > $NRF

conda deactivate