#!/bin/bash
. /opt/conda/etc/profile.d/conda.sh
conda activate python2

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="calculated FLD from dedup bam file"      # this is optional
source /opt/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--dedupbam',required=True, help='dedupbam file')
parser.add_argument('--fldout',required=True, help='output FLD file')
parser.add_argument('--scriptsfolder',required=False, default='/opt', help='folder where the scripts are... used for debuging without rebuilding the docker')
EOF

set -e -x -o pipefail
python ${SCRIPTSFOLDER}/ccbr_bam2FLD.py -i $DEDUPBAM -o $FLDOUT

conda deactivate
