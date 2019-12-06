#!/bin/bash
. /opt2/conda/etc/profile.d/conda.sh
conda activate python2

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="calculate width and chromosomal distributions from narrowPeak files"      # this is optional
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--narrowpeak',required=True, help='input narrowPeak file')
parser.add_argument('--widthout',required=True, help='output width dist file')
parser.add_argument('--chrout',required=True, help='output chr dist file')
parser.add_argument('--scriptsfolder',required=False, default='/opt2', help='folder where the scripts are... used for debuging without rebuilding the docker')
EOF

set -e -x -o pipefail
python ${SCRIPTSFOLDER}/ccbr_narrowPeak_distributions.py -i $NARROWPEAK -w $WIDTHOUT -c $CHROUT

conda deactivate
