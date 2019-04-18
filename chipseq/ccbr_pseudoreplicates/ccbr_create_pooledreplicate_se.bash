#!/bin/bash
set -e -x -o pipefail
cpus=`nproc`

# script_path='/home/kopardev/Desktop/Dockers/chipseq/ccbr_pseudoreplicates'
script_path='/opt'

ARGPARSE_DESCRIPTION="Create pooledreplicate for SE replicates"      # this is optional
source ${script_path}/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--tagAlignRep1',required=True, help='tagAlign.gz for Rep1')
parser.add_argument('--tagAlignRep2',required=True, help='tagAlign.gz for Rep2')
EOF


tagAlignRep1Base=`echo $TAGALIGNREP1|sed "s/.tagAlign.gz//g"`
tagAlignRep2Base=`echo $TAGALIGNREP2|sed "s/.tagAlign.gz//g"`

pooledRep0=${tagAlignRep1Base}__${tagAlignRep2Base}.pooled.tagAlign.gz

zcat $TAGALIGNREP1 > pooled
zcat $TAGALIGNREP2 >> pooled
pigz -p $cpus pooled
mv pooled.gz $pooledRep0

rm -f pooled.gz 