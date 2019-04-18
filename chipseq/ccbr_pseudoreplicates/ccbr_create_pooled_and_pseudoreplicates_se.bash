#!/bin/bash
set -e -x -o pipefail
cpus=`nproc`

# script_path='/home/kopardev/Desktop/Dockers/chipseq/ccbr_pseudoreplicates'
script_path='/opt'

ARGPARSE_DESCRIPTION="Create pooled and pseudoreplicates for SE replicates"      # this is optional
source ${script_path}/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--tagAlignRep1',required=True, help='tagAlign.gz for Rep1')
parser.add_argument('--tagAlignRep2',required=True, help='tagAlign.gz for Rep2')
EOF

tagalignRep1pr1=`echo $TAGALIGNREP1|sed "s/.tagAlign/.pr1.tagAlign/g"`
tagalignRep1pr2=`echo $TAGALIGNREP1|sed "s/.tagAlign/.pr2.tagAlign/g"`

tagalignRep2pr1=`echo $TAGALIGNREP2|sed "s/.tagAlign/.pr1.tagAlign/g"`
tagalignRep2pr2=`echo $TAGALIGNREP2|sed "s/.tagAlign/.pr2.tagAlign/g"`

bash ${script_path}/ccbr_create_pseudoreplicates_se.bash --tagAlign $TAGALIGNREP1
bash ${script_path}/ccbr_create_pseudoreplicates_se.bash --tagAlign $TAGALIGNREP2

tagAlignRep1Base=`echo $TAGALIGNREP1|sed "s/.tagAlign.gz//g"`
tagAlignRep2Base=`echo $TAGALIGNREP2|sed "s/.tagAlign.gz//g"`

pooledRep0=${tagAlignRep1Base}__${tagAlignRep2Base}.pooled.tagAlign.gz
pooledRep0pr1=${tagAlignRep1Base}__${tagAlignRep2Base}.pooled.pr1.tagAlign.gz
pooledRep0pr2=${tagAlignRep1Base}__${tagAlignRep2Base}.pooled.pr2.tagAlign.gz

zcat $TAGALIGNREP1 > pooled
zcat $TAGALIGNREP2 >> pooled
pigz -p $cpus pooled
mv pooled.gz $pooledRep0

zcat $tagalignRep1pr1 > ppr1
zcat $tagalignRep2pr1 >> ppr1
pigz -p $cpus ppr1
mv ppr1.gz $pooledRep0pr1

zcat $tagalignRep1pr2 > ppr2
zcat $tagalignRep2pr2 >> ppr2
pigz -p $cpus ppr2
mv ppr2.gz $pooledRep0pr2

rm -f pooled.gz ppr1.gz ppr2.gz