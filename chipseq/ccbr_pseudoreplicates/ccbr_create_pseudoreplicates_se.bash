#!/bin/bash
set -e -x -o pipefail
cpus=`nproc`

ARGPARSE_DESCRIPTION="Create pseudoreplicates for SE replicates"      # this is optional
source /opt/argparse.bash || exit 1

argparse "$@" <<EOF || exit 1
parser.add_argument('--tagAlign',required=True, help='tagAlign.gz')
EOF

tagalignpr1=`echo $TAGALIGN|sed "s/.tagAlign/.pr1.tagAlign/g"`
tagalignpr2=`echo $TAGALIGN|sed "s/.tagAlign/.pr2.tagAlign/g"`

# Get total number of read pairs
nlines=$( zcat ${TAGALIGN} | wc -l ) 
nlines=$(( (nlines + 1) / 2 ))
# Shuffle and split BED file into 2 equal parts
zcat ${TAGALIGN} | shuf | split -d -l ${nlines} - pseudoreplicate 
# Will produce ${PR_PREFIX}00 and ${PR_PREFIX}01
# Convert reads into standard tagAlign file 
pigz -p $cpus pseudoreplicate00 
mv pseudoreplicate00.gz $tagalignpr1
pigz -p $cpus pseudoreplicate01 
mv pseudoreplicate01.gz $tagalignpr2
rm -f pseudoreplicate0?