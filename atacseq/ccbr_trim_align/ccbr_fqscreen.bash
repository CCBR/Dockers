#!/bin/bash

. /opt2/conda/etc/profile.d/conda.sh
conda activate python3

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="run fastqscreen for contamination check"      # this is optional
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--infastq',required=True, help='fastq.gz file(s). Multiple fastqs should be space separated',action='append', nargs='+')
parser.add_argument('--confurl',required=True, help='fastqscreen conf file url eg. https://hpc.nih.gov/~CCBR/sbg_reference_bundle/fastq_screen.conf')
parser.add_argument('--dblisturl',required=True, help='url of file that contains the urls of databases (tarballs) to download eg. https://hpc.nih.gov/~CCBR/sbg_reference_bundle/fastq_screen_databases_list.txt')
parser.add_argument('--threads',required=True, help='number of threads')
EOF
# get rid of the python list box brackets, commas and single quotes
INFASTQ=$(echo $INFASTQ|sed "s/\[//g"|sed "s/\]//g"|sed "s/,//g"|sed "s/'//g")
#get the dblist file ... download all tarballs .. .untar them
curl -O $DBLISTURL
dblistfile=$(echo $DBLISTURL|awk -F"/" '{print $NF}')
while read a;do
b=$(echo $a|awk -F"/" '{print $NF}')
curl -O $a
tar xzvf $b
done < $dblistfile
### get conf file
conffile=$(echo $CONFURL|awk -F"/" '{print $NF}')
curl -O $CONFURL
### run fastq_screen
fastq_screen --threads $THREADS --subset 1000000 --conf fastq_screen.conf $INFASTQ