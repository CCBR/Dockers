#!/bin/bash
. /opt2/conda/etc/profile.d/conda.sh
conda activate python3

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="calculate TSS score from tagAlign file"     
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--tagaligngz',required=True, help='tagalign.gz file')
parser.add_argument('--genome',required=True, help='genome file')
parser.add_argument('--tsstxt',required=True, help='output TSS file')
parser.add_argument('--ncpus',required=False, default=2, help='Number of CPUs')
parser.add_argument('--scriptsfolder',required=False, default='/opt2', help='folder where the scripts are... used for debuging without rebuilding the docker')
EOF

set -e -x -o pipefail

# create Tn5 nicks file
nicksbed=$(echo $TAGALIGNGZ|sed "s/.tagAlign.gz/.tn5nicks.bed/g")
zcat $TAGALIGNGZ |awk -F"\t" -v OFS="\t" '{if ($6=="+") {print $1,$2,$2+1} else {print $1,$3,$3+1}}' > $nicksbed

# untar tss bed and create do file for gnu parallel
cp /${SCRIPTSFOLDER}/db/${GENOME}_tssbeds.tar.gz .
tar xzvf ${GENOME}_tssbeds.tar.gz
for f in `ls E*.bed|grep -v tn5nicks`;do echo "bedmap --echo --count $f $nicksbed|awk -F\"\t\" '{if (\$6!=\"+|0\" && \$6!=\"-|0\"){print \$4,\$6}}'|sed \"s/ /|/g\"|awk -F\"|\" -v OFS=\"\t\" '{print \$3,\$5}' > ${f}.counts";done > do_bedmap
head do_bedmap

# run bedmap, get counts and calculate density from all counts
# filter out TSS regions with <=20 nicks
parallel -j $NCPUS < do_bedmap
for f in $(ls E*bed.counts);do
	c=$(awk -F"\t" '{sum=sum+$2}END{print sum}' $f)
	if [ "$c" -gt "20" ];then
		echo $f
	fi
done > countfiles
while read a; do cat $a;done < countfiles | awk -F"\t" -v OFS="\t" '{if (NF==2){s[$1]+=$2}}END{for (var in s){print var,s[var]}}'|sort -k1,1n|awk -F"\t" -v OFS="\t" '{if ($2!=0){print $1,$2}}'|python ${SCRIPTSFOLDER}/ccbr_counts2density.py - > $TSSTXT
ntss=$(wc -l countfiles|awk '{print $1}')
echo "# TSS with 20 or more Tn5 nicking sites: $ntss" >> $TSSTXT


conda deactivate