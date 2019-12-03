#!/bin/bash
set -e -x -o pipefail
ARGPARSE_DESCRIPTION="use bedtools jaccard to calculated pairwise scores and then plot PCA plot"      # this is optional
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--inputfilelist',required=True, help='tab delimited file .. 1st column label, 2nd column filename')
parser.add_argument('--pairwise',required=True, help='output pairwise jaccard scores file')
parser.add_argument('--pcahtml',required=True, help='output PCA html filename')
parser.add_argument('--scriptsfolder',required=False, default='/opt2', help='folder where the scripts are... used for debuging without rebuilding the docker')
EOF

while read label1 grouplabel1 file1;do
  while read label2 grouplabel2 file2;do
    bedtools jaccard -a $file1 -b $file2 |tail -n1|awk -v a=$label1 -v b=$label2 '{printf("%s\t%s\t%s\n",a,b,$3)}'
  done < $INPUTFILELIST
done < $INPUTFILELIST > $PAIRWISE

if [ "$(wc -l $PAIRWISE|awk '{print $1}')" -eq "1" ];
then
	touch $PAIRWISE
	touch $PCAHTML
else
	cp ${SCRIPTSFOLDER}/pairwise.Rmd .
	Rscript -e "rmarkdown::render(\"pairwise.Rmd\",output_file=\"${PCAHTML}\", params = list(pairwise=\"${PAIRWISE}\",inputfilelist=\"${INPUTFILELIST}\"))"
fi