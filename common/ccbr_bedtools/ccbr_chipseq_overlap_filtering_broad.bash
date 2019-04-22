#!/bin/bash

get_max_min_ratio() {
min=$(($1>$2?$2:$1))
max=$(($1>$2?$1:$2))
r=`echo $max $min|awk '{printf("%.2f",$1/$2)}'`
echo "$r"
}


set -e -x -o pipefail
script_folder="/opt"
# script_folder="/home/kopardev/Desktop/Dockers/chipseq/ccbr_macs_peakcalling"
ARGPARSE_DESCRIPTION="Get Optimal/Conservative broadPeak calls using naive overlaps"      # this is optional
source ${script_folder}/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--PeaksRep1',required=True, help='broadPeak file')
parser.add_argument('--PeaksRep1Pr1',required=True, help='broadPeak file')
parser.add_argument('--PeaksRep1Pr2',required=True, help='broadPeak file')
parser.add_argument('--PeaksRep2',required=True, help='broadPeak file')
parser.add_argument('--PeaksRep2Pr1',required=True, help='broadPeak file')
parser.add_argument('--PeaksRep2Pr2',required=True, help='broadPeak file')
parser.add_argument('--PeaksRep0',required=True, help='pooled broadPeak file')
parser.add_argument('--PeaksRep0Pr1',required=True, help='pooled pr1 broadPeak file')
parser.add_argument('--PeaksRep0Pr2',required=True, help='pooled pr2 broadPeak file')
EOF

# ====================== # For broadPeak files
# ======================
# Find pooled peaks that overlap Rep1 and Rep2 where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5
intersectBed -wo -a $PEAKSREP0 -b $PEAKSREP1 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}' | cut -f 1-9 | sort | uniq | \
intersectBed -wo -a stdin -b $PEAKSREP2 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}' | cut -f 1-9 | sort | uniq > Nt.broadPeak
# Find pooled peaks that overlap PooledPseudoRep1 and PooledPseudoRep2 where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5
intersectBed -wo -a $PEAKSREP0 -b $PEAKSREP0PR1 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}' | cut -f 1-9 | sort | uniq | \
intersectBed -wo -a stdin -b $PEAKSREP0PR2 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}' | cut -f 1-9 | sort | uniq > Np.broadPeak

intersectBed -wo -a $PEAKSREP1 -b $PEAKSREP1PR1 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}' | cut -f 1-9 | sort | uniq | \
intersectBed -wo -a stdin -b $PEAKSREP1PR2 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}' | cut -f 1-9 | sort | uniq > N1.broadPeak

intersectBed -wo -a $PEAKSREP2 -b $PEAKSREP2PR1 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}' | cut -f 1-9 | sort | uniq | \
intersectBed -wo -a stdin -b $PEAKSREP2PR2 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}' | cut -f 1-9 | sort | uniq > N2.broadPeak

Nt=`cat Nt.broadPeak|wc -l`
Np=`cat Np.broadPeak|wc -l`
N1=`cat N1.broadPeak|wc -l`
N2=`cat N2.broadPeak|wc -l`

selfConsistencyRatio=`get_max_min_ratio $N1 $N2`
rescueRatio=`get_max_min_ratio $Nt $Np`

echo "Nt = $Nt" > results.broadPeak.txt
echo "Np = $Np" >> results.broadPeak.txt
echo "N1 = $N1" >> results.broadPeak.txt
echo "N2 = $N2" >> results.broadPeak.txt
echo "selfConsistencyRatio = $selfConsistencyRatio" >> results.broadPeak.txt
echo "rescueRatio = $rescueRatio" >> results.broadPeak.txt

if [ $Np -gt $Nt ]
then
	cp Np.broadPeak Optimal.broadPeak
else
	cp Nt.broadPeak Optimal.broadPeak
fi

cp Nt.broadPeak Conservative.broadPeak

awk -F"\t" -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3,$5,$6}' Optimal.broadPeak  | sort -k1,1 -k2,2n  > Optimal_broad_peaks.bed
awk -F"\t" -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3,$5,$6}' Conservative.broadPeak  | sort -k1,1 -k2,2n > Conservative_broad_peaks.bed