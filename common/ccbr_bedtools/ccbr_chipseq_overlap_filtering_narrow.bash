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
ARGPARSE_DESCRIPTION="Get Optimal/Conservative peak calls using naive overlaps"      # this is optional
source ${script_folder}/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--PeaksRep1',required=True, help='narrowPeak file')
parser.add_argument('--PeaksRep1Pr1',required=True, help='narrowPeak file')
parser.add_argument('--PeaksRep1Pr2',required=True, help='narrowPeak file')
parser.add_argument('--PeaksRep2',required=True, help='narrowPeak file')
parser.add_argument('--PeaksRep2Pr1',required=True, help='narrowPeak file')
parser.add_argument('--PeaksRep2Pr2',required=True, help='narrowPeak file')
parser.add_argument('--PeaksRep0',required=True, help='pooled narrowPeak file')
parser.add_argument('--PeaksRep0Pr1',required=True, help='pooled pr1 narrowPeak file')
parser.add_argument('--PeaksRep0Pr2',required=True, help='pooled pr2 narrowPeak file')
EOF

# ====================== # For narrowPeak files
# ======================
# Find pooled peaks that overlap Rep1 and Rep2 where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5
intersectBed -wo -a $PEAKSREP0 -b $PEAKSREP1 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | \
intersectBed -wo -a stdin -b $PEAKSREP2 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq > Nt.narrowPeak
# Find pooled peaks that overlap PooledPseudoRep1 and PooledPseudoRep2 where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5
intersectBed -wo -a $PEAKSREP0 -b $PEAKSREP0PR1 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | \
intersectBed -wo -a stdin -b $PEAKSREP0PR2 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq > Np.narrowPeak

intersectBed -wo -a $PEAKSREP1 -b $PEAKSREP1PR1 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | \
intersectBed -wo -a stdin -b $PEAKSREP1PR2 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq > N1.narrowPeak

intersectBed -wo -a $PEAKSREP2 -b $PEAKSREP2PR1 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | \
intersectBed -wo -a stdin -b $PEAKSREP2PR2 | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq > N2.narrowPeak

Nt=`cat Nt.narrowPeak|wc -l`
Np=`cat Np.narrowPeak|wc -l`
N1=`cat N1.narrowPeak|wc -l`
N2=`cat N2.narrowPeak|wc -l`

selfConsistencyRatio=`get_max_min_ratio $N1 $N2`
rescueRatio=`get_max_min_ratio $Nt $Np`

echo "Nt = $Nt" > results.narrowPeak.txt
echo "Np = $Np" >> results.narrowPeak.txt
echo "N1 = $N1" >> results.narrowPeak.txt
echo "N2 = $N2" >> results.narrowPeak.txt
echo "selfConsistencyRatio = $selfConsistencyRatio" >> results.narrowPeak.txt
echo "rescueRatio = $rescueRatio" >> results.narrowPeak.txt


if [ $Np -gt $Nt ]
then
	cp Np.narrowPeak Optimal.narrowPeak
else
	cp Nt.narrowPeak Optimal.narrowPeak
fi

cp Nt.narrowPeak Conservative.narrowPeak

awk -F"\t" -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3,$5,$6}' Optimal.narrowPeak | sort -k1,1 -k2,2n > Optimal_narrow_peaks.bed
awk -F"\t" -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3,$5,$6}' Conservative.narrowPeak | sort -k1,1 -k2,2n > Conservative_narrow_peaks.bed

# cat Nt.narrowPeak Np.narrowPeak | sort | uniq | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' | grep -P 'chr[0-9XY]+(?!_)' > FinalPeakList.narrowPeak