#!/bin/bash
. /opt2/conda/etc/profile.d/conda.sh
conda activate python3

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="call atac-seq peaks using macs2 ... two replicate wrapper script"
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--tagalign1',required=True, help='input tagAlign.gz file for replicate 1')
parser.add_argument('--tagalign2',required=True, help='input tagAlign.gz file for replicate 2')
parser.add_argument('--extsize',required=False, default=200, help='extsize')
parser.add_argument('--shiftsize',required=False, default=100, help='shiftsize')
parser.add_argument('--rep1name',required=True, help='samplename for replicate 1')
parser.add_argument('--rep2name',required=True, help='samplename for replicate 2')
parser.add_argument('--samplename',required=True, help='samplename')
parser.add_argument('--genomefilerep1',required=True, help='dedupbam based genome file for replicate 1 required by bedGraphToBigWig')
parser.add_argument('--genomefilerep2',required=True, help='dedupbam based genome file for replicate 2 required by bedGraphToBigWig')
parser.add_argument('--genomename',required=True, help='hg19/hg38/mm9/mm10')
# only required if filtering
parser.add_argument('--filterpeaks',required=False, default="True", help='filterpeaks by qvalue: True or False')
parser.add_argument('--qfilter',required=False, default=0.693147, help='default qfiltering value is 0.693147 (-log10 of 0.5) for q=0.5')
parser.add_argument('--scriptsfolder',required=False, default='/opt2', help='folder where the scripts are... used for debuging without rebuilding the docker')
EOF

extsize=$EXTSIZE
shiftsize=$SHIFTSIZE
genome=$GENOMENAME
samplename=$SAMPLENAME

rep1tagAlign=$TAGALIGN1
rep1name=$REP1NAME
prefix1=${rep1name}.macs2
genomefilerep1=$GENOMEFILEREP1

rep2tagAlign=$TAGALIGN2
rep2name=$REP2NAME
prefix2=${rep2name}.macs2
genomefilerep2=$GENOMEFILEREP2

rep1pr=$(echo $rep1tagAlign|awk -F"/" '{print $NF}'|awk -F".tagAlign" '{print $1".pr.tagAlign.gz"}')
rep2pr=$(echo $rep2tagAlign|awk -F"/" '{print $NF}'|awk -F".tagAlign" '{print $1".pr.tagAlign.gz"}')
pooled=$(echo ${rep1pr}_${rep2pr}|sed "s/.pr.tagAlign.gz//g"|awk '{print $1".pooled.tagAlign.gz"}')
pooledpr=$(echo $pooled|sed "s/pooled/pooled.pr/")

zcat $rep1tagAlign |awk -v v=50 'BEGIN {srand()} !/^$/ { if (rand() <= v/100) print $0}'|gzip -c - > $rep1pr
zcat $rep2tagAlign |awk -v v=50 'BEGIN {srand()} !/^$/ { if (rand() <= v/100) print $0}'|gzip -c - > $rep2pr
zcat $rep1tagAlign $rep2tagAlign > ${pooled}_tmp.tagAlign
bedSort ${pooled}_tmp.tagAlign ${pooled}_tmp.tagAlign
gzip ${pooled}_tmp.tagAlign
mv ${pooled}_tmp.tagAlign.gz $pooled
zcat $pooled |awk -v v=50 'BEGIN {srand()} !/^$/ { if (rand() <= v/100) print $0}'|gzip -c - > $pooledpr

bash ${SCRIPTSFOLDER}/ccbr_macs2_peak_calling.bash --tagalign $rep1tagAlign --samplename $rep1name --genomename $genome --filterpeaks False --savebigwig True --genomefile $genomefilerep1 --savetn5knicksbed True
bash ${SCRIPTSFOLDER}/ccbr_macs2_peak_calling.bash --tagalign $rep2tagAlign --samplename $rep2name --genomename $genome --filterpeaks False --savebigwig True --genomefile $genomefilerep2 --savetn5knicksbed True
bash ${SCRIPTSFOLDER}/ccbr_macs2_peak_calling.bash --tagalign $rep1pr --samplename ${rep1name}.pr --genomename $genome --filterpeaks False --savebigwig False
bash ${SCRIPTSFOLDER}/ccbr_macs2_peak_calling.bash --tagalign $rep2pr --samplename ${rep2name}.pr --genomename $genome --filterpeaks False --savebigwig False
bash ${SCRIPTSFOLDER}/ccbr_macs2_peak_calling.bash --tagalign $pooled --samplename ${samplename}.pooled --genomename $genome --filterpeaks False --savebigwig True --savetn5knicksbed True
bash ${SCRIPTSFOLDER}/ccbr_macs2_peak_calling.bash --tagalign $pooledpr --samplename ${samplename}.pooled.pr --genomename $genome --filterpeaks False --savebigwig False

rep1_peaks="${rep1name}.macs2.narrowPeak"
rep2_peaks="${rep2name}.macs2.narrowPeak"
rep1pr_peaks="${rep1name}.pr.macs2.narrowPeak"
rep2pr_peaks="${rep2name}.pr.macs2.narrowPeak"
pooled_peaks="${samplename}.pooled.macs2.narrowPeak"
pooledpr_peaks="${samplename}.pooled.pr.macs2.narrowPeak"

idr --samples $rep1_peaks $rep2_peaks --peak-list $pooled_peaks --output-file ${samplename}.idr.narrowPeak --rank p.value --soft-idr-threshold 0.1 --plot --use-best-multisummit-IDR --log-output-file ${samplename}.idr.log
# remove duplicate peak calls
mv ${samplename}.idr.narrowPeak ${samplename}.idr.narrowPeak.tmp
sort -k9,9gr ${samplename}.idr.narrowPeak.tmp|awk -F"\t" '!NF || !seen[$1":"$2"-"$3]++'|sort -k1,1 -k2,2n > ${samplename}.idr.narrowPeak
awk -F"\t" '{if ($5>=540) {print}}' ${samplename}.idr.narrowPeak > ${samplename}.idr.filt.narrowPeak

Ntr1=$(wc -l $rep1_peaks|awk '{print $1}')
Ntr2=$(wc -l $rep2_peaks|awk '{print $1}')

Nt=$(( $Ntr1 > $Ntr2 ? $Ntr1 : $Ntr2 ))

Np=$(wc -l $pooledpr_peaks|awk '{print $1}')
num=$(( $Nt > $Np ? $Nt : $Np ))
den=$(( $Nt < $Np ? $Nt : $Np ))

rescue_ratio=$(echo "scale=3;$num/$den"|bc)
N1=$(wc -l $rep1pr_peaks|awk '{print $1}')
N2=$(wc -l $rep2pr_peaks|awk '{print $1}')
num=$(( $N1 > $N2 ? $N1 : $N2 ))
den=$(( $N1 < $N2 ? $N1 : $N2 ))

self_consistency_ratio=$(echo "scale=3;$num/$den"|bc)

if [ "$Ntr1" -gt "$Ntr2" ]
then
conservative=$rep1_peaks
else
conservative=$rep2_peaks
fi

if [ "$Nt" -gt "$Np" ]
then
optimal=$conservative
else
optimal=$pooledpr_peaks
fi
if [ -f ${samplename}.conservative.narrowPeak ];then rm -f ${samplename}.conservative.narrowPeak; fi
if [ -f ${samplename}.optimal.narrowPeak ];then rm -f ${samplename}.optimal.narrowPeak; fi
cp $conservative ${samplename}.conservative.narrowPeak
cp $optimal ${samplename}.optimal.narrowPeak

echo -ne "$samplename\t" > ${samplename}.peakstats.txt
echo -ne "$rescue_ratio\t" >> ${samplename}.peakstats.txt
echo -ne "$self_consistency_ratio\t" >> ${samplename}.peakstats.txt
echo -ne "$Ntr1\t" >> ${samplename}.peakstats.txt
echo -ne "$Ntr2\t" >> ${samplename}.peakstats.txt

a=$(wc -l ${samplename}.conservative.narrowPeak|awk '{print $1}')
echo -ne "$a\t" >> ${samplename}.peakstats.txt
a=$(wc -l ${samplename}.optimal.narrowPeak|awk '{print $1}')
echo -ne "$a\t" >> ${samplename}.peakstats.txt
a=$(wc -l ${samplename}.idr.narrowPeak|awk '{print $1}')
echo -ne "$a\n" >> ${samplename}.peakstats.txt

conda deactivate

if [ $FILTERPEAKS == "True" ];then
  qvalue=$QFILTER
  for f in $rep1_peaks $rep2_peaks ${samplename}.idr.filt.narrowPeak ${samplename}.idr.narrowPeak;do
    filteredpeakfile=$(echo $f|sed "s/.narrowPeak/.qfilter.narrowPeak/g")
    awk -F"\t" -v q=$qvalue "{if (\$9>q){print}}" $f > $filteredpeakfile
    Rscript ${SCRIPTSFOLDER}/ccbr_annotate_peaks.R -n $f -a ${f}.annotated -g $GENOME -l ${f}.genelist -f ${f}.annotation_summary 
    Rscript ${SCRIPTSFOLDER}/ccbr_annotate_peaks.R -n $filteredpeakfile -a ${filteredpeakfile}.annotated -g $GENOME -l ${filteredpeakfile}.genelist -f ${filteredpeakfile}.annotation_summary
    cut -f1,2 ${f}.annotation_summary > ${f}.annotation_distribution
    cut -f1,2 ${filteredpeakfile}.annotation_summary > ${filteredpeakfile}.annotation_distribution
  done
fi