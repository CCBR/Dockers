#!/bin/bash
# setenv SINGULARITYENV_PYTHONNOUSERSITE "1"

set -e -x -o pipefail
. /opt2/conda/etc/profile.d/conda.sh
conda activate gimme
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

ARGPARSE_DESCRIPTION="use gimmemotifs to perform motif analysis a) scan and b) enrichment"      # this is optional
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--narrowpeak',required=True, help='narrowPeak input file')
parser.add_argument('--genome',required=True, help='hg19/hg38/mm9/mm10')
parser.add_argument('--outfolder',required=True, help='path to output folder')
parser.add_argument('--threads',required=False, default=2, type=int, help='number of threads')
EOF

ncpus=$THREADS
samplename=$(echo $NARROWPEAK|awk -F"/" '{print $NF}'|awk -F"." '{print $1}')
genomefa="/index/${GENOME}.fa"
outfolder_name=$(echo $OUTFOLDER|awk -F"/" '{print $NF}')
if [ ! -f /data2/$outfolder_name ];then mkdir -p /data2/$outfolder_name ;fi

if [ "$GENOME" == "hg19" ]; then motifdb="HOCOMOCOv11_HUMAN";fi
if [ "$GENOME" == "hg38" ]; then motifdb="HOCOMOCOv11_HUMAN";fi
if [ "$GENOME" == "mm9" ]; then motifdb="HOCOMOCOv11_MOUSE";fi
if [ "$GENOME" == "mm10" ]; then motifdb="HOCOMOCOv11_MOUSE";fi

# get peak width and number
wn=$(cat $NARROWPEAK|awk -F"\t" '{count=count+1;totalwidth=totalwidth+$3-$2}END{printf("%d\t%d\n",count,totalwidth/count)}')
npeaks=$(echo $wn|awk '{print $1}')
peakwidth=$(echo $wn|awk '{print $2}')
if [ "$npeaks" -lt "10000" ]; then npeaks=10000;fi


# create the background fasta file
head -n2 $NARROWPEAK
head -n2 $genomefa
cut -f1-3 $NARROWPEAK > /data2/${samplename}.bed
head -n2 /data2/${samplename}.bed

gimme background -i /data2/${samplename}.bed -s $peakwidth -n $npeaks -g $genomefa /data2/bg.fasta gc
bedtools getfasta -fi $genomefa -bed $NARROWPEAK -fo /data2/${samplename}.fa

# find hocomoco enriched motifs
# gimme motifs $NARROWPEAK /data2/$outfolder_name/enrichment --known -g $genomefa -p $motifdb -b /data2/bg.fasta -s 0 -f 0.5 -N $ncpus
gimme motifs /data2/${samplename}.fa /data2/${outfolder_name}/enrichment --known -g $genomefa -p $motifdb -b /data2/bg.fasta -s 0 -N $ncpus --noreport

# scan for motifs in all peaks
# if [ ! -f /data2/${outfolder_name}/scan ];then mkdir -p /data2/${outfolder_name}/scan;fi
# gimme scan /data2/${samplename}.fa -p $motifdb -f 0.01 -B /data2/bg.fasta -n 3 -b -z -N $ncpus > /data2/${outfolder_name}/scan/matches.bed
# gzip /data2/${outfolder_name}/scan/matches.bed

rm -f /data2/bg.fasta
rm -f /data2/${samplename}.fa
rm -f /data2/${samplename}.bed

conda deactivate
