#!/bin/bash

set -e -x -o pipefail
. /opt2/conda/etc/profile.d/conda.sh
conda activate python3
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

ARGPARSE_DESCRIPTION="use homer(findMotifsGenome.pl) and meme(ame) to perform motif enrichment analysis using HOCOMOCO_v11 motifs"      # this is optional
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--narrowpeak',required=True, help='narrowPeak input file')
parser.add_argument('--genome',required=True, help='hg19/hg38/mm9/mm10')
parser.add_argument('--ntop',required=False, default=50000, type=int, help='use the top x peaks ... default 50k')
parser.add_argument('--threads',required=False, default=2, type=int, help='number of threads')
EOF

ncpus=$THREADS
samplename=$(echo $NARROWPEAK|awk -F"/" '{print $NF}'|awk -F"." '{print $1}')
genomefa="/index/${GENOME}.fa"
outfolder_name=$(echo $OUTFOLDER|awk -F"/" '{print $NF}')
workdir=$(pwd)

if [ "$GENOME" == "hg19" ]; then species="HUMAN";fi
if [ "$GENOME" == "hg38" ]; then species="HUMAN";fi
if [ "$GENOME" == "mm9" ]; then species="MOUSE";fi
if [ "$GENOME" == "mm10" ]; then species="MOUSE";fi

homermotif="/opt2/db/HOCOMOCOv11_full_${species}_mono_homer_format_0.001.motif"
mememotiftar="/opt2/db/HOCOMOCOv11_full_${species}_mono_meme_format.tar.gz"

cp $homermotif $workdir
cp $mememotiftar $workdir
tar xzvf HOCOMOCOv11_full_${species}_mono_meme_format.tar.gz

sort -k9,9gr $NARROWPEAK|cut -f1-3 |awk -v n=$NTOP '{if (NR<=n) {print}}' > ${samplename}.input.bed
bedSort ${samplename}.input.bed ${samplename}.input.bed

findMotifsGenome.pl ${samplename}.input.bed $genomefa motif_enrichment \
-nomotif \
-size given \
-mknown HOCOMOCOv11_full_${species}_mono_homer_format_0.001.motif \
-N $(wc -l ${samplename}.input.bed|awk '{print $1*4}') \
-h \
-p $ncpus \
-dumpFasta \
-cpg \
-maxN 0.1 \
-len 10 \
-preparsedDir preparsed

conda deactivate

# deactivate conda to use ame version 5.1.0 else it will use 5.0.5 in python3 conda environment which does not work

# remove -nomotif from the above command if you want to look for de novo motifs as well..... this will take a long time.... as it is a serial homer2 job

targetfa="motif_enrichment/target.fa"
backgroundfa="motif_enrichment/background.fa"

ls *.meme |sort > memes
while read a;do echo "ame --o ${a}_ame_out --noseq --control $backgroundfa --seed 12345 --verbose 1 $targetfa $a";done < memes > do_memes
parallel -j $ncpus < do_memes
while read a;do cat ${a}_ame_out/ame.tsv;done < memes |grep "^rank"|sort|uniq > ame_results.txt
while read a;do cat ${a}_ame_out/ame.tsv;done < memes |grep -v "^#"|grep -v "^$"|grep -v "^rank" |sort -k6,6g|awk -F"\t" '{printf("%d\t%s\n",NR,$_)}'|cut -f1,3- >> ame_results.txt
if [ "$(grep CTCF_ ame_results.txt|wc -l)" -ne "0" ]; then
	ctcf_enrich=$(grep -m1 CTCF_ ame_results.txt |awk -F"\t" '{print $(NF-2)/$NF}')
else
	ctcf_enrich="NA"
fi
echo -ne "# CTCF_enrichment = $ctcf_enrich\n" >> ame_results.txt

mv ame_results.txt motif_enrichment
rm -f $targetfa $backgroundfa
