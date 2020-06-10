#!/bin/bash
# setenv SINGULARITYENV_PYTHONNOUSERSITE "1"

set -e -x -o pipefail
. /opt2/conda/etc/profile.d/conda.sh
conda activate python3
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

ARGPARSE_DESCRIPTION="Create all files required to generate a MultiQC report"
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--qcfolder',required=True, help='path to the qc folder')
parser.add_argument('--scriptsfolder',required=False, default='/opt2', help='folder where the scripts are... used for debuging without rebuilding the docker')
EOF

cd $QCFOLDER

#####
#make_nreads_table
#OUTPUT:Nreads_mqc.csv
#####
count=0
rm -f Nreads_mqc.csv
for f in `find . -name "*.nreads.txt"`;do
count=$((count+1))
if [ "$count" -eq "1" ]; then
echo -ne "SampleName\tTrimmed\tMitoRibo\tUnMapped\tDuplicates\tForPeakCalling\n" > Nreads_mqc.csv.header
fi
data=$(grep -v "^#" $f|sh ${SCRIPTSFOLDER}/transpose.sh|head -n1|awk -F"\t" '{printf("%d\t%d\t%d\t%d\t%d\n",$1-$2,$2-$3,$3-$4,$4-$5,$5)}')
samplename=$(basename $f|awk -F".nreads" '{print $1}')
echo -ne "$samplename\t$data\n" >> Nreads_mqc.csv
done
cat Nreads_mqc.csv.header > Nreads_mqc.csv.tmp
rm -f Nreads_mqc.csv.header
sort Nreads_mqc.csv >> Nreads_mqc.csv.tmp
mv Nreads_mqc.csv.tmp Nreads_mqc.csv
count=0
rm -f NRF_stats.txt

#####
#make_nrf_table.sh
#OUTPUT:NRF_stats.txt
#####
for f in `find . -name "*.nrf"`;do
count=$((count+1))
if [ "$count" -eq "1" ]; then
head=$(cat $f|sed "s/:/\t/g"|sh ${SCRIPTSFOLDER}/transpose.sh|head -n1)
echo -ne "SampleName\t$head\n" > NRF_stats.txt.header
fi
data=$(cat $f|sed "s/:/\t/g"|sh ${SCRIPTSFOLDER}/transpose.sh|tail -n1)
samplename=$(basename $f|awk -F".nrf" '{print $1}')
echo -ne "$samplename\t$data\n" >> NRF_stats.txt
done
cat NRF_stats.txt.header > NRF_stats.txt.tmp
rm -f NRF_stats.txt.header
sort NRF_stats.txt >> NRF_stats.txt.tmp
mv NRF_stats.txt.tmp NRF_stats.txt


#####
#annotated2peakwidthdensity.sh
#OUTPUT:.peak_width_density files in peak_annotation subfolder
#####
for f in `find . -name "*consensus*annotated"`;do 
echo $f
python ${SCRIPTSFOLDER}/qc_annotated2peakwidthdensity.py $f > ${f}.peak_width_density
gzip -n $f
done


#####
#make_tss_scatter_table.sh
#OUTPUT:data.tss_knicking_sites.txt
#####
rm -f data.tss_knicking_sites.txt
for f in `find . -name "*.tss.txt"`;do
samplename=$(basename $f|awk -F".tss" '{print $1}')
enrich=$(grep "enrichment" $f|awk '{print $NF}')
nsites=$(grep "20 or more Tn5" $f|awk '{print $NF}')
echo -ne "$samplename\t$nsites\t$enrich\n" >> data.tss_knicking_sites.txt
done


#####
#make_frip_table.sh
#OUTPUT:FRiP_stats.txt
#####
find . -name "*.frip" -not -name "*init*" -exec cat {} \; | python ${SCRIPTSFOLDER}/qc_get_frip_stats_table.py


#####
#make_fld_table.sh
#OUTPUT:FLD_stats.txt
#####
count=0
rm -f FLD_stats.txt
for f in `find . -name "*fld.txt"`;do
count=$((count+1))
if [ "$count" -eq "1" ]; then
head=$(grep "^#" $f|grep -v PEAKS|sed "s/# //g"|sed "s/://g"|sh ${SCRIPTSFOLDER}/transpose.sh|head -n1)
echo -ne "SampleName\t$head\n" > FLD_stats.txt.header
fi
data=$(grep "^#" $f|grep -v PEAKS|sed "s/# //g"|sed "s/://g"|sh ${SCRIPTSFOLDER}/transpose.sh|tail -n1)
samplename=$(basename $f|awk -F".fld" '{print $1}')
echo -ne "$samplename\t$data\n" >> FLD_stats.txt
done
cat FLD_stats.txt.header > FLD_stats.txt.tmp
rm -f FLD_stats.txt.header
sort FLD_stats.txt >> FLD_stats.txt.tmp
mv FLD_stats.txt.tmp FLD_stats.txt


#####
#make_macs2_annotation_table.sh
#OUTPUT:MACS2_Peak_Annotations_mqc.csv
#####
count=0
rm -f MACS2_Peak_Annotations_mqc.csv
for f in `find . -name "*.macs2.narrowPeak.annotation_distribution"|grep "_1.macs\|_2.macs\|_3.macs"`;do
count=$((count+1))
if [ "$count" -eq "1" ]; then
head=$(grep -v "^#" $f|sh ${SCRIPTSFOLDER}/transpose.sh|head -n1)
echo -ne "SampleName\t$head\n" > MACS2_Peak_Annotations_mqc.csv.header
fi
data=$(grep -v "^#" $f|sh ${SCRIPTSFOLDER}/transpose.sh|tail -n1)
samplename=$(basename $f|awk -F".macs2" '{print $1}')
echo -ne "$samplename\t$data\n" >> MACS2_Peak_Annotations_mqc.csv
done
cat MACS2_Peak_Annotations_mqc.csv.header > MACS2_Peak_Annotations_mqc.csv.tmp
rm -f MACS2_Peak_Annotations_mqc.csv.header
sort MACS2_Peak_Annotations_mqc.csv >> MACS2_Peak_Annotations_mqc.csv.tmp
mv MACS2_Peak_Annotations_mqc.csv.tmp MACS2_Peak_Annotations_mqc.csv


#####
#make_genrich_annotation_table.sh
#OUTPUT:Genrich_Peak_Annotations_mqc.csv
#####
count=0
rm -f Genrich_Peak_Annotations_mqc.csv
for f in `find . -name "*.genrich.narrowPeak.annotation_distribution"|grep "_1.gen\|_2.gen\|_3.gen"`;do
count=$((count+1))
if [ "$count" -eq "1" ]; then
head=$(grep -v "^#" $f|sh ${SCRIPTSFOLDER}/transpose.sh|head -n1)
echo -ne "SampleName\t$head\n" > Genrich_Peak_Annotations_mqc.csv.header
fi
data=$(grep -v "^#" $f|sh ${SCRIPTSFOLDER}/transpose.sh|tail -n1)
samplename=$(basename $f|awk -F".genrich" '{print $1}')
echo -ne "$samplename\t$data\n" >> Genrich_Peak_Annotations_mqc.csv
done
cat Genrich_Peak_Annotations_mqc.csv.header > Genrich_Peak_Annotations_mqc.csv.tmp
rm -f Genrich_Peak_Annotations_mqc.csv.header
sort Genrich_Peak_Annotations_mqc.csv >> Genrich_Peak_Annotations_mqc.csv.tmp
mv Genrich_Peak_Annotations_mqc.csv.tmp Genrich_Peak_Annotations_mqc.csv


#####
#Run multiqc
#OUTPUT:multiqc_report.html
#####
multiqc -c ${SCRIPTSFOLDER}/ccbr_atac_config.yaml -f --interactive -v .
