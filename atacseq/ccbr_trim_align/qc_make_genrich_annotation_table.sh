count=0
rm -f Genrich_Peak_Annotations_mqc.csv
for f in `find . -name "*.genrich.narrowPeak.annotation_distribution"|grep "_1.gen\|_2.gen\|_3.gen"`;do
count=$((count+1))
if [ "$count" -eq "1" ]; then
head=$(grep -v "^#" $f|sh ~/scripts/transpose.sh|head -n1) 
echo -ne "SampleName\t$head\n" > Genrich_Peak_Annotations_mqc.csv.header
fi
data=$(grep -v "^#" $f|sh ~/scripts/transpose.sh|tail -n1)
samplename=$(basename $f|awk -F".genrich" '{print $1}')
echo -ne "$samplename\t$data\n" >> Genrich_Peak_Annotations_mqc.csv
done
cat Genrich_Peak_Annotations_mqc.csv.header > Genrich_Peak_Annotations_mqc.csv.tmp
rm -f Genrich_Peak_Annotations_mqc.csv.header
sort Genrich_Peak_Annotations_mqc.csv >> Genrich_Peak_Annotations_mqc.csv.tmp
mv Genrich_Peak_Annotations_mqc.csv.tmp Genrich_Peak_Annotations_mqc.csv
