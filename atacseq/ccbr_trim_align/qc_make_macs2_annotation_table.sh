count=0
rm -f MACS2_Peak_Annotations_mqc.csv
for f in `find . -name "*.macs2.narrowPeak.annotation_distribution"|grep "_1.macs\|_2.macs\|_3.macs"`;do
count=$((count+1))
if [ "$count" -eq "1" ]; then
head=$(grep -v "^#" $f|sh ~/scripts/transpose.sh|head -n1) 
echo -ne "SampleName\t$head\n" > MACS2_Peak_Annotations_mqc.csv.header
fi
data=$(grep -v "^#" $f|sh ~/scripts/transpose.sh|tail -n1)
samplename=$(basename $f|awk -F".macs2" '{print $1}')
echo -ne "$samplename\t$data\n" >> MACS2_Peak_Annotations_mqc.csv
done
cat MACS2_Peak_Annotations_mqc.csv.header > MACS2_Peak_Annotations_mqc.csv.tmp
rm -f MACS2_Peak_Annotations_mqc.csv.header
sort MACS2_Peak_Annotations_mqc.csv >> MACS2_Peak_Annotations_mqc.csv.tmp
mv MACS2_Peak_Annotations_mqc.csv.tmp MACS2_Peak_Annotations_mqc.csv
