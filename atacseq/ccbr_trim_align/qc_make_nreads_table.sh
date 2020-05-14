count=0
rm -f Nreads_mqc.csv
for f in `find . -name "*.nreads.txt"`;do
count=$((count+1))
if [ "$count" -eq "1" ]; then
echo -ne "SampleName\tTrimmed\tMitoRibo\tUnMapped\tDuplicates\tForPeakCalling\n" > Nreads_mqc.csv.header
fi
data=$(grep -v "^#" $f|sh ~/scripts/transpose.sh|head -n1|awk -F"\t" '{printf("%d\t%d\t%d\t%d\t%d\n",$1-$2,$2-$3,$3-$4,$4-$5,$5)}')
samplename=$(basename $f|awk -F".nreads" '{print $1}')
echo -ne "$samplename\t$data\n" >> Nreads_mqc.csv
done
cat Nreads_mqc.csv.header > Nreads_mqc.csv.tmp
rm -f Nreads_mqc.csv.header
sort Nreads_mqc.csv >> Nreads_mqc.csv.tmp
mv Nreads_mqc.csv.tmp Nreads_mqc.csv
