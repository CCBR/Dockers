count=0
rm -f FLD_stats.txt
for f in `find . -name "*fld.txt"`;do
count=$((count+1))
if [ "$count" -eq "1" ]; then
head=$(grep "^#" $f|grep -v PEAKS|sed "s/# //g"|sed "s/://g"|sh ~/scripts/transpose.sh|head -n1)
echo -ne "SampleName\t$head\n" > FLD_stats.txt.header
fi
data=$(grep "^#" $f|grep -v PEAKS|sed "s/# //g"|sed "s/://g"|sh ~/scripts/transpose.sh|tail -n1)
samplename=$(basename $f|awk -F".fld" '{print $1}')
echo -ne "$samplename\t$data\n" >> FLD_stats.txt
done
cat FLD_stats.txt.header > FLD_stats.txt.tmp
rm -f FLD_stats.txt.header
sort FLD_stats.txt >> FLD_stats.txt.tmp
mv FLD_stats.txt.tmp FLD_stats.txt
