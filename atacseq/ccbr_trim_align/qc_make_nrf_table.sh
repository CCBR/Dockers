count=0
rm -f NRF_stats.txt
for f in `find . -name "*.nrf"`;do
count=$((count+1))
if [ "$count" -eq "1" ]; then
head=$(cat $f|sed "s/:/\t/g"|sh ~/scripts/transpose.sh|head -n1)
echo -ne "SampleName\t$head\n" > NRF_stats.txt.header
fi
data=$(cat $f|sed "s/:/\t/g"|sh ~/scripts/transpose.sh|tail -n1)
samplename=$(basename $f|awk -F".nrf" '{print $1}')
echo -ne "$samplename\t$data\n" >> NRF_stats.txt
done
cat NRF_stats.txt.header > NRF_stats.txt.tmp
rm -f NRF_stats.txt.header
sort NRF_stats.txt >> NRF_stats.txt.tmp
mv NRF_stats.txt.tmp NRF_stats.txt
