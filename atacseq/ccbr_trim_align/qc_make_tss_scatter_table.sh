rm -f data.tss_knicking_sites.txt
for f in `find . -name "*.tss.txt"`;do
samplename=$(basename $f|awk -F".tss" '{print $1}')
enrich=$(grep "enrichment" $f|awk '{print $NF}')
nsites=$(grep "20 or more Tn5" $f|awk '{print $NF}')
echo -ne "$samplename\t$nsites\t$enrich\n" >> data.tss_knicking_sites.txt
done
