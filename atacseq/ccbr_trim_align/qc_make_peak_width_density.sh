SCRIPTSFOLDER=$1
for f in `find . -name "*consensus*annotated"`;do 
echo $f
python ${SCRIPTSFOLDER}/annotated2peakwidthdensity.py $f > ${f}.peak_width_density
gzip -n $f
done
