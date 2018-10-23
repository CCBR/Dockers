#!/bin/bash
reftargz=$1
t_tagalignfile=$2
i_tagalignfile=$3
t_ppqt=$4

t_tagalign=`echo $t_tagalignfile|awk -F "/" '{print $NF}'`
i_tagalign=`echo $i_tagalignfile|awk -F "/" '{print $NF}'`

ref=`echo $reftargz|awk -F "/" '{print $NF}'|sed "s/.tar.gz//g"`
ref=`echo $ref|awk '{print substr($1,1,2)}'`
if [ "$ref" == "hg" ]; then
genome="hs"
elif [ "$ref" == "mm" ]; then
genome="mm"
fi

outname=${t_tagalign}_vs_${i_tagalign}
narrowPeak=${outname}_peaks.narrowPeak
xls=${outname}_peaks.xls
bed=${outname}_summits.bed

extsize=`cat $t_ppqt|awk -F"\t" '{print $3}'|awk -F"," '{print $1}'`
tagsize=`cat $t_ppqt|awk -F"\t" '{print $5}'|awk -F"," '{print $1}'`

macs2 callpeak -t $t_tagalignfile -c $i_tagalignfile -n $outname --nomodel --extsize $extsize --tsize $tagsize -q 0.01 -f BED -g $genome --keep-dup=auto
