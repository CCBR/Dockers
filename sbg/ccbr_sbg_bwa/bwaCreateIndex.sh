#!/bin/bash
fasta=$1
prefix=`echo ${fasta%.*}|awk -F "/" '{print $NF}'`
bwa index -p $prefix $fasta
rm -f $fasta
outtarfile=$prefix.tar.gz
tar -czvf $outtarfile $prefix.*

