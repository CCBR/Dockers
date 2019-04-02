#!/bin/bash
set -e -x -o pipefail
ARGPARSE_DESCRIPTION="Call peaks with macs2"      # this is optional
source /opt/argparse.bash || exit 1
# source /home/kopardev/Desktop/Dockers/chipseq/ccbr_macs_peakcalling/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--broad', action='store_true',
                    default=False, help='Provide this option to make broad peak calls')
parser.add_argument('--treatment',required=True, help='Treatment tagAlign file')
parser.add_argument('--input',required=True, help='Input tagAlign file')
parser.add_argument('--outprefix',required=True, help='Output name prefix')
parser.add_argument('--genomesize',required=True, help='effective genome size',type=int)
parser.add_argument('--treatmentppqt',required=True, help='ppqt output after running spp')
EOF

treatmentextsize=`cat $TREATMENTPPQT|awk -F"\t" '{print $3}'|awk -F"," '{print $1}'`
if [ $treatmentextsize -lt 0 ];then
    treatmentextsize=`cat $TREATMENTPPQT|awk -F"\t" '{print $3}'|awk -F"," '{print $2}'`
fi  
if [ $treatmentextsize -lt 0 ];then
    treatmentextsize="200"
fi 

# echo $BROAD
if [ $BROAD ]
then
	echo "Broad"
	macs2 callpeak -t $TREATMENT -c $CONTROL -g $GENOMESIZE -n $OUTPREFIX --broad --broad-cutoff 0.01 --keep-dup='all' --nomodel --extsize $treatmentextsize
else
	echo "Narrow"
	macs2 callpeak -t $TREATMENT -c $CONTROL -g $GENOMESIZE -n $OUTPREFIX -q 0.01 --keep-dup='all' --nomodel --extsize $treatmentextsize
fi
