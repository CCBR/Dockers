#!/bin/bash
set -e -x -o pipefail
ncpus=`nproc`
ARGPARSE_DESCRIPTION="Trim SE reads using cutadapt"      # this is optional
source /opt/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--infastq',required=True, help='input fastq.gz file')
parser.add_argument('--genome',required=True, help='mm9/10 or hg19/38')
parser.add_argument('--genomesize',required=True, help='effective genome size',type=int)
EOF

infastq=$INFASTQ
outfastq=`echo $infastq|sed "s/.fastq.gz/.trim.fastq.gz/g"`
samplename=`echo $infastq|sed "s/.fastq.gz//g"|sed "s/.R1//g"`

cutadapt \
--nextseq-trim=2 \
--trim-n \
-n 5 -O 5 \
-q 10,10 \
-m 35 \
-b file:/opt/TruSeq_and_nextera_adapters.consolidated.fa \
-j $ncpus \
-o $outfastq \
$infastq

nreads=`zcat $infastq|wc -l`
nreadstrim=`zcat $outfastq|wc -l`
echo "$nreads $nreadstrim"|awk '{printf("%d\tInput Nreads\n%d\tAfter trimming\n",$1/4,$2/4)}' > ${samplename}.nreads.txt

basefilename=`echo $outfastq|sed "s/.fastq.gz//g"`
bamfile=${basefilename}.bam
sortedbamfile=${basefilename}.sorted.bam
sortedbambaifile=${sortedbamfile}.bai
sortedq5bamfile=${basefilename}.sorted.Q5.bam
sortedq5bambaifile=${basefilename}.sorted.Q5.bam.bai
sortedbamflagstatfile=${sortedbamfile}.flagstat
sortedq5bamflagstatfile=${sortedq5bamfile}.flagstat

bwa mem -t $ncpus ${GENOME}_blacklist $outfastq | samtools view -@ $ncpus -f4 -b -o notAlignedToBlacklist.bam

bedtools bamtofastq -i notAlignedToBlacklist.bam -fq notAlignedToBlacklist.fastq

cat notAlignedToBlacklist.fastq|wc -l |awk '{printf("%d\tAfter removing reads aligning to blacklisted regions\n",$1/4)}' >> ${samplename}.nreads.txt

bwa mem -t $ncpus $GENOME notAlignedToBlacklist.fastq > $bamfile

samtools sort -@ $ncpus -o $sortedbamfile $bamfile

samtools index $sortedbamfile

samtools view -@ $ncpus -b -q 6 $sortedbamfile -o $sortedq5bamfile

samtools index $sortedq5bamfile

samtools flagstat $sortedbamfile > $sortedbamflagstatfile

samtools flagstat $sortedq5bamfile > $sortedq5bamflagstatfile

grep -m1 mapped $sortedbamflagstatfile |awk '{printf("%d\tMapped Reads\n",$1)}' >> ${samplename}.nreads.txt

grep -m1 mapped $sortedq5bamflagstatfile |awk '{printf("%d\tMAPQ>5 Reads\n",$1)}' >> ${samplename}.nreads.txt

tagalignfile=`echo $sortedq5bamfile |sed "s/.bam/.DD.tagAlign/g"`
outbamfile=`echo $sortedq5bamfile|sed "s/.bam/.DD.bam/g"`

macs2 filterdup -i $BAM -g $GENOMESIZE -o TmpTagAlign
awk -F"\t" -v OFS="\t" '{if ($2>0 && $3>0) {print}}' TmpTagAlign > TmpTagAlign2

samtools view -@ $ncpus -H $BAM|grep "^@SQ"|cut -f2,3|sed "s/SN://"|sed "s/LN://g" > GenomeFile
awk -F"\t" -v OFS="\t" '{print $1,1,$2}' GenomeFile |sort -k1,1 -k2,2n > GenomeFileBed
bedtools intersect -wa -f 1.0 -a TmpTagAlign2 -b GenomeFileBed > TagAlign
cat TagAlign|wc -l |awk '{printf("%d\tAfter Deduplication\n",$1)}' >> ${samplename}.nreads.txt

bedtools bedtobam -i TagAlign -g GenomeFile > OutBam

mv TagAlign $tagalignfile
mv OutBam $outbamfile

pigz -p $ncpus -f $tagalignfile
