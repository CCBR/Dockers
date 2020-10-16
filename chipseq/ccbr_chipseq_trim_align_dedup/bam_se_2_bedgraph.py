#!/opt/conda/bin/python
import pysam,os
import argparse

parser = argparse.ArgumentParser(description='convert single end BAM to bedgraph .. requires samtools,bedSort and bedtools')
parser.add_argument('-i',dest='inBam',required=True,help='Input Bam File')
parser.add_argument('-q',dest='mapQ',type=int,required=False,help='mapQ value ... default 6',default=6)
parser.add_argument('-e',dest='ext',type=int,required=False,help='extension value ... default 200',default=200)
args = parser.parse_args()
if not os.path.exists(args.inBam+".bai"):
	pysam.index(args.inBam)
ext=args.ext
samfile = pysam.AlignmentFile(args.inBam, "rb")

cmd="samtools view -H "+args.inBam+"|grep ^@SQ|cut -f2,3|sed \"s/SN://g\"|sed \"s/LN://g\" > "+args.inBam+".genome"
os.system(cmd)

lens=dict()
for i in open(args.inBam+".genome").readlines():
	i=i.split("\t")
	lens[i[0]]=int(i[1])

o=open(args.inBam+".bed",'w')
for read in samfile.fetch():
	if read.is_unmapped:
		continue
	if read.is_supplementary:
		continue
	if read.is_secondary:
		continue
	if read.is_duplicate:
		continue
	if read.mapping_quality < args.mapQ:
		continue
	chrom=read.reference_name
	if read.is_reverse:
		start=read.reference_start - ext
		if start < 1:
			start=1
		end=read.reference_end
	else:
		start=read.reference_start
		end=read.reference_end + ext
		if end > lens[chrom]:
			end = lens[chrom]
	o.write("%s\t%s\t%s\t.\t.\t.\n"%(chrom,start,end))

o.close()

cmd="bedSort "+args.inBam+".bed "+args.inBam+".bed"
os.system(cmd)

cmd="bedtools genomecov -i "+args.inBam+".bed -g "+args.inBam+".genome -bg > "+args.inBam+".bg"
os.system(cmd)

os.remove(args.inBam+".bed")
os.remove(args.inBam+".genome")
