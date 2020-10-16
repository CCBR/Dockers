#!/opt/conda/bin/python
import pysam,os
import argparse

parser = argparse.ArgumentParser(description='convert paired end BAM to bedgraph .. requires samtools,bedSort and bedtools')
parser.add_argument('-i',dest='inBam',required=True,help='Input Bam File')
parser.add_argument('-q',dest='mapQ',type=int,required=False,help='mapQ value ... default 6',default=6)
args = parser.parse_args()
if not os.path.exists(args.inBam+".bai"):
	pysam.index(args.inBam)
samfile = pysam.AlignmentFile(args.inBam, "rb")
mapq=dict()
for read in samfile.fetch():
	if read.is_unmapped:
		continue
	if read.is_supplementary:
		continue
	if read.is_secondary:
		continue
	if read.is_duplicate:
		continue
	if read.is_proper_pair:
		if read.mapping_quality < args.mapQ and read.query_name in mapq:
			del mapq[read.query_name]
		if read.mapping_quality >= args.mapQ:
			if not read.query_name in mapq:
				mapq[read.query_name]=dict()
				mapq[read.query_name]['chr']=read.reference_name
				mapq[read.query_name]['coords']=list()
		
			mapq[read.query_name]['coords'].append(read.reference_start)
			mapq[read.query_name]['coords'].append(read.reference_end)

o=open(args.inBam+".bed",'w')
for q in mapq.keys():
	chrom=mapq[q]['chr']
	start=min(mapq[q]['coords'])
	end=max(mapq[q]['coords'])
	o.write("%s\t%s\t%s\t.\t.\t.\n"%(chrom,start,end))
o.close()

cmd="bedSort "+args.inBam+".bed "+args.inBam+".bed"
os.system(cmd)

cmd="samtools view -H "+args.inBam+"|grep ^@SQ|cut -f2,3|sed \"s/SN://g\"|sed \"s/LN://g\" > "+args.inBam+".genome"
os.system(cmd)

cmd="bedtools genomecov -i "+args.inBam+".bed -g "+args.inBam+".genome -bg > "+args.inBam+".bg"
os.system(cmd)

os.remove(args.inBam+".bed")
os.remove(args.inBam+".genome")



