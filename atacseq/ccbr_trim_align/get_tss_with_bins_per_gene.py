import pybedtools
import sys
GTF=sys.argv[1]
gtf=pybedtools.BedTool(GTF)
count=0
include_chroms=[]
for i in range(1,23):
    include_chroms.append("chr"+str(i))
include_chroms.append("chrX")
include_chroms.append("chrY")
namecount=dict()
for i in gtf:
    if not i.chrom in include_chroms:
        continue
    if i.fields[2]==u'gene' and "coding" in str(i.attrs):
        tss_string=""
        count+=1
        genename=i.attrs['gene_id']+"|"+i.name
        bn=0
        for j in range(-2000,2000,10):
            s=int(i.start)+j
            e=s+10
            bn+=1
            realbin=bn
            if i.strand=="-":
               realbin=401-bn
            genenamebn=genename+"|"+str(realbin)
            tss_string+="%s\t%s\t%s\t%s\t%s\t%s\n"%(i.chrom,str(s),str(e),genenamebn,i.score,i.strand)
        tssfile=open(i.attrs['gene_id']+".bed",'w')
        tssfile.write(tss_string)
        tssfile.close()
