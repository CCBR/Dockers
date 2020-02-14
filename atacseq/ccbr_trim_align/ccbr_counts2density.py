import sys
n2bin=dict()
density=dict()
for i,b in enumerate(range(-1995,2000,10)):
	n2bin[i+1]=b
	density[b]=0.0
#for l in map(lambda x:x.strip().split("\t"),open(sys.argv[1]).readlines()):
for l in sys.stdin:
	l=l.strip().split("\t")
	density[n2bin[int(l[0])]]=int(l[1])
flanksum=0
for i in range(1,11):
	flanksum+=density[n2bin[i]]
for i in range(391,401):
	flanksum+=density[n2bin[i]]
flankavg=flanksum*1.0/20
alldensities=list()
for i in range(1,401):
	density[n2bin[i]]=density[n2bin[i]]*1.0/flankavg
	alldensities.append(density[n2bin[i]])
	print("%d\t%.4f"%(n2bin[i],density[n2bin[i]]))
print("# TSS enrichment: %.4f"%(max(alldensities)))