import sys
import pandas
import numpy
x=pandas.read_csv(sys.argv[1],header=0,sep="\t",usecols=["width"])
x.dropna(inplace=True)
a,b=numpy.histogram(x['width'],bins=400,range=(0,20000))
c=[]
for i in range(len(b)-1):
	c.append(0.5*(b[i]+b[i+1]))
d=a*100.0/sum(a)
for i in range(len(d)):
	print("%d\t%.5f"%(c[i],d[i]))
