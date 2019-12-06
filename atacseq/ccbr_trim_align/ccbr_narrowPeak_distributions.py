import sys,os,numpy,pandas
import argparse

import re 

def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

parser = argparse.ArgumentParser(description='calculate various QC distributions from narrowPeak file')
parser.add_argument('-i',dest='inNarrowPeak',required=True,help='Input NarrowPeak file')
parser.add_argument('-w',dest='widthdensity',required=False,help='Output peak width density distribution file')
parser.add_argument('-c',dest='chrdist',required=False,help='peaks per chromosome distribution file')
args = parser.parse_args()
d=pandas.read_csv(args.inNarrowPeak,sep="\t",header=None,names=["chr","start","end","peakid","score","strand","signal","pvalue","qvalue","summit"],usecols=range(10))
if args.widthdensity:
    width=d["end"]-d["start"]
    yvalues,xranges=numpy.histogram(width,density=True,bins=500)
    xvalues=[(xranges[i]+xranges[i+1])*0.5 for i in range(len(xranges)-1)]
    o=open(args.widthdensity,'w')
    for x,y in zip(xvalues,yvalues):
        if y*1000>0.01:
            o.write("%.2f\t%.6f\n"%(x,y*1000.0))
    o.close()
if args.chrdist:
    chrdist=dict()
    s=0
    for i in list(d["chr"]):
        if not i in chrdist:
            chrdist[i]=0
        chrdist[i]+=1
        s+=1
    o=open(args.chrdist,'w')
    o2=open(args.chrdist+".perc",'w')
    l=sorted_nicely(chrdist.keys())
    for x in l:
        y=chrdist[x]
        o.write("%s\t%d\n"%(x,y))
        o2.write("%s\t%.2f\n"%(x,y*100.0/s))
    o.close()
    o2.close()


