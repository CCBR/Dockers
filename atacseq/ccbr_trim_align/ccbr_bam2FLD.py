import pysam,sys,os,numpy
from scipy.signal import find_peaks_cwt
import argparse
def get_max_distance(l):
    mi=min(l)
    ma=max(l)
    return ma-mi

def is_peak_in_range(peaks,a,b):
    if len(numpy.where(numpy.logical_and(peaks>a,peaks<b))[0]) >= 1:
        return "PRESENT"
    else:
        return "ABSENT"

def get_sum_in_range(freq,a,b):
    sum=0
    for i in range(a,b,1):
        sum+=freq[i]
    return sum

parser = argparse.ArgumentParser(description='calculate normalized fragment length distribution from PE bam file')
parser.add_argument('-i',dest='inBam',required=True,help='Input Bam File')
parser.add_argument('-o',dest='out',required=True,help='Output tab-delimited File')
args = parser.parse_args()
if not os.path.exists(args.inBam+".bai"):
    pysam.index(args.inBam)
samfile = pysam.AlignmentFile(args.inBam, "rb")
refcoordsperread=dict()
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
        if not read.query_name in refcoordsperread:
            refcoordsperread[read.query_name]=list()
        refcoordsperread[read.query_name].append(int(read.reference_start))
        refcoordsperread[read.query_name].append(int(read.reference_end))
fraglens=list()
for k,v in refcoordsperread.items():
	fraglens.append(get_max_distance(refcoordsperread[k]))
o=open(args.out,'w')
s=len(fraglens)
freqs=dict()
for i in range(0,max(fraglens)+1,1):
    freqs[i]=0
for i in range(min(fraglens),max(fraglens)+1,1):
    c=fraglens.count(i)
    freqs[i]=c
    if c != 0:
        o.write("%d\t%.6f\n"%(i,c*1000.0/s))

peaks=find_peaks_cwt(freqs.values(),numpy.array([25]))

nfr_peak=is_peak_in_range(peaks,20,90)
mono_nuc_peak=is_peak_in_range(peaks,120,250)
di_nuc_peak=is_peak_in_range(peaks,300,500)

nfr_reads=get_sum_in_range(freqs,0,150)
mono_nuc_reads=get_sum_in_range(freqs,150,300)
non_nfr_reads=sum(freqs.values())-nfr_reads
nfr_to_mono_nuc_ratio=nfr_reads*1.0/mono_nuc_reads
nfr_to_other_ratio=nfr_reads*1.0/non_nfr_reads

o.write("# PEAKS : %s\n"%(peaks))
o.write("# NUCLEOSOME_FREE_PEAK_PRESENT : %s\n"%(nfr_peak))
o.write("# MONONUCLEOSOME_PEAK_PRESENT : %s\n"%(mono_nuc_peak))
o.write("# DINUCLEOSOME_PEAK_PRESENT : %s\n"%(di_nuc_peak))
o.write("# NUCLEOSOME_FREE_READ_FRACTION : %0.3f\n"%(nfr_reads*1.0/s))
o.write("# MONONUCLEOSOME_READ_FRACTION : %0.3f\n"%(mono_nuc_reads*1.0/s))
o.write("# NFR_TO_MONONUCLEOSOME_READS_RATIO : %0.3f\n"%(nfr_to_mono_nuc_ratio))
o.write("# NFR_TO_NON_NFR_READS_RATIO : %0.3f\n"%(nfr_to_other_ratio))

o.close()
