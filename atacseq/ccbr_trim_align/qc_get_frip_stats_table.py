import sys
import pandas
x=pandas.read_csv(sys.stdin,sep="\t",header=None,dtype="str")
x.columns=["SampleName","ReadsFile","Peakcaller","CallingMethod","PeaksFile","FRiP","FRiDHS","FRiPromoters","FRiEnhancers"]
x.drop(columns=["ReadsFile","PeaksFile"],inplace=True)
y=x[x['CallingMethod']=="individual_replicate"]
macs2=y[y['Peakcaller']=="macs2"]
genrich=y[y['Peakcaller']=="genrich"]
genrich=genrich[['SampleName','FRiP']]
genrich.set_index('SampleName',inplace=True)
genrich.columns=['FRiP_Genrich']
macs2=macs2[['SampleName','FRiP','FRiDHS','FRiPromoters','FRiEnhancers']]
macs2.set_index('SampleName',inplace=True)
macs2.columns=['FRiP_MACS2','FRiDHS','FRiPromoters','FRiEnhancers']
joined=genrich.join(macs2)
joined.to_csv("FRiP_stats.txt",sep="\t")
