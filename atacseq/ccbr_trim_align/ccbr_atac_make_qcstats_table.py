#!/usr/bin/env python
# coding: utf-8

import pandas


nreads=pandas.read_csv("Nreads_mqc.csv",sep="\t")
x=nreads.set_index("SampleName")
cols=list(x)
x['Nreads']=x.loc[:,cols].sum(axis=1)
for i in cols:
    j=i+"_perc"
    x[j]=x[i]*100.0/x['Nreads']
x.head()
cols2=['Nreads']
for i in cols:
    cols2.append(i)
    cols2.append(i+"_perc")
x=x.loc[:,cols2]
x.head()
tss_stats=pandas.read_csv("data.tss_knicking_sites.txt",header=None,sep="\t",names=["SampleName","N_TSS_gt20knicks","TSSscore"])
nrf=pandas.read_csv("NRF_stats.txt",sep="\t")
frip=pandas.read_csv("FRiP_stats.txt",sep="\t")
fld=pandas.read_csv("FLD_stats.txt",sep="\t")


allstats=x.join(tss_stats.set_index("SampleName"))
allstats=allstats.join(nrf.set_index("SampleName"))
allstats=allstats.join(frip.set_index("SampleName"))
allstats=allstats.join(fld.set_index("SampleName"))


macs=pandas.read_csv("MACS2_Peak_Annotations_mqc.csv",sep="\t")
x=macs.set_index("SampleName")
x.columns=list(map(lambda z:"macs2_"+z,list(x)))
cols=list(x)
x['macs2_Npeaks']=x.loc[:,cols].sum(axis=1)
for i in cols:
    j=i+"_perc"
    x[j]=x[i]*100.0/x['macs2_Npeaks']
cols2=['macs2_Npeaks']
for i in cols:
    cols2.append(i)
    cols2.append(i+"_perc")
x=x.loc[:,cols2]


allstats=allstats.join(x)


genrich=pandas.read_csv("Genrich_Peak_Annotations_mqc.csv",sep="\t")
x=genrich.set_index("SampleName")
x.columns=list(map(lambda z:"genrich_"+z,list(x)))
cols=list(x)
x['genrich_Npeaks']=x.loc[:,cols].sum(axis=1)
for i in cols:
    j=i+"_perc"
    x[j]=x[i]*100.0/x['genrich_Npeaks']
cols2=['genrich_Npeaks']
for i in cols:
    cols2.append(i)
    cols2.append(i+"_perc")
x=x.loc[:,cols2]


allstats=allstats.join(x)

allstats=allstats.round(3)
allstats.to_csv("QCStats.txt",sep="\t")
