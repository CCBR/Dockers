#!/bin

# genome           blacklistversion
# mm10             v1
# hg19             v3
# hg38             v3

# bed files are downloaded from https://sites.google.com/site/anshulkundaje/projects/blacklists
# fasta files are created using bedtools getfasta on biowulf and then transferred here
# chrM and chr_rDNA are added to blacklist fasta to create
# 1.mm10.blacklist.chrM.chr_rDNA.fa
# 2.hg19.blacklist_v3.chrM.chr_rDNA.fa          
# 3.hg38.blacklist_v3.chrM.chr_rDNA.fa          

# create index
bwa index -p mm10_blacklist mm10.blacklist.chrM.chr_rDNA.fa
bwa index -p hg19_blacklist hg19.blacklist_v3.chrM.chr_rDNA.fa
bwa index -p hg38_blacklist hg38.blacklist_v3.chrM.chr_rDNA.fa

# gzip files
gzip mm10.blacklist.chrM.chr_rDNA.fa
gzip hg19.blacklist_v3.chrM.chr_rDNA.fa
gzip hg38.blacklist_v3.chrM.chr_rDNA.fa
