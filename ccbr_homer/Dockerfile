FROM kopardev/ccbr_samtools_1.9

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y 

RUN apt-get install -y r-base-core

COPY install_R_packages.R /data

RUN cd /data && Rscript install_R_packages.R

RUN mkdir -p /data/homer

COPY configureHomer.pl /data/homer

RUN apt-get install -y wget

RUN cd /data/homer && perl configureHomer.pl -install

ENV PATH="/data/homer/bin":$PATH

RUN perl configureHomer.pl -install mouse 
RUN perl configureHomer.pl -install human 
RUN perl configureHomer.pl -install mm9 
RUN perl configureHomer.pl -install mm10 
RUN perl configureHomer.pl -install hg19
RUN perl configureHomer.pl -install hg38
 
WORKDIR /data

MAINTAINER vishal.koparde@nih.gov
