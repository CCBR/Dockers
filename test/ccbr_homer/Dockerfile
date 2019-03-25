FROM	bitnami/minideb:jessie

RUN	apt-get update \
&&	apt-get install -y r-base-core wget samtools

RUN 	mkdir -p /data/homer \
&&	mkdir -p /opt

COPY install_R_packages.R /data
WORKDIR /data
RUN Rscript install_R_packages.R

COPY configureHomer.pl /data/homer
WORKDIR /data/homer
RUN perl configureHomer.pl -install
RUN perl configureHomer.pl -install mouse 
RUN perl configureHomer.pl -install human 
RUN perl configureHomer.pl -install mm9 
RUN perl configureHomer.pl -install mm10 
RUN perl configureHomer.pl -install hg19
RUN perl configureHomer.pl -install hg38

ENV PATH="/data/homer/bin":$PATH
WORKDIR /data

COPY Dockerfile /opt

MAINTAINER vishal.koparde@nih.gov
