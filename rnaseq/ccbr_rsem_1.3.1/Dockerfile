FROM ubuntu:16.04

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y 

RUN apt-get install -y zlibc zlib1g zlib1g-dev make gcc g++ wget
RUN apt-get install -y r-base-core libboost-dev
RUN apt-get install -y libncurses-dev libbz2-dev liblzma-dev
RUN apt-get install -y bowtie bowtie2

WORKDIR /opt
COPY 2.6.0a.tar.gz /opt
RUN tar -xzf 2.6.0a.tar.gz
WORKDIR /opt/STAR-2.6.0a/source
RUN make STAR
ENV PATH="/opt/STAR-2.6.0a/source":$PATH

WORKDIR /opt
COPY v1.3.1.tar.gz /opt
RUN tar -xzf v1.3.1.tar.gz
WORKDIR /opt/RSEM-1.3.1
RUN make
RUN make ebseq

ENV PATH="/opt/RSEM-1.3.1":$PATH

WORKDIR /data
COPY Dockerfile /opt
RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
