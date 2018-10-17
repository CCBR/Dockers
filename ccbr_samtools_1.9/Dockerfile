FROM ubuntu:16.04

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y 

RUN apt-get install -y zlibc zlib1g zlib1g-dev bzip2 make gcc 
RUN apt-get install -y libncurses-dev libbz2-dev liblzma-dev

COPY samtools-1.9.tar.bz2 /opt
COPY Dockerfile /opt

WORKDIR /opt
RUN tar -xjvf samtools-1.9.tar.bz2
WORKDIR /opt/samtools-1.9
RUN ls /opt/samtools-1.9

RUN make

ENV PATH="/opt/samtools-1.9":$PATH

WORKDIR /data

MAINTAINER vishal.koparde@nih.gov
