FROM ubuntu:16.04

RUN mkdir -p /data2
RUN mkdir -p /opt2

RUN apt-get update && \
    apt-get upgrade -y 

RUN apt-get install -y --fix-missing zlibc zlib1g zlib1g-dev make gcc g++ wget git

WORKDIR /opt2
RUN git clone --recursive https://github.com/bwa-mem2/bwa-mem2 && cd bwa-mem2 && make
ENV PATH="/opt2/bwa-mem2":$PATH

WORKDIR /data2
COPY Dockerfile /opt2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
