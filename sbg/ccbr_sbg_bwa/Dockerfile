FROM ubuntu:16.04

RUN mkdir -p /data
RUN mkdir -p /opt
RUN mkdir -p /scripts

COPY bwaCreateIndex.sh /scripts
RUN chmod a+x /scripts/bwaCreateIndex.sh
ENV PATH=$PATH:/scripts

RUN apt-get update && \
    apt-get upgrade -y 

RUN apt-get install -y zlibc zlib1g zlib1g-dev bzip2 make gcc

COPY bwa-0.7.17.tar.bz2 /opt

WORKDIR /opt
RUN tar -xjvf bwa-0.7.17.tar.bz2 && rm -f bwa-0.7.17.tar.bz2
WORKDIR /opt/bwa-0.7.17
RUN ls /opt/bwa-0.7.17

RUN make

ENV PATH="/opt/bwa-0.7.17":$PATH

WORKDIR /data

COPY Dockerfile /opt
RUN apt-get clean
MAINTAINER vishal.koparde@nih.gov
