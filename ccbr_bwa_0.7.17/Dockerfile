FROM ubuntu:16.04

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y 

RUN apt-get install -y zlibc zlib1g zlib1g-dev bzip2 make gcc

COPY bwa-0.7.17.tar.bz2 /opt
COPY Dockerfile /opt

WORKDIR /opt
RUN tar -xjvf bwa-0.7.17.tar.bz2
WORKDIR /opt/bwa-0.7.17
RUN ls /opt/bwa-0.7.17

RUN make

ENV PATH="/opt/bwa-0.7.17":$PATH

WORKDIR /data

MAINTAINER vishal.koparde@nih.gov
