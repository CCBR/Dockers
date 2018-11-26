FROM ubuntu:16.04

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y 

RUN apt-get install -y zlibc zlib1g zlib1g-dev make gcc g++ wget
RUN apt-get install -y bowtie bowtie2

WORKDIR /data
COPY Dockerfile /opt

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
