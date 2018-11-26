FROM ubuntu:16.04

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y 

RUN apt-get install -y zlibc zlib1g zlib1g-dev make gcc g++ wget

WORKDIR /opt
RUN wget https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz
RUN tar -xzf 2.6.0a.tar.gz
RUN apt-get install -y g++
WORKDIR /opt/STAR-2.6.0a/source
RUN make STAR

ENV PATH="/opt/STAR-2.6.0a/source":$PATH

WORKDIR /data
COPY Dockerfile /opt

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
