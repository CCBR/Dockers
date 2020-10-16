FROM ubuntu:16.04

RUN mkdir -p /data2
RUN mkdir -p /opt2

RUN apt-get update && \
    apt-get upgrade -y 

RUN apt-get install -y --fix-missing zlibc zlib1g zlib1g-dev make gcc g++ wget

WORKDIR /opt2
RUN wget https://github.com/alexdobin/STAR/archive/2.7.6a.tar.gz
RUN tar -xzf 2.7.6a.tar.gz
RUN apt-get install -y g++
WORKDIR /opt2/STAR-2.7.6a/source
RUN make STAR

ENV PATH="/opt2/STAR-2.7.6a/source/":$PATH

WORKDIR /data2
COPY Dockerfile /opt2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
