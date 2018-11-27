FROM ubuntu:16.04

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y pigz

RUN apt-get install -y bzip2 make gcc

RUN apt-get install -y python3 && \
    apt-get install -y python3-pip

RUN apt-get install -y zlibc zlib1g zlib1g-dev

RUN pip3 install --upgrade pip

RUN pip3 install RSeQC

COPY Dockerfile /opt

WORKDIR /data

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
