FROM ubuntu:16.04

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y 

RUN apt-get install -y zlibc zlib1g zlib1g-dev

RUN apt-get install -y python && \
    apt-get install -y python-pip

RUN pip install --upgrade pip

RUN pip install deeptools

COPY Dockerfile /opt

WORKDIR /data

MAINTAINER vishal.koparde@nih.gov
