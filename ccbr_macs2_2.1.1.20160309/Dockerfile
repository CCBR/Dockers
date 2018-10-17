FROM ubuntu:16.04

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y parallel

RUN apt-get install -y python && \
    apt-get install -y python-pip

RUN apt-get install -y python-numpy

RUN pip install --upgrade pip

RUN pip install MACS2

COPY Dockerfile /opt

WORKDIR /data

MAINTAINER vishal.koparde@nih.gov
