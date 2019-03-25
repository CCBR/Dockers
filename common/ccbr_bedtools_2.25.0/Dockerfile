FROM ubuntu:16.04

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y 

RUN apt-get install -y bedtools

RUN apt-get clean

COPY Dockerfile /opt

WORKDIR /data

MAINTAINER vishal.koparde@nih.gov
