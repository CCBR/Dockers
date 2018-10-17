FROM ubuntu:16.04

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y pigz

RUN apt-get install -y python3 && \
    apt-get install -y python3-pip

RUN pip3 install --upgrade pip

RUN pip3 install cutadapt==1.18

COPY TruSeq_and_nextera_adapters.consolidated.fa /opt
COPY Dockerfile /opt

WORKDIR /data

RUN ls /opt

RUN head /opt/TruSeq_and_nextera_adapters.consolidated.fa

RUN ls /data

MAINTAINER vishal.koparde@nih.gov
