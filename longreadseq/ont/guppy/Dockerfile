FROM ubuntu:18.04

# installing cpu version of guppy basecaller for ONT

RUN mkdir -p /data2
RUN mkdir -p /opt2

RUN apt-get update && \
    apt-get upgrade -y 

RUN apt-get install -y ont-guppy

COPY Dockerfile /opt2

RUN chmod a+r /opt2/Dockerfile

WORKDIR /data2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
