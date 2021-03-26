FROM ubuntu:16.04

# installing cpu version of guppy basecaller for ONT

RUN mkdir -p /data2
RUN mkdir -p /opt2

# Need to add cloud repo to install R v3.6 on ubuntu v16.04
# else default is R v3.2.3 which does not have required libraries like RCurl
RUN apt-get update && \
    apt-get upgrade -y 

# ENV PLATFORM=$(lsb_release -cs)
ENV PLATFORM="xenial"

RUN apt-get install -y wget lsb-release apt-transport-https
RUN wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add -
RUN echo "deb http://mirror.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | tee /etc/apt/sources.list.d/nanoporetech.sources.list 

RUN apt-get update && \
  apt install -y ont-guppy-cpu

COPY Dockerfile /opt2

RUN chmod a+r /opt2/Dockerfile

WORKDIR /data2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
