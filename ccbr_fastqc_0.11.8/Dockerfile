FROM ubuntu:16.04

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y

RUN apt-get install -y cpanminus

RUN cpanm FindBin 

RUN apt-get install -y  software-properties-common && \
    add-apt-repository ppa:webupd8team/java -y && \
    apt-get update && \
    echo oracle-java7-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections && \
    apt-get install -y oracle-java8-installer

RUN apt-get install -y unzip

COPY fastqc_v0.11.8.zip /opt
COPY Dockerfile /opt

WORKDIR /opt
RUN unzip fastqc_v0.11.8.zip
WORKDIR /opt/FastQC
RUN chmod 755 fastqc
RUN ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc

WORKDIR /data

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
