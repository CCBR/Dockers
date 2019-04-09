FROM ubuntu:18.04
MAINTAINER vishal.koparde@nih.gov


RUN 	mkdir -p /data \
&&		mkdir -p /opt \
&&		apt-get update \
&&		apt-get upgrade -y 

RUN 	apt-get install -y \
		samtools \
		bedtools \
		zlibc \
		zlib1g \
		zlib1g-dev

RUN 	apt-get install -y python \
		python-pip

RUN pip install --upgrade pip

RUN pip install deeptools==3.2.0

COPY Dockerfile /opt

WORKDIR /data

