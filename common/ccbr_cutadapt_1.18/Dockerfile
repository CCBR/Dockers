FROM ubuntu:18.04

RUN	mkdir -p /data \
&&	mkdir -p /opt \
&&	apt-get update \
&&	apt-get install -y \
	pigz \
	python3 \
	python3-pip \
&&	pip3 install --upgrade pip
RUN	pip3 install cutadapt==1.18
RUN	rm -rf /var/lib/apt/lists/* \
&& 	apt-get clean

COPY TruSeq_and_nextera_adapters.consolidated.fa /opt
COPY Dockerfile /opt

WORKDIR /data

MAINTAINER vishal.koparde@nih.gov

RUN	apt-get remove --yes --purge build-essential \
&&	apt autoremove --yes
