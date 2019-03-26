FROM ubuntu:18.04

RUN	mkdir -p /data \
&&	mkdir -p /opt \
&&	apt-get update \
&&	apt-get install -y \
	pigz \
	python3 \
	python3-pip \
&&	pip3 install --upgrade pip
RUN	apt-get install -y \
	libcurl4 \
	libcurl4-openssl-dev \
	zlibc zlib1g zlib1g-dev
RUN	pip3 install deeptools==3.2.0
RUN	rm -rf /var/lib/apt/lists/* \
&& 	apt-get clean

COPY Dockerfile /opt

WORKDIR /data

MAINTAINER vishal.koparde@nih.gov

RUN	apt autoremove --yes
