FROM ubuntu:18.04

RUN	mkdir -p /data \
&&	mkdir -p /opt \
&&	apt-get update \
&&	apt-get install -y \
	samtools \
&&	rm -rf /var/lib/apt/lists/* \
&& 	apt-get clean \
&&	apt autoremove --yes

WORKDIR /data

MAINTAINER vishal.koparde@nih.gov
