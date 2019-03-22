FROM 	bitnami/minideb:jessie

RUN	mkdir -p /data \
&&	mkdir -p /opt \
&&	apt-get update \
&&	apt-get install -y bowtie bowtie2 \
&&	apt-get clean

COPY Dockerfile /opt

WORKDIR /data

MAINTAINER vishal.koparde@nih.gov
