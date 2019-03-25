#tools
#bedtools 	v2.26.0
#samtools 	1.7
#bowtie		1.2.2
#bowtie2	2.3.4-1
#bwa		0.7.17-r1188
#sambamba	0.6.7-pre1
#
#
#also included ...
#python		2.7.15rc1


#FROM 	bitnami/minideb:jessie
FROM 	ubuntu:18.04

RUN	mkdir /data \
&&	mkdir /opt \
&&	apt-get update \
&&	apt-get install -y \
	bedtools \
	samtools \
	bwa \
	bowtie \
	bowtie2 \
	libsys-hostname-long-perl \
	sambamba \
&&	rm -rf /var/lib/apt/lists/* \
&& 	apt-get clean

WORKDIR /data
COPY Dockerfile /opt

MAINTAINER vishal.koparde@nih.gov

