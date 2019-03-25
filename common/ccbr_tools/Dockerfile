#tools
#bedtools 	v2.21.0
#samtools 	0.1.19-96b5f2294a
#bowtie		1.1.1-2
#bowtie2	2.2.4-1
#bwa		0.7.10-1
#
#
#also included ...
#python		2.7.9-1


FROM 	bitnami/minideb:jessie

RUN	apt-get update \
&&	apt-get install -y \
	bedtools \
	samtools \
	bwa \
	bowtie \
	bowtie2 \
&&	rm -rf /var/lib/apt/lists/* \
&& 	apt-get clean


