FROM ubuntu:18.04
MAINTAINER vishal.koparde@nih.gov

RUN	mkdir -p /data \
&&	mkdir -p /opt \
&&	apt-get update \
&&	apt-get install -y \
	bwa \
	samtools \
	bedtools \
&&	rm -rf /var/lib/apt/lists/* \
&& 	apt-get clean \
&&	apt autoremove --yes


COPY Dockerfile /opt
COPY ccbr_bwa_chipseq_align_se.bash /opt
COPY argparse.bash /opt

RUN chmod a+x /opt/*.bash
ENV PATH=/opt:$PATH
WORKDIR /data
