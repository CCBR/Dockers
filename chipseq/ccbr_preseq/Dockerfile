FROM 	ubuntu:18.04
MAINTAINER vishal.koparde@nih.gov

RUN	mkdir -p /data \
&&	mkdir -p /opt \
&&	apt-get update \
&&	apt-get upgrade -y

#RUN		apt-cache policy samtools

RUN apt-get install -y \
	samtools=1.7-1 \
	build-essential \
	zlibc zlib1g zlib1g-dev libgsl-dev

WORKDIR		/opt
COPY		preseq_v2.0.3.tar.bz2 /opt
RUN			tar xjvf preseq_v2.0.3.tar.bz2
WORKDIR		/opt/preseq
RUN			make all
ENV 		PATH=/opt/preseq:$PATH

COPY 		argparse.bash /opt
COPY 		ccbr_chipseq_preseq.bash /opt
COPY		ccbr_nrf.py /opt
COPY 		Dockerfile /opt

RUN 		chmod a+x /opt/*.bash
ENV 		PATH=/opt:$PATH
WORKDIR 	/data