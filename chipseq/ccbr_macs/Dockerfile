FROM ubuntu:18.04
MAINTAINER vishal.koparde@nih.gov

RUN	mkdir -p /data \
&&	mkdir -p /opt \
&&	apt-get update \
&&	apt-get upgrade -y

RUN apt-get install -y \
	python \
	python-pip \
	python-numpy

RUN pip install --upgrade pip

RUN pip install MACS2

COPY Dockerfile /opt
COPY ccbr_macs_callpeak_se.bash /opt
COPY ccbr_macs_callpeak_pe.bash /opt
COPY argparse.bash /opt

RUN chmod a+x /opt/*.bash
ENV PATH=/opt:$PATH

WORKDIR /data
