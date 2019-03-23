FROM 	bitnami/minideb:jessie

RUN	mkdir -p /data \
&&	mkdir -p /opt \
&&	apt-get update \
&&	apt-get install -y \
	python3-dev \
	python3-numpy \
	python3-setuptools \
	python3-matplotlib \
&&	apt-get clean

COPY Dockerfile /opt
COPY 2.0.4.2.tar.gz /opt
WORKDIR /opt
RUN apt-get install -y build-essential

RUN apt-get install -y gfortran

RUN apt-get install -y \
	cython \
	libopenblas-dev \
	liblapack-dev \
	libatlas-base-dev
RUN tar xzvf 2.0.4.2.tar.gz \
&&	cd idr-2.0.4.2 \
&&	python3 setup.py install

WORKDIR /data

MAINTAINER vishal.koparde@nih.gov
