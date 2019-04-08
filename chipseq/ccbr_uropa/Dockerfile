#r-base-core
#FROM 	bitnami/minideb:jessie
FROM 	ubuntu:18.04
MAINTAINER vishal.koparde@nih.gov


RUN	mkdir -p /data \
&& 	mkdir -p /opt \
&&	apt-get update \
&&	apt-get install -y \
	python3 python3-pip
RUN	apt-get install -y \
	zlib1g-dev \
	libbz2-dev \
	liblzma-dev

WORKDIR /opt
COPY	3.1.0.tar.gz /opt
RUN	tar xzvf 3.1.0.tar.gz
WORKDIR	/opt/UROPA-3.1.0
RUN	python3 setup.py install
#RUN	pip3 install uropa

ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt update && apt-get install -y r-base-core
RUN	apt-get -y install libcurl4-openssl-dev libssl-dev libboost-dev
RUN	rm -rf /var/lib/apt/lists/* \
&& 	apt-get clean

COPY	install_R_packages.R /opt
COPY 	Dockerfile /opt

WORKDIR	/opt

RUN	apt-get install python3.6 \
&&	ln -s /usr/bin/python3.6 /usr/bin/python
RUN	Rscript install_R_packages.R

COPY argparse.bash /opt
COPY ccbr_uropa_annotate_tf.bash /opt

COPY tf.json /data

RUN chmod a+x /opt/*.bash
ENV PATH=/opt:$PATH
WORKDIR /data