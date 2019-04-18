FROM ubuntu:18.04
MAINTAINER vishal.koparde@nih.gov

RUN	mkdir -p /data \
&&	mkdir -p /opt \
&&	apt-get update \
&&	apt-get install -y \
	pigz \
	python3 \
	python3-pip \
&&	pip3 install --upgrade pip
RUN	pip3 install cutadapt==1.18
RUN	rm -rf /var/lib/apt/lists/* \
&& 	apt-get clean \
&&	apt autoremove --yes

RUN ln -s /usr/bin/python3 /usr/bin/python

COPY TruSeq_and_nextera_adapters.consolidated.fa /opt
COPY ccbr_cutadapt_se.bash /opt
COPY argparse.bash /opt
COPY Dockerfile /opt

RUN chmod a+x /opt/*.bash
ENV PATH="/opt:${PATH}"

WORKDIR /data
