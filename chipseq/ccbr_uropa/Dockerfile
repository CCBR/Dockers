#r-base-core
FROM 	bitnami/minideb:jessie

RUN	mkdir -p /data \
&& 	mkdir -p /opt \
&&	apt-get update \
&&	apt-get install -y \
	python3 python3-pip
RUN	apt-get install -y \
	zlib1g-dev \
	libbz2-dev \
	liblzma-dev
RUN	pip3 install uropa
RUN	rm -rf /var/lib/apt/lists/* \
&& 	apt-get clean

COPY	install_R_packages.R /opt
COPY 	Dockerfile /opt

WORKDIR	/opt

RUN	Rscript install_R_packages.R




