FROM 	ubuntu:18.04
MAINTAINER vishal.koparde@nih.gov


RUN	mkdir -p /data \
&& 	mkdir -p /opt \
&&	apt-get update 

ENV 	TZ=America/New_York
RUN 	ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN 	apt update && apt-get install -y r-base-core
RUN	apt-get -y install libcurl4-openssl-dev libssl-dev libboost-dev

RUN	apt-get -y install libxml2-dev
RUN	apt-get -y install libmariadbclient-dev libmariadb-client-lgpl-dev

RUN	rm -rf /var/lib/apt/lists/* \
&& 	apt-get clean

COPY	install_R_packages.R /opt
COPY 	Dockerfile /opt

WORKDIR	/opt
RUN	Rscript install_R_packages.R

COPY argparse.bash /opt

RUN chmod a+x /opt/*.bash
ENV PATH=/opt:$PATH
WORKDIR /data

