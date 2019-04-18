FROM 		ubuntu:18.04
MAINTAINER	vishal.koparde@nih.gov

RUN 		mkdir -p /data \
&&			mkdir -p /opt

RUN 		apt-get update \
&&			apt-get upgrade -y \ 
&&			apt-get install -y\
			pigz \
			python3

RUN 		ln -s /usr/bin/python3 /usr/bin/python

RUN 		apt-get install -y bedtools=2.26.0+dfsg-5

RUN			rm -rf /var/lib/apt/lists/* \
&& 			apt-get clean \
&&			apt autoremove --yes

COPY 		Dockerfile /opt

WORKDIR 	/data
