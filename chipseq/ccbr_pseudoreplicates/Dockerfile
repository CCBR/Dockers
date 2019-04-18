FROM 		ubuntu:18.04
MAINTAINER	vishal.koparde@nih.gov

RUN			mkdir -p /data \
&&			mkdir -p /opt \
&&			apt-get update \
&&			apt-get install -y \
			pigz \
			python3

RUN			rm -rf /var/lib/apt/lists/* \
&& 			apt-get clean \
&&			apt autoremove --yes

RUN 		ln -s /usr/bin/python3 /usr/bin/python

COPY 		Dockerfile /opt
COPY 		argparse.bash /opt
COPY		ccbr_create_pooled_and_pseudoreplicates_se.bash /opt
COPY		ccbr_create_pseudoreplicates_se.bash /opt
COPY		ccbr_create_pooledreplicate_se.bash /opt

RUN 		chmod a+x /opt/*.bash
ENV 		PATH=/opt:$PATH
WORKDIR		/data
