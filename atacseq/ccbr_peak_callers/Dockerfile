FROM			nciccbr/ccbr_macs:v0.0.3
MAINTAINER		vishal.koparde@nih.gov

RUN				apt-get install -y zlib1g zlib1g-dev

COPY			Genrich_v0.5.tar.gz /opt

RUN				cd /opt \
&&				tar xzvf Genrich_v0.5.tar.gz \
&&				cd Genrich-0.5 \
&&				make

ENV				PATH=/opt/Genrich-0.5:$PATH

