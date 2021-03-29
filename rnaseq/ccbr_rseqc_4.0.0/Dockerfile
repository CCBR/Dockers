FROM nciccbr/ccbr_ubuntu_base_20.04:v1.0

RUN apt-get install -y zlibc zlib1g zlib1g-dev

RUN apt-get install -y r-base-core libboost-dev

RUN pip3 install --upgrade pip && pip3 install RSeQC

COPY Dockerfile /opt2

WORKDIR /data2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
