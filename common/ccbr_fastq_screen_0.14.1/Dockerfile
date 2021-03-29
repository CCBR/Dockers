FROM nciccbr/ccbr_ubuntu_base_20.04:v1.0

RUN apt-get update 

RUN apt-get install -y zlibc zlib1g zlib1g-dev make gcc g++ wget \
   bowtie2 libgd-perl pkg-config libgd-dev

RUN cpanm GD::Graph::bars

ADD https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v0.14.1.tar.gz /opt2
COPY v0.14.1.tar.gz /opt2
WORKDIR /opt2
RUN tar xzvf v0.14.1.tar.gz
WORKDIR /opt2/FastQ-Screen-0.14.1
RUN chmod 755 fastq_screen
ENV PATH="/opt2/FastQ-Screen-0.14.1":$PATH

COPY Dockerfile /opt2

WORKDIR /data2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
