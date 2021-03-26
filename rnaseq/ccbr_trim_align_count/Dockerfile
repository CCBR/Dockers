FROM ubuntu:16.04

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y

RUN apt-get install -y default-jre \
bowtie bedtools bowtie2 samtools bwa \
pigz python3 python3-pip r-base-core \
cpanminus unzip zlibc make gcc g++ \
zlib1g-dev \
libboost-dev \
libncurses-dev \
libbz2-dev \
liblzma-dev

RUN	pip3 install --upgrade pip
RUN	pip3 install cutadapt==1.18

RUN cpanm FindBin 

COPY fastqc_v0.11.8.zip /opt

WORKDIR /opt
RUN unzip fastqc_v0.11.8.zip
WORKDIR /opt/FastQC
RUN chmod 755 fastqc
RUN ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc

WORKDIR /opt
COPY 2.6.0a.tar.gz /opt
RUN tar -xzf 2.6.0a.tar.gz
WORKDIR /opt/STAR-2.6.0a/source
RUN make STAR

ENV PATH="/opt/STAR-2.6.0a/source":$PATH

WORKDIR /opt
COPY v1.3.1.tar.gz /opt
RUN tar -xzf v1.3.1.tar.gz
WORKDIR /opt/RSEM-1.3.1
RUN make
RUN make ebseq

ENV PATH="/opt/RSEM-1.3.1":$PATH

COPY TruSeq_and_nextera_adapters.consolidated.fa /opt
COPY Dockerfile /opt
COPY picardcloud.jar /opt
COPY argparse.bash /opt

RUN chmod a+x /opt/*.bash

ENV PATH="/opt:${PATH}"

WORKDIR /data

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
