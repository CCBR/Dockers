FROM ubuntu:16.04

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y

RUN apt-get install -y default-jre
RUN apt-get install -y bedtools bowtie2 samtools bwa

RUN	apt-get install -y \
	pigz \
	python3 \
	python3-pip \
&&	pip3 install --upgrade pip
RUN	pip3 install cutadapt==1.18

RUN apt-get install -y cpanminus unzip

RUN cpanm FindBin 

COPY fastqc_v0.11.8.zip /opt

WORKDIR /opt
RUN unzip fastqc_v0.11.8.zip
WORKDIR /opt/FastQC
RUN chmod 755 fastqc
RUN ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc

COPY TruSeq_and_nextera_adapters.consolidated.fa /opt
COPY Dockerfile /opt
COPY picardcloud.jar /opt
COPY bedSort /opt
COPY argparse.bash /opt
COPY ccbr_cutadapt_pe.bash /opt
COPY ccbr_remove_blacklisted_reads_pe.bash /opt
COPY ccbr_atac_trim_align_pe.bash /opt
COPY ccbr_bowtie2_align_pe.bash /opt
COPY atac_assign_multimappers.py /opt
COPY the_blacklists.fa /opt
COPY the_blacklists.fa.amb /opt
COPY the_blacklists.fa.ann /opt
COPY the_blacklists.fa.bwt /opt
COPY the_blacklists.fa.pac /opt
COPY the_blacklists.fa.sa /opt

RUN chmod a+x /opt/bedSort
RUN chmod a+x /opt/*.bash
RUN chmod a+x /opt/*.py

ENV PATH="/opt:${PATH}"

WORKDIR /data

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
