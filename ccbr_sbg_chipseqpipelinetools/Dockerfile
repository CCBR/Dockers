FROM ubuntu:16.04

RUN mkdir -p /data
RUN mkdir -p /opt
RUN mkdir -p /scripts

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y pigz

RUN apt-get install -y python3 && \
    apt-get install -y python3-pip && \
    apt-get install -y python3-numpy 

RUN apt-get install -y python3-venv
WORKDIR /opt
RUN python3.5 -m venv py3
ENV PATH=$PATH:/opt/py3/bin
SHELL ["/bin/bash","-c"]


# install cutadapt
RUN source /opt/py3/bin/activate && \
    pip3 install --upgrade pip && \
    pip3 install cutadapt==1.18

COPY TruSeq_and_nextera_adapters.consolidated.fa /opt

# install bwa
RUN apt-get install -y zlibc zlib1g zlib1g-dev bzip2 make gcc

COPY bwa-0.7.17.tar.bz2 /opt

WORKDIR /opt
RUN tar -xjvf bwa-0.7.17.tar.bz2 && rm -f bwa-0.7.17.tar.bz2
WORKDIR /opt/bwa-0.7.17

RUN make

ENV PATH="/opt/bwa-0.7.17":$PATH

# install macs2
RUN apt-get install -y parallel && \
    apt-get install -y python python-pip python-numpy

RUN pip install --upgrade pip

RUN pip install MACS2

# install samtools
RUN apt-get install -y zlibc zlib1g zlib1g-dev bzip2 make gcc 
RUN apt-get install -y libncurses-dev libbz2-dev liblzma-dev

COPY samtools-1.9.tar.bz2 /opt

WORKDIR /opt
RUN tar -xjvf samtools-1.9.tar.bz2 && rm -f samtools-1.9.tar.bz2
WORKDIR /opt/samtools-1.9
RUN ls /opt/samtools-1.9

RUN make

ENV PATH="/opt/samtools-1.9":$PATH

# install bedtools
RUN apt-get install -y bedtools

# install deeptools
RUN pip install deeptools

#install spp
RUN apt-get install -y r-base-core libboost-dev

COPY spp_1.14.tar.gz /opt

COPY install_R_packages.R /opt

RUN cd /opt && Rscript install_R_packages.R && rm -f spp_1.14.tar.gz

COPY run_spp.R /opt
RUN chmod a+rx /opt/run_spp.R

ENV PATH="/opt":$PATH

#scripts
COPY ccbr_bwaAlignSE.sh /scripts
COPY ccbr_macsFilterdup.sh /scripts
COPY ccbr_macsSEPeakCalling.sh /scripts
COPY ccbr_fingerprintPlot.sh /scripts
COPY ccbr_bam2bw.sh /scripts

RUN chmod a+rx /scripts/*.sh
ENV PATH=$PATH:"/scripts"


#cleanup
RUN apt-get clean
COPY tool_list.txt /opt
COPY Dockerfile /opt

WORKDIR /data

MAINTAINER vishal.koparde@nih.gov
