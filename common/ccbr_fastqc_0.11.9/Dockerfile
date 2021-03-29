FROM nciccbr/ccbr_ubuntu_base_20.04:v1.0


RUN apt-get update 

RUN apt install -y openjdk-8-jdk unzip

# COPY fastqc_v0.11.9.zip /opt2
ADD https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip /opt2/
COPY Dockerfile /opt2

WORKDIR /opt2
RUN unzip fastqc_v0.11.9.zip
WORKDIR /opt2/FastQC
RUN chmod 755 fastqc
RUN ln -s /opt2/FastQC/fastqc /usr/local/bin/fastqc

WORKDIR /data2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
