FROM ubuntu:16.04

RUN mkdir -p /data2
RUN mkdir -p /opt2

# Need to add cloud repo to install R v3.6 on ubuntu v16.04
# else default is R v3.2.3 which does not have required libraries like RCurl
RUN apt-get update && \
    apt-get upgrade -y

RUN apt-get install -y libxml2-dev software-properties-common apt-transport-https
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" >> /etc/apt/sources.list
# RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -s -c) main"
RUN apt-get update && \
    apt-get upgrade -y --allow-unauthenticated

RUN apt-get install -y cpanminus unzip

RUN cpanm FindBin Term::ReadLine

RUN apt-get install -y \
    bc \
    bedtools \
    bowtie2 \
    bwa \
    cmake \
    default-jre \
    g++ \
    gcc \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libssl-dev \
    make \
    pandoc \
    parallel \
    pigz \
    pkg-config \
    samtools \
    wget \
    zlib1g-dev

# INSTALL R packages

RUN apt-get install -y --allow-unauthenticated r-base-core
COPY install_R_packages.R /opt2
WORKDIR /opt2
RUN Rscript install_R_packages.R

# INSTALL ChIPseeker
COPY install_R_packages_chipseeker.R /opt2
WORKDIR /opt2
RUN Rscript install_R_packages_chipseeker.R

