FROM ubuntu:18.04

MAINTAINER Skyler Kuhn <kuhnsa@nih.gov>

RUN mkdir -p /data2
RUN mkdir -p /opt2

# Apt-get packages
# Install python (3.6), bowtie2=2.3.4.1-1 (ubuntu:18.04 default)
RUN apt-get update && apt-get -y upgrade
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
    build-essential \
    apt-utils \
    git-all \
    python3 \ 
    python3-pip \
    bowtie2 \
    wget

WORKDIR /opt2

# Build Samtools 1.9, Telescope requires htslib=1.9 (must install SAMtools=1.9)
# SAMtools installation information: https://github.com/samtools/samtools/blob/develop/INSTALL
# HTSlib installation information: https://github.com/samtools/htslib/blob/1.9/INSTALL
# Apt-get remaining dependencies  
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
    gcc \
    make \
    perl \
    bzip2 \
    zlibc \
    libssl-dev \
    libbz2-dev \
    zlib1g-dev \
    libncurses5-dev \ 
    libncursesw5-dev \
    libcurl4-gnutls-dev \
    liblzma-dev && \
apt-get clean && apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Build SAMtools 1.9
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xjvf samtools-1.9.tar.bz2 && \
    rm samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    ./configure --prefix $(pwd) && \
    make

# Build HTSlib 1.9 (required by telescope)
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -vxjf htslib-1.9.tar.bz2 && \
    rm htslib-1.9.tar.bz2 && \
    cd htslib-1.9 && \
    ./configure --prefix $(pwd) && \
    make

# Add SAMtools and HTSlib to PATH
ENV PATH=${PATH}:/opt2/samtools-1.9
ENV PATH=${PATH}:/opt2/htslib-1.9
ENV HTSLIB_INCLUDE_DIR="/opt2/htslib-1.9"

# pip install: Cutadapt, Telescope dependencies, and then Telescope
RUN pip3 install --upgrade pip
RUN pip3 install cutadapt==2.10 

# Python requirements from github page some are need for compiling telescope, installing now
RUN pip3 install future pyyaml cython==0.29.7 numpy==1.16.3 scipy==1.2.1 pysam==0.15.2 intervaltree==3.0.2
RUN pip3 install git+git://github.com/mlbendall/telescope.git 

# Adpater sequences for cutadapt and HERV reference files
RUN mkdir -p /opt2/refs
COPY refs/trimmonatic_TruSeqv3_adapters.fa /opt2/refs
COPY refs/HERV_rmsk.hg38.v2.genes.gtf /opt2/refs
COPY refs/HERV_rmsk.hg38.v2.transcripts.gtf /opt2/refs

# Set environment variable(s)
ENV LANG en_US.UTF-8

# Copy the Dockerfile used to create image in /opt2
COPY Dockerfile /opt2
WORKDIR /data2
