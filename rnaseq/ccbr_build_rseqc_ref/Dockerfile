# ccbr_build_rseqc_ref Dockerfile (v0.0.1)
FROM ubuntu:18.04

MAINTAINER <vishal.koparde@nih.gov>, <kuhnsa@nih.gov>

RUN mkdir -p /data2
RUN mkdir -p /opt2/ucsc/

# Update apt-get before downloading packages 
RUN apt-get update && \
    apt-get upgrade -y 

# Install dependencies: python, wget to download UCSC binaries   
# samtools to install libkrb5.so.3: open shared object file (needed bby UCSC tools)
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
	wget \
	python \
        bedtools \
        samtools \
	locales && \
apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /opt2/ucsc

# Get a controlled release of UCSC tool from my mirror: https://github.com/skchronicles/ucsc_mirror/ 
RUN wget https://github.com/skchronicles/ucsc_mirror/raw/main/bedSort && \
	wget https://github.com/skchronicles/ucsc_mirror/raw/main/genePredToBed && \
	wget https://github.com/skchronicles/ucsc_mirror/raw/main/gtfToGenePred && \
	chmod a+rx /opt2/ucsc/*

# Set environment variable(s)
# Configure "locale", see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8

# Add UCSC tools to PATH
ENV PATH="/opt2/ucsc":$PATH

COPY Dockerfile /opt2/
WORKDIR /data2
