# TRUST4/v0.1.0 Dockerfile (tag: v0.0.1)

# Overview of TRUST4 Dependencies
# 1. pthreads (bundled with OS)
# 2. build-essential
# 3. samtools
# Installation TLDR
# 1. git clone https://github.com/liulab-dfci/TRUST4.git
# 2. Run make in the repo directory
# 3. Add the directory of TRUST4 to the environment variable $PATH

FROM ubuntu:18.04

MAINTAINER Skyler Kuhn <kuhnsa@nih.gov>

RUN mkdir -p /data2
RUN mkdir -p /opt2

WORKDIR /opt2

# Update apt-get before downloading packages
RUN apt-get update && \
    apt-get upgrade -y

# Download packages
# Install Build and runtime dependencies for Samtools and HTSlib
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
    build-essential \
    git-all \
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
    liblzma-dev \
    locales && \
apt-get clean && apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install TRUST4
RUN git clone https://github.com/liulab-dfci/TRUST4.git && \
    cd TRUST4/ && \
    make

# Set environment variable(s)
# Configure "locale", see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8

# Add run-trust4 to PATH
ENV PATH="/opt2/TRUST4":$PATH

# Clean-up Image
RUN apt-get clean && apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Copy the Dockerfile used to create image in /opt2
COPY Dockerfile /opt2
WORKDIR /data2
