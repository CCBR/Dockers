# TRUST4 Dockerfile
# Define release with --build-arg TRUST4_VERSION=$RELEASE_VERSION
# See README for an Example

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

# Setting build time arg for defining a release
# Default = 1.0.2-beta   
ARG TRUST4_VERSION=1.0.2-beta

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
    wget \
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

# Install Defined Release of TRUST4
# Release Archive also contains the example datasets 
RUN wget https://github.com/liulab-dfci/TRUST4/archive/v${TRUST4_VERSION}.tar.gz && \
    tar -xzf v${TRUST4_VERSION}.tar.gz && \
    rm v${TRUST4_VERSION}.tar.gz && \
    mv TRUST4-${TRUST4_VERSION} TRUST4 && \
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
