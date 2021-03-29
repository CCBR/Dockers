FROM ubuntu:16.04

# Tools installed in this docker image:
# gimmemotifs     0.14       conda:python3

RUN mkdir -p /data2
RUN mkdir -p /opt2

RUN apt-get update && \
    apt-get upgrade -y

RUN apt-get install -y \
  autoconf \
  automake \
  build-essential \
  cython \
  gfortran \
  ghostscript \
  ghostscript \
  libatlas-base-dev \
  libatlas-base-dev \
  libbz2-dev \
  libexpat1-dev \
  libexpat1-dev \
  libfreetype6-dev \
  libgd-dev \
  libgs-dev \
  libgsl0-dev \
  libhtml-template-compiled-perl \
  liblapack-dev \
  libopenmpi-dev \
  libtool \
  libxml-libxml-debugging-perl \
  libxml-opml-simplegen-perl \
  libxml2-dev \
  libxslt-dev \
  openmpi-bin \
  openssh-server \
  sudo \
  libpng12-dev \
  zlib1g-dev



# INSTALL miniconda

WORKDIR /opt2
COPY Miniconda3-4.7.12.1-Linux-x86_64.sh /opt2/miniconda.sh
# ADD https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh /opt2
RUN /bin/bash /opt2/miniconda.sh -b -p /opt2/conda && \
    rm /opt2/miniconda.sh && \
    ln -s /opt2/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt2/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# CREATE MINICONDA PYTHON3 ENVIRONMENT

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda config --add channels defaults && \
  conda config --add channels bioconda && \
  conda config --add channels conda-forge && \
  conda create -n gimme python=3.6 gimmemotifs 

COPY Dockerfile /opt2
COPY argparse.bash /opt2

RUN chmod a+rX /opt2/Dockerfile
RUN chmod a+rX /opt2/argparse.bash

WORKDIR /data2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
