FROM ubuntu:16.04

RUN mkdir -p /data2
RUN mkdir -p /opt2

RUN apt-get update && \
    apt-get upgrade -y

RUN apt-get install -y \
build-essential \
inkscape \
wget \
freetype2-demos \
libpng-dev \
libfreetype6-dev \
libx11-dev

RUN ln -s /usr/include/freetype2/ft2build.h /usr/include/

WORKDIR /opt2
COPY Miniconda3-4.7.12.1-Linux-x86_64.sh /opt2/miniconda.sh
# downloaded from https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh
RUN /bin/bash /opt2/miniconda.sh -b -p /opt2/conda && \
    rm /opt2/miniconda.sh && \
    ln -s /opt2/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt2/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# CREATE MINICONDA PYTHON2 ENVIRONMENT

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda create -n python2 python=2.7

# INSTALL packages in python2
# bam2fld scripts requires older version of pysam ... newer version gives errors ...hence installing pysam in python2 and python3
# metaseq required for TSS score calculation is also installed in python2 conda environment
RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python2 && \
  pip install --upgrade pip && \
  pip install pkgconfig && \
  pip install freetype-py && \
  pip install scipy==0.16.0 && \
  pip install numpy==1.10.1 && \
  pip install scikit-learn==0.17 && \
  pip install matplotlib==1.4.3

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python2 && \
  pip install blargs

RUN apt-get install -y git
WORKDIR /opt2
RUN git clone https://github.com/phbradley/tcr-dist.git
WORKDIR /opt2/tcr-dist
RUN ls -alrth && pwd
RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python2 && \
  echo "Y" | python setup.py
RUN chmod -R a+rX /opt2/tcr-dist



# WORKDIR /opt2
# ADD https://bootstrap.pypa.io/get-pip.py /opt2/get-pip.py
# RUN python get-pip.py

# RUN pip install \
# pkgconfig \
# freetype-py \
# scipy==0.16.0 \
# scikit-learn==0.17 \
# matplotlib==1.4.3 \
# numpy==1.10.1

# WORKDIR /opt2
# RUN git clone https://github.com/phbradley/pubtcrs.git
# WORKDIR /opt2/pubtcrs
# RUN make && mv /opt2/pubtcrs/bin/* /opt2 && chmod a+x /opt2/*

# ENV PATH="/opt2/:${PATH}"

WORKDIR /data2
MAINTAINER vishal.koparde@nih.gov
