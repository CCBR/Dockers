FROM ubuntu:16.04

RUN mkdir -p /data2
RUN mkdir -p /opt2

RUN apt-get update && \
    apt-get upgrade -y

RUN apt-get install -y \
build-essential 

WORKDIR /opt2
COPY Miniconda3-4.7.12.1-Linux-x86_64.sh /opt2/miniconda.sh
# downloaded from https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh
RUN /bin/bash /opt2/miniconda.sh -b -p /opt2/conda && \
    rm /opt2/miniconda.sh && \
    ln -s /opt2/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt2/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc


RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda create -n python3 python=3.7 && \
  conda install -c anaconda pip 

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  pip install --upgrade pip && \
  pip install DeepTCR 

WORKDIR /data2
MAINTAINER vishal.koparde@nih.gov
