FROM ubuntu:16.04

# installing cpu version of guppy basecaller for ONT

RUN mkdir -p /data2
RUN mkdir -p /opt2

RUN apt-get update && \
    apt-get upgrade -y 


# install python 3.7
RUN apt-get install -y software-properties-common && add-apt-repository ppa:deadsnakes/ppa
RUN apt update && \
  apt install -y python3.7 python3.7-distutils

WORKDIR /opt2
RUN apt-get install -y curl && \
  curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
  python3.7 get-pip.py && \
  rm -f get-pip.py 
RUN apt-get install -y build-essential zlibc zlib1g zlib1g-dev bzip2 libbz2-dev libncurses-dev libcurl4-openssl-dev liblzma-dev python3.7-dev

RUN pip install pycoQC

COPY Dockerfile /opt2

RUN chmod a+r /opt2/Dockerfile

WORKDIR /data2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
