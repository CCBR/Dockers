FROM ubuntu:16.04

# installing cpu version of guppy basecaller for ONT

RUN mkdir -p /data2
RUN mkdir -p /opt2

RUN apt-get update && \
    apt-get upgrade -y 


# install python 3.8
RUN apt-get install -y software-properties-common && add-apt-repository ppa:deadsnakes/ppa
RUN apt update && \
  apt install -y python3.8 python3.8-distutils

WORKDIR /opt2
RUN apt-get install -y curl && \
  curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
  python3.8 get-pip.py && \
  rm -f get-pip.py && \
  pip install qcat 


COPY Dockerfile /opt2

RUN chmod a+r /opt2/Dockerfile

WORKDIR /data2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
