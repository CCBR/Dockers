FROM ubuntu:18.04
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update 

# install python 3.7
RUN apt-get install -y software-properties-common && add-apt-repository ppa:deadsnakes/ppa
RUN apt update && \
  apt install -y python3.7 python3.7-distutils

WORKDIR /opt2
RUN apt-get install -y curl && \
  curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
  python3.7 get-pip.py && \
  rm -f get-pip.py 

RUN pip install HTSeq

RUN mkdir -p /opt2 && mkdir -p /data2
RUN if [ -f /usr/bin/python ];then rm -f /usr/bin/python;fi && ln -s /usr/bin/python3.7 /usr/bin/python
RUN if [ -f /usr/bin/python3 ];then rm -f /usr/bin/python3;fi && ln -s /usr/bin/python3.7 /usr/bin/python3

ADD Dockerfile /opt2
RUN chmod -R a+rX /opt2
WORKDIR /data2

MAINTAINER vishal.koparde@nih.gov
