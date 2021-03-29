FROM ubuntu:18.04
RUN mkdir -p /opt2 && mkdir -p /data2
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update
RUN apt-get install -y cpanminus unzip zlibc zlib1g zlib1g-dev bzip2 make gcc g++ wget git 
RUN cpanm FindBin Term::ReadLine
# Set the locale
RUN apt-get install -y locales && localedef -i en_US -f UTF-8 en_US.UTF-8
#ENV LANG en_US.UTF-8  
#ENV LANGUAGE en_US:en  
#ENV LC_ALL en_US.UTF-8
#ENV LC_CTYPE en_US.UTF-8


# install python 3.7
RUN apt-get install -y software-properties-common && add-apt-repository ppa:deadsnakes/ppa
RUN apt update && \
  apt install -y python3.7 python3.7-distutils python3.7-dev
RUN if [ -f /usr/bin/python ];then rm -f /usr/bin/python;fi && ln -s /usr/bin/python3.7 /usr/bin/python
RUN if [ -f /usr/bin/python3 ];then rm -f /usr/bin/python3;fi && ln -s /usr/bin/python3.7 /usr/bin/python3
# install pip
WORKDIR /opt2
RUN apt-get install -y curl && \
  curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
  python3.7 get-pip.py && \
  rm -f get-pip.py 


ADD Dockerfile /opt2/Dockerfile.baseimage
RUN chmod -R a+rX /opt2
WORKDIR /data2

MAINTAINER vishal.koparde@nih.gov
