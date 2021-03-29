FROM ubuntu:20.04
RUN mkdir -p /opt2 && mkdir -p /data2
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update
RUN apt-get install -y cpanminus unzip zlibc zlib1g zlib1g-dev bzip2 make gcc g++ wget git 
RUN cpanm FindBin Term::ReadLine
# Set the locale
RUN apt-get install -y locales && localedef -i en_US -f UTF-8 en_US.UTF-8

#install pip
RUN apt-get install -y python3-pip

ADD Dockerfile /opt2/Dockerfile.baseimage
RUN chmod -R a+rX /opt2
WORKDIR /data2

RUN apt-get clean
MAINTAINER vishal.koparde@nih.gov
