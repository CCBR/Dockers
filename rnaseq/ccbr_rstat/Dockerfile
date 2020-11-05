# rstat Dockerfile (v0.0.1)
FROM ubuntu:18.04

MAINTAINER kuhnsa@nih.gov

RUN mkdir -p /data2
RUN mkdir -p /opt2/picard/

# Update apt-get before downloading packages 
RUN apt-get update && \
    apt-get upgrade -y 

# Install Standard Java Runtime, python/3.6.5-3, pip3, samtools/1.7-1
RUN apt-get install -y \
	openjdk-8-jre \
	samtools \
	python3 \
	python3-pip \
	locales \
	wget && \
apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install pysam from pypi
RUN python3.6 -m pip install --upgrade pip && \
	python3.6 -m pip install pysam

# Picard Cloud
WORKDIR /opt2/picard

# Add Cloud Picard JAR release
COPY resources/picardcloud.jar /opt2/picard/picard.jar
RUN chmod a+x /opt2/picard/*.jar
ENV PICARDJARPATH=/opt2/picard


# Set environment variable(s)
# Configure "locale", see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
  && locale-gen en_US.utf8 \
  && /usr/sbin/update-locale LANG=en_US.UTF-8


COPY Dockerfile /opt2/
WORKDIR /data2
