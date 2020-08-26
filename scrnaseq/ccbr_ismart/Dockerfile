FROM ubuntu:18.04

MAINTAINER Skyler Kuhn <kuhnsa@nih.gov>

RUN mkdir -p /data2
RUN mkdir -p /opt2

# Install python3
RUN apt-get update && apt-get -y upgrade
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
    build-essential \
    apt-utils \
    git-all \
    python \ 
    python-pip \
    wget && \
apt-get clean && apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# pip install: numpy, pandas, matplotlib, argparse
RUN pip install --upgrade pip
RUN pip install numpy pandas 
RUN pip install biopython==1.76

# Clone iSMART repo from Github
WORKDIR /opt2
RUN git clone https://github.com/s175573/iSMART.git
RUN chmod a+x /opt2/iSMART/*.py
ENV PATH=${PATH}:/opt2/iSMART

# Patch she-bang
RUN sed -i 's@^#!usr/bin/python@#!/usr/bin/env python@' /opt2/iSMART/iSMARTv3.py

# Set environment variable(s)
ENV LANG en_US.UTF-8

COPY Dockerfile /opt2
WORKDIR /data2
