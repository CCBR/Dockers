FROM continuumio/miniconda3

RUN mkdir -p /data2
RUN mkdir -p /opt2

RUN apt-get update && apt-get install -y build-essential zlib1g-dev zlib1g

# INSTALL TOBIAS 
RUN pip install tobias

# # INSTALL CLEAR
WORKDIR /opt2

RUN chmod -R a+rX /opt2
WORKDIR /data2

MAINTAINER vishal.koparde@nih.gov
