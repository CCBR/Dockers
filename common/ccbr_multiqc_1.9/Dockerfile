# MultiQC/1.9 Dockerfile (v0.0.1)
FROM ubuntu:18.04

MAINTAINER Skyler Kuhn <kuhnsa@nih.gov>

RUN mkdir -p /data2
RUN mkdir -p /opt2/

# Update apt-get before downloading packages 
RUN apt-get update && \
    apt-get upgrade -y 

# Install Standard python/3.6.5-3, pip3, python3-pandas (pyparser dependency)
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
	python3 \
	python3-pip \
	python3-pandas \
	locales \
	wget && \
apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /opt2/

# Install MultiQC from pypi
RUN python3.6 -m pip install --upgrade pip && \
	python3.6 -m pip install multiqc==1.9

# Set environment variable(s)
# Configure "locale", see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
  && locale-gen en_US.utf8 \
  && /usr/sbin/update-locale LANG=en_US.UTF-8

# Resolves the a click Runtime (click is a dependency of MultiQC)
# Note: In Python 3.7 and later you will no longer get a RuntimeError in many cases
# thanks to PEP 538 and PEP 540, which changed the default assumption in unconfigured environments.
# RuntimeError: Click will abort further execution because Python 3 
# was configured to use ASCII as encoding for the environment.
# Consult https://click.palletsprojects.com/python3/ for mitigation steps.
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

COPY Dockerfile /opt2/
WORKDIR /data2
