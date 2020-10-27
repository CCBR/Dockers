# cDNA_Cupcake Dockerfile
FROM ubuntu:18.04

MAINTAINER Skyler Kuhn <kuhnsa@nih.gov>

RUN mkdir -p /data2
RUN mkdir -p /opt2

# Apt-get packages to build python (3.8) from source
RUN apt-get update && apt-get -y upgrade
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
	build-essential \
	zlib1g-dev \
	libncurses5-dev \
	libgdbm-dev \
	libnss3-dev \
	libssl-dev \
	libreadline-dev \
	libffi-dev \
	libsqlite3-dev \
	wget \
	libbz2-dev \
	locales \
	git && \
apt-get clean && apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /opt2

# Cupcake needs python>=3.7
# Build python 3.8 from source
RUN wget https://www.python.org/ftp/python/3.8.0/Python-3.8.0.tgz && \
	tar -xvf Python-3.8.0.tgz && \
	rm Python-3.8.0.tgz && \
	cd Python-3.8.0/ && \
	./configure && \
	make && \
	make altinstall

WORKDIR /opt2

# Install cDNA_Cupcake from Github
RUN git clone https://github.com/Magdoll/cDNA_Cupcake.git && \
	cd cDNA_Cupcake/ && \
	python3.8 -m pip install --upgrade pip && \
	python3.8 -m pip install -r requirements.txt && \
	python3.8 -m pip install Cython

# Create symlinks to Python interpreters
RUN ln -s /usr/local/bin/python3.8 /usr/local/bin/python3
RUN ln -s /usr/local/bin/python3.8 /usr/local/bin/python

# Build cDNA_Cupcake and add as site package
WORKDIR /opt2/cDNA_Cupcake/
RUN python3.8 setup.py build
RUN python3.8 setup.py install
ENV PATH=$PATH:/opt2/cDNA_Cupcake/sequence/

# Set environment variable(s)
# Configure "locale", see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
  && locale-gen en_US.utf8 \
  && /usr/sbin/update-locale LANG=en_US.UTF-8

# Copy the Dockerfile used to create image in /opt2
COPY Dockerfile /opt2
WORKDIR /data2
