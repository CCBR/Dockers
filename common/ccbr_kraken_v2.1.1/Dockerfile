# ccbr_kraken_v2.1.1 Dockerfile (v0.0.1)
FROM ubuntu:18.04

MAINTAINER Skyler Kuhn <kuhnsa@nih.gov>

RUN mkdir -p /data2
RUN mkdir -p /opt2

# Update apt-get before downloading packages 
RUN apt-get update && \
    apt-get upgrade -y 

# Install dependencies for latest versions of Kraken (v2.1.1) and Krona (v.2.7.1)
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
	build-essential \
	wget \
	curl \
	zlib1g-dev \
	make \
	python \
	perl \
	g++ \
	locales && \
apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /opt2

# Get latest kraken release (v2.1.1): https://github.com/DerrickWood/kraken2/releases 
RUN wget https://github.com/DerrickWood/kraken2/archive/v2.1.1.tar.gz && \
	tar -xzf v2.1.1.tar.gz && \
	rm v2.1.1.tar.gz && \
	cd kraken2-2.1.1 && \
 	./install_kraken2.sh . && \
	chmod a+rX /opt2/kraken2-2.1.1/*

# Add kraken2 to PATH
ENV PATH="/opt2/kraken2-2.1.1":$PATH

WORKDIR /opt2

# Get latest krona release (v2.7.1): https://github.com/marbl/Krona/releases
# Build Krona Taxonomy Databases
# ktImportTaxonomy relies on NCBI taxonomy: https://github.com/marbl/Krona/wiki/Installing 
RUN wget https://github.com/marbl/Krona/archive/v2.7.1.tar.gz && \
	tar -xzf v2.7.1.tar.gz && \
	rm v2.7.1.tar.gz && \ 
	cd Krona-2.7.1/KronaTools && \
	./install.pl  && \
	mkdir -p /opt2/Krona-2.7.1/KronaTools/taxonomy && \
	./updateTaxonomy.sh

# Set environment variable(s)
# Configure "locale", see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8

# Dockerfile to build image located at /opt2/Dockerfile
COPY Dockerfile /opt2/
WORKDIR /data2
