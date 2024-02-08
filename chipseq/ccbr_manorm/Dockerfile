FROM ubuntu:20.04

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

RUN mkdir -p /opt2 && mkdir -p /data2
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update && apt-get -y upgrade
# Set the locale
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
		locales build-essential cmake cpanminus && \
	localedef -i en_US -f UTF-8 en_US.UTF-8 && \
	cpanm FindBin Term::ReadLine

# install basic dependencies with apt-get
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential \
	figlet \
	g++ \
	gcc \
	gfortran \
	git \
	libatlas-base-dev \
	libblas-dev \
	libboost-dev \
	libbz2-dev \
	libcurl4-openssl-dev \
	libexpat1-dev \
	libfreetype6-dev \
	libgd-dev \
	libgd-perl \
	libglib2.0-dev \
    libgpgme11-dev \
	libgs-dev \
	libgsl-dev \
	libgsl0-dev \
	libhtml-template-compiled-perl \
	libicu-dev \
	libjudy-dev \
	liblapack-dev \
	liblzma-dev \
	libmysqlclient-dev \
	libncurses-dev \
	libopenmpi-dev \
	libpng-dev \
	librtmp-dev \
    libseccomp-dev \
	libssl-dev \
	libtool \
	libxml-libxml-debugging-perl \
	libxml-opml-simplegen-perl \
	libxml2-dev \
	libxslt-dev \
	make \
	manpages-dev \
	openjdk-17-jre-headless \
	parallel \
	pigz \
    pkg-config \
	python3-pip \
    python3-dev \
	rsync \
    squashfs-tools \
	unzip \
    uuid-dev \
	wget \
	zlib1g \
	zlib1g-dev \
	zlibc

# Install conda and give write permissions to conda folder
RUN echo 'export PATH=/opt2/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" -O ~/miniforge3.sh && \
    /bin/bash ~/miniforge3.sh -b -p /opt2/conda && \
    rm ~/miniforge3.sh && chmod 777 -R /opt2/conda/
ENV PATH="/opt2/conda/bin:$PATH"

# install packages from conda
COPY environment.txt /data2/
RUN mamba install -c bioconda -c conda-forge --file /data2/environment.txt

# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
