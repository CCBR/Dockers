# rNA/v0.1.0 Dockerfile (tag: v0.0.1)

# Overview of rNA Dependencies
# Python >= 3
#    Python Packages <optional>
#       argparse, numpy==1.18.5, pandas==0.25.3, python-dateutil==2.8.1, 
#       pytz==2020.1, six==1.15.0, xlrd==1.2.0
# R >= 3.5
#   R Packages <required>
#       CRAN: 'plyr', 'plotly', 'ggplot2', 'RColorBrewer', 'gridExtra', 
#             'crosstalk', 'DT', 'reshape2', 'circlize', 'flexdashboard'
#             'knitr', 'rmarkdown', 'argparse', 'viridis'
#       Bioconductor: 'limma', 'edgeR', 'ComplexHeatmap'  

FROM ubuntu:18.04

MAINTAINER Skyler Kuhn <kuhnsa@nih.gov>

RUN mkdir -p /data2
RUN mkdir -p /opt2

WORKDIR /opt2

# Update apt-get before downloading packages
RUN apt-get update && \
    apt-get upgrade -y

# Download packages
# Install Python 3.6.5-3, build and runtime dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
    python3 \
    python3-pip \
    build-essential \
    make \
    git \
    gcc \
    g++ \
    ca-certificates \
    libcurl4-openssl-dev \
    libxml2-dev \
    wget \
    zlibc \
    zlib1g \
    zlib1g-dev \
    libssl-dev \
    locales \
    pandoc \
    software-properties-common

# Install Python Packages
RUN git clone https://github.com/CCBR/rNA.git && \
    cd rNA/ && \
    python3 -m pip install --upgrade pip && \
    python3 -m pip install -r requirements.txt


# Install R (3.6) -- default through apt-get is 3.4.4 (edgeR needs 3.6) -- and R packages
# For more information, check out: https://cran.r-project.org/bin/linux/ubuntu/
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 51716619E084DAB9
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/"
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt install --yes r-base

# Install Required R packages
RUN Rscript -e 'install.packages(c("argparse", "knitr", "plyr", "plotly", "ggplot2", "RColorBrewer"), repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages(c("shiny", "gridExtra", "flexdashboard", "rmarkdown", "crosstalk", "DT", "reshape2", "circlize", "viridis"), repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("limma", "edgeR", "ComplexHeatmap"))'

# PLEASE NOTE: The bioconductor Package, Repitools, expects a cache directory hold cache files
# It is convenient to use  '/root/.cache/R/R.cache' because it follows the standard
# on your operating system. If not, a temporary directory '/tmp/RtmpyUY1Zs/.Rcache'
# that is specific to this R session will be used. This may cause problems on with
# docker/singularity as the container filesystem is read-only. To over-come any potential
# issues the /tmp directory may need to be added to the container's bind PATH.

# Set environment variable(s)
# Configure "locale", see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8

# Add rNA to PATH
ENV PATH="/opt2/rNA":$PATH

# Clean-up Image
RUN apt-get clean && apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Copy the Dockerfile used to create image in /opt2
COPY Dockerfile /opt2
WORKDIR /data2
