FROM ubuntu:16.04


RUN mkdir -p /data2
RUN mkdir -p /opt2

# Need to add cloud repo to install R v3.6 on ubuntu v16.04
# else default is R v3.2.3 which does not have required libraries like RCurl
RUN apt-get update && \
    apt-get upgrade -y

RUN apt-get install -y libxml2-dev software-properties-common apt-transport-https
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" >> /etc/apt/sources.list
# RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -s -c) main"

RUN apt-get update && \
    apt-get upgrade -y --allow-unauthenticated

RUN apt-get install --fix-missing -y --allow-unauthenticated \
    cmake \
    default-jre \
    g++ \
    gcc \
    make \
    r-base

RUN apt-get install -y --allow-unauthenticated \
	libcurl4-gnutls-dev \
	librtmp-dev \
	libssl-dev

# INSTALL R packages

RUN apt-get install -y --allow-unauthenticated r-base-core
COPY install_R_packages.R /opt2
WORKDIR /opt2
RUN Rscript install_R_packages.R

# INSTALL miniconda

# WORKDIR /opt2
# COPY Miniconda3-4.7.12.1-Linux-x86_64.sh /opt2/miniconda.sh

# downloaded from https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh

# RUN /bin/bash /opt2/miniconda.sh -b -p /opt2/conda && \
#     rm /opt2/miniconda.sh && \
#     ln -s /opt2/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
#     echo ". /opt2/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
#     echo "conda activate base" >> ~/.bashrc

# CREATE MINICONDA PYTHON3 ENVIRONMENT

# RUN . /opt2/conda/etc/profile.d/conda.sh && \
#   conda create -n python3 python=3.6 





# SCRIPTS AND DB files

COPY Dockerfile /opt2
COPY argparse.bash /opt2
RUN mkdir -p /opt2/db
COPY BlueprintEncodeData.rds /opt2/db
COPY DatabaseImmuneCellExpressionData.rds /opt2/db
COPY HumanPrimaryCellAtlasData.rds /opt2/db
COPY ImmGenData.rds /opt2/db
COPY MonacoImmuneData.rds /opt2/db
COPY MouseRNAseqData.rds /opt2/db
COPY newNovershternHematopoieticData.rds /opt2/db


# python scripts

# bash scripts

# R scripts

COPY singleR_standalone.0.0.3.R /opt2/singleR_standalone.R

RUN chmod a+r /opt2/*.*
RUN chmod a+r /opt2/db/*.rds
RUN chmod a+r Dockerfile

ENV PATH="/opt2:${PATH}"

WORKDIR /data2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
