# QualiMap/v2.2.1 Dockerfile (tag: v0.0.1)
FROM ubuntu:18.04

MAINTAINER Skyler Kuhn <kuhnsa@nih.gov>

# Here is QualiMap's installation documentation needs: http://gensoft.pasteur.fr/docs/qualimap/v2.2.1/QualimapManual.pdf
# 1. JAVA runtime > 6
# 2. R enviroment > 3.1
# 3. R Packages: optparse (CRAN), NOISeq (Bioconductor), Repitools (Bioconductor), Rsamtools (Bioconductor), GenomicFeatures (Bioconductor), rtracklayer (Bioconductor)

RUN mkdir -p /data2
RUN mkdir -p /opt2

WORKDIR /opt2

# Update apt-get before downloading packages
RUN apt-get update && \
    apt-get upgrade -y

# Download packages
# Install Java, build and runtime dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
      wget \
      unzip \
      openjdk-8-jre \
      locales \
      software-properties-common \
      build-essential \
      libcurl4-openssl-dev \
      libxml2-dev \
      libssl-dev && \
    apt-get clean && apt-get purge && \
      rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Qualimap
RUN wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip && \
    unzip qualimap_v2.2.1.zip && \
    rm qualimap_v2.2.1.zip

# Install R (3.6) -- default through apt-get is 3.4.4 -- and R packages
# For more information, check out: https://cran.r-project.org/bin/linux/ubuntu/
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 51716619E084DAB9
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/"
RUN apt update && DEBIAN_FRONTEND=noninteractive apt install --yes r-base
RUN Rscript -e "install.packages('optparse')"
RUN Rscript -e "install.packages('XML')"
RUN Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install(c('Repitools', 'GenomicFeatures'))"
RUN Rscript -e "BiocManager::install(c('NOISeq', 'Rsamtools', 'rtracklayer'))"

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

# Add QualiMap to PATH
ENV PATH="/opt2/qualimap_v2.2.1":$PATH

# Copy the Dockerfile used to create image in /opt2
COPY Dockerfile /opt2
WORKDIR /data2
