# BBtools/38.87 Dockerfile (v0.0.1)
FROM ubuntu:18.04

MAINTAINER Skyler Kuhn <kuhnsa@nih.gov>

# BBTools requires `java >= 7` to be installed on the target system. 
# In addition, it is recommended install the following programs: 
#   `samtools` (sam to bam conversion)
#   `pigz` to speed up compression
#   `unzip` to support zip files
#   `bzip2` to support bzip files  


RUN mkdir -p /data2
RUN mkdir -p /opt2/bin

WORKDIR /opt2

# Update apt-get before downloading packages
RUN apt-get update && \
    apt-get upgrade -y

# Download packages
# Install Java, Samtools/1.7-1, build and runtime dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
      groff \
      less \
      man \
      pigz \
      wget \
      unzip \ 
      bzip2 \
      samtools \
      openjdk-8-jre \
      locales && \
    apt-get clean && apt-get purge && \
      rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install BBtools
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_38.87.tar.gz && \
    tar -xzf BBMap_38.87.tar.gz && \
    rm BBMap_38.87.tar.gz

# Set environment variable(s)
# Configure "locale", see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8

# Add BBtools to PATH
ENV PATH="/opt2/bbmap":$PATH

# Copy the Dockerfile used to create image in /opt2
COPY Dockerfile /opt2
COPY bin/bbtools /opt2/bin/
RUN chmod a+rx /opt2/bin/bbtools
ENV PATH="/opt2/bin":$PATH
WORKDIR /data2
