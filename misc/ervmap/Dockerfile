FROM ubuntu:16.04

# Tools included:
#
# bbmap           38.86
# bedtools        2.25.0
# bowtie2         2.2.6
# btrim           0.3.0 (from http://graphics.med.yale.edu/trim/)
# bwa             0.7.12-r1039
# cufflinks       2.2.1
# cutadapt        1.18
# samtools        0.1.19-96b5f2294a
# tophat          2.1.0
#

RUN mkdir -p /data2
RUN mkdir -p /opt2

# Need to add cloud repo to install R v3.6 on ubuntu v16.04
# else default is R v3.2.3 which does not have required libraries like RCurl
RUN apt-get update && \
    apt-get upgrade -y

# add btrim
# ADD http://graphics.med.yale.edu/trim/btrim64-static /opt2/btrim
COPY btrim64-static /opt2/btrim
RUN chmod a+x /opt2/btrim

RUN apt-get install -y bedtools samtools tophat bwa python

# add cuffdiff
COPY cufflinks-2.2.1.Linux_x86_64.tar.gz /opt2
WORKDIR /opt2
RUN tar xzvf cufflinks-2.2.1.Linux_x86_64.tar.gz && cd cufflinks-2.2.1.Linux_x86_64 && chmod a+x *
ENV PATH="/opt2/cufflinks-2.2.1.Linux_x86_64/:${PATH}"

# add ERVmap
# ERVmap requires HTSeq so installing that first
RUN apt-get install -y build-essential python2.7-dev python-numpy python-matplotlib python-pysam python-htseq
WORKDIR /opt2
RUN apt-get install -y git cpanminus && \
    git clone https://github.com/mtokuyama/ERVmap.git && \
    cd ERVmap && \
    chmod a+x *.pl && \
    cpanm File::Type
ENV PATH="/opt2/ERVmap/:${PATH}"


# add bbmap ... will be used for interleaving fastq files
RUN apt-get install -y default-jre
# curl -L "https://sourceforge.net/projects/bbmap/files/BBMap_38.86.tar.gz/download" > BBMap.tar.gz
COPY BBMap.tar.gz /opt2
WORKDIR /opt2
RUN tar xzvf BBMap.tar.gz
ENV PATH="/opt2/bbmap:${PATH}"


# INSTALL cutadapt
# INSTALL miniconda

WORKDIR /opt2
# COPY Miniconda3-4.7.12.1-Linux-x86_64.sh /opt2/miniconda.sh
ADD https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh /opt2/miniconda.sh
RUN /bin/bash /opt2/miniconda.sh -b -p /opt2/conda && \
    rm /opt2/miniconda.sh && \
    ln -s /opt2/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt2/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc
RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda create -n python3 python=3.6
RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  pip install --upgrade pip && \
  pip install cutadapt==1.18 && \
  pip install numpy && \
  pip install matplotlib===3.1.1



COPY Dockerfile /opt2
COPY TruSeq_and_nextera_adapters.consolidated.fa /opt2
COPY ccbr_trim_interleave_pe.bash /opt2
COPY argparse.bash /opt2
RUN apt-get clean && \
chmod a+r /opt2/TruSeq_and_nextera_adapters.consolidated.fa && \
chmod a+r /opt2/Dockerfile && \
chmod a+rx /opt2/argparse.bash && \
chmod a+rx /opt2/ccbr_trim_interleave_pe.bash 

ENV PATH="/opt2:${PATH}"

WORKDIR /data2


MAINTAINER vishal.koparde@nih.gov
