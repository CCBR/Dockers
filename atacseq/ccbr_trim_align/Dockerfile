FROM ubuntu:16.04

# Tools installed in this docker image:
# bedGraphToBigWig
# bedSort
# bedtools     2.25.0
# bowtie2      2.2.6
# bwa          0.7.12
# cutadapt     1.18       conda:python3
# fastqc       0.11.8
# idr          2.0.4.2    conda:python3
# jre          1.8.0_222
# macs2        2.2.5      conda:python3
# miniconda3   4.7.12.1
# picardcloud
# samtools     0.1.19

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y

RUN apt-get install -y cpanminus unzip

RUN cpanm FindBin Term::ReadLine

RUN apt-get install -y default-jre \
    bedtools \
    bowtie2 \
    samtools \
    bwa \
    wget \
    pigz \
    gcc \
    g++ \
    make \
    cmake \
    bc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev

# INSTALL miniconda

WORKDIR /opt
COPY Miniconda3-4.7.12.1-Linux-x86_64.sh /opt/miniconda.sh
# downloaded from https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh
RUN /bin/bash /opt/miniconda.sh -b -p /opt/conda && \
    rm /opt/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# CREATE MINICONDA PYTHON3 ENVIRONMENT

RUN . /opt/conda/etc/profile.d/conda.sh && \
  conda create -n python3 python=3.6 && \
  conda create -n python2 python=2.7

# INSTALL packages in python2
# bam2fld scripts requires older version of pysam ... newer version gives errors ...hence installing pysam in python2 and python3
RUN . /opt/conda/etc/profile.d/conda.sh && \
  conda activate python2 && \
  pip install --upgrade pip && \
  pip install argparse && \
  pip install numpy && \
  pip install pysam && \
  pip install scipy


# INSTALL cutadapt

RUN . /opt/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  pip install --upgrade pip && \
  pip install cutadapt==1.18 && \
  pip install numpy && \
  pip install matplotlib

# INSTALL idr

WORKDIR /opt
# conda config --add channels conda-forge && \
RUN . /opt/conda/etc/profile.d/conda.sh && \
    conda activate python3 && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda install idr && \
    conda update idr

# INSTALL macs2

RUN . /opt/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  pip install macs2==2.2.5

# ADD pysam

RUN . /opt/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  conda config --add channels r && \
  conda install pysam

# INSTALL fastqc

WORKDIR /opt
COPY fastqc_v0.11.8.zip /opt
# downloaded from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
RUN unzip fastqc_v0.11.8.zip && \
  cd /opt/FastQC && \
  chmod 755 fastqc && \
  ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc


COPY TruSeq_and_nextera_adapters.consolidated.fa /opt
COPY Dockerfile /opt
COPY picardcloud.jar /opt
COPY bedSort /opt
COPY bedGraphToBigWig /opt
COPY argparse.bash /opt
COPY the_blacklists.fa /opt
COPY the_blacklists.fa.amb /opt
COPY the_blacklists.fa.ann /opt
COPY the_blacklists.fa.bwt /opt
COPY the_blacklists.fa.pac /opt
COPY the_blacklists.fa.sa /opt

COPY atac_assign_multimappers.py /opt
COPY ccbr_atac_trim_align_pe.bash /opt
COPY ccbr_bam_filter_by_mapq.py /opt
COPY ccbr_bowtie2_align_pe.bash /opt
COPY ccbr_cutadapt_pe.bash /opt
COPY ccbr_macs2_peak_calling.bash /opt
COPY ccbr_remove_blacklisted_reads_pe.bash /opt

RUN chmod a+x /opt/bedSort
RUN chmod a+x /opt/bedGraphToBigWig
RUN chmod a+x /opt/*.bash
RUN chmod a+x /opt/*.py
RUN chmod a+r /opt/*.*

ENV PATH="/opt:${PATH}"

WORKDIR /data

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
