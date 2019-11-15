FROM ubuntu:16.04

# Tools installed in this docker image:
# bedSort
# bedtools     2.25.0
# bowtie2      2.2.6
# bwa          0.7.12
# cutadapt     1.18     conda:python3
# fastqc       0.11.8
# idr          2.0.2    conda:python3
# jre          1.8.0_222
# macs2        2.2.5    conda:python3
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
    cmake

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
  conda create -n python3 python=3.6

# INSTALL cutadapt

RUN . /opt/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  pip install --upgrade pip && \
  pip install cutadapt==1.18 && \
  pip install numpy && \
  pip install matplotlib

# INSTALL idr

WORKDIR /opt
COPY 2.0.2.zip /opt
# downloaded from https://github.com/nboley/idr/archive/2.0.2.zip
RUN . /opt/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  unzip 2.0.2.zip && \
  cd idr-2.0.2 && \
  python setup.py install && \
  export PATH=$PATH:/opt/idr-2.0.2/bin/ && \
  rm -f /opt/2.0.2.zip

# INSTALL macs2

RUN . /opt/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  pip install macs2==2.2.5

# ADD pysam

RUN . /opt/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  conda config --add channels r && \
  conda config --add channels bioconda && \
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
COPY argparse.bash /opt
COPY ccbr_cutadapt_pe.bash /opt
COPY ccbr_remove_blacklisted_reads_pe.bash /opt
COPY ccbr_atac_trim_align_pe.bash /opt
COPY ccbr_bowtie2_align_pe.bash /opt
COPY atac_assign_multimappers.py /opt
COPY ccbr_bam_filter_by_mapq.py /opt
COPY the_blacklists.fa /opt
COPY the_blacklists.fa.amb /opt
COPY the_blacklists.fa.ann /opt
COPY the_blacklists.fa.bwt /opt
COPY the_blacklists.fa.pac /opt
COPY the_blacklists.fa.sa /opt

RUN chmod a+x /opt/bedSort
RUN chmod a+x /opt/*.bash
RUN chmod a+x /opt/*.py
RUN chmod a+r /opt/*.*

ENV PATH="/opt:${PATH}"

WORKDIR /data

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
