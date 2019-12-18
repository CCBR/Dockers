FROM ubuntu:16.04

# Tools installed in this docker image:
# bedGraphToBigWig
# bedSort
# bedtools     2.25.0
# bowtie2      2.2.6
# bwa          0.7.12
# cutadapt     1.18       conda:python3
# fastqc       0.11.8
# genrich      0.6
# idr          2.0.4.2    conda:python3
# jre          1.8.0_222
# macs2        2.2.5      conda:python3
# miniconda3   4.7.12.1
# multiqc      1.8        conda:python3
# picardcloud
# samtools     0.1.19

RUN mkdir -p /data2
RUN mkdir -p /opt2

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
    liblzma-dev \
    pandoc \
    libcurl4-openssl-dev \
    libssl-dev

# INSTALL R

RUN apt-get install -y r-base-core
COPY install_R_packages.R /opt2
WORKDIR /opt2
RUN Rscript install_R_packages.R

# INSTALL miniconda

WORKDIR /opt2
COPY Miniconda3-4.7.12.1-Linux-x86_64.sh /opt2/miniconda.sh
# downloaded from https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh
RUN /bin/bash /opt2/miniconda.sh -b -p /opt2/conda && \
    rm /opt2/miniconda.sh && \
    ln -s /opt2/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt2/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# CREATE MINICONDA PYTHON3 ENVIRONMENT

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda create -n python3 python=3.6 && \
  conda create -n python2 python=2.7

# INSTALL packages in python2
# bam2fld scripts requires older version of pysam ... newer version gives errors ...hence installing pysam in python2 and python3
RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python2 && \
  pip install --upgrade pip && \
  pip install argparse && \
  pip install numpy && \
  pip install pysam && \
  pip install scipy


# INSTALL cutadapt

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  pip install --upgrade pip && \
  pip install cutadapt==1.18 && \
  pip install numpy && \
  pip install matplotlib

# INSTALL idr

WORKDIR /opt2
# conda config --add channels conda-forge && \
RUN . /opt2/conda/etc/profile.d/conda.sh && \
    conda activate python3 && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda install idr && \
    conda update idr

# INSTALL macs2

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  pip install macs2==2.2.5

# ADD pysam

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  conda config --add channels r && \
  conda install pysam

# INSTALL fastqc

WORKDIR /opt2
COPY fastqc_v0.11.8.zip /opt2
# downloaded from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
RUN unzip fastqc_v0.11.8.zip && \
  cd /opt2/FastQC && \
  chmod 755 fastqc && \
  ln -s /opt2/FastQC/fastqc /usr/local/bin/fastqc

# INSTALL Genrich v0.6

COPY Genrich-master.zip /opt2
WORKDIR /opt2
RUN unzip Genrich-master.zip && \
    cd Genrich-master && \
    make
ENV PATH="/opt2/Genrich-master:${PATH}"

# INSTALL multiqc

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  pip install multiqc

# INSTALL meme and homer2

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  conda install meme homer parallel blat

COPY TruSeq_and_nextera_adapters.consolidated.fa /opt2
COPY Dockerfile /opt2
COPY picardcloud.jar /opt2
COPY bedSort /opt2
COPY bedGraphToBigWig /opt2
COPY argparse.bash /opt2
COPY the_blacklists.fa /opt2
COPY the_blacklists.fa.amb /opt2
COPY the_blacklists.fa.ann /opt2
COPY the_blacklists.fa.bwt /opt2
COPY the_blacklists.fa.pac /opt2
COPY the_blacklists.fa.sa /opt2

COPY atac_assign_multimappers.py /opt2
COPY ccbr_bam2FLD.py /opt2
COPY ccbr_bam_filter_by_mapq.py /opt2
COPY ccbr_narrowPeak_distributions.py /opt2
  
COPY ccbr_atac_trim_align_pe.bash /opt2
COPY ccbr_bam2FLD.bash /opt2
COPY ccbr_bowtie2_align_pe.bash /opt2
COPY ccbr_cutadapt_pe.bash /opt2
COPY ccbr_genrich_peak_calling.bash /opt2
COPY ccbr_genrich_peak_calling_two_replicates.bash /opt2
COPY ccbr_jaccard_pca.bash /opt2
COPY ccbr_macs2_peak_calling.bash /opt2
COPY ccbr_macs2_peak_calling_two_replicates.bash /opt2
COPY ccbr_narrowPeak_distribution.bash /opt2
COPY ccbr_remove_blacklisted_reads_pe.bash /opt2

COPY pairwise.Rmd /opt2

COPY db /opt2/db 

RUN chmod a+x /opt2/bedSort
RUN chmod a+x /opt2/bedGraphToBigWig
RUN chmod a+x /opt2/*.bash
RUN chmod a+x /opt2/*.py
RUN chmod a+r /opt2/*.*
RUN chmod -R a+rX /opt2/db

ENV PATH="/opt2:${PATH}"

WORKDIR /data2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
