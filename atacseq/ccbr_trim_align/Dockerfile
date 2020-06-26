FROM ubuntu:16.04

# Tools installed in this docker image:
# bedGraphToBigWig
# bedops       2.4.37     conda:python3
# bedSort
# bedtools     2.25.0
# bowtie2      2.2.6
# bwa          0.7.12
# cutadapt     1.18       conda:python3
# fastqc       0.11.8
# genrich      0.6
# homer 			 4.11       conda:python3
# idr          2.0.4.2    conda:python3
# jre          1.8.0_222
# meme				 5.1.0          
# macs2        2.2.5      conda:python3
# miniconda3   4.7.12.1
# multiqc      1.8        conda:python3
# picardcloud
# samtools     0.1.19

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

RUN apt-get install -y cpanminus unzip

RUN cpanm FindBin Term::ReadLine

RUN apt-get install -y \
    bc=1.06.95-9build1 \
    bedtools=2.25.0-1 \
    bowtie2=2.2.6-2 \
    bwa=0.7.12-5 \
    cmake=3.5.1-1ubuntu3 \
    default-jre=2:1.8-56ubuntu2 \
    g++=4:5.3.1-1ubuntu1 \
    gcc=4:5.3.1-1ubuntu1 \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libssl-dev \
    make=4.1-6 \
    pandoc=1.16.0.2~dfsg-1 \
    parallel=20161222-1~ubuntu0.16.04.1 \
    pigz=2.3.1-2 \
    pkg-config \
    samtools=0.1.19-1ubuntu1 \
    wget \
    zlib1g-dev

# INSTALL R packages


RUN apt-get install -y --allow-unauthenticated r-base-core libblas-dev liblapack-dev

# latest Rcpp version 1.0.4 causes errors while installing GoSemSim --> enrichplot --> ChIPseeker .. hence install previous version
COPY Rcpp_1.0.3.tar.gz /opt2
WORKDIR /opt2
RUN R CMD INSTALL Rcpp_1.0.3.tar.gz

COPY install_R_packages.R /opt2
WORKDIR /opt2
RUN Rscript install_R_packages.R

# INSTALL ChIPseeker
COPY install_R_packages_chipseeker.R /opt2
WORKDIR /opt2
RUN Rscript install_R_packages_chipseeker.R


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
# metaseq required for TSS score calculation is also installed in python2 conda environment
RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python2 && \
  pip install --upgrade pip && \
  pip install argparse && \
  pip install numpy && \
  pip install pysam && \
  pip install scipy && \
  pip install pybedtools
#  pip install metaseq


# INSTALL cutadapt

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  pip install --upgrade pip && \
  pip install cutadapt==1.18 && \
  pip install numpy && \
  pip install matplotlib===3.1.1

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
  pip install multiqc==1.8

# INSTALL homer2
# conda meme (ame) does not work ... hence installed separately

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  conda install homer blat

# INSTALL meme

WORKDIR /opt2
# ADD http://meme-suite.org/meme-software/5.1.0/meme-5.1.0.tar.gz meme.tar.gz
COPY meme-5.1.0.tar.gz meme.tar.gz
RUN tar xzvf meme.tar.gz && \
  cd meme-5.1.0 && \
  ./configure --prefix=/opt2/meme --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt && \
  make && \
  make install
ENV PATH=$PATH:/opt2/meme/bin:/opt2/meme/libexec/meme-5.1.0

# INSTALL bedops

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  conda install bedops

# INSTALL deeptools

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  conda install -c bioconda deeptools

# INSTALL pandas

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  pip install --upgrade pip && \
  pip install pandas
  
# INSTALL preseq
# download from https://github.com/smithlabcode/preseq/releases/download/v2.0.1/preseq_v2.0.1.tar.bz2
RUN apt-get install -y --allow-unauthenticated libgsl2 libgsl-dev
COPY preseq_v2.0.1.tar.bz2 /opt2
WORKDIR /opt2
RUN tar xjvf preseq_v2.0.1.tar.bz2 && \
	cd preseq && \
	make
ENV PATH="/opt2/preseq:${PATH}"


# SCRIPTS AND DB files

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
COPY transpose.sh /opt2
COPY ccbr_atac_config.yaml /opt2

# python scripts

COPY atac_assign_multimappers.py /opt2
COPY ccbr_atac_make_qcstats_table.py /opt2
COPY ccbr_bam2FLD.py /opt2
COPY ccbr_bam_filter_by_mapq.py /opt2
COPY ccbr_counts2density.py /opt2
COPY ccbr_get_consensus_peaks.py /opt2
COPY ccbr_narrowPeak_distributions.py /opt2
COPY get_tss_with_bins_per_gene.py /opt2
COPY nrf.py /opt2
COPY qc_annotated2peakwidthdensity.py /opt2
COPY qc_get_frip_stats_table.py /opt2

# bash scripts

COPY ccbr_atac_motif_enrichment.bash /opt2  
COPY ccbr_atac_trim_align_pe.bash /opt2
COPY ccbr_atac_trim_align_pe_with_reference_tarball.bash /opt2
COPY ccbr_atac_qc.bash /opt2
COPY ccbr_bam2FLD.bash /opt2
COPY ccbr_bam2nrf.bash /opt2
COPY ccbr_bowtie2_align_pe.bash /opt2
COPY ccbr_cutadapt_pe.bash /opt2
COPY ccbr_frip.bash /opt2
COPY ccbr_genrich_peak_calling.bash /opt2
COPY ccbr_jaccard_pca.bash /opt2
COPY ccbr_macs2_peak_calling.bash /opt2
COPY ccbr_narrowPeak_distribution.bash /opt2
COPY ccbr_remove_blacklisted_reads_pe.bash /opt2
COPY ccbr_tagAlign2TSS.bash /opt2

# R scripts

COPY pairwise.Rmd /opt2
COPY ccbr_annotate_peaks.R /opt2
COPY ccbr_annotate_bed.R /opt2
COPY db /opt2/db 

RUN chmod a+x /opt2/bedSort
RUN chmod a+x /opt2/bedGraphToBigWig
RUN chmod a+x /opt2/*.bash
RUN chmod a+x /opt2/*.py
RUN chmod a+r /opt2/*.*
RUN chmod a+r Dockerfile
RUN chmod -R a+rX /opt2/db

ENV PATH="/opt2:${PATH}"

WORKDIR /data2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
