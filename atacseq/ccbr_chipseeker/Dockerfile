FROM ubuntu:16.04

RUN mkdir -p /data2
RUN mkdir -p /opt2

# Need to add cloud repo to install R v3.6 on ubuntu v16.04
# else default is R v3.2.3 which does not have required libraries like RCurl
RUN apt-get update && \
    apt-get upgrade -y

RUN apt-get install -y cpanminus unzip

RUN cpanm FindBin Term::ReadLine

# 
WORKDIR /opt2
COPY Miniconda3-4.7.12.1-Linux-x86_64.sh /opt2/miniconda.sh
# # downloaded from https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh
RUN /bin/bash /opt2/miniconda.sh -b -p /opt2/conda && \
    rm /opt2/miniconda.sh && \
    ln -s /opt2/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt2/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# # CREATE MINICONDA PYTHON3 ENVIRONMENT
# 
RUN . /opt2/conda/etc/profile.d/conda.sh && \
   conda create -n python3 python=3.6 

# 
RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python3 && \
  pip install --upgrade pip && \
  conda install -c r r-base


#  conda install -c bioconda bioconductor-chipseeker
#RUN . /opt2/conda/etc/profile.d/conda.sh && \
#  conda activate python3 && \
#  conda install -c bioconda bioconductor-org.hs.eg.db bioconductor-org.mm.eg.db bioconductor-txdb.mmusculus.ucsc.mm9.knowngene bioconductor-txdb.mmusculus.ucsc.mm10.knowngene bioconductor-txdb.hsapiens.ucsc.hg19.knowngene bioconductor-txdb.hsapiens.ucsc.hg38.knowngene
