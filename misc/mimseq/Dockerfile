FROM nciccbr/ccbr_ubuntu_base_20.04:v1.4

# According to https://github.com/nedialkova-lab/mim-tRNAseq/issues/33#issuecomment-1168411432
# mimseq may work best with python/3.7
# python3.8 gave error described in https://github.com/nedialkova-lab/mim-tRNAseq/issues/33
# default version is 3.8 hence install 3.7 and setting that as default
# distutils dev packages are required to run pip and install packages via pip for python3.7
RUN add-apt-repository -y ppa:deadsnakes/ppa && \
  apt-get update && \
  apt-get install -y python3.7 python3.7-distutils python3.7-dev && \
  update-alternatives --install /usr/bin/python python /usr/bin/python3.8 1 && \
  update-alternatives --install /usr/bin/python python /usr/bin/python3.7 2 && \
  update-alternatives --set python /usr/bin/python3.7

# install usearch
WORKDIR /opt2
RUN wget https://drive5.com/downloads/usearch10.0.240_i86linux32.gz && \
	gunzip usearch10.0.240_i86linux32.gz && \
	chmod a+x usearch10.0.240_i86linux32 && \
	mv usearch10.0.240_i86linux32 usearch && \
	rm -f usearch10.0.240_i86linux32.gz
ENV PATH=$PATH:/opt2

# install mimseq
# RUN pip3 install mimseq
# the above command will install mimseq with the default python version .. python3.8
# use the following command to install mimseq with the new python3.7 version
RUN python -m pip install mimseq

# install cmalign
RUN apt-get install -y infernal

# install gmap gsnap
WORKDIR /opt2
RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2019-02-26.tar.gz && \
 tar xvzf gmap-gsnap-2019-02-26.tar.gz && \
 rm -f gmap-gsnap-2019-02-26.tar.gz && \
 mkdir -p /opt2/gmap
WORKDIR /opt2/gmap-2019-02-26
RUN ./configure --prefix=/opt2/gmap && \
  make && \
  make install
ENV PATH=$PATH:/opt2/gmap/bin

# install R>=4.0
RUN apt-get update && \
  apt install -y dirmngr gnupg apt-transport-https ca-certificates software-properties-common --fix-missing && \
  apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
  add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
  apt install -y r-base
# install R packages

RUN apt-get install -y libcurl4-openssl-dev libssl-dev
#RUN apt-get install -y r-cran-devtools
RUN R -e "install.packages('devtools')"
RUN R -e "install.packages('RColorBrewer')"
RUN R -e "install.packages('pheatmap')"
RUN R -e "install.packages('calibrate')"
RUN R -e "install.packages('gridExtra')"
RUN R -e "install.packages('plyr')"
RUN R -e "install.packages('dplyr')"
RUN R -e "install.packages('reshape2')"
RUN R -e "install.packages('circlize')"
RUN R -e "install.packages('tidyverse')"
RUN R -e "install.packages('ggplot2')"
RUN R -e "install.packages('ggpol')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('DESeq2',update=TRUE,ask=FALSE)"
RUN R -e "BiocManager::install('ComplexHeatmap', update=TRUE, ask=FALSE)"

# install bedtools
WORKDIR /opt2
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && \
  tar -zxvf bedtools-2.30.0.tar.gz && \
  cd bedtools2 && \
  make
ENV PATH=/opt2/bedtools2/bin:$PATH

# install samtools
WORKDIR /opt2
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
   tar xjvf samtools-1.11.tar.bz2  && \
   cd samtools-1.11/ && \
   make
ENV PATH=/opt2/samtools-1.11:$PATH

# install blastn
WORKDIR /opt2
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/ncbi-blast-2.10.1+-x64-linux.tar.gz && \
  tar xzvf ncbi-blast-2.10.1+-x64-linux.tar.gz
ENV PATH=/opt2/ncbi-blast-2.10.1+/bin:$PATH

# cleanup
WORKDIR /data2
RUN apt-get clean && \
  rm -rf /opt2/gmap-2019-02-26

