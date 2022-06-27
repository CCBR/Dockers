FROM nciccbr/ccbr_ubuntu_base_20.04:v1.4

# install usearch
WORKDIR /opt2
RUN wget https://drive5.com/downloads/usearch10.0.240_i86linux32.gz && \
	gunzip usearch10.0.240_i86linux32.gz && \
	chmod a+x usearch10.0.240_i86linux32 && \
	mv usearch10.0.240_i86linux32 usearch && \
	rm -f usearch10.0.240_i86linux32.gz
ENV PATH=$PATH:/opt2

# install mimseq
RUN pip3 install mimseq

# install cmalign
RUN apt-get install -y infernal

# install blastn
RUN apt-get install -y ncbi-blast+

# install gmap gsnap
WORKDIR /opt2
RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2021-12-17.tar.gz && \
 tar xvzf gmap-gsnap-2021-12-17.tar.gz && \
 rm -f gmap-gsnap-2021-12-17.tar.gz && \
 mkdir -p /opt2/gmap
WORKDIR /opt2/gmap-2021-12-17
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

# cleanup etc
COPY Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2
WORKDIR /data2
RUN apt-get clean && \
  rm -rf /opt2/gmap-2021-12-17

