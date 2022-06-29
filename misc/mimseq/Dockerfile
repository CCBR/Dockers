FROM nciccbr/ccbr_ubuntu_base_20.04:v1.4

# install bedtools
WORKDIR /opt2
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && \
  tar -zxvf bedtools-2.30.0.tar.gz && \
  rm -f bedtools-2.30.0.tar.gz && \
  cd bedtools2 && \
  make
ENV PATH=/opt2/bedtools2/bin:$PATH

# install samtools
WORKDIR /opt2
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
   tar xjvf samtools-1.11.tar.bz2  && \
   rm -f samtools-1.11.tar.bz2  && \
   cd samtools-1.11/ && \
   make
ENV PATH=/opt2/samtools-1.11:$PATH

RUN apt-get update && apt-get install -y libffi-dev libdeflate-dev libsqlite3-dev

# According to https://github.com/nedialkova-lab/mim-tRNAseq/issues/33#issuecomment-1168411432
# mimseq may work best with python/3.7
# python3.8 gave error described in https://github.com/nedialkova-lab/mim-tRNAseq/issues/33
# default version is 3.8 hence install 3.7 and setting that as default
# distutils dev packages are required to run pip and install packages via pip for python3.7
# Build python 3.8 from source
WORKDIR /opt2
RUN wget https://www.python.org/ftp/python/3.7.1/Python-3.7.1.tgz && \
	tar -xvf Python-3.7.1.tgz && \
	rm Python-3.7.1.tgz
WORKDIR /opt2/Python-3.7.1/ 
RUN ./configure --enable-loadable-sqlite-extensions --enable-optimizations && \
  make && \
  make altinstall

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.8 1 && \
  update-alternatives --install /usr/bin/python python /usr/local/bin/python3.7 2 && \
  update-alternatives --set python /usr/local/bin/python3.7 && \
  python -m pip install --upgrade pip

RUN python -m pip install numpy==1.21.1
RUN python -m pip install pandas==1.3.1
RUN python -m pip install Biopython==1.79
RUN python -m pip install pyfiglet==0.8.post1
RUN python -m pip install pysam==0.16.0.1
RUN python -m pip install seaborn==0.11.1
RUN python -m pip install pybedtools==0.8.2
RUN python -m pip install requests==2.26.0

WORKDIR /opt2
RUN git clone https://github.com/nedialkova-lab/mim-tRNAseq.git && \
  cd mim-tRNAseq && \
  python setup.py build && \
  python setup.py install

# install usearch
WORKDIR /opt2
RUN wget https://drive5.com/downloads/usearch10.0.240_i86linux32.gz && \
        gunzip usearch10.0.240_i86linux32.gz && \
        chmod a+x usearch10.0.240_i86linux32 && \
        mkdir -p /opt2/usearch && \
        mv usearch10.0.240_i86linux32 /opt2/usearch/usearch && \
        rm -f usearch10.0.240_i86linux32.gz
ENV PATH=/opt2/usearch:$PATH

# install blast
WORKDIR /opt2
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/ncbi-blast-2.10.1+-x64-linux.tar.gz && \
  tar xzvf ncbi-blast-2.10.1+-x64-linux.tar.gz && rm -f ncbi-blast-2.10.1+-x64-linux.tar.gz
ENV PATH=/opt2/ncbi-blast-2.10.1+/bin:$PATH

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

# install cmalign
RUN apt-get install -y infernal

# cleanup
RUN apt-get clean && apt-get purge && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
  rm -rf /opt2/gmap-2019-02-26

COPY Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2

WORKDIR /data2
