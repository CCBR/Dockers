FROM nciccbr/ccbr_ubuntu_base_20.04:v1.4

RUN apt-get update && apt-get install -y libffi-dev libdeflate-dev libsqlite3-dev libcurl4-openssl-dev

# install miniconda
ENV CONDA_DIR /opt2/conda
WORKDIR /opt2
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt2/miniconda.sh && \
  /bin/bash /opt2/miniconda.sh -b -p /opt2/conda
ENV PATH=$CONDA_DIR/bin:$PATH

# install mamba
RUN conda create -n mimseq -y python=3.7
RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate mimseq && \
  conda config --add channels conda-forge && \
  conda install -c conda-forge -y mamba
RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate mimseq && \
  mamba install -c bioconda -y mimseq
WORKDIR /opt2

# create wrapper for mimseq conda env
COPY mimseq /opt2
RUN chmod a+x /opt2/mimseq
ENV PATH=/opt2:$PATH

# install usearch
WORKDIR /opt2
RUN wget https://drive5.com/downloads/usearch10.0.240_i86linux32.gz && \
        gunzip usearch10.0.240_i86linux32.gz && \
        chmod a+x usearch10.0.240_i86linux32 && \
        mkdir -p /opt2/usearch && \
        mv usearch10.0.240_i86linux32 /opt2/usearch/usearch && \
        rm -f usearch10.0.240_i86linux32.gz
ENV PATH=/opt2/usearch:$PATH

# cleanup
RUN apt-get clean && apt-get purge && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 

COPY Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2

WORKDIR /data2
