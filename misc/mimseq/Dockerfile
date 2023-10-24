FROM nciccbr/ccbr_ubuntu_base_20.04:v6

RUN apt-get update && apt-get install -y libffi-dev libdeflate-dev libsqlite3-dev libcurl4-openssl-dev

# install usearch
WORKDIR /opt2
RUN wget https://drive5.com/downloads/usearch10.0.240_i86linux32.gz && \
        gunzip usearch10.0.240_i86linux32.gz && \
        chmod a+x usearch10.0.240_i86linux32 && \
        mkdir -p /opt2/usearch && \
        mv usearch10.0.240_i86linux32 /opt2/usearch/usearch && \
        rm -f usearch10.0.240_i86linux32.gz
ENV PATH=/opt2/usearch:$PATH

# install mimseq
RUN mamba install -c conda-forge mimseq=1.3.7

# cleanup
RUN apt-get clean && apt-get purge && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 

COPY Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2/Dockerfile

WORKDIR /data2
