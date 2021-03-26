FROM continuumio/miniconda3

RUN mkdir -p /data2
RUN mkdir -p /opt2

RUN apt-get update && \
    apt-get upgrade -y 

COPY fastqc.yaml /opt2
WORKDIR /opt2
RUN conda-env create -f /opt2/fastqc.yaml
RUN echo "source activate fastqc_env" > ~/.bashrc
ENV PATH /opt/conda/envs/fastqc_env/bin:$PATH

COPY Dockerfile /opt2

RUN chmod a+r /opt2/Dockerfile

WORKDIR /data2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
