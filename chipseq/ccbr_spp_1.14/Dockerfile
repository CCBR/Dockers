FROM kopardev/ccbr_samtools_1.9

RUN mkdir -p /data
RUN mkdir -p /opt

RUN apt-get update && \
    apt-get upgrade -y 

RUN apt-get install -y r-base-core libboost-dev

COPY spp_1.14.tar.gz /opt

COPY install_R_packages.R /opt

RUN cd /opt && Rscript install_R_packages.R

COPY run_spp.R /opt
RUN chmod a+rx /opt/run_spp.R

ENV PATH="/opt":$PATH

WORKDIR /data

MAINTAINER vishal.koparde@nih.gov
