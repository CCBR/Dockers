FROM ubuntu:20.04
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN mkdir -p /opt2 && mkdir -p /data2
RUN apt-get update && apt-get install -y dirmngr gnupg apt-transport-https ca-certificates software-properties-common

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
&& add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

RUN apt-get -y install r-base r-base-core r-recommended r-base-dev
RUN apt-get -y install libcurl4-openssl-dev libssl-dev libboost-dev libxml2-dev
WORKDIR /opt2
COPY install_packages.R /opt2
RUN Rscript install_packages.R
COPY install_packages2.R /opt2
RUN Rscript install_packages2.R
COPY install_packages3.R /opt2
RUN Rscript install_packages3.R
RUN apt-get -y install libmagick++-dev
COPY install_packages4.R /opt2
RUN Rscript install_packages4.R
COPY install_packages5.R /opt2
RUN Rscript install_packages5.R
RUN apt-get -y install pandoc

COPY Dockerfile /opt2
WORKDIR /data2
RUN chmod -R a+rX /opt2 && apt-get clean

MAINTAINER vishal.koparde@nih.gov
