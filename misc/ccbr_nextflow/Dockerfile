FROM ubuntu:18.04

RUN mkdir -p /data2
RUN mkdir -p /opt2

RUN apt-get update && \
    apt-get upgrade -y 

#RUN apt-get install -y  software-properties-common && \
#    add-apt-repository ppa:webupd8team/java -y && \
#    apt-get update && \
#    echo oracle-java7-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections && \
#    apt-get install -y oracle-java8-installer

RUN apt-get install -y default-jre wget
WORKDIR /opt2 
RUN mkdir -p nextflow 
WORKDIR /opt2/nextflow 
RUN wget -qO- https://get.nextflow.io | bash
ENV PATH="/opt2/nextflow":$PATH

COPY Dockerfile /opt2

RUN chmod a+r /opt2/Dockerfile

WORKDIR /data2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
