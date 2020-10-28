# Picard/2.11.0 Dockerfile (v0.0.1)
FROM ubuntu:18.04

MAINTAINER kuhnsa@nih.gov

RUN mkdir -p /data2
RUN mkdir -p /opt2/picard/

# Update apt-get before downloading packages 
RUN apt-get update && \
    apt-get upgrade -y 

# Install Standard Java Runtime
RUN apt-get install -y \
	default-jre \
	wget && \
apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /opt2/picard

# Get specified release of Picard
ENV PICARD_VERSION="2.11.0"
RUN wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar
RUN chmod a+x /opt2/picard/*.jar
ENV PICARDJARPATH=/opt2/picard

COPY Dockerfile /opt2/
WORKDIR /data2
