FROM ubuntu:16.04

RUN mkdir -p /data2
RUN mkdir -p /opt2

RUN apt-get update && \
    apt-get upgrade -y

RUN apt-get update && apt-get install -y \
    g++ \
    libboost-dev \
    make git

WORKDIR /opt2
RUN git clone https://github.com/phbradley/pubtcrs.git
WORKDIR /opt2/pubtcrs
RUN make && mv /opt2/pubtcrs/bin/* /opt2 && chmod a+x /opt2/*

ENV PATH="/opt2/:${PATH}"

WORKDIR /data2
MAINTAINER vishal.koparde@nih.gov
