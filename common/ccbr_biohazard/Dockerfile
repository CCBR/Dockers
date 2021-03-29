FROM ubuntu:18.04

RUN mkdir -p /opt2 && mkdir -p /data2

RUN apt-get update && apt-get install -y git libjudy-dev haskell-platform
RUN cabal update
WORKDIR /opt2
RUN apt-get install -y wget
RUN wget https://get.haskellstack.org/stable/linux-x86_64.tar.gz && tar xzvf linux-x86_64.tar.gz
ENV PATH="/opt2/stack-2.5.1-linux-x86_64:${PATH}"
RUN apt-get install -y xz-utils make
RUN stack ghc -- --version
#RUN git clone https://ustenzel@bitbucket.org/ustenzel/biohazard-tools.git
#RUN cd /opt2 && cabal install biohazard-tools/
