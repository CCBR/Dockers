FROM nciccbr/ccbr_ubuntu_base_20.04:v1.4

# install METRO
WORKDIR /opt2
RUN wget https://github.com/CCBR/METRO/archive/refs/tags/v1.4.tar.gz && \
tar xzvf v1.4.tar.gz
WORKDIR /opt2/METRO-1.4
RUN pip install -r requirements.txt 
ENV PATH="/opt2/METRO-1.4":$PATH

# cleanup etc
COPY Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2
WORKDIR /data2
RUN apt-get clean
