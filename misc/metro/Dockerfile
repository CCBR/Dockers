FROM nciccbr/ccbr_ubuntu_base_20.04:v1.4
LABEL maintainer "Vishal Koparde *(kopardev on GitHub)*"

# install cufflinks
WORKDIR /opt2
RUN wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz && \
    tar -zxvf cufflinks-2.2.1.Linux_x86_64.tar.gz && \
    rm -rf cufflinks-2.2.1.Linux_x86_64.tar.gz
ENV PATH="/opt2/cufflinks-2.2.1.Linux_x86_64:$PATH"

# install METRO
WORKDIR /opt2
RUN wget https://github.com/CCBR/METRO/archive/refs/tags/v1.4.tar.gz && \
tar xzvf v1.4.tar.gz
WORKDIR /opt2/METRO-1.4
RUN pip install -r requirements.txt 
ENV PATH="/opt2/METRO-1.4":$PATH

# cleanup etc
COPY Dockerfile /opt2/Dockerfile
RUN chmod a+rX /opt2/Dockerfile
WORKDIR /data2
RUN apt-get clean && rm -rf /var/lib/apt/lists/*
