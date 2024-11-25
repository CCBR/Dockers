FROM nciccbr/ccbr_ubuntu_base_20.04:v7

RUN pip3 install --upgrade pip && \
    pip3 install RSeQC==4.0.0 && \
    tin.py -h

COPY Dockerfile /opt2

WORKDIR /data2

RUN apt-get clean
