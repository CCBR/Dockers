FROM nciccbr/ccbr_ubuntu_base_18.04:v1.0
RUN apt-get update 
RUN apt-get install -y bedtools

WORKDIR /opt2
#ADD v5.0.tar.gz /opt2
#RUN pip install pybedtools==0.8.1
RUN git clone https://github.com/kopardev/TALON.git
WORKDIR /opt2/TALON
RUN pip install .

ADD Dockerfile /opt2
RUN chmod -R a+rX /opt2
WORKDIR /data2

MAINTAINER vishal.koparde@nih.gov
