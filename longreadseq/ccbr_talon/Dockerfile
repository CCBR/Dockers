FROM nciccbr/ccbr_ubuntu_base_18.04:v1.0
RUN apt-get update 

WORKDIR /opt2
#ADD v5.0.tar.gz /opt2
#RUN pip install pybedtools==0.8.1
ADD https://github.com/mortazavilab/TALON/archive/refs/tags/v5.0.tar.gz /opt2
RUN tar xzvf v5.0.tar.gz && rm -f v5.0.tar.gz
WORKDIR /opt2/TALON-5.0 
RUN pip install .

ADD Dockerfile /opt2
RUN chmod -R a+rX /opt2
WORKDIR /data2

MAINTAINER vishal.koparde@nih.gov
