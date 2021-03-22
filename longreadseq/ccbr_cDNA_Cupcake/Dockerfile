FROM nciccbr/ccbr_ubuntu_base_18.04:v1.0
RUN apt-get update 

WORKDIR /opt2
ADD https://github.com/Magdoll/cDNA_Cupcake/archive/refs/tags/v22.0.0.tar.gz /opt2
RUN tar xzvf v22.0.0.tar.gz && rm -f v22.0.0.tar.gz
WORKDIR /opt2/cDNA_Cupcake-22.0.0
RUN pip install numpy Cython
#RUN python setup.py build && python setup.py install
RUN pip install .

ADD Dockerfile /opt2
RUN chmod -R a+rX /opt2
WORKDIR /data2

MAINTAINER vishal.koparde@nih.gov
