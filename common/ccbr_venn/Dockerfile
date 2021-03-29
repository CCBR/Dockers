FROM nciccbr/ccbr_baser

RUN apt-get install -y libudunits2-dev libgdal-dev

WORKDIR /opt2
# COPY install_packages.R /opt2
# RUN Rscript install_packages.R

COPY Dockerfile /opt2
COPY 2set_venn.R /opt2
RUN chmod a+rX /opt2/2set_venn.R
ENV PATH=/opt2:$PATH
WORKDIR /data2
RUN chmod -R a+rX /opt2 && apt-get clean

MAINTAINER vishal.koparde@nih.gov
