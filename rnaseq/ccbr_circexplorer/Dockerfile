FROM nciccbr/ccbr_ubuntu_base_20.04:v2.1
ENV TAG="v1.0"

# install circxplorer2
WORKDIR /opt2
RUN pip install circexplorer2

# add dependencies
WORKDIR /usr/local/bin
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf && chmod a+x genePredToGtf
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred && chmod a+x gtfToGenePred
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && chmod a+x bedGraphToBigWig
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed && chmod a+x bedToBigBed

COPY Dockerfile /opt2/Dockerfile.${TAG}
RUN if [ -f /opt2/Dockerfile ];then rm -f /opt2/Dockerfile;fi && \
    chmod a+r /opt2/Dockerfile.${TAG} && \
    ln -s Dockerfile.${TAG} /opt2/Dockerfile 
WORKDIR /data2

# clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

# details
LABEL maintainer="vishal.koparde@nih.gov"
LABEL github_handle="kopardev"
ARG BUILD_DATE
LABEL build_date=$BUILD_DATE
