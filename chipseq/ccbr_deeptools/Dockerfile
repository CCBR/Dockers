FROM nciccbr/ccbr_deeptools:v0.0.2
# base image has samtools bedtools deeptools==3.2.0
MAINTAINER vishal.koparde@nih.gov

COPY ccbr_treatmentbam2bw.bash /opt
COPY Dockerfile /opt

RUN chmod a+x /opt/*.bash
ENV PATH=/opt:$PATH
WORKDIR /data
