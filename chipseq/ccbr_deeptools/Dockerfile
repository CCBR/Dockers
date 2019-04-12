FROM nciccbr/ccbr_deeptools:v0.0.1
# base image has samtools bedtools deeptools==3.2.0
MAINTAINER vishal.koparde@nih.gov

COPY argparse.bash /opt
COPY ccbr_bam2bw.bash /opt
COPY Dockerfile /opt

RUN chmod a+x /opt/*.bash
ENV PATH=/opt:$PATH
WORKDIR /data
