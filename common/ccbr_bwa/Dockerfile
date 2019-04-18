FROM nciccbr/ccbr_bwa:v0.0.1
MAINTAINER vishal.koparde@nih.gov

COPY Dockerfile /opt
COPY ccbr_bwa_chipseq_align_se.bash /opt

RUN chmod a+x /opt/*.bash
ENV PATH=/opt:$PATH
WORKDIR /data
