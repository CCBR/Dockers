FROM nciccbr/ccbr_macs:v0.0.2
MAINTAINER vishal.koparde@nih.gov

COPY Dockerfile /opt
COPY ccbr_macs_filterdup_se.bash /opt
COPY ccbr_macs_callpeak_se.bash /opt

RUN chmod a+x /opt/*.bash
ENV PATH=/opt:$PATH

WORKDIR /data
