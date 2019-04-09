FROM 	nciccbr/ccbr_uropa:v0.0.1
MAINTAINER vishal.koparde@nih.gov

COPY argparse.bash /opt
COPY ccbr_uropa_annotate_tf.bash /opt
COPY tf.json /opt

RUN chmod a+x /opt/*.bash
ENV PATH=/opt:$PATH
WORKDIR /data