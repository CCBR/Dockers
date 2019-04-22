FROM 		nciccbr/ccbr_bedtools:v0.0.1
MAINTAINER	vishal.koparde@nih.gov

COPY 		Dockerfile /opt
COPY		ccbr_chipseq_overlap_filtering_narrow.bash /opt
COPY		ccbr_chipseq_overlap_filtering_broad.bash /opt
COPY		argparse.bash /opt

RUN 		chmod a+x /opt/*.bash
ENV 		PATH=/opt:$PATH
WORKDIR 	/data