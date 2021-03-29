FROM	quay.io/biocontainers/gimmemotifs:0.13.1--py36h516909a_3
MAINTAINER	vishal.koparde@nih.gov

COPY	seqlogo /usr/local/lib/python3.6/site-packages/gimmemotifs/included_tools/weblogo/
RUN 	chmod a+rx /usr/local/lib/python3.6/site-packages/gimmemotifs/included_tools/weblogo/seqlogo

RUN	mkdir -p /data
WORKDIR	/data

COPY	Dockerfile /opt




