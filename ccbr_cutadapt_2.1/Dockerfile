FROM python:3.5-alpine3.8

RUN 	mkdir -p /data \
	&& mkdir -p /opt \
	&& apk add --no-cache gcc musl-dev \
	&& pip install cutadapt==2.1

COPY TruSeq_and_nextera_adapters.consolidated.fa /opt
COPY Dockerfile /opt
WORKDIR /data

MAINTAINER vishal.koparde@nih.gov
