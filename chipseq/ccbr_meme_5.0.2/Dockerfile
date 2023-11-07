FROM mambaorg/micromamba:1.5-bookworm

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

RUN mkdir -p /opt2 && mkdir -p /data2

# meme error workaround
# https://stackoverflow.com/a/66319448
# https://www.biostars.org/p/401084/
RUN mamba install -c bioconda \
  meme=5.0.2 icu=58.2 parallel

# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
