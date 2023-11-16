FROM nciccbr/ccbr_ubuntu_base_20.04:v6.1

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

# uninstall R from base image
RUN DEBIAN_FRONTEND=noninteractive apt-get --purge remove \
    r-base r-base-dev r-recommended r-cran-* && \
    apt-get autoremove && \
    apt-get update

# install R packages with conda
COPY environment.txt /data2/
RUN mamba install -c bioconda -c conda-forge --file /data2/environment.txt
ENV R_LIBS_USER=/opt2/conda/lib/R/library/

# cleanup etc
# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
ENV PATH="/opt2/:$PATH"
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
WORKDIR /data2