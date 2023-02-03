FROM nciccbr/ccbr_baser:latest
ENV TAG="v0.1"
# install FRASER
RUN Rscript -e 'BiocManager::install("FRASER")'
# cleanup etc
COPY Dockerfile /opt2/Dockerfile.${TAG}
RUN chmod a+rX /opt2/Dockerfile.${TAG}
WORKDIR /data2
RUN apt-get clean
# details
LABEL maintainer="kopardevn@nih.gov"
LABEL github_handle="kopardev"
ARG BUILD_DATE
LABEL build_date=$BUILD_DATE
