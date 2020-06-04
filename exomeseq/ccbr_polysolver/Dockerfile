# Dockerfile for polysolver (v4): 
FROM sachet/polysolver:v4

RUN mkdir -p /data2/
RUN mkdir -p /opt2/

#Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

# Set environment variables
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV DEBIAN_FRONTEND noninteractive


# Package information for polysolver: https://software.broadinstitute.org/cancer/cga/polysolver 
# Install r-base 
RUN apt-get -qq update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends r-base && \
    apt-get clean && \
    apt-get purge && \
  rm -rf /var/lib/apt/lists/*


WORKDIR /data2/

# Copy Dockerfile and README into /opt2
COPY Dockerfile /opt2
COPY README.md /opt2

# PATCH 1
# Singularity containers have read-only filesystems
# Cannot write to /home/polysolver directory which author specifies as a TMP_DIR
# Do not trust authors source code to allow for write!
RUN sed -i 's@TMP_DIR=/home/polysolver@TMP_DIR=$outDir@g' /home/polysolver/scripts/config.bash /home/polysolver/scripts/shell_call_hla_type /home/polysolver/scripts/shell_call_hla_type_test 

# PATCH 2
# Found a bug in the `shell_call_hla_type` script
# iFile is an uninitialized variable if `insertCalc` argument is set to 1
# Causes several problems downstream in `first_allele_calculations_fork.pl` and other child processes
RUN sed -i 's/# calculate insert size distribution/# calculate insert size distribution\niFile=1/' /home/polysolver/scripts/shell_call_hla_type /home/polysolver/scripts/shell_call_hla_type_test

