# Arriba/v2.0.0 Dockerfile (v0.0.1)

# Arriba/v2.0.0 requires `STAR >= 2.7.6a` to be installed on the target system.
# In addition, the build includes install the following programs:
#   `r-base` arriba relies on R for generating figures  <required>
#   `wget` pull in releases from Github                 <build_dependency>
#   `python` script to find optimal read lengths        <optional>
#   `samtools` sam to bam conversion                    <optional>
#   `picard` mark dups and collect stats                <optional>
#   `samblaster` marking dups in bam file               <optional>
#        - dependencies `git` and `g++`

FROM ubuntu:18.04

MAINTAINER Skyler Kuhn <kuhnsa@nih.gov>

RUN mkdir -p /data2
RUN mkdir -p /opt2/bin

WORKDIR /opt2

# Update apt-get before downloading packages
RUN apt-get update && \
    apt-get upgrade -y

# Download packages
# Install Java 8, R 3.4.4-1, Samtools/1.7-1, build and runtime dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
      build-essential \
      make \
      gcc \
      g++ \
      samtools \
      r-base \
      ca-certificates \
      libcurl4-openssl-dev \
      libxml2-dev \
      wget \
      samtools \
      openjdk-8-jre \
      zlibc \
      zlib1g \
      zlib1g-dev \
      locales

# Install Required Arriba R packages
RUN Rscript -e 'install.packages("circlize", repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("GenomicRanges", "GenomicAlignments"))'

# Install and build Arriba/2.0.0
RUN wget 'https://github.com/suhrig/arriba/releases/download/v2.0.0/arriba_v2.0.0.tar.gz' && \
    tar -xzf arriba_v2.0.0.tar.gz && \
    rm arriba_v2.0.0.tar.gz && \
    cd /opt2/arriba_v2.0.0 && \
    make

# Install and build STAR/2.7.6a
RUN wget 'https://github.com/alexdobin/STAR/archive/2.7.6a.tar.gz' && \
    tar -xzf 2.7.6a.tar.gz && \
    rm 2.7.6a.tar.gz && \
    cd /opt2/STAR-2.7.6a/source && \
    make STAR

# Install and build samblaster/v.0.1.26
RUN wget 'https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.26/samblaster-v.0.1.26.tar.gz' && \
    tar -xzf samblaster-v.0.1.26.tar.gz && \
    rm samblaster-v.0.1.26.tar.gz && \
    cd /opt2/samblaster-v.0.1.26 && \
    make

# Adding Picard Cloud JAR
RUN mkdir -p /opt2/picard/
COPY resources/picardcloud.jar /opt2/picard/picard.jar
RUN chmod a+rx /opt2/picard/*.jar

# Arriba Wrapper Scripts
WORKDIR /opt2/bin

# download_references.sh wrapper
RUN echo '#!/bin/bash\n\
cd /references\n\
/opt2/arriba_v2.0.0/download_references.sh $1 && \\\n\
ASSEMBLY=$(sed -e "s/viral+.*//" -e "s/+.*//" <<<"$1") && \\\n\
cp /opt2/arriba_v2.0.0/database/*$ASSEMBLY* /references' > /usr/local/bin/download_references.sh && \
chmod a+x /usr/local/bin/download_references.sh

# run_arriba.sh wrapper
RUN echo '#!/bin/bash\n\
cd /opt2/arriba_v2.0.0/test/ && \\\n\
/opt2/arriba_v2.0.0/run_arriba.sh /references/STAR_index_* /references/*.gtf /references/*.fa /references/blacklist_*.tsv.gz /references/known_fusions_*.tsv.gz /references/protein_domains_*.gff3 ${THREADS-8} read1.fastq.gz $(ls read2.fastq.gz 2> /dev/null)' > /usr/local/bin/arriba.sh && \
chmod a+x /usr/local/bin/arriba.sh

# draw_fusions.R wrapper
RUN echo '#!/bin/bash\n\
Rscript /opt2/arriba_v2.0.0/draw_fusions.R --annotation=$(ls /references/*.gtf) --fusions=/opt2/arriba_v2.0.0/test/fusions.tsv --output=fusions.pdf --proteinDomains=$(ls /references/protein_domains_*.gff3) --alignments=Aligned.sortedByCoord.out.bam --cytobands=$(ls /references/cytobands_*.tsv)' > /usr/local/bin/draw_fusions.sh && \
chmod a+x /usr/local/bin/draw_fusions.sh

# Configure "locale", see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8

# Add Tools to PATH
ENV PICARDJARPATH=/opt2/picard
ENV PATH="/opt2/STAR-2.7.6a/source":$PATH
ENV PATH="/opt2/arriba_v2.0.0":$PATH
ENV PATH="/opt2/samblaster-v.0.1.26":$PATH
ENV PATH="/opt2/bin":$PATH

# Clean-up Image
RUN apt-get clean && apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Copy the Dockerfile used to create image in /opt2
COPY Dockerfile /opt2
WORKDIR /data2
