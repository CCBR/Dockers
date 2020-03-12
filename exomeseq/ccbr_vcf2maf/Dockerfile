FROM 		nciccbr/ccbr_vep_99:v0.0.1
MAINTAINER	vishal.koparde@nih.gov

USER		root
RUN mkdir -p /data2

RUN		apt update && apt install -y zlib1g-dev libncurses5-dev libncursesw5-dev wget liblzma-dev libbz2-dev 

# INSTALL htslib samtools bcftools

RUN	cd /usr/bin && \
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
tar -vxjf htslib-1.9.tar.bz2 && \
cd htslib-1.9 && \
make 

RUN	cd /usr/bin && \
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
tar -vxjf samtools-1.9.tar.bz2 && \
cd samtools-1.9 && \
make

RUN	cd /usr/bin && \
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
tar -vxjf bcftools-1.9.tar.bz2 && \
cd bcftools-1.9 && \
make

ENV	PATH=$PATH:/usr/bin/bcftools-1.9:/usr/bin/samtools-1.9:/usr/bin/htslib-1.9

# INSTALL vcf2maf

COPY		v1.6.17.tar.gz /opt2

WORKDIR		/opt2
RUN		tar xzvf v1.6.17.tar.gz
ENV		PATH=$PATH:/opt2/vcf2maf-1.6.17
RUN		chmod a+x /opt2/vcf2maf-1.6.17/*.pl

# INSTALL MINICONDA

# ADD https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh /opt2/miniconda.sh
COPY Miniconda3-4.7.12.1-Linux-x86_64.sh /opt2/miniconda.sh
RUN /bin/bash /opt2/miniconda.sh -b -p /opt2/conda && \
    rm /opt2/miniconda.sh && \
    ln -s /opt2/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt2/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# CREATE MINICONDA PYTHON2 ENVIRONMENT

RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda create -n python2 python=2.7

WORKDIR		/data2

COPY 		Dockerfile /opt2
COPY        argparse.bash /opt2
COPY        vcf2maf_wrapper.bash /opt2
RUN 		chmod a+r /opt2/Dockerfile && chmod a+r /opt2/argparse.bash && chmod a+r /opt2/vcf2maf_wrapper.bash

RUN         apt-get clean