FROM 		nciccbr/ccbr_vep_97.3:v0.0.2
MAINTAINER	vishal.koparde@nih.gov

USER		root
RUN		apt update && apt install -y libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev wget

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

COPY		v1.6.17.tar.gz /opt

WORKDIR		/opt
RUN		tar xzvf v1.6.17.tar.gz
ENV		PATH=$PATH:/opt/vcf2maf-1.6.17
RUN		chmod a+x /opt/vcf2maf-1.6.17/*.pl

WORKDIR		/data

COPY 		Dockerfile.0.0.3 /opt
