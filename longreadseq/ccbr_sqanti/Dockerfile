FROM continuumio/miniconda3
RUN apt-get update && mkdir /opt2 && mkdir /data2

#add minimap2
WORKDIR /opt2
RUN apt-get install -y curl 
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 -O && \
    tar xjvf minimap2-2.17_x64-linux.tar.bz2 && rm -rf minimap2-2.17_x64-linux.tar.bz2
ENV PATH=$PATH:/opt2/minimap2-2.17_x64-linux

#install SQANTI3
WORKDIR /opt2
RUN apt-get install -y git && \
    git clone https://github.com/ConesaLab/SQANTI3.git
WORKDIR /opt2/SQANTI3/ 
RUN conda env update -n base --file SQANTI3.conda_env.yml
RUN echo "source activate base" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH

WORKDIR /opt2
RUN echo "options(repos = c(CRAN = \"https://cran.rstudio.com\"))" >.Rprofile
COPY install_packages.R /opt2
RUN Rscript install_packages.R

#install cDNA Cupcake
RUN apt-get install -y zlibc zlib1g zlib1g-dev bzip2 make gcc g++ wget
#ADD https://github.com/Magdoll/cDNA_Cupcake/archive/refs/tags/v22.0.0.tar.gz /opt2
COPY v22.0.0.tar.gz /opt2/
RUN tar xzvf v22.0.0.tar.gz && rm -f v22.0.0.tar.gz
WORKDIR /opt2/cDNA_Cupcake-22.0.0
RUN pip install numpy Cython && pip install .
ENV PYTHONPATH="${PYTHONPATH}:/opt2/cDNA_Cupcake-22.0.0/sequence"
ENV PATH=$PATH:/opt2/SQANTI3
WORKDIR /opt2/SQANTI3/utilities
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred && \
    chmod a+x /opt2/SQANTI3/utilities/gtfToGenePred
RUN chmod a+x /opt2/SQANTI3/sqanti3_qc.py && \
    mv /opt2/SQANTI3/sqanti3_RulesFilter.py /opt2/SQANTI3/sqanti3_RulesFilter.py.tmp && \
    echo "#!/usr/bin/env python" > /opt2/SQANTI3/sqanti3_RulesFilter.py && \
    cat /opt2/SQANTI3/sqanti3_RulesFilter.py.tmp >> /opt2/SQANTI3/sqanti3_RulesFilter.py && \
    rm -f /opt2/SQANTI3/sqanti3_RulesFilter.py.tmp && \
    chmod a+x /opt2/SQANTI3/sqanti3_RulesFilter.py


ADD Dockerfile /opt2
RUN chmod -R a+rX /opt2
WORKDIR /data2

MAINTAINER vishal.koparde@nih.gov
