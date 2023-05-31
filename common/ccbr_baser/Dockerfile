FROM nciccbr/ccbr_ubuntu_base_20.04:v3.1
RUN apt-get update && apt-get install -y \
	software-properties-common \
	dirmngr gnupg apt-transport-https \
	ca-certificates \
	libcurl4-openssl-dev \
	libssl-dev \
	libboost-dev \
	libxml2-dev \
	libharfbuzz-dev \
	libfribidi-dev \
	libmagick++-dev \
	libgit2-dev 
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
&& add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN apt-get update && apt-get -y install r-base r-base-dev r-recommended r-base-dev
WORKDIR /opt2
COPY Dockerfile /opt2
RUN mkdir -p /repo
WORKDIR /repo
COPY . .
RUN Rscript cran_tools.1.R 
RUN apt-get install -y libudunits2-dev
RUN Rscript bc_tools.1.R
COPY cran_tools.2.R .
RUN Rscript cran_tools.2.R
COPY bc_tools.2.R .
RUN Rscript bc_tools.2.R
RUN apt-get clean
MAINTAINER vishal.koparde@nih.gov
#WORKDIR /opt2
#COPY install_packages.R /opt2
#RUN Rscript install_packages.R
