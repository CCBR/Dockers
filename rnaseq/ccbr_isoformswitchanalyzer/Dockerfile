FROM 		nciccbr/ccbr_rbase:3.6

RUN 		apt-get update && \
   			apt-get upgrade -y 

RUN			apt-get install -y libcurl4-openssl-dev \
			libxml2-dev

RUN			apt-get install -y libssl-dev

COPY 		install_R_packages.R /opt

RUN 		cd /opt && Rscript install_R_packages.R
