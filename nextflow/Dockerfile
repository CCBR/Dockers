FROM 		ubuntu:16.04
MAINTAINER	vishal.koparde@nih.gov

# USER		root
RUN 		mkdir -p /data2 && \
			mkdir -p /opt2 && \
			apt-get update && \
    		apt-get upgrade -y

# INSTALL JAVA
RUN			apt-get install -y default-jre curl

# INSTALL NEXTFLOW
RUN			curl -s https://get.nextflow.io | bash && mv /nextflow /opt2/nextflow

WORKDIR		/data2
COPY 		Dockerfile /opt2
RUN 		chmod a+r /opt2/Dockerfile && chmod a+rx /opt2/nextflow
ENV			PATH=$PATH:/opt2
RUN         apt-get clean