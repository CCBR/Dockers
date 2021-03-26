FROM continuumio/miniconda3

RUN mkdir -p /data2
RUN mkdir -p /opt2

RUN apt-get update && \
    apt-get upgrade -y 

COPY flair.yaml /opt2
WORKDIR /opt2
RUN conda-env create -f /opt2/flair.yaml
RUN echo "source activate flair_env" > ~/.bashrc
RUN git clone https://github.com/BrooksLabUCSC/flair.git && \
cd flair && mv flair.py flair.py.old && \
echo "#!//opt/conda/envs/flair_env/bin/python" > flair.py && \
cat flair.py.old >> flair.py && \
chmod a+x flair.py
ENV PATH /opt/conda/envs/flair_env/bin:/opt2/flair:$PATH

COPY Dockerfile /opt2

RUN chmod a+r /opt2/Dockerfile

WORKDIR /data2

RUN apt-get clean

MAINTAINER vishal.koparde@nih.gov
