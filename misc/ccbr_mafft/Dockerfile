# Dockerfile for MAFFT multiple sequence alignment tool (version 7.453)
FROM ubuntu:18.04

RUN mkdir -p /data/
RUN mkdir -p /opt


# Package information for MAFFT: https://mafft.cbrc.jp/alignment/software/merge.html
# Install MAFFT
RUN apt-get update && \
    apt-get install -y build-essential curl && \
    curl -L http://mafft.cbrc.jp/alignment/software/mafft-7.453-with-extensions-src.tgz > /tmp/mafft-7.453-with-extensions-src.tgz && \
    cd /tmp && tar zxvf mafft-7.453-with-extensions-src.tgz && \
    cd /tmp/mafft-7.453-with-extensions/core && \
    sed -e "s/^PREFIX = \/usr\/local/PREFIX = \/opt/" Makefile > Makefile.tmp && \
    mv Makefile.tmp Makefile && \
    make clean && make && make install && \
    cd /tmp && rm -rf mafft-7.453-with-extensions && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# Add /opt/bin to PATH
ENV PATH /opt/bin:$PATH

WORKDIR /data

# Download latest covid-19 data for MSA: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049
RUN mkdir -p /opt/example/
RUN curl --http1.1 --retry 5 --verbose -L 'https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/?q=*:*&fq=%7B!tag=SeqType_s%7DSeqType_s:(%22Nucleotide%22)&fq=VirusLineageId_ss:(2697049)&cmd=download&sort=SourceDB_s%20desc,CreateDate_dt%20desc&dlfmt=fasta&fl=id,Definition_s,Nucleotide_seq' > /opt/example/sequences.fa 

COPY Dockerfile /opt
COPY README.md /opt
COPY covid /opt/
