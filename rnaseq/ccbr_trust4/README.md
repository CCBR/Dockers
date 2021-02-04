## TRUST4
Tcr Receptor Utilities for Solid Tissue (TRUST) is a computational tool to analyze TCR and BCR sequences using unselected RNA sequencing data, profiled from solid tissues, including tumors. TRUST4 performs de novo assembly on V, J, C genes including the hypervariable complementarity-determining region 3 (CDR3) and reports consensus of BCR/TCR sequences. TRUST4 then realigns the contigs to IMGT reference gene sequences to report the corresponding information. TRUST4 supports both single-end and paired-end sequencing data with any read length.

For more information or to contact the authors, please visit [TRUST4's GitHub repository](https://github.com/liulab-dfci/TRUST4).

### Dependencies
TRUST4 depends on [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads) and samtools depends on [zlib](http://en.wikipedia.org/wiki/Zlib). For MacOS, TRUST4 has been successfully compiled with gcc_darwin17.7.0 and gcc_9.2.0 installed by Homebrew.

### Building Docker Image for TRUST4 (v1.0.2-beta) 

Directly below are instructions for building an image using the provided Dockerfile:
```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile using defined release
RELEASE_VERSION=1.0.2-beta 
docker build --build-arg TRUST4_VERSION=${RELEASE_VERSION} --no-cache -f Dockerfile --tag=ccbr_trust4:${RELEASE_VERSION} .

# Testing
docker run -ti ccbr_trust4:${RELEASE_VERSION} /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag ccbr_trust4:${RELEASE_VERSION} skchronicles/ccbr_trust4:v${RELEASE_VERSION}
docker tag ccbr_trust4:${RELEASE_VERSION} skchronicles/ccbr_trust4         # latest
docker tag ccbr_trust4:${RELEASE_VERSION} nciccbr/ccbr_trust4:v${RELEASE_VERSION}
docker tag ccbr_trust4:${RELEASE_VERSION} nciccbr/ccbr_trust4              # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_trust4:v${RELEASE_VERSION}
docker push skchronicles/ccbr_trust4:latest
docker push nciccbr/ccbr_trust4:v${RELEASE_VERSION}
docker push nciccbr/ccbr_trust4:latest 
```

> **Please Note**: Any references to `skchronicles` or `nciccbr` should be replaced with your username or org account.

### Run with Singularity

Run TRUST4 using the provided example dataset:
```
# Assumes Singularity in $PATH
# Pull Image from DockerHub
SINGULARITY_CACHEDIR=$PWD singularity pull -F docker://nciccbr/ccbr_trust4:latest
# Run with example BAM file as Input
singularity exec -B $PWD:/data2 ccbr_trust4_latest.sif run-trust4 -b /opt2/TRUST4/example/example.bam -f /opt2/TRUST4/hg38_bcrtcr.fa --ref /opt2/TRUST4/human_IMGT+C.fa -o /data2/test1
# Run with example FastQ file as Input
singularity exec -B $PWD:/data2 ccbr_trust4_latest.sif run-trust4 -f /opt2/TRUST4/human_IMGT+C.fa --ref /opt2/TRUST4/human_IMGT+C.fa -1 /opt2/TRUST4/example/example_1.fq -2 /opt2/TRUST4/example/example_2.fq  -o /data2/test2
```
