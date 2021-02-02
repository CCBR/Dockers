### Building Docker Image for TRUST4 (v0.1.0) 

Tcr Receptor Utilities for Solid Tissue (TRUST) is a computational tool to analyze TCR and BCR sequences using unselected RNA sequencing data, profiled from solid tissues, including tumors. TRUST4 performs de novo assembly on V, J, C genes including the hypervariable complementarity-determining region 3 (CDR3) and reports consensus of BCR/TCR sequences. TRUST4 then realigns the contigs to IMGT reference gene sequences to report the corresponding information. TRUST4 supports both single-end and paired-end sequencing data with any read length.

### Install

1. Clone the [GitHub repo](https://github.com/liulab-dfci/TRUST4), e.g. with `git clone https://github.com/liulab-dfci/TRUST4.git`
2. Run `make` in the repo directory

You will find the executable files in the downloaded directory. If you want to run TRUST4 without specifying the directory, you can either add the directory of TRUST4 to the environment variable PATH or create a soft link ("ln -s") of the file "run-trust4" to a directory in PATH.

TRUST4 depends on [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads) and samtools depends on [zlib](http://en.wikipedia.org/wiki/Zlib). For MacOS, TRUST4 has been successfully compiled with gcc_darwin17.7.0 and gcc_9.2.0 installed by Homebrew.


Directly below are instructions for building an image using the provided Dockerfile:
```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=ccbr_trust4:v0.0.1 .

# Testing
docker run -ti ccbr_trust4:v0.0.1 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag ccbr_trust4:v0.0.1 skchronicles/ccbr_trust4:v0.0.1
docker tag ccbr_trust4:v0.0.1 skchronicles/ccbr_trust4         # latest
docker tag ccbr_trust4:v0.0.1 nciccbr/ccbr_trust4:v0.0.1
docker tag ccbr_trust4:v0.0.1 nciccbr/ccbr_trust4              # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_trust4:v0.0.1
docker push skchronicles/ccbr_trust4:latest
docker push nciccbr/ccbr_trust4:v0.0.1
docker push nciccbr/ccbr_trust4:latest 
```

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a non-org account.
