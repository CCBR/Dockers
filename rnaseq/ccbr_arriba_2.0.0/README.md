### Building Docker Image for Arriba (v2.0.0.0) 

Arriba requires the latest version of STAR, i.e. >= 2.7.6a, to be installed on the target system. 

In addition, it is recommended install the following programs:
 - `python` script to find optimal read lengths 
 - `samtools` (sam to bam conversion)
 - `r-base` arriba relies on R for generating figures
 - `wget` pull in releases from Github
 - `samblaster` optional tool for marking dups in bam file
  - dependencies `git` and `g++`

Directly below are instructions for building an image using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=ccbr_arriba_2.0.0:v0.0.1 .

# Testing
docker run -ti ccbr_arriba_2.0.0:v0.0.1 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag ccbr_arriba_2.0.0:v0.0.1 skchronicles/ccbr_arriba_2.0.0:v0.0.1
docker tag ccbr_arriba_2.0.0:v0.0.1 skchronicles/ccbr_arriba_2.0.0         # latest
docker tag ccbr_arriba_2.0.0:v0.0.1 nciccbr/ccbr_arriba_2.0.0:v0.0.1
docker tag ccbr_arriba_2.0.0:v0.0.1 nciccbr/ccbr_arriba_2.0.0              # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_arriba_2.0.0:v0.0.1
docker push skchronicles/ccbr_arriba_2.0.0:latest
docker push nciccbr/ccbr_arriba_2.0.0:v0.0.1
docker push nciccbr/ccbr_arriba_2.0.0:latest 
```

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a non-org account.
