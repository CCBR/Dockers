### Building Docker Image for building RSeQC reference files

This Docker image is used to build reference files for [RSeQC tin.py](http://rseqc.sourceforge.net/#tin-py).

Here is more information about the about TIN:
- https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0922-z

Directly below are instructions for building an image using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=ccbr_build_rseqc_ref:v0.0.1 .

# Testing
docker run -ti ccbr_build_rseqc_ref:v0.0.1 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag ccbr_build_rseqc_ref:v0.0.1 skchronicles/ccbr_build_rseqc_ref:v0.0.1
docker tag ccbr_build_rseqc_ref:v0.0.1 skchronicles/ccbr_build_rseqc_ref         # latest
docker tag ccbr_build_rseqc_ref:v0.0.1 nciccbr/ccbr_build_rseqc_ref:v0.0.1
docker tag ccbr_build_rseqc_ref:v0.0.1 nciccbr/ccbr_build_rseqc_ref              # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_build_rseqc_ref:v0.0.1
docker push skchronicles/ccbr_build_rseqc_ref:latest
docker push nciccbr/ccbr_build_rseqc_ref:v0.0.1
docker push nciccbr/ccbr_build_rseqc_ref:latest 
```

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a non-org account.

