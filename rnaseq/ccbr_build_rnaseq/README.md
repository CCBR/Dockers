### Building Docker Image for building RNA-seq References

This is a general purpose Docker image is used to build reference files for RNA-seek pipeline. It includes dependencies to build RSeQC tin.py reference file along with any other python package dependencies to build other resources. RNA-seek build contains a few python scripts to help build any require reference files.

Here is a current list of bundled python packages:
 - argparse
 - Bio
 - HTSeq
 - json
 - numpy
 - scipy
 - pysam

Here is more information about the about TIN:
- https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0922-z

Directly below are instructions for building an image using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=ccbr_build_rnaseq:v0.0.1 .

# Testing
docker run -ti ccbr_build_rnaseq:v0.0.1 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag ccbr_build_rnaseq:v0.0.1 skchronicles/ccbr_build_rnaseq:v0.0.1
docker tag ccbr_build_rnaseq:v0.0.1 skchronicles/ccbr_build_rnaseq         # latest
docker tag ccbr_build_rnaseq:v0.0.1 nciccbr/ccbr_build_rnaseq:v0.0.1
docker tag ccbr_build_rnaseq:v0.0.1 nciccbr/ccbr_build_rnaseq              # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_build_rnaseq:v0.0.1
docker push skchronicles/ccbr_build_rnaseq:latest
docker push nciccbr/ccbr_build_rnaseq:v0.0.1
docker push nciccbr/ccbr_build_rnaseq:latest 
```

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a non-org account.

