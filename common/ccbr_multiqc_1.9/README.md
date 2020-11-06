### Building Docker Image for MultiQC/1.9

Dependencies:
 - python3 
   - multiqc (current latest release == 1.9)

The latest version of multiqc no longer supports python2 (or compatability is no longer guaranteed!).

Directly below are instructions for building an image using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=ccbr_multiqc_1.9:v0.0.1 .

# Testing
docker run -ti ccbr_multiqc_1.9:v0.0.1 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag ccbr_multiqc_1.9:v0.0.1 skchronicles/ccbr_multiqc_1.9:v0.0.1
docker tag ccbr_multiqc_1.9:v0.0.1 skchronicles/ccbr_multiqc_1.9         # latest
docker tag ccbr_multiqc_1.9:v0.0.1 nciccbr/ccbr_multiqc_1.9:v0.0.1
docker tag ccbr_multiqc_1.9:v0.0.1 nciccbr/ccbr_multiqc_1.9              # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_multiqc_1.9:v0.0.1
docker push skchronicles/ccbr_multiqc_1.9:latest
docker push nciccbr/ccbr_multiqc_1.9:v0.0.1
docker push nciccbr/ccbr_multiqc_1.9:latest 
```
