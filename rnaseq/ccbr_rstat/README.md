### Building Docker Image for rule stat

Depedencies:
 - picard
 - samtools
 - python3 

Directly below are instructions for building an image using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=ccbr_rstat:v0.0.1 .

# Testing
docker run -ti ccbr_rstat:v0.0.1 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag ccbr_rstat:v0.0.1 skchronicles/ccbr_rstat:v0.0.1
docker tag ccbr_rstat:v0.0.1 skchronicles/ccbr_rstat         # latest
docker tag ccbr_rstat:v0.0.1 nciccbr/ccbr_rstat:v0.0.1
docker tag ccbr_rstat:v0.0.1 nciccbr/ccbr_rstat              # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_rstat:v0.0.1
docker push skchronicles/ccbr_rstat:latest
docker push nciccbr/ccbr_rstat:v0.0.1
docker push nciccbr/ccbr_rstat:latest 
```
