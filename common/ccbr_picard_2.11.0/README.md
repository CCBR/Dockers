### Building Docker Image for Picard (2.11) 

Directly below are instructions for building an image using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=ccbr_picard:v0.0.1 .

# Testing
docker run -ti ccbr_picard:v0.0.1 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag ccbr_picard:v0.0.1 skchronicles/ccbr_picard:v0.0.1
docker tag ccbr_picard:v0.0.1 skchronicles/ccbr_picard         # latest
docker tag ccbr_picard:v0.0.1 nciccbr/ccbr_picard:v0.0.1
docker tag ccbr_picard:v0.0.1 nciccbr/ccbr_picard              # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_picard:v0.0.1
docker push skchronicles/ccbr_picard:latest
docker push nciccbr/ccbr_picard:v0.0.1
docker push nciccbr/ccbr_picard:latest 
```
