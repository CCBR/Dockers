### Building Docker Image for BBTools (38.87) 

BBTools requires `java >= 7` to be installed on the target system. In addition, it is recommended install the following programs: 
 - `samtools` (sam to bam conversion)
 - `pigz` to speed up compression
 - `unzip` to support zip files
 - `bzip2` to support bzip files  

Directly below are instructions for building an image using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=ccbr_bbtools_38.87:v0.0.1 .

# Testing
docker run -ti ccbr_bbtools_38.87:v0.0.1 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag ccbr_bbtools_38.87:v0.0.1 skchronicles/ccbr_bbtools_38.87:v0.0.1
docker tag ccbr_bbtools_38.87:v0.0.1 skchronicles/ccbr_bbtools_38.87         # latest
docker tag ccbr_bbtools_38.87:v0.0.1 nciccbr/ccbr_bbtools_38.87:v0.0.1
docker tag ccbr_bbtools_38.87:v0.0.1 nciccbr/ccbr_bbtools_38.87              # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_bbtools_38.87:v0.0.1
docker push skchronicles/ccbr_bbtools_38.87:latest
docker push nciccbr/ccbr_bbtools_38.87:v0.0.1
docker push nciccbr/ccbr_bbtools_38.87:latest 
```

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a non-org account.
