### Building Docker Image for bam2strandedbw 

Directly below are instructions for building an image using the provided Dockerfile:

If using GNU parallel, please use the `--will-cite` flag to ignore interactive prompt (i.e. `parallel --will-cite < mycmds.sh`).

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=ccbr_bam2strandedbw:v0.0.1 .

# Testing
docker run -ti ccbr_bam2strandedbw:v0.0.1 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag ccbr_bam2strandedbw:v0.0.1 skchronicles/ccbr_bam2strandedbw:v0.0.1
docker tag ccbr_bam2strandedbw:v0.0.1 skchronicles/ccbr_bam2strandedbw         # latest
docker tag ccbr_bam2strandedbw:v0.0.1 nciccbr/ccbr_bam2strandedbw:v0.0.1
docker tag ccbr_bam2strandedbw:v0.0.1 nciccbr/ccbr_bam2strandedbw              # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_bam2strandedbw:v0.0.1
docker push skchronicles/ccbr_bam2strandedbw:latest
docker push nciccbr/ccbr_bam2strandedbw:v0.0.1
docker push nciccbr/ccbr_bam2strandedbw:latest 
```

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a non-org account.

