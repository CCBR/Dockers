### Building Docker Image for Kraken (2.1.1) and Krona (2.7.1)

Directly below are instructions for building an image using the provided Dockerfile:

A pre-built krona's taxonomy database can be found within the following location in the container's filesystem:
 - `/opt2/Krona-2.7.1/KronaTools/taxonomy` 

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=ccbr_kraken_v2.1.1:v0.0.1 .

# Testing
docker run -ti ccbr_kraken_v2.1.1:v0.0.1 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag ccbr_kraken_v2.1.1:v0.0.1 skchronicles/ccbr_kraken_v2.1.1:v0.0.1
docker tag ccbr_kraken_v2.1.1:v0.0.1 skchronicles/ccbr_kraken_v2.1.1         # latest
docker tag ccbr_kraken_v2.1.1:v0.0.1 nciccbr/ccbr_kraken_v2.1.1:v0.0.1
docker tag ccbr_kraken_v2.1.1:v0.0.1 nciccbr/ccbr_kraken_v2.1.1              # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_kraken_v2.1.1:v0.0.1
docker push skchronicles/ccbr_kraken_v2.1.1:latest
docker push nciccbr/ccbr_kraken_v2.1.1:v0.0.1
docker push nciccbr/ccbr_kraken_v2.1.1:latest 
```

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a non-org account.
