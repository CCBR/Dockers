### Building Docker Image for STAR (2.7) v0.0.1

Directly below are instructions for building an image using the following Dockerfile:
- Dockerfile

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --tag=ccbr_star_2.7.0f:v0.0.1 .

# Testing
docker run -ti ccbr_star_2.7.0f:v0.0.1 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag ccbr_star_2.7.0f:v0.0.1 skchronicles/ccbr_star_2.7.0f:v0.0.1
docker tag ccbr_star_2.7.0f:v0.0.1 skchronicles/ccbr_star_2.7.0f         # latest
docker tag ccbr_star_2.7.0f:v0.0.1 nciccbr/ccbr_star_2.7.0f:v0.0.1
docker tag ccbr_star_2.7.0f:v0.0.1 nciccbr/ccbr_star_2.7.0f              # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_star_2.7.0f:v0.0.1
docker push skchronicles/ccbr_star_2.7.0f:latest
docker push nciccbr/ccbr_star_2.7.0f:v0.0.1
docker push nciccbr/ccbr_star_2.7.0f:latest 
```


### Building Docker Image for STAR (2.7) v0.0.2

Directly below are instructions for building an image using the following Dockerfile:
- Dockerfile.v0.0.2

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile.v0.0.2 --tag=ccbr_star_2.7.0f:v0.0.2 .

# Testing
docker run -ti ccbr_star_2.7.0f:v0.0.2 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag ccbr_star_2.7.0f:v0.0.2 skchronicles/ccbr_star_2.7.0f:v0.0.2
docker tag ccbr_star_2.7.0f:v0.0.2 skchronicles/ccbr_star_2.7.0f         # latest
docker tag ccbr_star_2.7.0f:v0.0.2 nciccbr/ccbr_star_2.7.0f:v0.0.2
docker tag ccbr_star_2.7.0f:v0.0.2 nciccbr/ccbr_star_2.7.0f              # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_star_2.7.0f:v0.0.2
docker push skchronicles/ccbr_star_2.7.0f:latest
docker push nciccbr/ccbr_star_2.7.0f:v0.0.2
docker push nciccbr/ccbr_star_2.7.0f:latest 
```
