### Build from Dockerfile
```bash
# See listing of images on computer
docker image ls

# Build
docker build --tag=ccbr_ismart:v0.0.1 .

# Updating tag(s) before pushing to DockerHub
docker tag ccbr_ismart:v0.0.1 skchronicles/ccbr_ismart:v0.0.1
docker tag ccbr_ismart:v0.0.1 skchronicles/ccbr_ismart        # latest
docker tag ccbr_ismart:v0.0.1 nciccbr/ccbr_ismart:v0.0.1
docker tag ccbr_ismart:v0.0.1 nciccbr/ccbr_ismart             # latest

# Check out new tag(s)
docker image ls

# Peak around the container: verify things run correctly
docker run -ti ccbr_ismart:v0.0.1 /bin/bash

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_ismart:v0.0.1
docker push skchronicles/ccbr_ismart:latest
docker push nciccbr/ccbr_ismart:v0.0.1
docker push nciccbr/ccbr_ismart:latest
```

### Run using Singularity
```bash
# Pull Image from Dockerhub
module load singularity
singularity pull -F ccbr_ismart_v0.0.1.sif docker://nciccbr/ccbr_ismart:v0.0.1

# Display Usage
singularity exec -B "$PWD:/data2/" ccbr_ismart_v0.0.1.sif iSMARTv3.py --help
```
