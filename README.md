# :anchor: Dockers :anchor:

This repository contains recipes to Dockers created by [CCBR](https://bioinformatics.ccr.cancer.gov/ccbr/) team members to be used :

- in our on-prem pipelines orchestrated using Snakemake and Nextflow. These pipelines run on [Biowulf](https://hpc.nih.gov/).
- in pipelines running on DNAnexus (and SBG).
- in pipelines running on AWS using AWS [Genomics CLI](https://aws.amazon.com/genomics-cli/)
  
The docker images were pushed to dockerhub and are available [here](https://hub.docker.com/u/nciccbr).

## General conventions:

 * Each docker image has the following folders:
   * `/data2` &rarr; Default working directory
   * `/opt2` &rarr; Tools/software installed inside the docker container go here

> NOTE: The suffix `2` ensures that there is no conflict with the hosts' `/data` and `/opt` folders.

 * Original recipe `Dockerfile` is copied into the docker image itself in the `/opt2` folder.
 
 * Most docker images are built using our own base image `nciccbr/ccbr_ubuntu_20.04:latest`

> NOTE: Some of the older docker images may use one of the following base image:
>   * ```ubuntu:16.04```
>   * ```ubuntu:18.04```
>   * ```bitnami/minideb:jessie``` : Docker images built using this base image tend to have a smaller digital footprint.

 * Docker images should have the following environmental variables:
   * `BUILD_DATE`
   * `BUILD_TAG`
   * `REPONAME`

## How the base image `nciccbr/ccbr_ubuntu_20.04` is built:

- uses `ubuntu:20.04` LTS as its base image
- includes `/data2` and `/opt2` folder
- workdir is always `/data2`
- contains most commonly used bioinformatics tools like
  - bwa
  - bowtie2
  - samtools
  - bedtools
  - etc.

## How to build your own docker image:

- use `nciccbr/ccbr_ubuntu_20.04:SOMETAG` as base image
- add "layers" on top of the base image in order to add tools of interest

Use these lines to kick start your recipe:
```bash
FROM nciccbr/ccbr_ubuntu_20.04:SOMETAG

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}
```

- Once your recipe is ready, `build` and `push` scripts are from the `scripts` folder can be used to:
  - automatically build and push built docker images to the `nciccbr` dockerhub account
  - both scripts take only 1 argument, ie., TAG
  - `build` expects the folder to contain a file named `Dockerfile.TAG` which is symlinked to `Dockerfile` prior to initiating a build
  - `build` extracts the name of the repo from the current folder name
  - `push` pushes the built docker over to dockerhub

eg.: Running the following:

```bash
> /path/to/scripts/build v2
```

inside the folder `/some/path/ccbr_xyz` will 
- look for a file `Dockerfile.v2` in the folder
- symlink it to `Dockerfile`
- build a docker with the tag `nciccbr/ccbr_xyz:v2`

Following this up with:

```bash
> /path/to/scripts/path v2
```

will 
- push `nciccbr/ccbr_xyz:v2` dockerhub
- replace `nciccbr/ccbr_xyz:latest` with `nciccbr/ccbr_xyz:v2` on dockerhub

