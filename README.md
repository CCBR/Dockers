# Dockers
Dockers created by [CCBR](https://ccbr.ccr.cancer.gov) to be used on DNAnexus (and SBG) in the following pipelines:

 * RNASeq
 * ChIPSeq
 * ExomeSeq

The docker images were pushed to dockerhub and are available [here](https://hub.docker.com/u/nciccbr).

## General convention followed in all Docker images:

 * All docker images have the following folders:
   * ```/data2``` --> Default working directory
   * ```/opt2``` --> Tools/software installed inside the docker container goes here
 * The suffix ```2``` ensures that there is no conflict with the host ```/data``` and ```/opt``` folders, which is likely if the docker images is run as a converted singularity container on HPC systems.
 * Original Dockerfile is copied into the docker image as ```/opt2/Dockerfile```. One can view this directly by running

  ```docker run nciccbr/ccbr_xxx_yyy cat /opt2/Dockerfile```
 
 * All docker images are built using one of the following base images:
   * ```ubuntu:16.04```
   * ```bitnami/minideb:jessie``` : Docker images built using this base image tend to have a smaller digital footprint.
