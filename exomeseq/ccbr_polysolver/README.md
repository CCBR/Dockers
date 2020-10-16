# polysolver

`POLYSOLVER` is a tool developed by the Broad for HLA typing based on whole exome sequencing (WES) data; it infers alleles for the three major MHC class I genes: HLA-A, -B, -C.

The latest docker image of polysolver does not have R bundled, and this causes an issue when running the positional argument `-insertCalc` is set to 1. This argument indicates whether empirical insert size distribution should be used in the model (0 or 1). Under the hood, it internally calls `picard CollectInsertSizeMetrics` which needs R to be installed.

This patch includes the latest version of [polysolver](https://software.broadinstitute.org/cancer/cga/polysolver) as the base image with R installed, and it includes a few bugs fixes that were found when testing the program.
It should noted that there are still a few bugs in `polsolver` when setting `-insertCalc` to 1; specically, in the author perl program `/home/polysolver/scripts/first_allele_calculations_fork.pl` which also effects its child processes. I have reached out to the original author about this issue. With that being said, it would aviod setting this argument to 1.

## Building the Docker Image for polysolver
``` bash
# See listing of images on computer
docker image ls

# Build Dockerfile
docker build --tag=ccbr_polysolver:v0.0.1 .

# Peak around the container: verify things run correctly
docker run -ti ccbr_polysolver:v0.0.1 /bin/bash

# Updating tag(s) before pushing to DockerHub
docker tag ccbr_polysolver:v0.0.1 skchronicles/ccbr_polysolver:v0.0.1-beta
docker tag ccbr_polysolver:v0.0.1 skchronicles/ccbr_polysolver         # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_polysolver:v0.0.1-beta
docker push skchronicles/ccbr_polysolver:latest
```

### Run `polysolver` with included test data
``` bash
# Docker
docker run -v $PWD:/data2 ccbr_polysolver:v0.0.1 bash /home/polysolver/scripts/shell_call_hla_type_test /home/polysolver/test/test.bam Unknown 1 hg19 STDFQ 1 test_docker_dir 

# Singularity
singularity pull docker://skchronicles/ccbr_polysolver:v0.0.1-beta
singularity exec -B "$PWD:/data2/" ccbr_polysolver_v0.0.1-beta.sif bash /home/polysolver/scripts/shell_call_hla_type_test /home/polysolver/test/test.bam Unknown 0 hg19 STDFQ 0 test_singularity_dir
```
