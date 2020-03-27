# `covid`
Automatically download the latest SARS-CoV-2 sequences from GenBank and RefSeq using NCBI's Virus portal, and run multiple sequence alignment using MAFFT.

## Building the Docker Container for MAFFT
``` bash
# See listing of images on computer
docker image ls

# Build Dockerfile
docker build --tag=ccbr_mafft:v0.0.1 .

# Peak around the container: verify things run correctly
docker run -ti ccbr_mafft:v0.0.1 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag ccbr_mafft:v0.0.1 skchronicles/ccbr_mafft:v0.0.1
docker tag ccbr_mafft:v0.0.1 skchronicles/ccbr_mafft         # latest
docker tag ccbr_mafft:v0.0.1 nciccbr/ccbr_mafft:v0.0.1
docker tag ccbr_mafft:v0.0.1 nciccbr/ccbr_mafft              # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_mafft:v0.0.1
docker push skchronicles/ccbr_mafft:latest
docker push nciccbr/ccbr_mafft:v0.0.1
docker push nciccbr/ccbr_mafft:latest 
```

## Run `covid`
``` bash
docker run -v $PWD:/data ccbr_mafft:v0.0.1 /opt/covid help
docker run -v $PWD:/data ccbr_mafft:v0.0.1 /opt/covid list
docker run -v $PWD:/data ccbr_mafft:v0.0.1 /opt/covid view
docker run -v $PWD:/data ccbr_mafft:v0.0.1 /opt/covid download 
docker run -v $PWD:/data ccbr_mafft:v0.0.1 /opt/covid msa input_seq [out_msa.fa]
```

## USAGE
```bash
USAGE:
    covid <view|download|msa|list|help>  input_sequences.fa  [output_msa_results.fa]

Help Information:
    view       Displays an example input fasta file
    download   Downloads the latest SARS-CoV-2 sequences from RefSeq and GenBank
    msa        Run Multiple Sequence Alignment using MAFFT
    list       Lists all available sub-commands
    help       Displays usage and this help page

Examples:
    covid list
    covid view
    covid download
    covid msa mysequences.fa
    covid help
```

