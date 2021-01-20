### Building Docker Image for rNA (v0.1.0) 

rNA requires R and Python to be installed on the target system.

**Overview of rNA Dependencies**
```
Python >= 3
    Python Packages <optional>
       argparse, numpy==1.18.5, pandas==0.25.3, python-dateutil==2.8.1, 
       pytz==2020.1, six==1.15.0, xlrd==1.2.0
R >= 3.6
   R Packages <required>
       CRAN: 'plyr', 'plotly', 'ggplot2', 'RColorBrewer', 'gridExtra', 
             'crosstalk', 'DT', 'reshape2', 'circlize'
       Bioconductor: 'limma', 'edgeR', 'ComplexHeatmap'  
```


Directly below are instructions for building an image using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=ccbr_rna:v0.0.1 .

# Testing
docker run -ti ccbr_rna:v0.0.1 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag ccbr_rna:v0.0.1 skchronicles/ccbr_rna:v0.0.1
docker tag ccbr_rna:v0.0.1 skchronicles/ccbr_rna         # latest
docker tag ccbr_rna:v0.0.1 nciccbr/ccbr_rna:v0.0.1
docker tag ccbr_rna:v0.0.1 nciccbr/ccbr_rna              # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_rna:v0.0.1
docker push skchronicles/ccbr_rna:latest
docker push nciccbr/ccbr_rna:v0.0.1
docker push nciccbr/ccbr_rna:latest 
```

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a non-org account.
