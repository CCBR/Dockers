#!/bin/bash

. /opt/conda/etc/profile.d/conda.sh
conda activate python3
python --version
cutadapt --version
idr -h
macs2 --version
echo $PATH
python /opt/ccbr_bam_filter_by_mapq.py --help
conda deactivate
fastqc --version

