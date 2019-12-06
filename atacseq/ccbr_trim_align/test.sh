#!/bin/bash

. /opt2/conda/etc/profile.d/conda.sh
conda activate python3
export LC_ALL=C.UTF-8
export LANG=C.UTF-8
#python --version
#cutadapt --version
#idr -h
#macs2 --version
#echo $PATH
python /opt2/ccbr_bam_filter_by_mapq.py --help
multiqc --version
conda deactivate
#fastqc --version

