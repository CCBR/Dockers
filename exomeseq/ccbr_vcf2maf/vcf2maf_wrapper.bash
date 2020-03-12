#!/bin/bash
# module load singularity
. /opt2/conda/etc/profile.d/conda.sh && \
  conda activate python2

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="vcf2maf wrapper script... dont forget to bind /data/CCBR_Pipeliner/db/PipeDB/lib/vcf2maf_resources at /vepresourcebundlepath"
source /opt2/argparse.bash || exit 1
# source /data/CCBR_Pipeliner/db/PipeDB/lib/vcf2maf_resources/scripts/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--vcf',required=True, help='input vcf file')
parser.add_argument('--maf',required=True, help='output maf file')
parser.add_argument('--genome',required=True, help='hg19/hg38/mm10')
parser.add_argument('--tid',required=True, help='TumorID ... sample id of tumor sample')
parser.add_argument('--nid',required=False, default='', help='NormalID ... sample id of normal sample ... do not provide in case of tumor-only')
EOF

# This assumes that the VEP resource bundle has been bound to the docker/singularity container at "/vepresourcebundlepath"
# On biowulf it should be '/data/CCBR_Pipeliner/db/PipeDB/lib/vcf2maf_resources'

ncbi_build="unknown"
species="unknown"
if [ $GENOME == "hg38" ]; then
        ncbi_build="GRCh38"
        species="homo_sapiens"
fi

if [ $GENOME == "hg19" ]; then
        ncbi_build="GRCh37"
        species="homo_sapiens"
fi

if [ $GENOME == "mm10" ]; then
        ncbi_build="GRCm38"
        species="mus_musculus"
fi

fa="/vepresourcebundlepath/fastas/${species}/${ncbi_build}.fa"
dotvep="/vepresourcebundlepath/VEP_tarballs/.vep"
filtervcf="/vepresourcebundlepath/filtervcf/${species}/${ncbi_build}.filter.vcf.gz"

if [ "$NID" == "" ]; then
vcf2maf.pl \
--vep-forks 1 \
--input-vcf $VCF \
--output-maf $MAF \
--tumor-id $TID \
--vep-path /opt/vep/src/ensembl-vep \
--vep-data $dotvep \
--filter-vcf $filtervcf \
--ncbi-build $ncbi_build \
--species $species \
--ref-fasta $fa
else
vcf2maf.pl \
--vep-forks 1 \
--input-vcf $VCF \
--output-maf $MAF \
--tumor-id $TID \
--normal-id $NID \
--vep-path /opt/vep/src/ensembl-vep \
--vep-data $dotvep \
--filter-vcf $filtervcf \
--ncbi-build $ncbi_build \
--species $species \
--ref-fasta $fa
fi

