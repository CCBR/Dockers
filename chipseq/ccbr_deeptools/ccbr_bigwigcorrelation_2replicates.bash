#!/bin/bash
set -e -x -o pipefail
scripts_folder="/opt"
# scripts_folder="/home/kopardev/Desktop/Dockers/chipseq/ccbr_deeptools"
ARGPARSE_DESCRIPTION="Produce bigwig_correlation_heatmap and bigwig_correlation_scatterplot plots for 2 replicate ChiPSeq data"      # this is optional
source /${scripts_folder}/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--treatmentbigwigrep1',required=True, help='treatment bigwig file for Replicate 1')
parser.add_argument('--treatmentbigwigrep2',required=True, help='treatment bigwig file for Replicate 2')
parser.add_argument('--inputbigwigrep1',required=True, help='Input bigwig file for Replicate 1')
parser.add_argument('--inputbigwigrep2',required=True, help='Input bigwig file for Relicate 2')
parser.add_argument('--ncpus',required=False, type=int, default=1, help='nproc')
EOF

multiBigwigSummary bins \
-b $TREATMENTBIGWIGREP1 $TREATMENTBIGWIGREP2 $INPUTBIGWIGREP1 $INPUTBIGWIGREP2 \
-o results.npz \
--smartLabels \
--numberOfProcessors $NCPUS

plotCorrelation \
--corData results.npz \
--plotFile bigwig_correlation_heatmap.pdf \
--corMethod pearson \
--whatToPlot heatmap \
--skipZeros 

plotCorrelation \
--corData results.npz \
--plotFile bigwig_correlation_scatterplot.pdf \
--corMethod pearson \
--whatToPlot scatterplot \
--skipZeros 

rm -f results.npz
