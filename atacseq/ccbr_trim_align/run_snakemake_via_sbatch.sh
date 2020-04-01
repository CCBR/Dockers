#!/bin/bash

module load python/3.5
module load snakemake/5.1.3
#cd $SLURM_SUBMIT_DIR
#R="$SLURM_SUBMIT_DIR"
R=$(pwd)
echo $R
if [ ! -d Reports ]; then mkdir Reports;fi
if [ ! -d slurmfiles ]; then mkdir slurmfiles;fi

snakemake $@ -s $R/atacseqpipeline.snakemake \
-d $R --printshellcmds --latency-wait 120 --cluster-config $R/cluster.json \
--cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem}" \
-j 500 --rerun-incomplete --keep-going --restart-times 0 --stats $R/Reports/atacseq.stats -T 2>&1|tee -a $R/Reports/snakemake.log

mv slurm*out slurmfiles
