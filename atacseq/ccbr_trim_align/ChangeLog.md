# Change Log

- v0.1.1
  - adding *threads* argument to all scripts in TAD ... nproc was always 1 when running with sbatch... hence the change
  - all output files written to /lscratch/$SLURM_JOB and the required file then moved to appropriate folder
  - *qc* will have subfolders eg. *fastqc*
  - *tagAlign.gz* files will be moved to its own folder *tagAlign*
- v0.1.2
  - intermediate alignment file (containing multimapped reads, prior to MAPQ-based filtering, sorted by name) is saved as input for Genrich
  - *--paired-end* argument added to *atac_assign_multimappers* script
- v0.1.3
  - Added *bedGraphToBigwig* from *http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig*
  - *ccbr_macs2_peak_calling* added
    - Macs2 bdg files can be saved as bigwig
    - peaks can be filtered using qvalue (default 0.05)
  - *bc* installed ... required for macs2 two replicate wrapper script
  - *ccbr_macs2_peak_calling_two_replicates.bash* added
  - some peaks called by macs have the exact same start and end coordinates but are label a, b, etc. with different p-/q-values.... we now sort them by q-value and only pick the first peak call. These filtered peakcalls are fed into idr. The idr output is also filtered to find peaks with unique **chr:start-end** string
  - 