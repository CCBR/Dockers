Irreproducible Discovery Rate (IDR)
===

<p align="justify">The IDR (Irreproducible Discovery Rate) framework is a uniﬁed approach to measure the reproducibility of ﬁndings identiﬁed from replicate experiments and provide highly stable thresholds based on reproducibility. Unlike the usual scalar measures of reproducibility, the IDR approach creates a curve, which quantitatively assesses when the ﬁndings are no longer consistent across replicates. In layman's terms, the IDR method compares a pair of ranked lists of identifications (such as ChIP-seq peaks). These ranked lists should not be pre-thresholded i.e. they should provide identifications across the entire spectrum of high confidence/enrichment (signal) and low confidence/enrichment (noise). The IDR method then fits the bivariate rank distributions over the replicates in order to separate signal from noise based on a defined confidence of rank consistency and reproducibility of identifications i.e the IDR threshold.</p>

<p align="justify">The method was developed by <a href="http://www.personal.psu.edu/users/q/u/qul12/index.html">Qunhua Li</a> and <a href="http://www.stat.berkeley.edu/~bickel/">Peter Bickel</a>'s group and is extensively used by the ENCODE and modENCODE  projects and is part of their ChIP-seq guidelines and standards.</p>


Installation
------------

* Get the current repo
```
wget https://github.com/kundajelab/idr/archive/2.0.4.zip
```

* Install the dependencies
- python3
- python3 headers 
- numpy
- setuptools
- matplotlib (only required for plotting the results)

In Ubuntu 14.04+ one can run: 
(sudo) apt-get install python3-dev python3-numpy python3-setuptools python3-matplotlib

In a shared environment, the dependencies and idr package may need to be installed locally. [Anaconda](http://continuum.io/downloads#py34) largely automates this process. To install anaconda, which includes all the neessary dependencies:

```
Download [Anaconda3-2.2.0-Linux-x86_64.sh](http://continuum.io/downloads#py34) 
bash Anaconda3-2.2.0-Linux-x86_64.sh
```

* Download and unzip the idr code
```
wget https://github.com/kundajelab/idr/archive/2.0.4.zip
unzip 2.0.4.zip
cd 2.0.4/
```

* Then install idr 
```
python3 setup.py install
```

Usage
-----

List all the options
 
```
idr -h
```

Sample idr run using test peak files in the repo

```
idr --samples ../idr/test/data/peak1 ../idr/test/data/peak2
```

Run idr using an oracle peak list (e.g. peaks called from merged replicates):

```
idr --samples ../idr/test/data/peak1 ../idr/test/data/peak2 --peak-list ../idr/test/data/merged_peaks
```

### Peak matching

The method in which peaks are matched can significantly affect the output. We have chosen defaults that we believe are reaosnable in the vast majoroity of cases, but it may be worth exploring the various options for your data set.

* --peak-list is *not* provided

Peaks are grouped by overlap and then merged. The merged peak aggregate value is determined by --peak-merge-method. 

Peaks that don't overlap another peak in every other replicate are not included unless --use-nonoverlapping-peaks is set. 

* --peak-list *is* provided 

For each oracle peak a single peak from each replicate is chosen that overlaps the oracle peak. If there are multiple peaks that overlap the oracle, then ties are broken by applying the following criteria in order: 1) choose the replicate peak with a summit closest to the oracle peak's summit 2) choose the replicate peak that has the largest overlap with the oracle peak 3) choose the replicate peak with the highest score

Output
------

### Output file format
The output format mimics the input file type, with some additional fields. 

We provide an example for narrow peak files - note that the first 6 columns
are a standard bed6, the first 10 columns are a standard narrowPeak. Also, 
for columns 7-10, only the score that the IDR code used for rankign 
will be set - the remaining two columns will be set to -1.

Broad peak 
output files are the same *except* that they do not include the the summit 
columns (e.g. columns 10, 18, and 22 for samples with 2 replicates)

1.  chrom             string  
Name of the chromosome for common peaks

2.  chromStart        int     
The starting position of the feature in the chromosome or scaffold for common 
peaks, shifted based on offset. The first base in a chromosome is numbered 0.

3.  chromEnd          int     
The ending position of the feature in the chromosome or scaffold for common 
peaks. The chromEnd base is not included in the display of the feature.

4.  name              string  
Name given to a region (preferably unique) for common peaks. Use '.' 
if no name is assigned.

5.  score             int     
Contains the scaled IDR value, min(int(log2(-125*IDR), 1000). e.g. peaks with 
an IDR of 0 have a score of 1000, idr 0.05 have a score of int(-125*log2(0.05))
= 540, and idr 1.0 has a score of 0.

6.  strand         [+-.]   Use '.' if no strand is assigned.

7.  signalValue       float   
Measurement of enrichment for the region for merged peaks. When a peak list is provided this is the value from the peak list.

8.  p-value           float   
Merged peak p-value. When a peak list is provided this is the value from the peak list.

9.  q-value           float   
Merged peak q-value. When a peak list is provided this is the value from the peak list.

10. summit            int     
Merged peak summit

11. localIDR          float 
-log10(Local IDR value)

12. globalIDR         float 
-log10(Global IDR value)

13. rep1_chromStart   int     
The starting position of the feature in the chromosome or scaffold for common 
replicate 1 peaks, shifted based on offset. The first base in a chromosome is 
numbered 0.

14. rep1_chromEnd     int     
The ending position of the feature in the chromosome or scaffold for common 
replicate 1 peaks. The chromEnd base is not included in the display of the 
feature.

15. rep1_signalValue  float   
Signal measure from replicate 1. Note that this is determined by the --rank 
option. e.g. if --rank is set to signal.value, this corresponds to the 7th 
column of the narrowPeak, whereas if it is set to p.value it corresponds to
the 8th column. 

16. rep1_summit       int     
The summit of this peak in replicate 1. 

[rep 2 data]

...

[rep N data]


### Plot output

Upper Left: 
Replicate 1 peak ranks versus replicate 2 peak ranks - peaks that do not pass the specified idr threshold are colered red.

Upper Right: 
Replicate 1 log10 peak scores versus replicate 2 log10 peak scores - peaks that do not pass the specified idr threshold are colered red.

Bottom Row: 
Peaks rank versus idr scores are plotted in black. The overlayed boxplots display the distribution of idr values in each 5% quantile. The idr values are thresholded at the optimization precision - 1e-6 bny default.

Command Line Arguments
----------------------
`````
optional arguments:
  -h, --help            show this help message and exit
  --samples SAMPLES SAMPLES, -s SAMPLES SAMPLES
                        Files containing peaks and scores.
  --peak-list PEAK_LIST, -p PEAK_LIST
                        If provided, all peaks will be taken from this file.
  --input-file-type {narrowPeak,broadPeak,bed,gff}
                        File type of --samples and --peak-list.
  --rank RANK           Which column to use to rank peaks.  
                        Options: signal.value p.value q.value columnIndex
                        Defaults:
                          narrowPeak/broadPeak: signal.value
                          bed: score
  --output-file OUTPUT_FILE, -o OUTPUT_FILE
                        File to write output to.
                        Default: idrValues.txt
  --output-file-type {narrowPeak,broadPeak,bed}
                        Output file type. Defaults to input file type when available, otherwise bed.
  --log-output-file LOG_OUTPUT_FILE, -l LOG_OUTPUT_FILE
                        File to write output to. Default: stderr
  --idr-threshold IDR_THRESHOLD, -i IDR_THRESHOLD
                        Only return peaks with a global idr threshold below this value.
                        Default: report all peaks
  --soft-idr-threshold SOFT_IDR_THRESHOLD
                        Report statistics for peaks with a global idr below this value but return all peaks with an idr below --idr.
                        Default: 0.05
  --use-old-output-format
                        Use old output format.
  --plot                Plot the results to [OFNAME].png
  --use-nonoverlapping-peaks
                        Use peaks without an overlapping match and set the value to 0.
  --peak-merge-method {sum,avg,min,max}
                        Which method to use for merging peaks.
                          Default: 'sum' for signal/score/column indexes, 'min' for p/q-value.
  --initial-mu INITIAL_MU
                        Initial value of mu. Default: 0.10
  --initial-sigma INITIAL_SIGMA
                        Initial value of sigma. Default: 1.00
  --initial-rho INITIAL_RHO
                        Initial value of rho. Default: 0.20
  --initial-mix-param INITIAL_MIX_PARAM
                        Initial value of the mixture params. Default: 0.50
  --fix-mu              Fix mu to the starting point and do not let it vary.
  --fix-sigma           Fix sigma to the starting point and do not let it vary.
  --dont-filter-peaks-below-noise-mean
                        Allow signal points that are below the noise mean (should only be used if you know what you are doing).
  --use-best-multisummit-IDR
                        Set the IDR value for a group of multi summit peaks (a group of peaks with the same chr/start/stop but different summits) to the best value across all of these peaks. This is a work around for peak callers that don't do a good job splitting scores across multi summit peaks (e.g. MACS). If set in conjunction with --plot two plots will be created - one with alternate summits and one without.  Use this option with care.
  --allow-negative-scores
                        Allow negative values for scores. (should only be used if you know what you are doing)
  --random-seed RANDOM_SEED
                        The random seed value (sor braking ties). Default: 0
  --max-iter MAX_ITER   The maximum number of optimization iterations. Default: 3000
  --convergence-eps CONVERGENCE_EPS
                        The maximum change in parameter value changes for convergence. Default: 1.00e-06
  --only-merge-peaks    Only return the merged peak list.
  --verbose             Print out additional debug information
  --quiet               Don't print any status messages
  --version             show program's version number and exit
`````


Contributors
------------

The main contributors of IDR code:

  * Nathan Boley        - Kundaje Lab, Dept. of Genetics, Stanford University
  * Anshul Kundaje      - Assistant Professor, Dept. of Genetics, Stanford University
  * Peter J. Bickel     - Professor, Dept. of Statistics, University of California at Berkeley
  * Jin Lee             - Software Developer, Dept. of Genetics, Stanford University

References
----------
"Measuring reproducibility of high-throughput experiments" (2011), Annals of Applied Statistics, Vol. 5, No. 3, 1752-1779, by Li, Brown, Huang, and Bickel

Issues
------

If you notice any problem with the code, please file an issue over [here](https://github.com/kundajelab/idr/issues)
