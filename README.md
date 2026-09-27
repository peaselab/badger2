# badger2
BADGER2: Bivariate Analysis of Differential Gene Expression Reactions 2

**Authors:** Ellen I. Weinheimer, James B. Pease

**Description:** BADGER2 is used for detection of differentially reacting genes across multi-treatment gene expression time courses. The pipeline can be applied to time courses of any length or number of treatments, but requires at least two of each. 

**Python Non-Standard Library Requirements:**  numpy, scipy

Basic Usage: BADGER requires as input (1) a csv count table of normalized RNAseq reads (--input); (2) a csv file summarizing the experimental design factors with one sample per line and columns indicating treatment, timepoint, individual, and sample ID (--factors), and (3) the desired control condition label as seen in the treatment column of the experimental design factors file (--control). See files in test_files/ and command in ./test.sh as an example.

Questions? Please contact James Pease at [http://u.osu.edu/peaselab]
