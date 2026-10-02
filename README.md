# badger2
BADGER2: Bivariate Analysis of Differential Gene Expression Reactions 2

**Authors:** Ellen I. Weinheimer, James B. Pease

**Description:** BADGER2 is used for detection of differentially reacting genes across multi-treatment gene expression time courses. The pipeline can be applied to time courses of any length or number of treatments, but requires at least two of each. 
BADGER2 adds an unpaired mode

**Python Non-Standard Library Requirements:**  numpy, scipy

Basic Usage: BADGER requires as input (1) a csv count table of normalized RNAseq reads (`--input`); (2) a csv file summarizing the experimental design factors with one sample per line and columns indicating treatment, timepoint, individual, and sample ID (`--metadata`), and (3) the contrasts for your metadata header variables as `--contrast-a` and `--contrast-b`.  

For further details run the script with `--help` or see the example files in `test_files`.

If you use this software: 
- Please include the URL [https://github.com/peaselab/badger2].
- Please cite Weinheimer, et al. 2025 *The Plant Journal* [https://onlinelibrary.wiley.com/doi/10.1111/tpj.70385]

Questions? Please contact James Pease at [http://u.osu.edu/peaselab]
