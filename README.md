# fishpond <img id="fishpond_logo" src="man/figures/fishpond.png" align="right" width="125"/>

[![R build status](https://github.com/mikelove/fishpond/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/mikelove/fishpond/actions/workflows/check-bioc.yml)

## Fishpond: downstream methods and tools for expression data

Fishpond contains a method, `swish()`, for differential transcript and
gene expression analysis of RNA-seq data using inferential replicates.
Also the package contains utilities for working with *Salmon*,
*alevin*, and *alevin-fry* quantification data, including
`loadFry()`.

## Quick start

The following paradigm is used for running a Swish analysis:

```
y <- tximeta(coldata) # reads in counts and inf reps
y <- scaleInfReps(y) # scales counts
y <- labelKeep(y) # labels features to keep
set.seed(1) # for reproducibility
y <- swish(y, x="condition") # simplest Swish case
```

## How does Swish work

Swish accounts for inferential uncertainty in expression estimates
by averaging test statistics over a number of *inferential replicate*
datasets, either posterior samples or bootstrap samples. This is
inspired by a method called 
[SAMseq](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4605138/), 
hence we named our method *Swish*, for "SAMseq With Inferential
Samples Helps". Averaging over inferential replicates produces a
different test statistic than what one would obtain using only point
estimates for expression level.

For example, one of the tests possible with `swish()` is a correlation
test of expression level over a condition variable. We can visualize
the distribution of inferential replicates with `plotInfReps()`:

![](man/figures/plotInfReps.png)

The test statistic is formed by averaging over these sets of data:

![](man/figures/swish.gif)

p-values and q-values are computed through permutation of samples (see
vignette for details on permutation schemes). 

The *Swish* method is described in the following publication:

> Zhu, A., Srivastava, A., Ibrahim, J.G., Patro, R., Love, M.I. 
> "Nonparametric expression analysis using inferential replicate counts" 
> *Nucleic Acids Research* (2019) 47(18):e105
> [PMC6765120](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6765120/)

The *SEESAW* method for allelic expression analysis is described in
the following preprint:

> Euphy Wu, Noor P. Singh, Kwangbom Choi, Mohsen Zakeri, Matthew
> Vincent, Gary A. Churchill, Cheryl L. Ackert-Bicknell, Rob Patro,
> Michael I. Love.
> "Detecting isoform-level allelic imbalance accounting for
> inferential uncertainty" *bioRxiv* (2022)
> [doi: 10.1101/2022.08.12.503785](https://doi.org/10.1101/2022.08.12.503785)

## Installation

This package can be installed via Bioconductor:

```
BiocManager::install("fishpond")
```

## Funding

This work was funded by NIH NHGRI R01-HG009937.
