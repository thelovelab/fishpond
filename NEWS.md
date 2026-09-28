# fishpond 2.8.0

* Corrected a bug in the two-group interaction (without pairing)
  functionality, when the groups were imbalanced, as identified by
  Samuel Chen. Fixes GitHub issue #35.

# fishpond 2.7.1

* Corrected a bug in the two-group interaction (without pairing)
  functionality, when the groups were imbalanced, as identified by
  Samuel Chen. Fixes GitHub issue #35.

# fishpond 2.5.4

* As CellRanger 7 includes both spliced and unspliced counts in their
  count matrix, we want to mimic this behavior by adding more
  pre-defined output formats in the loadFry function. We added "all"
  and "U+S+A" to include all counts in the count matrix. Moreover, now
  the "scRNA" output format has an "unspliced" field, which contains
  the unspliced count matrix.

# fishpond 2.5.3

* Fix bug where salmonEC did not correct equivalence class
  names for going from 0-indexing to 1-indexing internally. In prior  
  versions, to correctly link equivalence classes to gene names,
  users would have needed to manually add a value of 1 to the 
  equivalence class names, which was erroneously not mentioned
  in the man files. After this bug fix, if the equivalence class 
  identifier reads 1|2|8, then the equivalence class is immediatly 
  compatible with the transcripts and their respective genes in rows
  1, 2 and 8 of 'tx2gene_matched', without any further user 
  intervention.

# fishpond 2.5.1

* Fix plotAllelicGene() so that when samples have an allele
  with no expression at the gene level, it doesn't throw
  an error trying to divide by 0.

# fishpond 2.4.0

* For simple paired swish analysis, adding a `fast=1` method 
  which uses a one-sample z-score on the paired LFCs
  (averaged over samples, then median over inferential 
  replicates). The permutation is computed by changing the
  signs of the LFC matrix and recomputing z-scores.
  Testing on the vignette example, but using all the transcripts,
  the one-sample z-score method takes < 20 seconds while the
  signed rank method takes > 200 seconds (12x speedup), 
  while they have a high rate of agreement on the detected set 
  (30:1 in common vs discordant).
* Removed the fast=0 methods that were previously implemented
  where ranks could optionally be recomputed for every
  permutation. This was much slower and didn't have any
  appreciable benefit.
* readEDS() has moved to the `eds` package, such that
  fishpond no longer requires *Rcpp* and C++ code compilation.
* Fix bug identified by GitHub user @JosephLalli, where
  importAllelicCounts would find a1 and a2 strings
  internal to gene IDs, instead of at the suffix.

# fishpond 2.3.22

* Fix bug identified by GitHub user @JosephLalli, where
  importAllelicCounts would find a1 and a2 strings
  internal to gene IDs, instead of at the suffix.

# fishpond 2.3.14

* For simple paired swish analysis, adding a `fast=1` method 
  which uses a one-sample z-score on the paired LFCs
  (averaged over samples, then median over inferential 
  replicates). The permutation is computed by changing the
  signs of the LFC matrix and recomputing z-scores.
  Testing on the vignette example, but using all the transcripts,
  the one-sample z-score method takes <20 seconds while the
  signed rank method takes >200 seconds (12x speedup), 
  while they have a high rate of agreement on the detected set (
  30:1 in common vs discordant).
* Removed the fast=0 methods that were previously implemented
  where ranks could optionally be recomputed for every
  permutation. This was much slower and didn't have any
  appreciable benefit.

# fishpond 2.3.7

* readEDS() has moved to the `eds` package.

# fishpond 2.2.0

* New vignette demonstrating allelic analysis at isoform, 
  TSS, or gene-level. See more details below.
* New import functions for equivalence class analysis of Salmon
  or alevin data, written by Jeroen Gilis. See salmonEC() and
  alevinEC() man pages.
* New plotAllelicGene() and plotAllelicHeatmap() functions
  for plotting results from allelic expression analysis.
* New makeTx2Tss() helper function for allelic analysis.
* Now importFromAllelicCounts() can take a GRanges
  object as the `tx2gene` argument, so that ranges will
  be distributed to the rowRanges of the output 
  SummarizedExperiment.
* Adding `shiftX` argument to plotInfReps() for numeric
  x variable, to help with overplotting.
* Re-organized package for new pkgdown homepage:
  https://mikelove.github.io/fishpond

# fishpond 2.0.0

* New loadFry() function, written by Dongze He with
  contributions from Steve Lianoglou and Wes Wilson.
  loadFry() helps users to import and process
  alevin-fry quantification results. Can process
  spliced, unspliced and ambiguous counts separately 
  and flexibly. Has specific output formats designed
  for downstream use with scVelo or velocity analysis.
  See ?loadFry for more details.
* Adding correlation tests: Spearman or Pearson
  correlations of a numeric covariate with the
  log counts, or with the log fold changes across
  pairs. The Spearman correlation test with counts
  was already implemented in the original SAMseq
  method as response type = "Quantitative".
  For new functionality see 'cor' argument in the
  ?swish man page.
* Adding importAllelicCounts() to facilitate importing
  Salmon quantification data against a diploid
  transcriptome. Can import either as a 'wide'
  format or as 'assays'. Leverages tximeta().
  For gene-level summarization, importAllelicCounts()
  can create an appropriate tx2gene table
  with the necessary a1 and a2 suffices,
  and it will automatically set txOut=FALSE, see
  ?importAllelicCounts for more details.
* Added a 'q' argument to plotInfReps to change the
  intervals when making point and line plots.
* Switched the legend of plotInfReps so that
  reference levels will now be on the bottom,
  and non-reference (e.g. treatment) on top.

# fishpond 1.99.18

* Added helper functionality to importAllelicCounts,
  so it will create an appropriate tx2gene table
  with the necessary a1 and a2 suffices,
  and it will automatically set txOut=FALSE.
* Added a 'q' argument to plotInfReps to change the
  intervals when making point and line plots.
* Switched the legend of plotInfReps so that
  reference levels will now be on the bottom,
  and non-reference (e.g. treatment) on top.
* Added loadFry() to process alevin-fry 
  quantification result. Can process spliced, 
  unspliced and ambiguous counts separately 
  and flexibly.

# fishpond 1.99.15

* Adding correlation tests: Spearman or Pearson
  correlations of a numeric covariate with the
  log counts, or with the log fold changes across
  pairs. The Spearman correlation test with counts
  was already implemented in the original SAMseq
  method as response type = "Quantitative".
  For new functionality see 'cor' argument in the
  ?swish man page.
* Adding importAllelicCounts() to facilitate importing
  Salmon quantification data against a diploid
  transcriptome. Can import either as a 'wide'
  format or as 'assays'. Leverages tximeta().

# fishpond 1.9.6

* Specifying ties.method in matrixStats::rowRanks.

# fishpond 1.9.1

* Added importAllelicCounts() with options for importing
  Salmon quantification on diploid transcriptomes.

# fishpond 1.8.0

* Added note in vignette about how to deal with estimated
  batch factors, e.g. from RUVSeq or SVA. Two strategies are
  outlined: either discretizing the estimate batch factors
  and performing stratified analysis, or regressing out the
  batch-associated variation using limma's removeBatchEffect.
  Demonstation code is included.

# fishpond 1.6.0

* Added makeInfReps() to create pseudo-inferential replicates
  via negative binomial simulation from mean and variance
  matrices. Note: the mean and the variance provide the
  _inferential_ distribution per element of the count matrix.
  See preprint for details, doi: 10.1101/2020.07.06.189639.
* Added splitSwish() and addStatsFromCSV(), which can be used
  to distribute running of Swish across a number of jobs
  managed by `Snakemake`. See vignette for description of
  a suggested workflow. For a large single-cell dataset
  with mean and variance summaries of inferential uncertainty,
  splitSwish() avoids generating the inferential replicate
  counts until the data has been split into smaller pieces and
  sent to different jobs, then only the necessary summary
  statistics are gathered and q-values computed by
  addStatsFromCSV().
* plotInfReps() gains many new features to facilitate plotting of
  inferential count distributions for single cells, as quantified
  with alevin and imported with tximport. E.g. allow for numeric
  `x` argument plus grouping with `cov` for showing
  counts over pseudotime across groups of cells. Also added
  `applySF` argument which can be used to divide out a
  size factor, and the `reorder` argument which will re-order
  the samples/cells within groups by the count. plotInfReps()
  will draw boxplots with progressively thinner visual features
  as the number of cells grows to make the plots still legible.

# fishpond 1.5.2

* First version of makeInfReps(), to create pseudo-infReps
  via negative binomial simulation from set of mean and
  variance matrices in the assays of the SummarizedExperiment.

# fishpond 1.4.0

* Added isoformProportions(), which can be run after
  scaleInfReps() and optionally after filtering out
  transcripts using labelKeep(). Running swish() after
  isoformProportions() will produce differential transcript
  usage (DTU) results, instead of differential transcript
  expression (DTE) results. Example in vignette.
* Default number of permutations increased from 30 to 100.
  It was observed that there was too much fluctuation in the
  DE called set for nperms=30 across different seeds, and
  setting to 100 helped to stabilize results across seeds,
  without increasing running time too much. For further reduced
  dependence on the seed, even higher values of nperms
  (e.g. 200, 300) can be used.

# fishpond 1.3.8

* Added isoformProportions(), which can be run after
  scaleInfReps() and optionally after filtering out
  transcripts using labelKeep(). Running swish() after
  isoformProportions() will produce differential transcript
  usage (DTU) results, instead of differential transcript
  expression (DTE) results. Example in vignette.

# fishpond 1.3.4

* Default number of permutations increased from 30 to 100.
  It was observed that there was too much fluctuation in the
  DE called set for nperms=30 across different seeds, and
  setting to 100 helped to stabilize results across seeds,
  without increasing running time too much. For further reduced
  dependence on the seed, even higher values of nperms
  (e.g. 200, 300) can be used.

# fishpond 1.2.0

* Switching to a faster version of Swish which only
  computes the ranks of the data once, and then re-uses
  this for the permutation distribution. This bypasses
  the addition of uniform noise per permutation and
  is 10x faster. Two designs which still require
  re-computation of ranks per permutation are the
  paired analysis and the general interaction analysis.
  Two-group, stratified two-group, and the paired
  interaction analysis now default to the new fast
  method, but the original, slower method can be used
  by setting fast=0 in the call to swish().
* Adding Rcpp-based function readEDS() written by
  Avi Srivastava which imports the sparse counts stored
  in Alevin's Efficient Data Storage (EDS) format.
* Changed the vignette so that it (will) use a linkedTxome,
  as sometime the build would break if the Bioc build
  machine couldn't access ftp.ebi.ac.uk.
* Add 'computeInfRV' function. InfRV is not used in the
  Swish methods, only for visualization purposes in the
  Swish paper.
* removed 'samr' from Imports, as it required source
  installation, moved to Suggests, for optional qvalue
  calculation

# fishpond 0.99.30

* added two interaction tests, described in ?swish
* incorporate qvalue package for pvalue, locfdr and qvalue
* added plotMASwish() to facilitate plotting
* wilcoxP is removed, and the mean is used instead

# fishpond 0.99.0

* fishpond getting ready for submission to Bioc
