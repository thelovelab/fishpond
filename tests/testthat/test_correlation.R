context("correlation")
library(SummarizedExperiment)
library(fishpond)

test_that("swish can detect correlations", {

  set.seed(1)

  y <- makeSimSwishData(null=TRUE)
  y <- scaleInfReps(y, quiet=TRUE)
  y <- labelKeep(y)

  y <- swish(y, x="condition")
  plot(-log10(mcols(y)$pvalue[1:30]))
  
})
