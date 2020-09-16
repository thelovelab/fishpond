context("ase")
library(SummarizedExperiment)
library(fishpond)

test_that("swish for ASE works", {

  set.seed(1)
  
  y <- makeSimSwishData(n=20)
  y <- scaleInfReps(y, quiet=TRUE)
  y <- labelKeep(y)
  y$pair <- rep(1:10,2)

  y <- swish(y, x="condition", pair="pair")
  hist(mcols(y)$pvalue)
  z <- mcols(y)$stat
  hist(z[abs(z)<20], breaks=100)
  
})
