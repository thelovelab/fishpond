context("ase")
library(SummarizedExperiment)
library(fishpond)

test_that("swish for ASE works", {

  set.seed(1)

  par(mfrow=c(3,1))
  for (np in c(8,10,12)) {
    y <- makeSimSwishData(n=np*2)
    y <- scaleInfReps(y, quiet=TRUE)
    y <- labelKeep(y)
    y$pair <- rep(1:np,2)
    y <- swish(y, x="condition", pair="pair")
    z <- mcols(y)$stat
    hist(z[abs(z)<20], breaks=100, xlim=c(-30,30))
  }
  
})
