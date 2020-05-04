context("compressing uncertainty")
library(SummarizedExperiment)
library(fishpond)

test_that("compressing uncertainty works", {

  set.seed(1)
  y <- makeSimSwishData()

  a <- assays(y)
  infRepIdx <- grep("infRep",names(a))
  infReps <- a[infRepIdx]
  infRepsCube <- abind::abind(as.list(infReps), along=3)
  a[["mean"]] <- apply(infRepsCube, 1:2, mean)
  a[["variance"]] <- apply(infRepsCube, 1:2, var)
  assays(y) <- a

  # remove inf reps
  infRepIdx <- grep("infRep",assayNames(y))
  assays(y) <- assays(y)[-infRepIdx]
  
  y <- makeInfReps(y, numReps=50)
  y <- scaleInfReps(y)
  y <- labelKeep(y)
  y <- swish(y, x="condition")
  #hist(mcols(y)$stat, breaks=40, col="grey")
  #cols = rep(c("blue","purple","red"),each=2)
  ## for (i in 1:6) {
  ##   arrows(mcols(y)$stat[i], 20,
  ##          mcols(y)$stat[i], 10,
  ##          col=cols[i], length=.1, lwd=2)
  ## }
  #plotInfReps(y, 1, x="condition")
  #plotInfReps(y, 3, x="condition")
  #plotInfReps(y, 5, x="condition")
  
})
