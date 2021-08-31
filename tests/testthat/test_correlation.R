context("correlation")
library(SummarizedExperiment)
library(fishpond)

test_that("swish can detect correlations", {

  set.seed(1)

  y <- makeSimSwishData(n=20, null=TRUE)
  nms <- c("counts",paste0("infRep",1:20))
  lambda1 <- seq(40, 80, length.out=20)
  lambda2 <- rev(lambda1)
  for (a in nms) {
    assays(y)[[a]][1,] <- rpois(20,lambda1)
    assays(y)[[a]][2,] <- rpois(20,lambda2)
  }
  y$condition <- sort(runif(20))
  y <- scaleInfReps(y, quiet=TRUE)
  y <- labelKeep(y)

  y <- swish(y, x="condition", cor="spearman")
  hist(mcols(y)$pvalue)
  plot(-log10(mcols(y)$pvalue[1:30]))
  
})
