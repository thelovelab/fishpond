context("interaction")
library(SummarizedExperiment)
library(fishpond)

test_that("matched samples interactions work", {

  set.seed(1)
  
  # two group matched samples interaction
  y <- makeSimSwishData(m=200, n=20, null=TRUE)
  nms <- c("counts",paste0("infRep",1:20))
  lambda1 <- rep(c(40,80),length=20)
  lambda2 <- c(rep(c(40,80),length=10),rep(c(40,160),length=10))
  for (a in nms) {
    assays(y)[[a]][1,] <- rpois(20,lambda1)
    assays(y)[[a]][2,] <- rpois(20,lambda2)
  }
  y$condition <- factor(rep(1:2,length=20))
  y$pair <- factor(rep(1:10,each=2))
  y$group <- factor(rep(1:2,each=10))

  y <- scaleInfReps(y, quiet=TRUE)
  y <- labelKeep(y)
  y <- swish(y, x="condition", cov="group", pair="pair", interaction=TRUE, quiet=TRUE)

  expect_true(mcols(y)$pvalue[2] < .01)

})

test_that("two group interactions work", {

  set.seed(1)
  
  # two group interaction
  y <- makeSimSwishData(m=200, n=20, null=TRUE)
  nms <- c("counts",paste0("infRep",1:20))
  lambda1 <- rep(c(40,80,40,80),each=5)
  lambda2 <- rep(c(40,80,40,160),each=5)
  for (a in nms) {
    assays(y)[[a]][1,] <- rpois(20,lambda1)
    assays(y)[[a]][2,] <- rpois(20,lambda2)
  }
  y$condition <- factor(rep(c(1,2,1,2),each=5))
  y$group <- factor(rep(1:2,each=10))

  y <- scaleInfReps(y, quiet=TRUE)
  y <- labelKeep(y)
  y <- swish(y, x="condition", cov="group", interaction=TRUE, quiet=TRUE)

  expect_true(mcols(y)$pvalue[2] < .01)

  # two groups with imbalanced sample sizes
  y <- makeSimSwishData(m=200, n=20, null=TRUE)
  nms <- c("counts",paste0("infRep",1:20))
  lambda1 <- rep(c(40,80,40,80),c(4,6,6,4))
  lambda2 <- rep(c(40,80,40,160),c(4,6,6,4))
  for (a in nms) {
    assays(y)[[a]][1,] <- rpois(20,lambda1)
    assays(y)[[a]][2,] <- rpois(20,lambda2)
  }
  y$condition <- factor(rep(c(1,2,1,2),c(4,6,6,4)))
  y$group <- factor(rep(1:2,each=10))
  table(y$condition, y$group)

  y <- scaleInfReps(y, quiet=TRUE)
  y <- labelKeep(y)
  y <- swish(y, x="condition", cov="group", interaction=TRUE, quiet=TRUE)

  expect_true(mcols(y)$pvalue[2] < .01)
  
})
