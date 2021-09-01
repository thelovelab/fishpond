context("swish")
library(SummarizedExperiment)
library(fishpond)

test_that("basic variable errors thrown", {

  set.seed(1)
  
  y <- makeSimSwishData()
  y <- scaleInfReps(y, quiet=TRUE)
  y <- labelKeep(y)

  # too many levels of condition
  y2 <- y
  y2$condition <- gl(3,3,10)
  expect_error(swish(y2, "condition"))

  # batch and pair together
  y2 <- y
  y2$batch <- rep(1:2,5)
  y2$pair <- rep(1:5,2)
  expect_error(swish(y2, "condition", "batch", "pair"))

  # wrong number of pairs
  y2 <- y
  y2$pair <- c(1,2,3,4,5,1,1,2,2,3)
  expect_error(swish(y2, "condition", pair="pair"), "single sample for both levels")
  
  # no inferential replicates
  y <- makeSimSwishData()
  assays(y) <- assays(y)[c("counts","abundance","length")]
  expect_error(scaleInfReps(y), "no inferential")
  expect_error(swish(y, "condition"), "no inferential")

  # too many permutations requested
  y <- makeSimSwishData(m=100, n=4)
  y <- scaleInfReps(y, quiet=TRUE)
  y <- labelKeep(y)
  # there are 4! = 24 permutations
  expect_message(swish(y, x="condition", nperms=25, quiet=TRUE), "less permutations")

  # too few samples
  y <- makeSimSwishData(m=100, n=2)
  y <- scaleInfReps(y, quiet=TRUE)
  y2 <- labelKeep(y)
  expect_error(swish(y2, x="condition"), "all rows")
  y <- labelKeep(y, minN=2)
  expect_error(swish(y, x="condition"), "too few samples")
  y <- makeSimSwishData(m=100, n=6)
  y$batch <- factor(c(1:3,1:3))
  y <- scaleInfReps(y, quiet=TRUE)
  y <- labelKeep(y)
  expect_error(swish(y, x="condition", cov="batch"), "too few samples")
  
})

test_that("basic swish analyses", {

  set.seed(1)
  
  # two group
  y <- makeSimSwishData(m=200)
  y <- scaleInfReps(y, quiet=TRUE)
  y <- labelKeep(y)
  y <- swish(y, x="condition", quiet=TRUE)
  expect_true("qvalue" %in% colnames(mcols(y)))

  plotInfReps(y, 1, "condition")
  dev.off()

  # differential transcript usage
  # requires a gene ID column
  mcols(y)$gene_id <- as.character(rep(1:(nrow(y)/5), each=5))
  iso <- isoformProportions(y, quiet=TRUE)
  iso <- swish(iso, x="condition", quiet=TRUE)
  
  # try the old, slow method
  y <- swish(y, x="condition", fast=0, quiet=TRUE)
  
  # estimate pi0
  y <- swish(y, x="condition", estPi0=TRUE, qvaluePkg="qvalue", quiet=TRUE)
  y <- swish(y, x="condition", estPi0=TRUE, qvaluePkg="samr", quiet=TRUE)

  # use samr for qvalue
  y <- swish(y, x="condition", qvaluePkg="samr", quiet=TRUE)

  # two group with batch covariate
  y <- makeSimSwishData(m=200, n=20)
  y$batch <- factor(rep(c(1,2,1,2),each=5))
  y <- scaleInfReps(y, quiet=TRUE)
  y <- labelKeep(y)
  y <- swish(y, x="condition", cov="batch", quiet=TRUE)
  plotInfReps(y, 1, "condition", "batch")
  dev.off()

  # try the old, slow method
  y <- swish(y, x="condition", cov="batch", fast=0, quiet=TRUE)
  
  # two group, matched samples
  y <- makeSimSwishData(m=200)
  y$pair <- rep(1:5,2)
  y <- scaleInfReps(y, quiet=TRUE)
  y <- labelKeep(y)
  y <- swish(y, x="condition", pair="pair", quiet=TRUE)

  # alternative scaling

  y <- makeSimSwishData(m=200)
  y <- scaleInfReps(y, lengthCorrect=FALSE, quiet=TRUE)
  y <- makeSimSwishData()
  y <- scaleInfReps(y, sfFun=function(x) colSums(x)/mean(colSums(x)), quiet=TRUE)
  
})

test_that("infRV calculation and plotting", {

  y <- makeSimSwishData(m=200)
  y <- computeInfRV(y)
  #mcols(y)$meanCts <- rowMeans(assays(y)[["counts"]])
  #with(mcols(y), plot(meanCts, meanInfRV))
  
})

test_that("alevin to fishpond", {

  if (FALSE) {

    dir <- system.file("extdata", package="tximportData")
    files <- file.path(dir,"alevin/neurons_900_v014/alevin/quants_mat.gz")
    file.exists(files)
    library(tximeta)
    se <- tximeta(files, type="alevin") # this makes a ~1 Gb SE
    se <- se[1:1000,1:100] # 12 Mb
    y <- se

    y <- labelKeep(y, minCount=3, minN=10)
    y <- y[mcols(y)$keep,]

    assays(y) <- lapply(assays(y), as.matrix)
    y <- scaleInfReps(y, lengthCorrect=FALSE)
    y$condition <- factor(rep(1:2, each=50))
    y <- swish(y, x="condition")
    hist(mcols(y)$pvalue)

  }
  
})

test_that("basic deswish analyses", {

  # y <- makeSimSwishData()
  # y <- labelKeep(y)
  # y <- deswish(y, ~condition, "condition_2_vs_1")
  
})

