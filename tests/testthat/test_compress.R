context("compressing uncertainty")
library(SummarizedExperiment)
library(fishpond)

test_that("compressing uncertainty works", {

  set.seed(1)
  y <- makeSimSwishData(m=200, n=100)
  y$batch <- factor(rep(3:1,c(30,40,30)))
  plotInfReps(y, 5, x="condition", reorder=TRUE)
  dev.off()
  plotInfReps(y, 5, x="condition", cov="batch", reorder=TRUE)
  dev.off()

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

  plotInfReps(y, idx=5, x="condition", reorder=TRUE)
  dev.off()
  plotInfReps(y, idx=5, x="condition", cov="batch", reorder=TRUE)
  dev.off()
  
  y2 <- y
  y2 <- makeInfReps(y2, numReps=50)
  y2 <- scaleInfReps(y2)
  y2 <- labelKeep(y2)
  y2 <- swish(y2, x="condition")
  #hist(mcols(y2)$stat, breaks=40, col="grey")
  #cols = rep(c("blue","purple","red"),each=2)
  ## for (i in 1:6) {
  ##   arrows(mcols(y2)$stat[i], 20,
  ##          mcols(y2)$stat[i], 10,
  ##          col=cols[i], length=.1, lwd=2)
  ## }
  #plotInfReps(y2, 1, x="condition")
  #plotInfReps(y2, 3, x="condition")
  #plotInfReps(y2, 5, x="condition")

  sfFun <- function(m) {
    DESeq2::estimateSizeFactorsForMatrix(
      m, geoMeans=exp(rowSums(log(m) * as.numeric(m > 0))/ncol(m))
    )
  }
  sf <- sfFun(assays(y)[["counts"]])
  colData(y)$sizeFactor <- sf

  y <- labelKeep(y)
  y <- y[mcols(y)$keep,]

  #library(Matrix)
  #assays(y)[["mean"]] <- as(assays(y)[["mean"]], "sparseMatrix")
  #assays(y)[["variance"]] <- as(assays(y)[["variance"]], "sparseMatrix")
  #splitSwish(y, 4, prefix="foo/swish", snakefile="foo/Snakefile")
  #miniSwish("foo/swish1.rds", "foo/swish1.csv", x="condition")
  #y <- addStatsFromCSV(y, "foo/summary.csv")

  if (FALSE) {

    # test plotting of neurons data
    dir <- system.file("extdata", package="tximportData")
    files <- file.path(dir,"alevin/neurons_900_v014/alevin/quants_mat.gz")
    file.exists(files)
    coldata <- data.frame(files, names="neurons")
    library(SummarizedExperiment)
    library(tximeta)
    se <- tximeta(coldata, type="alevin", dropInfReps=TRUE, skipMeta=TRUE)
    se$condition <- factor(sample(rep(1:2,length=ncol(se))))
    se <- labelKeep(se, minCount=10, minN=10)
    table(mcols(se)$keep)
    se <- se[mcols(se)$keep,]
    plotInfReps(se, 1, x="condition")

    sfFun <- function(m) {
      DESeq2::estimateSizeFactorsForMatrix(
        m, geoMeans=exp(rowSums(log(m) * as.numeric(m > 0))/ncol(m))
      )
    }
    sf <- sfFun(as.matrix(assays(se)[["mean"]]))
    se$sizeFactor <- pmax(sf, .25)

    plotInfReps(se, 1, x="condition", applySF=TRUE)
    plotInfReps(se, 1, x="condition", reorder=FALSE)
    plotInfReps(se, 1, x="condition", applySF=TRUE, reorder=FALSE)
    
  }
  
})
