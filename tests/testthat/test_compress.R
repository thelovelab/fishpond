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

  set.seed(1)
  y <- makeSimSwishData(m=200, n=100, meanVariance=TRUE)
  y$batch <- factor(rep(3:1,c(30,40,30)))
  plotInfReps(y, idx=5, x="condition", reorder=TRUE)
  dev.off()
  plotInfReps(y, idx=5, x="condition", cov="batch", reorder=TRUE)
  dev.off()
  
  y2 <- y
  y2 <- makeInfReps(y2, numReps=50)
  y2 <- scaleInfReps(y2, quiet=TRUE)
  y2 <- labelKeep(y2)
  y2 <- swish(y2, x="condition", quiet=TRUE)
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

  # poscounts not recommended actually... use scran::calculateSumFactors()
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

    # convert to SCE
    library(SingleCellExperiment)
    sce <- as(se, "SingleCellExperiment")
    
    sce$condition <- factor(sample(rep(1:2,length=ncol(sce))))
    time <- rnorm(ncol(sce))
    time[sce$condition == 1] <- sort(time[sce$condition == 1])
    time[sce$condition == 2] <- sort(time[sce$condition == 2], decreasing=TRUE)
    sce$pseudotime <- time
    sce <- labelKeep(sce, minCount=10, minN=10)
    table(mcols(sce)$keep)
    sce <- sce[mcols(sce)$keep,]
    mcols(sce)$name <- rep(LETTERS, length=nrow(sce))
    plotInfReps(sce, 1, x="condition")
    plotInfReps(sce, 1, x="pseudotime")
    plotInfReps(sce, 1, x="pseudotime", cov="condition")
    plotInfReps(sce, 1, x="condition", mainCol="name")

    # DESeq2 "poscounts" normalization
    sfFun <- function(m) {
      geoMeanNZ <- function(x) {
        if (all(x == 0)) { 0 } else { exp( sum(log(x[x > 0])) / length(x) ) }
      }
      geoMeans <- apply(m, 1, geoMeanNZ)
      DESeq2::estimateSizeFactorsForMatrix(m, geoMeans=geoMeans)
    }
    sf <- sfFun(as.matrix(assays(sce)[["mean"]]))

    # set size factors using SingleCellExperiment setter
    sizeFactors(sce) <- pmax(sf, .25)

    plotInfReps(sce, 1, x="condition", applySF=TRUE)
    plotInfReps(sce, 1, x="condition", reorder=FALSE)
    plotInfReps(sce, 1, x="condition", applySF=TRUE, reorder=FALSE)
    
  }
  
})
