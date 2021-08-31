context("ase")
library(SummarizedExperiment)
library(fishpond)

test_that("swish for ASE works", {

  set.seed(1)

  # looking at how the null distribution varies by sample size
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

  # test reading in wide format for ASE
  if (FALSE) {
    # this needs to be added as test data somewhere...
    dir <- "../../../ase_quants"
    names <- list.files(dir)
    files <- file.path(dir, names, "quant.sf")
    coldata <- data.frame(files, names, condition=factor(c("A","A","B","B")))
    suppressPackageStartupMessages(library(SummarizedExperiment))
    library(tximeta)
    se <- importAllelicCounts(coldata, a1="P", a2="M", format="wide")
    colData(se)
    assayNames(se)
    metadata(se)$alleles
    se <- importAllelicCounts(coldata, a1="P", a2="M", format="assays")
    colData(se)
    assayNames(se)
    metadata(se)$alleles
  }
    
})
