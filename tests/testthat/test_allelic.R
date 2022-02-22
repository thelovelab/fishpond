context("ase")
library(SummarizedExperiment)
library(fishpond)

test_that("swish for ASE works", {

  set.seed(1)

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
