context("ase")
library(SummarizedExperiment)
library(fishpond)

test_that("fishpond code for ASE works", {

  set.seed(1)

  # test reading in wide format for ASE
  if (FALSE) {
    dir <- "../../../../ase_quants"
    names <- list.files(dir, pattern="sample")
    files <- file.path(dir, names, "quant.sf")
    coldata <- data.frame(files, names, condition=factor(c("A","B")))
    suppressPackageStartupMessages(library(SummarizedExperiment))
    library(tximeta)
    se <- importAllelicCounts(coldata, a1="P", a2="M", format="wide")
    colData(se)
    assayNames(se)
    metadata(se)$alleles

    # check TSS aggregation
    library(ensembldb)
    library(AnnotationHub)
    ah <- AnnotationHub()
    query(ah, c("Drosophila melanogaster", "release-100"))
    Gtf <- ah["AH80008"]
    DbFile <- ensDbFromAH(Gtf)
    edb <- EnsDb(DbFile)
    t2t <- makeTx2Tss(edb)

    se <- importAllelicCounts(coldata,
                              a1="P", a2="M",
                              format="wide",
                              tx2gene=t2t,
                              ignoreAfterBar=TRUE)
    gr <- rowRanges(se)
    gr[2,]
    t2t[ gr$tx_id[[2]] ]

    se <- importAllelicCounts(coldata, a1="P", a2="M", format="assays")
    colData(se)
    assayNames(se)
    metadata(se)$alleles
  }
    
})
