context("readEDS")
library(fishpond)
test_that("Reading in Alevin EDS format works", {

  dir <- system.file("extdata", package="tximportData")
  files <- file.path(dir,"alevin/neurons_900_v014/alevin/quants_mat.gz")
  file.exists(files)
  dir <- sub("/alevin$","",dirname(files))  
  barcode.file <- file.path(dir, "alevin/quants_mat_rows.txt")
  gene.file <- file.path(dir, "alevin/quants_mat_cols.txt")
  cell.names <- readLines(barcode.file)
  gene.names <- readLines(gene.file)
  num.cells <- length(cell.names)
  num.genes <- length(gene.names)
  out <- readEDS(files, num.genes, num.cells)
  
})
