context("loadFry")
library(Matrix)
library(SingleCellExperiment)

test_that("Reading in Alevin-fry USA count matrix works", {
  
  dir <- system.file("extdata", "alevin", "test_loadFry", package="fishpond")
  fry.dir <- file.path(dir)
  file.exists(fry.dir)
  expect_equal(file.exists(fry.dir), TRUE)
  
  # reading in quants
  sce <- loadFry(fry.dir)

  expect_equal(nrow(sce), 2)
  expect_equal(ncol(sce), 3)
  cts <- counts(sce)
  m = Matrix(nrow = 2, ncol = 3, data = c(1,2,2,2,2,1), sparse = TRUE,
             dimnames = list(c("gene1", "gene2"),
                             c("bc1", "bc2", "bc3"))
             )
  m = as(m, "dgCMatrix")
  expect_equal(cts , m)
})


##############################################
# following are the code to create the test data
# library(Matrix)
# library(jsonlite)
# library(tibble)
# # create mtx
# m <- Matrix(nrow = 3, ncol = 6, data = 1, sparse = TRUE)
# m <- as(m, "dgCMatrix") # by default, Matrix() returns dgCMatrix
# m[1,1] <- 0
# m[2,3] <- 0
# m[3,6] <- 0
# m
# 
# Assuming the current work dir is the place where `test_loadFry.R` locates
# i.e., `fishpond/tests/testthat`
# fry.dir = file.path("../../inst/extdata/alevin/test_loadFry")
# dir.create(file.path(fry.dir, "alevin"),recursive = TRUE,showWarnings = FALSE)
# writeMM(m,file.path(fry.dir, "alevin", "quants_mat.mtx"))
# m = readMM(file.path(fry.dir, "alevin", "quants_mat.mtx"))
# # create feature names
# 
# write.table(c("gene1", "gene2"), file =file.path(fry.dir, "alevin", "quants_mat_cols.txt"),quote = FALSE,col.names = FALSE, row.names = FALSE, sep = "\t")
# 
# # create cellbarcodes
# 
# write.table(c("bc1", "bc2", "bc3"), file =file.path(fry.dir, "alevin", "quants_mat_rows.txt"),quote = FALSE,col.names = FALSE, row.names = FALSE, sep = "\t")
# 
# meta_info = list()
# meta_info[["alt_resolved_cell_numbers"]] = list()
# meta_info[["cmd"]] = ""
# meta_info[["dump_eq"]] = FALSE
# meta_info[["num_genes"]] = 6
# meta_info[["num_quantified_cells"]] = 3
# meta_info[["resolution_strategy"]] = "CellRangerLike"
# meta_info[["usa_mode"]] = TRUE
# 
# write(toJSON(meta_info, pretty=TRUE), file=file.path(fry.dir, "meta_info.json"))
# 
# fromJSON(file.path(fry.dir, "meta_info.json"))



