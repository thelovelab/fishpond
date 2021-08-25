context("loadFry")
library(Matrix)
library(SingleCellExperiment)

test_that("Reading in Alevin-fry USA count matrix works", {
  dir <- system.file("extdata", "alevin", "test_loadFry", package="fishpond")
  fry.dir <- file.path(dir)
  file.exists(fry.dir)
  expect_equal(file.exists(fry.dir), TRUE)

  # reading in quants with no velocity data
  sce <- loadFry(fry.dir)
  expect_equal(nrow(sce), 2)
  expect_equal(ncol(sce), 3)
  expect_equal(assayNames(sce), c("counts"))

  cts <- counts(sce)

  # load the matrix manually so we don't hard code expected values
  quants <- .loadFry(fry.dir, ncol(sce))
  m <- t(mraw[, idx$S, drop = FALSE]  + mraw[, idx$A, drop = FALSE])
  ng <- nrow(sce)
  idx <- list(
    S = seq(1, ng),
    U = seq(ng + 1, 2 * ng),
    A = seq(2 * ng + 1, 3 * ng))
  mraw <- Matrix::readMM(file.path(fry.dir, "alevin", "quants_mat.mtx"))

  dimnames(m) <- list(
    readLines(file.path(fry.dir, "alevin", "quants_mat_cols.txt")),
    readLines(file.path(fry.dir, "alevin", "quants_mat_rows.txt")))

  expect_equal(cts, m)
})

test_that("Main gene-level quantiation is same when velocity = TRUE or FALSE", {
  fdir <- system.file("extdata", "alevin", "test_loadFry", package = "fishpond")

  # read in counts using default S,A counting
  sce <- loadFry(fdir, which_counts = c("S", "A"))

  # return counts and unspliced estimates
  scev <- loadFry(fdir, which_counts = c("S", "A"), velocity = TRUE)
  expect_equal(assayNames(scev), c("counts", "unspliced"))

  expect_equal(assay(scev, "counts"), assay(sce, "counts"))
  # gene-level quantitation should be the same w/ and w/o velocity
})

# Test data creation receipt ---------------------------------------------------
# Create a test data set that implicitly does not get converted to a logical
# (dgcTMatrix) because it's all 0s and 1s. This will allow us to properly
# test the `velocity=TRUE` parameter.
if (FALSE) {
  # Assume the current work dir is root of Rpkg
  stopifnot(
    basename(getwd()) == "fishpond",
    dir.exists(file.path(getwd(), "inst/extdata/alevin")))

  library(Matrix)
  library(jsonlite)
  # create mtx
  m <- Matrix(nrow = 3, ncol = 6, data = 2, sparse = TRUE)
  m <- as(m, "dgCMatrix") # by default, Matrix() returns dgCMatrix
  m[1,1] <- 0
  m[2,3] <- 0
  m[2,5] <- 0
  m[3,6] <- 0
  m

  fry.dir = file.path("inst/extdata/alevin/test_loadFry")
  dir.create(file.path(fry.dir, "alevin"), recursive = TRUE, showWarnings = FALSE)
  writeMM(m,file.path(fry.dir, "alevin", "quants_mat.mtx"))
  m = readMM(file.path(fry.dir, "alevin", "quants_mat.mtx"))

  # create feature names
  write.table(c("gene1", "gene2"), file =file.path(fry.dir, "alevin", "quants_mat_cols.txt"),quote = FALSE,col.names = FALSE, row.names = FALSE, sep = "\t")

  # create cellbarcodes
  write.table(c("bc1", "bc2", "bc3"), file =file.path(fry.dir, "alevin", "quants_mat_rows.txt"),quote = FALSE,col.names = FALSE, row.names = FALSE, sep = "\t")

  # create metadata
  meta_info = list()
  meta_info[["alt_resolved_cell_numbers"]] = list()
  meta_info[["cmd"]] = ""
  meta_info[["dump_eq"]] = FALSE
  meta_info[["num_genes"]] = 6
  meta_info[["num_quantified_cells"]] = 3
  meta_info[["resolution_strategy"]] = "CellRangerLike"
  meta_info[["usa_mode"]] = TRUE

  write(toJSON(meta_info, pretty=TRUE), file=file.path(fry.dir, "meta_info.json"))
}

