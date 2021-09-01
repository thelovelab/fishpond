context("loadFry")

# Notes on how to generate test data sets for use here can be round in the
# R/test-helpers.R file

test_that("Reading in Alevin-fry USA count matrix works", {
  dat <- fishpond:::readExampleFryData("fry-usa-basic")
  expect_true(dir.exists(dat$parent_dir))

  # read default quantification (S + A)
  sce <- loadFry(dat$parent_dir, which_counts = c('S', 'A'))
  expect_equal(nrow(sce), length(dat$genes))
  expect_equal(ncol(sce), length(dat$barcodes))
  expect_equal(SummarizedExperiment::assayNames(sce), "counts")

  cts <- SummarizedExperiment::assay(sce, "counts")
  expect_s4_class(cts, "dgCMatrix")

  # Add the spliced and ambiguous reads manualy and convert into a matrix
  M <- local({
    m <- dat$matrix[, dat$usa$S, drop = FALSE] + dat$matrix[, dat$usa$A, drop = FALSE]
    dimnames(m) <- list(dat$barcodes, dat$genes)
    Matrix::Matrix(t(m), sparse = TRUE)
  })

  expect_equal(cts, M)
})

test_that("Main gene-level quantiation is same when velocity = TRUE or FALSE", {
  dat <- fishpond:::readExampleFryData("fry-usa-basic")

  # read in counts using default S,A counting
  sce <- loadFry(dat$parent_dir, which_counts = c("S", "A"))

  # return counts and unspliced estimates
  scev <- loadFry(dat$parent_dir, which_counts = c("S", "A"), velocity = TRUE)
  expect_equal(SummarizedExperiment::assayNames(scev), c("counts", "unspliced"))

  # gene-level quantitation should be the same w/ and w/o velocity
  expect_equal(assay(scev, "counts"), assay(sce, "counts"))

  # ensure spliced counts are same as manually assembled ones
  unspliced <- assay(scev, "unspliced")
  expect_s4_class(unspliced, "dgCMatrix")

  U <- local({
    u <- dat$matrix[, dat$usa$U, drop = FALSE]
    dimnames(u) <- list(dat$barcodes, dat$genes)
    Matrix::Matrix(t(u), sparse = TRUE)
  })

  expect_equal(unspliced, U)
})
