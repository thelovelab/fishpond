context("loadFry")

# Notes on how to generate test data sets for use here can be round in the
# R/test-helpers.R file

test_that("Reading in Alevin-fry USA count matrix works", {
  dat <- fishpond:::readExampleFryData("fry-usa-basic")
  expect_true(dir.exists(dat$parent_dir))

  # read default quantification (S + A)
  sce <- loadFry(dat$parent_dir)
  expect_equal(nrow(sce), length(dat$genes))
  expect_equal(ncol(sce), length(dat$barcodes))
  expect_equal(sort(SummarizedExperiment::assayNames(sce)), sort(c("counts", "unspliced")))

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
