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

# The paths of the test directories in inst/extdata/alevin/example-quants
fryExample <- function(...) {
  system.file("extdata", "alevin", "example-quants", ..., package = "fishpond",
              mustWork = TRUE)
}

# The base R content of a SingleCellExperiment, for a comparison with saved values
snapSCE <- function(sce) {
  list(assays = lapply(SummarizedExperiment::assays(sce), as.matrix),
       rowData = as.data.frame(SummarizedExperiment::rowData(sce)),
       colData = as.data.frame(SummarizedExperiment::colData(sce)),
       dimnames = dimnames(sce), metadata = S4Vectors::metadata(sce))
}

test_that("loadFry gives the same output as before for directories that are not multiplex", {
  # loadFry_regression.rds was made with the loadFry of fishpond 2.19.2
  # (commit 111f6c0), with the same directories and arguments as below.
  expected <- readRDS(test_path("loadFry_regression.rds"))
  dirs <- c(fry_usa_basic = fryExample("fry-usa-basic"),
            legacy_usa = fryExample("flex", "legacy_usa"))
  cases <- list(
    scRNA = list(outputFormat = "scRNA"), snRNA = list(outputFormat = "snRNA"),
    all = list(outputFormat = "all"), `S+A` = list(outputFormat = "S+A"),
    `U+S+A` = list(outputFormat = "U+S+A"), velocity = list(outputFormat = "velocity"),
    scVelo = list(outputFormat = "scVelo"), raw = list(outputFormat = "raw"),
    custom = list(outputFormat = list(spliced = c("S", "A"), unspliced = "U")),
    nonzero = list(outputFormat = "scRNA", nonzero = TRUE))
  for (d in names(dirs)) {
    for (k in names(cases)) {
      sce <- do.call(loadFry, c(list(fryDir = dirs[[d]], quiet = TRUE), cases[[k]]))
      expect_identical(snapSCE(sce), expected[[d]][[k]], info = paste(d, k))
    }
  }
  # a single-sample run: featureDump.txt has no sample_name column
  single <- fryExample("flex", "flex_single")
  expect_identical(snapSCE(loadFry(single, quiet = TRUE)), expected$flex_single$default)
  expect_identical(snapSCE(loadFry(single, nonzero = TRUE, quiet = TRUE)),
                   expected$flex_single$nonzero)
})
