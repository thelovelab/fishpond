context("swish")
library(SummarizedExperiment)
library(fishpond)

test_that("basic variable errors thrown", {

  y <- makeSimSwishData()
  y <- scaleInfReps(y)
  y <- labelKeep(y)

  # too many levels of condition
  y2 <- y
  y2$condition <- gl(3,3,10)
  expect_error(swish(y2, "condition"))

  # batch and pair together
  y2 <- y
  y2$batch <- rep(1:2,5)
  y2$pair <- rep(1:5,2)
  expect_error(swish(y2, "condition", "batch", "pair"))

  # wrong number of pairs
  y2 <- y
  y2$pair <- c(1,2,3,4,5,1,1,2,2,3)
  expect_error(swish(y2, "condition", pair="pair"), "single sample for both levels")
  
  # no inferential replicates
  y <- makeSimSwishData()
  assays(y) <- assays(y)[c("counts","abundance","length")]
  expect_error(scaleInfReps(y), "no inferential")
  expect_error(labelKeep(y), "no inferential")
  expect_error(swish(y, "condition"), "no inferential")
  
})
