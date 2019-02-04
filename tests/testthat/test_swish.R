context("swish")
library(SummarizedExperiment)
library(swish)

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
  
})
