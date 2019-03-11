swish.interx.pair <- function(infRepsArray, condition, covariate, pair,
                              nperms=30, wilcoxP, pc=5, quiet=FALSE) {
  stopifnot(is.numeric(pair) | is.character(pair) | is.factor(pair)) 
  pair <- as.integer(factor(pair))
  if (!all(table(pair, condition) == 1))
    stop("'pair' should have a single sample for both levels of condition")
  lfc.mat <- getLog2FCPair(infRepsArray, condition, pair, pc, mat=TRUE)
  NULL
}
