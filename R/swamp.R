#' swamp: Stratified Wilcoxon accounts for multiple predictors
#'
#' @param y a SummarizedExperiment containing the inferential replicate
#' matrices of median-ratio-scaled TPM as assays
#' @param x the name of the condition variable
#' @param cov the name of the covariate for adjustment
#' @param nperms the number of permutations of x
#'
#' @return a SummarizedExperiment with metadata columns added
#'
#' @references
#'
#' The \code{SAMseq} function in the \code{samr} package.
#' 
#' Jun Li and Robert Tibshirani "Finding consistent patterns:
#' A nonparametric approach for identifying differential expression
#' in RNA-Seq data" Stat Methods Med Res (2013).
#' 
#' @export
swamp <- function(y, x, cov, nperms=5) {
  if (is.null(metadata(y)$preprocessed) || !metadata(y)$preprocessed) {
    y <- preprocess(y)
    metadata(y)$preprocessed <- TRUE
  }
  ys <- y[mcols(y)$keep,]
  # rename 'y' to make it more clear
  resamp <- abind::abind(as.list(assays(ys)), along=3)
  # rename 'x' to make it more clear
  condition <- colData(y)[[x]]
  covariate <- colData(y)[[cov]]

  ngroups <- nlevels(covariate)
  #for (g in seq_len(ngroups)) {
  #...
  
  df <- data.frame(stat, qvalue)
  y <- postprocess(y, df)
  y
}

