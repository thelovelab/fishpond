#' Scale inferential replicate counts
#'
#' A helper function to scale the inferential replicates
#' to the mean sequencing depth. The scaling takes into account
#' a robust estimator of size factor (median ratio method is used).
#' First, counts are converted to TPM, then scaled to
#' the mean sequence depth, finally the scaledTPM are
#' adjusted up or down by the median ratio size factor to
#' minimize systematic differences across samples.
#'
#' @param se a SummarizedExperiment with: \code{infReps} a list of
#' inferential replicate count matrices, \code{counts} the
#' estimated counts matrix, and \code{length} the effective
#' lengths matrix
#' @param meanDepth (optional) user can
#' specify a different mean sequencing depth. By default
#' the geometric mean sequencing depth is computed
#' @param sfFun (optional) size factors function. An
#' alternative to the median ratio can be provided here to adjust
#' the scaledTPM so as to remove remaining library size differences
#' @param minCount for internal filtering, the minimum count 
#' @param minN for internal filtering, the minimum sample size at \code{minCount}
#'
#' @return a SummarizedExperiment wiht the inferential replicates
#' as scaledTPM with library size already corrected (no need for further
#' normalization)
#'
#' @export
scaleInfRep <- function(se, meanDepth=NULL, sfFun=NULL, minCount=10, minN=3) {
  infReps <- assays(se)[grep("infRep",assayNames(se))]
  counts <- assays(se)[["counts"]]
  length <- assays(se)[["length"]]
  nreps <- length(infReps)
  if (is.null(meanDepth)) {
    meanDepth <- exp(mean(log(colSums(counts))))
  }
  for (k in seq_len(nreps)) {
    cat(k,"")
    # we don't technically create TPM, but something proportion to
    tpm <- infReps[[k]] / length
    # divide out the column sum, then set all to the meanDepth
    tpm <- t(t(tpm) / colSums(tpm)) * meanDepth
    # filtering for calculting median ratio size factors
    use <- rowSums(infReps[[k]] >= minCount) >= minN
    if (is.null(sfFun)) {
      loggeomeans <- rowMeans(log(tpm[use,]))
      sf <- apply(tpm[use,], 2, function(s) {
        exp(median((log(s) - loggeomeans)[is.finite(loggeomeans)]))
      })
    } else {
      sf <- sfFun(tpm)
    }
    infReps[[k]] <- t( t(tpm)/sf )
  }
  cat("\n")
  SummarizedExperiment(assays=infReps, colData=colData(se), rowRanges=rowRanges(se))
}

#' Label rows to keep based on minimal count
#'
#' Adds a column \code{keep} to \code{mcols(y)} that specifies
#' which rows of the SummarizedExperiment will be included
#' in statistical testing. Rows are not removed, just marked
#' with the logical \code{keep}.
#'
#' @param y a SummarizedExperiment
#' @param minCount the minimum count
#' @param minN the minimum sample size at \code{minCount}
#'
#' @return a SummarizedExperiment with a new column \code{keep}
#' in \code{mcols(y)}
#'
#' @export
labelKeep <- function(y, minCount=10, minN=3) {
  rep.idx <- seq_along(assayNames(y))
  nreps <- length(assayNames(y))
  keep.mat <- matrix(nrow=nrow(y), ncol=nreps)
  for (k in seq_len(nreps)) {
    cat(k,"")
    keep.mat[,k] <- rowSums(assays(y)[[k]] >= minCount) >= minN
  }
  keep <- apply(keep.mat, 1, all)
  mcols(y)$keep <- keep
  metadata(y)$preprocessed <- TRUE
  cat("\n")
  y
}

postprocess <- function(y, df) {
  for (stat in names(df)) {
    mcols(y)[[stat]] <- numeric(nrow(y))
    mcols(y)[[stat]][mcols(y)$keep] <- df[[stat]]
    mcols(y)[[stat]][!mcols(y)$keep] <- NA
  }
  y
}
