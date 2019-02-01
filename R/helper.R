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
#' @param y a SummarizedExperiment with: \code{infReps} a list of
#' inferential replicate count matrices, \code{counts} the
#' estimated counts matrix, and \code{length} the effective
#' lengths matrix
#' @param lengthCorrect whether to use effective length correction
#' (default is TRUE)
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
scaleInfReps <- function(y, lengthCorrect=TRUE, meanDepth=NULL, sfFun=NULL, minCount=10, minN=3) {
  infReps <- assays(y)[grep("infRep",assayNames(y))]
  counts <- assays(y)[["counts"]]
  length <- assays(y)[["length"]]
  nreps <- length(infReps)
  if (is.null(meanDepth)) {
    meanDepth <- exp(mean(log(colSums(counts))))
  }
  for (k in seq_len(nreps)) {
    cat(k,"")
    # we don't technically create TPM, but something proportion to
    if (lengthCorrect) {
      tpm <- infReps[[k]] / length
    } else {
      # for 3' tagged scRNA-seq for example, don't length correct
      tpm <- infReps[[k]]
    }
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
  assays(y) <- infReps
  y
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

#' Make simulated data for swish for examples/testing
#'
#' Makes a small swish dataset for examples and testing.
#'
#' @param m number of genes/txps
#' @param n number of samples
#' @param numReps how many inferential replicates
#'
#' @return a SummarizedExperiment
#'
#' @export
makeSwishData <- function(m=400, n=10, numReps=20) {
  stopifnot(m>8)
  stopifnot(n %% 2 == 0)
  cts <- matrix(rpois(m*n, lambda=100), ncol=n)
  cts[1:4,1:(n/2)] <- rpois(2*n, lambda=200)
  cts[5:6,1:(n/2)] <- rpois(n, lambda=150)
  cts[7:8,] <- 0
  length <- matrix(1000, nrow=m, ncol=n)
  abundance <- t(t(cts)/colSums(cts))*1e6
  infReps <- lapply(seq_len(numReps), function(i) {
    m <- matrix(rpois(m*n, lambda=100), ncol=n)
    m[1:4,1:(n/2)] <- rpois(2*n, lambda=200)
    m[5:6,1:(n/2)] <- rpois(n, lambda=150)
    m[3:4,] <- round(m[3:4,] * sample(c(.5,1.5),2*n,TRUE))
    m[5:6,1:(n/2)] <- round(m[5:6,1:(n/2)] * sample(c(.33,1.66),n,TRUE))
    m[7:8,] <- 0
    m
  })
  names(infReps) <- paste0("infRep", seq_len(numReps))
  assays <- list(counts=cts, abundance=abundance, length=length)
  assays <- c(assays, infReps)
  se <- SummarizedExperiment::SummarizedExperiment(assays=assays)
  metadata(se) <- list(countsFromAbundance="no")
  colData(se) <- DataFrame(condition=gl(2,n/2))
  se
}

postprocess <- function(y, df) {
  for (stat in names(df)) {
    mcols(y)[[stat]] <- numeric(nrow(y))
    mcols(y)[[stat]][mcols(y)$keep] <- df[[stat]]
    mcols(y)[[stat]][!mcols(y)$keep] <- NA
  }
  y
}
