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
#' @param minN for internal filtering, the minimum sample size
#' at \code{minCount}
#' @param quiet display no messages
#'
#' @return a SummarizedExperiment wiht the inferential replicates
#' as scaledTPM with library size already corrected (no need for further
#' normalization)
#'
#' @examples
#'
#' y <- makeSimSwishData()
#' y <- scaleInfReps(y)
#' 
#' @export
scaleInfReps <- function(y, lengthCorrect=TRUE,
                         meanDepth=NULL, sfFun=NULL,
                         minCount=10, minN=3, quiet=FALSE) {
  infRepIdx <- grep("infRep",assayNames(y))
  infRepError(infRepIdx)
  infReps <- assays(y)[infRepIdx]
  counts <- assays(y)[["counts"]]
  length <- assays(y)[["length"]]
  nreps <- length(infReps)
  if (is.null(meanDepth)) {
    meanDepth <- exp(mean(log(colSums(counts))))
  }
  for (k in seq_len(nreps)) {
    if (!quiet) progress(k, max.value=nreps, init=(k==1))
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
  if (!quiet) message("")
  assays(y)[grep("infRep",assayNames(y))] <- infReps
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
#' @param quiet display no messages
#'
#' @return a SummarizedExperiment with a new column \code{keep}
#' in \code{mcols(y)}
#'
#' @examples
#' 
#' y <- makeSimSwishData()
#' y <- scaleInfReps(y)
#' y <- labelKeep(y)
#' 
#' @export
labelKeep <- function(y, minCount=10, minN=3, quiet=FALSE) {
  infRepIdx <- grep("infRep",assayNames(y))
  infRepError(infRepIdx)
  infReps <- assays(y)[infRepIdx]
  nreps <- length(infReps)
  keep.mat <- matrix(nrow=nrow(y), ncol=nreps)
  for (k in seq_len(nreps)) {
    if (!quiet) progress(k, max.value=nreps, init=(k==1))
    keep.mat[,k] <- rowSums(infReps[[k]] >= minCount) >= minN
  }
  keep <- apply(keep.mat, 1, all)
  mcols(y)$keep <- keep
  metadata(y)$preprocessed <- TRUE
  if (!quiet) message("")
  y
}

#' Make simulated data for swish for examples/testing
#'
#' Makes a small swish dataset for examples and testing.
#' The first six genes have some differential expression
#' evidence in the counts, with varying degree of inferential
#' variance across inferential replicates (1-2: minor,
#' 3-4: some, 5-6: substantial). The 7th and 8th
#' genes have all zeros to demonstrate \code{labelKeep}.
#' 
#' @param m number of genes
#' @param n number of samples
#' @param numReps how many inferential replicates
#'
#' @return a SummarizedExperiment
#'
#' @examples
#'
#' library(SummarizedExperiment)
#' y <- makeSimSwishData()
#' assayNames(y)
#' 
#' @export
makeSimSwishData <- function(m=1000, n=10, numReps=20) {
  stopifnot(m>8)
  stopifnot(n %% 2 == 0)
  cts <- matrix(rpois(m*n, lambda=80), ncol=n)
  grp2 <- (n/2+1):n
  cts[1:6,grp2] <- rpois(3*n, lambda=120)
  cts[7:8,] <- 0
  length <- matrix(1000, nrow=m, ncol=n)
  abundance <- t(t(cts)/colSums(cts))*1e6
  infReps <- lapply(seq_len(numReps), function(i) {
    m <- matrix(rpois(m*n, lambda=80), ncol=n)
    # these row numbers are fixed for the demo dataset
    m[1:6,grp2] <- rpois(3*n, lambda=120)
    m[3:4,] <- round(m[3:4,] * runif(2*n,.5,1.5))
    m[5:6,grp2] <- round(pmax(m[5:6,grp2] + runif(n,-120,80),0))
    m[7:8,] <- 0
    m
  })
  names(infReps) <- paste0("infRep", seq_len(numReps))
  assays <- list(counts=cts, abundance=abundance, length=length)
  assays <- c(assays, infReps)
  se <- SummarizedExperiment(assays=assays)
  rownames(se) <- paste0("gene-",seq_len(nrow(se)))
  metadata(se) <- list(countsFromAbundance="no")
  colData(se) <- DataFrame(condition=gl(2,n/2))
  se
}

#' Plot inferential replicates for a gene or transcript
#'
#' @param y a SummarizedExperiment (see \code{swish})
#' @param idx the name or row number of the gene or transcript
#' @param x the name of the condition variable
#' @param cov the name of the covariate for adjustment
#' @param cols.drk dark colors for the lines of the boxes
#' @param cols.lgt light colors for the inside of the boxes
#' @param xaxis logical, whether to label the sample numbers
#'
#' @return nothing, a plot is displayed
#' 
#' @examples
#'
#' y <- makeSimSwishData()
#' plotInfReps(y, 3, "condition")
#'
#' y <- makeSimSwishData(n=40)
#' y$batch <- factor(rep(c(1,2,3,1,2,3),c(5,10,5,5,10,5)))
#' plotInfReps(y, 3, "condition", "batch", xaxis=FALSE)
#' 
#' @export
plotInfReps <- function(y, idx, x, cov=NULL,
                        cols.drk=c("dodgerblue","goldenrod4"),
                        cols.lgt=c("lightblue1","goldenrod1"),
                        xaxis=TRUE) {
  infReps <- assays(y[idx,])[grep("infRep",assayNames(y))]
  condition <- colData(y)[[x]]
  if (is.null(cov)) {
    cts <- unlist(infReps)[,order(condition)]
    samp.nums <- as.vector(sapply(table(condition), seq_len))
    cols <- rep(cols.drk, table(condition))
    cols.in <- rep(cols.lgt, table(condition))
  } else {
    covariate <- factor(colData(y)[[cov]])
    ngrp <- nlevels(covariate)
    cts <- unlist(infReps)[,order(covariate, condition)]
    vec.tab <- as.vector(table(condition, covariate))
    samp.nums <- unlist(lapply(vec.tab, seq_len))
    cols <- rep(rep(cols.drk, ngrp), vec.tab)
    cols.in <- rep(rep(cols.lgt, ngrp), vec.tab)
  }
  main <- if (is.null(rownames(y))) {
            ""
          } else {
            if (is.character(idx)) idx else rownames(y)[idx]
          }
  boxplot(cts,range=0,border=cols,col=cols.in,xaxt="n",
          ylim=c(0,max(cts)),
          xlab="samples",ylab="scaled counts",
          main=main)
  if (xaxis) axis(1, seq_along(condition), samp.nums)
  if (!is.null(cov)) {
    cuts <- cumsum(table(covariate))
    segments(c(1,cuts[-ngrp]+1),0,cuts,0,lwd=3,
             col=rep(c("black","grey"),length=ngrp))
  }
}

postprocess <- function(y, df) {
  for (stat in names(df)) {
    mcols(y)[[stat]] <- numeric(nrow(y))
    mcols(y)[[stat]][mcols(y)$keep] <- df[[stat]]
    mcols(y)[[stat]][!mcols(y)$keep] <- NA
  }
  y
}

infRepError <- function(infRepIdx) {
  if (length(infRepIdx) == 0) {
    stop("there are no inferential replicates in the assays of 'y'")
  }
}
