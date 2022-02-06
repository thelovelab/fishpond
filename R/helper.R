#' Scale inferential replicate counts
#'
#' A helper function to scale the inferential replicates
#' to the mean sequencing depth. The scaling takes into account
#' a robust estimator of size factor (median ratio method is used).
#' First, counts are corrected per row using the effective lengths
#' (for gene counts, the average transcript lengths), then scaled
#' per column to the geometric mean sequence depth, and finally are
#' adjusted per-column up or down by the median ratio size factor to
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
#' the scaledTPM so as to remove remaining library size differences.
#' Alternatively, one can provide a numeric vector of size factors
#' @param minCount for internal filtering, the minimum count 
#' @param minN for internal filtering, the minimum sample size
#' at \code{minCount}
#' @param saveMeanScaled store the mean of scaled inferential
#' replicates as an assay 'meanScaled'
#' @param quiet display no messages
#'
#' @return a SummarizedExperiment with the inferential replicates
#' as scaledTPM with library size already corrected (no need for further
#' normalization). A column \code{log10mean} is also added which is the
#' log10 of the mean of scaled counts across all samples and all inferential
#' replicates.
#'
#' @examples
#'
#' y <- makeSimSwishData()
#' y <- scaleInfReps(y)
#'
#' @import Rcpp
#' @export
scaleInfReps <- function(y, lengthCorrect=TRUE,
                         meanDepth=NULL, sfFun=NULL,
                         minCount=10, minN=3,
                         saveMeanScaled=FALSE,
                         quiet=FALSE) {
  if (!interactive()) {
    quiet <- TRUE
  }
  if (!is.null(metadata(y)$infRepsScaled)) {
    if (metadata(y)$infRepsScaled) stop("inferential replicates already scaled")
  }
  infRepIdx <- grep("infRep",assayNames(y))
  infRepError(infRepIdx)
  infReps <- assays(y)[infRepIdx]
  counts <- assays(y)[["counts"]]
  length <- assays(y)[["length"]]
  nreps <- length(infReps)
  if (is.null(meanDepth) & !is(sfFun,"numeric")) {
    meanDepth <- exp(mean(log(colSums(counts))))
  }
  means <- matrix(nrow=nrow(y), ncol=nreps)
  if (is.null(length)) {
    if (lengthCorrect) {
      if (!quiet) message("not correcting for feature length (lengthCorrect=FALSE)")
    }
    lengthCorrect <- FALSE
  }
  for (k in seq_len(nreps)) {
    if (!quiet) svMisc::progress(k, max.value=nreps, init=(k==1), gui=FALSE)
    if (lengthCorrect) {
      # new length bias correction matrix centered on 1
      length <- length / exp(rowMeans(log(length)))
      # a temporary matrix 'cts' which will store
      # the inferential replicate counts
      cts <- infReps[[k]] / length
    } else {
      # for 3' tagged scRNA-seq for example, don't length correct
      cts <- infReps[[k]]
    }
    # if size factors (numeric) were _not_ provided...
    if (!is(sfFun, "numeric")) {
      # divide out the column sum, then set all to the meanDepth
      cts <- t(t(cts) / colSums(cts)) * meanDepth
      # filtering for calculting median ratio size factors
      use <- rowSums(infReps[[k]] >= minCount) >= minN
      # calculate size factors
      if (is.null(sfFun)) {
        loggeomeans <- rowMeans(log(cts[use,]))
        sf <- apply(cts[use,], 2, function(s) {
          exp(median((log(s) - loggeomeans)[is.finite(loggeomeans)]))
        })
      } else if (is(sfFun, "function")) {
        sf <- sfFun(cts)
      }
      # ...otherwise we just divide counts by provided size factors
    } else {
      sf <- sfFun
    }
    infReps[[k]] <- t( t(cts)/sf )
    means[,k] <- rowMeans(infReps[[k]])
  }
  if (!quiet) message("")
  assays(y)[grep("infRep",assayNames(y))] <- infReps
  mcols(y)$log10mean <- log10(rowMeans(means) + 1)
  metadata(y)$infRepsScaled <- TRUE
  if (saveMeanScaled) {
    infRepsArray <- abind::abind(as.list(infReps), along=3)
    meanScaled <- apply(infRepsArray, 1:2, mean)
    assays(y)[["meanScaled"]] <- meanScaled
  }
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
#' @param x the name of the condition variable, will
#' use the smaller of the two groups to set \code{minN}.
#' Similar to edgeR's \code{filterByExpr}, as the smaller group
#' grows past 10, \code{minN} grows only by 0.7 increments
#' of sample size
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
labelKeep <- function(y, minCount=10, minN=3, x) {
  if (!missing(x)) {
    stopifnot(x %in% names(colData(y)))
    minN <- min(table(colData(y)[[x]]))
    # this modeled after edgeR::filterByExpr()
    if (minN > 10) {
      minN <- 10 + (minN - 10) * 0.7
    }
  }
  cts <- assays(y)[["counts"]]
  if (is(cts, "dgCMatrix")) {
    keep <- Matrix::rowSums(cts >= minCount) >= minN
  } else {
    keep <- rowSums(cts >= minCount) >= minN
  }
  mcols(y)$keep <- keep
  metadata(y)$preprocessed <- TRUE
  if (!"infRepScaled" %in% names(metadata(y))) {
    metadata(y)$infRepsScaled <- FALSE
  }
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
#' @param numReps how many inferential replicates to generate
#' @param null logical, whether to make an all null dataset
#' @param meanVariance logical, whether to output only mean and
#' variance of inferential replicates
#' @param allelic logical, whether to make an allelic sim dataset
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
makeSimSwishData <- function(m=1000, n=10, numReps=20,
                             null=FALSE, meanVariance=FALSE,
                             allelic=FALSE) {
  stopifnot(m > 8)
  if (allelic) {
    n <- 2 * n
  }
  stopifnot(n %% 2 == 0)
  cts <- matrix(rpois(m*n, lambda=80), ncol=n)
  if (!null) {
    grp2 <- (n/2+1):n
    cts[1:6,grp2] <- rpois(3*n, lambda=120)
    cts[7:8,] <- 0
  }
  length <- matrix(1000, nrow=m, ncol=n)
  abundance <- t(t(cts)/colSums(cts))*1e6
  infReps <- lapply(seq_len(numReps), function(i) {
    m <- matrix(rpois(m*n, lambda=80), ncol=n)
    if (!null) {
      # these row numbers are fixed for the demo dataset
      m[1:6,grp2] <- rpois(3*n, lambda=120)
      m[3:4,] <- round(m[3:4,] * runif(2*n,.5,1.5))
      m[5:6,grp2] <- round(pmax(m[5:6,grp2] + runif(n,-120,80),0))
      m[7:8,] <- 0
    }
    m
  })
  names(infReps) <- paste0("infRep", seq_len(numReps))
  if (meanVariance) {
    infRepsCube <- abind::abind(infReps, along=3)
    mu <- apply(infRepsCube, 1:2, mean)
    variance <- apply(infRepsCube, 1:2, var)
    assays <- list(counts=cts, abundance=abundance, length=length,
                   mean=mu, variance=variance)
  } else {
    assays <- list(counts=cts, abundance=abundance, length=length)
    assays <- c(assays, infReps)
  }
  se <- SummarizedExperiment(assays=assays)
  rownames(se) <- paste0("gene-",seq_len(nrow(se)))
  colnames(se) <- paste0("s",seq_len(n))
  metadata(se) <- list(countsFromAbundance="no")
  if (allelic) {
    als <- c("a2","a1")
    coldata <- DataFrame(allele=factor(rep(als, each=n/2), levels=als),
                         sample=factor(paste0("sample",rep(1:(n/2),2))),
                         row.names=colnames(se))
  } else {
    coldata <- DataFrame(condition=gl(2,n/2),
                         row.names=colnames(se))
  }
  colData(se) <- coldata
  se
}

#' Compute inferential relative variance (InfRV)
#'
#' \code{InfRV} is used the Swish publication for visualization.
#' This function provides computation of the mean InfRV, a simple
#' statistic that measures inferential uncertainty.
#' It also computes and adds the mean and variance of inferential
#' replicates, which can be useful ahead of \code{\link{plotInfReps}}.
#' Note that InfRV is not used in the \code{swish}
#' statistical method at all, it is just for visualization.
#' See function code for details.
#'
#' @param y a SummarizedExperiment
#' @param pc a pseudocount parameter for the denominator
#' @param shift a final shift parameter
#' @param meanVariance logical, use pre-computed inferential mean
#' and variance assays instead of \code{counts} and
#' computed variance from \code{infReps}. If missing,
#' will use pre-computed mean and variance when present
#' @param useCounts logical, whether to use the MLE
#' count matrix for the mean instead of mean of inferential replicates.
#' this argument is for backwards compatability, as previous
#' versions used counts. Default is FALSE
#'
#' @return a SummarizedExperiment with \code{meanInfRV} in the metadata columns
#'
#' @export
computeInfRV <- function(y, pc=5, shift=.01, meanVariance, useCounts=FALSE) {
  if (missing(meanVariance)) {
    meanVariance <- all(c("mean","variance") %in% assayNames(y))
  }
  if (meanVariance) {
    stopifnot(all(c("mean","variance") %in% assayNames(y)))
    infVar <- assays(y)[["variance"]]
    infMean <- assays(y)[["mean"]]
  } else {
    infReps <- assays(y)[grep("infRep",assayNames(y))]
    infReps <- abind::abind(as.list(infReps), along=3)
    infMean <- apply(infReps, 1:2, mean)
    infVar <- apply(infReps, 1:2, var)
    assays(y)[["mean"]] <- infMean
    assays(y)[["variance"]] <- infVar
  }
  if (useCounts) {
    infMean <- assays(y)[["counts"]]
  }
  # the InfRV computation:
  InfRV <- pmax(infVar - infMean, 0)/(infMean + pc) + shift
  mcols(y)$meanInfRV <- rowMeans(InfRV)
  y
}

#' Obtain a trace of inferential replicates for a sample
#'
#' Simple helper function to obtain a trace (e.g. MCMC trace)
#' of the ordered inferential replicates for one samples.
#' Supports either multiple features, \code{idx}, or multiple
#' samples, \code{samp_idx} (not both). Returns a tidy
#' data.frame for easy plotting.
#'
#' @param y a SummarizedExperiment with inferential replicates
#' as assays \code{infRep1} etc.
#' @param idx the names or row numbers
#' of the gene or transcript to plot
#' @param samp_idx the names or column numbers
#' of the samples to plot
#'
#' @return a data.frame with the counts along the interential
#' replicates, possible with additional columns specifying
#' feature or sample
#'
#' @examples
#'
#' y <- makeSimSwishData()
#' getTrace(y, "gene-1", "s1")
#' 
#' @export 
getTrace <- function(y, idx, samp_idx) {
  stopifnot(length(idx) == 1 | samp_idx == 1)
  stopifnot(is(idx, "character") | is(idx, "numeric"))
  stopifnot(is(samp_idx, "character") | is(samp_idx, "numeric"))
  infRepIdx <- grep("infRep",assayNames(y))
  nrep <- length(infRepIdx)
  if (length(idx) == 1 & length(samp_idx) == 1) {
    count <- sapply(infRepIdx, function(k) assay(y, i=k)[idx,samp_idx])
    data.frame(infRep=seq_along(infRepIdx), count)
  } else if (length(idx) == 1) {
    out <- lapply(samp_idx, function(j) {
      sapply(infRepIdx, function(k) assay(y, i=k)[idx,j])
    })
    data.frame(infRep=seq_along(infRepIdx), 
               count=do.call(c, out),
               sample=rep(samp_idx, each=nrep))
  } else if (length(samp_idx) == 1) {
    out <- lapply(idx, function(i) {
      sapply(infRepIdx, function(k) assay(y, i=k)[i,samp_idx])
    })
    data.frame(infRep=seq_along(infRepIdx),
               count=do.call(c, out),
               feature=rep(idx, each=nrep))
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
    stop("there are no inferential replicates in the assays of 'y';
see Quick Start in the swish vignette for details")
  }
}
