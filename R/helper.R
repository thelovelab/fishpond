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
#' the scaledTPM so as to remove remaining library size differences
#' @param minCount for internal filtering, the minimum count 
#' @param minN for internal filtering, the minimum sample size
#' at \code{minCount}
#' @param quiet display no messages
#'
#' @return a SummarizedExperiment with the inferential replicates
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
  if (is.null(meanDepth)) {
    meanDepth <- exp(mean(log(colSums(counts))))
  }
  means <- matrix(nrow=nrow(y), ncol=nreps)
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
    # divide out the column sum, then set all to the meanDepth
    cts <- t(t(cts) / colSums(cts)) * meanDepth
    # filtering for calculting median ratio size factors
    use <- rowSums(infReps[[k]] >= minCount) >= minN
    if (is.null(sfFun)) {
      loggeomeans <- rowMeans(log(cts[use,]))
      sf <- apply(cts[use,], 2, function(s) {
        exp(median((log(s) - loggeomeans)[is.finite(loggeomeans)]))
      })
    } else {
      sf <- sfFun(cts)
    }
    infReps[[k]] <- t( t(cts)/sf )
    means[,k] <- rowMeans(infReps[[k]])
  }
  if (!quiet) message("")
  assays(y)[grep("infRep",assayNames(y))] <- infReps
  mcols(y)$log10mean <- log10(rowMeans(means) + 1)
  metadata(y)$infRepsScaled <- TRUE
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
  keep <- rowSums(assays(y)[["counts"]] >= minCount) >= minN
  mcols(y)$keep <- keep
  metadata(y)$preprocessed <- TRUE
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
#' @param null logical, whether to make an all null dataset
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
makeSimSwishData <- function(m=1000, n=10, numReps=20, null=FALSE) {
  stopifnot(m>8)
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
#' @param xaxis logical, whether to label the sample numbers.
#' default is \code{TRUE} if there are less than 30 samples
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
#' plotInfReps(y, 3, "condition", "batch")
#' 
#' @export
plotInfReps <- function(y, idx, x, cov=NULL,
                        cols.drk=c("dodgerblue","goldenrod4"),
                        cols.lgt=c("lightblue1","goldenrod1"),
                        xaxis) {
  infReps <- assays(y[idx,])[grep("infRep",assayNames(y))]
  stopifnot(x %in% names(colData(y)))
  condition <- colData(y)[[x]]
  if (missing(xaxis)) {
    xaxis <- ncol(y) < 30
  }
  if (is.null(cov)) {
    cts <- unlist(infReps)[,order(condition)]
    samp.nums <- unlist(lapply(table(condition), seq_len))
    col <- rep(cols.drk, table(condition))
    col.in <- rep(cols.lgt, table(condition))
  } else {
    stopifnot(cov %in% names(colData(y)))
    covariate <- factor(colData(y)[[cov]])
    ngrp <- nlevels(covariate)
    cts <- unlist(infReps)[,order(covariate, condition)]
    vec.tab <- as.vector(table(condition, covariate))
    samp.nums <- unlist(lapply(vec.tab, seq_len))
    col <- rep(rep(cols.drk, ngrp), vec.tab)
    col.in <- rep(rep(cols.lgt, ngrp), vec.tab)
  }
  main <- if (is.null(rownames(y))) {
            ""
          } else {
            if (is.character(idx)) idx else rownames(y)[idx]
          }
  ymax <- max(cts)
  ymin <- if (is.null(cov)) 0 else -0.02 * ymax
  boxplot2(cts, col=col, col.in=col.in, ylim=c(ymin,ymax),
           xlab="samples", ylab="scaled counts", main=main)
  if (xaxis) axis(1, seq_along(condition), samp.nums)
  if (!is.null(cov)) {
    cuts <- cumsum(table(covariate))
    segments(c(1,cuts[-ngrp]+1),ymin,cuts,ymin,lwd=3,
             col=rep(c("black","grey60"),length=ngrp))
  }
}

#' MA plot
#'
#' @param y a SummarizedExperiment (see \code{swish})
#' @param alpha the FDR threshold for coloring points
#' @param sigcolor the color for the significant points
#' @param ... passed to plot
#'
#' @return nothing, a plot is displayed
#' 
#' @examples
#'
#' y <- makeSimSwishData()
#' y <- scaleInfReps(y)
#' y <- labelKeep(y)
#' y <- swish(y, x="condition")
#' plotMASwish(y)
#' 
#' @export
plotMASwish <- function(y, alpha=.05, sigcolor="blue", ...) {
  dat <- data.frame(log10mean=mcols(y)$log10mean,
                    log2FC=mcols(y)$log2FC,
                    sig=mcols(y)$qvalue < alpha)
  dat <- dat[order(dat$sig),]
  with(dat, plot(log10mean, log2FC, pch=20, cex=.4,
                 col=ifelse(sig, sigcolor, "grey60"), ...))
  abline(h=0, col="grey40")
}

#' Compute inferential relative variance (InfRV)
#'
#' \code{InfRV} is used the Swish publication for visualization.
#' This function provides computation of the mean InfRV, a simple
#' statistic that measures inferential uncertainty.
#' Note that InfRV is not used in the \code{swish}
#' statistical method at all, it is just for visualization.
#' See function code for details.
#'
#' @param y a SummarizedExperiment
#' @param pc a pseudocount parameter for the denominator
#' @param shift a final shift parameter
#'
#' @return a SummarizedExperiment with \code{meanInfRV} in the metadata columns
#'
#' @export
computeInfRV <- function(y, pc=5, shift=.01) {
  infReps <- assays(y)[grep("infRep",assayNames(y))]
  infReps <- abind::abind(as.list(infReps), along=3)
  infVar <- apply(infReps, 1:2, var)
  mu <- assays(y)[["counts"]]
  # the InfRV computation:
  InfRV <- pmax(infVar - mu, 0)/(mu + pc) + shift
  mcols(y)$meanInfRV <- rowMeans(InfRV)
  y
}

boxplot2 <- function(x, w=.4, ylim, col, col.in, xlab="", ylab="", main="") {
  qs <- matrixStats::rowQuantiles(t(x), probs=0:4/4)
  if (missing(ylim)) {
    ylim <- c(min(x),max(x))
  }
  plot(qs[,3], type="n", xlim=c(0.5,ncol(x)+.5), xaxt="n",
       xlab=xlab, ylab=ylab, main=main, ylim=ylim)
  s <- seq_len(ncol(x))
  rect(s-w,qs[,2],s+w,qs[,4], col=col.in, border=col)
  segments(s-w, qs[,3], s+w, qs[,3], col=col, lwd=3, lend=1)
  segments(s, qs[,2], s, qs[,1], col=col, lty=2, lend=1)
  segments(s, qs[,4], s, qs[,5], col=col, lty=2, lend=1)
  segments(s-w/2, qs[,1], s+w/2, qs[,1], col=col)
  segments(s-w/2, qs[,5], s+w/2, qs[,5], col=col)
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
