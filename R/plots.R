#' Plot inferential replicates for a gene or transcript
#'
#' For datasets with inferential replicates, boxplots are
#' drawn for the two groups and potentially grouped by
#' covariates. For datasets with only mean and variance,
#' points and intervals (95% intervals using Normal
#' approximation) are drawn.
#' 
#' @param y a SummarizedExperiment (see \code{swish})
#' @param idx the name or row number of the gene or transcript
#' @param x the name of the condition variable
#' @param cov the name of the covariate for adjustment
#' @param colsDrk dark colors for the lines of the boxes
#' @param colsLgt light colors for the inside of the boxes
#' @param xaxis logical, whether to label the sample numbers.
#' default is \code{TRUE} if there are less than 30 samples
#' @param xlab the x-axis label
#' @param ylim y limits
#' @param main title
#' @param mainCol name of metadata column to use for title
#' (instead of rowname)
#' @param useMean logical, when inferential replicates
#' are not present, use the \code{mean} assay or the
#' \code{counts} assay for plotting
#' @param applySF logical, when inferential replicates are
#' not present, should \code{y$sizeFactor} be divided out
#' from the mean and interval plots (default FALSE)
#' @param reorder logical, should points within a group
#' defined by condition and covariate be re-ordered by
#' their count value (default is FALSE, except for alevin data)
#' @param thin integer, should the mean and interval lines
#' be drawn thin (the default switches from 0 [not thin]
#' to 1 [thinner] at n=150 cells, and from 1 [thinner]
#' to 2 [thinnest] at n=400 cells)
#' @param legend logical, show simple legend (default FALSE)
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
                        colsDrk=c("dodgerblue","goldenrod4","green4",
                                   "red3","purple4","royalblue4"),
                        colsLgt=c("lightblue1","goldenrod1","lightgreen",
                                   "salmon1","orchid1","royalblue1"),
                        xaxis, xlab, ylim,
                        main, mainCol,
                        legend=FALSE,
                        useMean=TRUE,
                        applySF=FALSE,
                        reorder,
                        thin) {
  
  hasInfReps <- any(grepl("infRep", assayNames(y)))
  if (!hasInfReps) {
    if (!("variance" %in% assayNames(y)))
      stop("if inferential replicates not present, requires 'variance' assay")
  }
  stopifnot(x %in% names(colData(y)))
  condition <- colData(y)[[x]]
  stopifnot(is(condition, "factor"))
  ncond <- nlevels(condition)
  stopifnot(length(colsDrk) == length(colsLgt))
  stopifnot(ncond <= length(colsDrk))
  colsDrk <- colsDrk[seq_len(ncond)]
  colsLgt <- colsLgt[seq_len(ncond)]
  if (missing(xaxis)) {
    xaxis <- ncol(y) < 30
  }
  if (missing(thin)) {
    thin <- if (ncol(y) >= 400) 2 else if (ncol(y) >= 150) 1 else 0
  } else {
    stopifnot(thin >= 0 & thin <= 2)
  }
  if (!is.null(cov)) {
    stopifnot(cov %in% names(colData(y)))
    covariate <- factor(colData(y)[[cov]])
    stopifnot(is(covariate, "factor"))
    ngrp <- nlevels(covariate)
  }
  infRepsScaled <- FALSE
  if (!is.null(metadata(y)$infRepsScaled)) {
    infRepsScaled <- metadata(y)$infRepsScaled
  }
  # single cell?
  sc <- FALSE 
  if (!is.null(metadata(y)$tximetaInfo$type)) {
    if (metadata(y)$tximetaInfo$type == "alevin") {
      sc <- TRUE
    }
  }
  if (missing(xlab)) {
    xlab <- if (sc) "cells" else "samples"
  }
  if (missing(reorder)) {
    reorder <- sc
  }
  ylab <- if (infRepsScaled) "scaled counts" else "counts"
  # this is a dummy variable used when making the plot()
  # if we don't put x-axis ticks, then we will move
  # the label up higher using title()
  xlabel <- if (xaxis) xlab else ""
  if (missing(main)) {
    if (missing(mainCol)) {
      if (is.character(idx)) {
        main <- idx
      } else {
        main <- rownames(y)[idx]
      }
    } else {
      stopifnot(mainCol %in% names(mcols(y)))
      main <- mcols(y)[idx,mainCol]
    }
  }
  if (reorder) {
    if (hasInfReps) {
      infReps <- assays(y[idx,])[grep("infRep",assayNames(y))]
      value <- colMeans(unlist(infReps))
    } else {
      if (useMean) stopifnot("mean" %in% assayNames(y))
      which.assay <- if (useMean) "mean" else "counts"
      value <- assays(y)[[which.assay]][idx,]
      if (applySF & !is.null(y$sizeFactor)) {
        value <- value/y$sizeFactor
      }
    }
    if (is.null(cov)) {
      o <- order(condition, value)
    } else {
      o <- order(covariate, condition, value)
    }    
  } else {
    if (is.null(cov)) {
      o <- order(condition)
    } else {
      o <- order(covariate, condition)
    }
  }
  
  if (is.null(cov)) {
    samp.nums <- unlist(lapply(table(condition), seq_len))
    col <- rep(colsDrk, table(condition))
    col.in <- rep(colsLgt, table(condition))
  } else {
    vec.tab <- as.vector(table(condition, covariate))
    samp.nums <- unlist(lapply(vec.tab, seq_len))
    col <- rep(rep(colsDrk, ngrp), vec.tab)
    col.in <- rep(rep(colsLgt, ngrp), vec.tab)
  }
  ### boxplot ###
  if (hasInfReps) {
    infReps <- assays(y[idx,])[grep("infRep",assayNames(y))]
    cts <- unlist(infReps)[,o]
    ymax <- max(cts)
    ymin <- if (is.null(cov)) 0 else -0.02 * ymax
    if (missing(ylim)) {
      ylim <- c(ymin,ymax)
    } else {
      stopifnot(length(ylim) == 2)
    }
    boxplot2(cts, col=col, col.in=col.in, ylim=ylim,
             xlab=xlabel, ylab=ylab, main=main)
    ### point and line plot ###
  } else {
    which.assay <- if (useMean) "mean" else "counts"
    cts <- assays(y)[[which.assay]][idx,o]
    sds <- sqrt(assays(y)[["variance"]][idx,o])
    Q <- qnorm(.975)
    if (applySF & !is.null(y$sizeFactor)) {
      cts <- cts / y$sizeFactor[o]
      sds <- sds / y$sizeFactor[o]
      ylab <- "scaled counts"
    }
    ymax <- max(cts + Q*sds)
    ymin <- if (is.null(cov)) 0 else -0.02 * ymax
    if (missing(ylim)) {
      ylim <- c(ymin, ymax)
    } else {
      stopifnot(length(ylim) == 2)
    }
    plot(cts, type="n", main=main,
         xaxt="n", ylim=ylim,
         xlab=xlabel, ylab=ylab)
    seg.lwd <- if (thin == 0) 2 else if (thin == 1) 1 else 3
    seg.col <- if (thin < 2) col else col.in
    segments(seq_along(cts), pmax(cts - Q*sds, 0),
             seq_along(cts), cts + Q*sds,
             col=seg.col, lwd=seg.lwd)
    pts.pch <- if (thin == 0) 22 else 15
    pts.lwd <- if (thin == 0) 1. else 1
    pts.cex <- if (thin == 0) 1 else 0.5
    points(cts, col=col, pch=pts.pch, bg=col.in,
           cex=pts.cex, lwd=pts.lwd)
  }
  if (xaxis) axis(1, seq_along(condition), samp.nums)
  if (!xaxis) {
    title(xlab=xlab, mgp=c(1,1,0))
  }
  if (!is.null(cov)) {
    cuts <- cumsum(table(covariate))
    segments(c(1,cuts[-ngrp]+1),ymin,cuts,ymin,lwd=3,
             col=rep(c("black","grey60"),length=ngrp))
  }
  if (legend) {
    legend("topleft", legend=levels(condition),
           col=colsDrk, pt.bg=colsLgt, pch=22,
           cex=.8, bg="white", box.col=NA, inset=.01)
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
