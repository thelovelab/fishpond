#' soggy: SAMseq on Gibbs-generated y
#'
#' @param y a SummarizedExperiment containing the inferential replicate
#' matrices of median-ratio-scaled TPM as assays
#' @param x the name of the condition variable
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
soggy <- function(y, x, nperms=5) {
  if (is.null(metadata(y)$preprocessed) || !metadata(y)$preprocessed) {
    y <- preprocess(y)
    metadata(y)$preprocessed <- TRUE
  }
  ys <- y[mcols(y)$keep,]
  # rename 'y' to make it more clear
  resamp <- abind::abind(as.list(assays(ys)), along=3)
  # rename 'x' to make it more clear
  condition <- colData(y)[[x]] 
  stat <- getSamStat(resamp, condition)
  perms <- samr:::getperms(condition, nperms)
  nulls <- matrix(nrow=nrow(ys), ncol=nperms)
  for (p in seq_len(nperms)) {
    cat(p, "")
    nulls[,p] <- getSamStat(resamp, condition[perms$perms[p,]])
  }
  nulls.vec <- as.vector(nulls)
  pi0 <- estimatePi0(stat, nulls.vec)
  samr.const.twoclass.unpaired.response <- "Two class unpaired"
  samr.obj <- list(
    resp.type=samr.const.twoclass.unpaired.response,
    tt=stat,
    ttstar0=nulls,
    foldchange.star=sign(stat),
    evo=sign(stat),
    pi0=pi0,
    assay.type="seq")
  delta.table <- samr:::samr.compute.delta.table(samr.obj)
  sig <- list(pup=which(stat >= 0), plo=which(stat < 0))
  qlist <- samr:::qvalue.func(samr.obj, sig, delta.table)
  qvalue <- numeric(length(stat))
  qvalue[stat >= 0] <- qlist$qvalue.up/100
  qvalue[stat < 0] <- qlist$qvalue.lo/100
  qvalue <- pmin(qvalue, 1)
  #par(mar=c(5,5,2,5))
  #hist(stat, breaks=100, col="grey", freq=FALSE)
  #d <- density(as.vector(nulls))
  #lines(d$x, pi0*d$y, col="blue", lwd=3)
  #lines(sort(stat), qvalues[order(stat)] * .03, col="red", lwd=3)
  #axis(4, c(0,.015,.03), c(0,.5,1))
  df <- data.frame(stat, qvalue)
  y <- postprocess(y, df)
  y
}

getSamStat <- function(resamp, condition) {
  dims <- dim(resamp)
  for (k in seq_len(dims[3])) {
    # modified from samr:::resample
    resamp[,,k] <- t(samr:::rankcol(t(resamp[,,k] + 0.1 * runif(dims[1]*dims[2]))))
  }
  fit <- samr:::wilcoxon.unpaired.seq.func(xresamp=resamp, y=condition)
  fit$tt
}

estimatePi0 <- function(stat, nulls.vec) {
  # modified from samr::samr
  qq <- quantile(nulls.vec, c(0.25, 0.75))
  sum(stat > qq[1] & stat < qq[2])/(0.5 * length(stat))
}
