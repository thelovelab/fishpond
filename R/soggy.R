#' soggy: SAMseq on Gibbs-generated y
#'
#' @param y a SummarizedExperiment containing the inferential replicate
#' matrices of median-ratio-scaled TPM as assays
#' @param x the condition
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
  y <- preprocess(y)
  ys <- y[mcols(y)$keep,]
  # rename 'y' to make it more clear
  resamp <- abind::abind(as.list(assays(ys)), along=3)
  # rename 'x' to make it more clear
  condition <- x 
  stat <- getSamStat(resamp, condition)
  perms <- samr:::getperms(condition, nperms)
  nulls <- matrix(nrow=nrow(ys), ncol=nperms)
  for (p in seq_len(nperms)) {
    cat(p, "")
    nulls[,p] <- getSamStat(resamp, condition[perms$perms[p,]])
  }
  ## not done, still need to incorporate perms as local FDR as in SAMseq
  #hist(as.vector(nulls), breaks=100, col="grey")
  #hist(stat, breaks=100, col="grey", freq=FALSE)
  #lines(density(as.vector(nulls)), col="blue", lwd=3)
  y <- postprocess(y, stat)
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
