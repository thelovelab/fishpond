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
  }
  ys <- y[mcols(y)$keep,]
  # rename 'y' to make it more clear
  infRepsArray <- abind::abind(as.list(assays(ys)), along=3)
  # rename 'x' to make it more clear
  condition <- colData(y)[[x]]
  covariate <- colData(y)[[cov]]

  ngroups <- nlevels(covariate)
  groups <- levels(covariate)
  stats <- matrix(nrow=nrow(ys), ncol=ngroups)
  nulls.big <- array(dim=list(nrow(ys), nperms, ngroups))
  for (i in seq_len(ngroups)) {
    g <- groups[i]
    infRepsArray.sub <- infRepsArray[,covariate == g,]
    cond.sub <- condition[covariate == g]
    stats[,i] <- getSamStat(infRepsArray.sub, cond.sub)
    perms <- samr:::getperms(cond.sub, nperms)
    for (p in seq_len(nperms)) {
      cat(p, "")
      nulls.big[,p,i] <- getSamStat(infRepsArray.sub, cond.sub[perms$perms[p,]])
    }
  }
  ns <- unname(table(covariate))
  wts <- 1/(ns + 1)
  stat <- as.vector(stats %*% wts)
  nulls <- matrix(nrow=nrow(ys), ncol=nperms)
  for (p in seq_len(nperms)) {
    nulls[,p] <- nulls.big[,p,] %*% wts
  }
  nulls.vec <- as.vector(nulls)
  pi0 <- estimatePi0(stat, nulls.vec)
  locfdr <- makeLocFDR(stat, nulls, pi0)

  ## par(mar=c(5,5,2,5))
  ## hist(stat, breaks=100, col="grey", freq=FALSE)
  ## d <- density(nulls.vec)
  ## lines(d$x, pi0*d$y, col="blue", lwd=3)
  ## lines(sort(stat), qvalue[order(stat)] * 1, col="red", lwd=3)
  ## axis(4, c(0,.5,1), c(0,.5,1))
  
  df <- data.frame(stat, locfdr)
  y <- postprocess(y, df)
  y
}
