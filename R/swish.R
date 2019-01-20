#' swish: SAMseq With Inferential Samples Helps
#'
#' @param y a SummarizedExperiment containing the inferential replicate
#' matrices of median-ratio-scaled TPM as assays
#' @param x the name of the condition variable
#' @param cov the name of the covariate for adjustment.
#' if provided a stratified Wilcoxon in performed
#' @param nperms the number of permutations
#' @param estPi0 logical, whether to estimate pi0
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
swish <- function(y, x, cov=NULL, nperms=30, estPi0=FALSE) {
  if (is.null(metadata(y)$preprocessed) || !metadata(y)$preprocessed) {
    y <- preprocess(y)
  }
  ys <- y[mcols(y)$keep,]
  # rename 'y' to make it more clear
  infRepsArray <- abind::abind(as.list(assays(ys)), along=3)
  # rename 'x' to make it more clear
  condition <- colData(y)[[x]] 

  if (is.null(cov)) {
    stat <- getSamStat(infRepsArray, condition)
    perms <- samr:::getperms(condition, nperms)
    nulls <- matrix(nrow=nrow(ys), ncol=nperms)
    for (p in seq_len(nperms)) {
      cat(p, "")
      nulls[,p] <- getSamStat(infRepsArray, condition[perms$perms[p,]])
    }
  } else {
    covariate <- colData(y)[[cov]]
    out <- swish.strat(infRepsArray, condition, covariate, nperms=nperms)
    stat <- out$stat
    nulls <- out$nulls
  }
  nulls.vec <- as.vector(nulls)
  if (estPi0) {
    pi0 <- estimatePi0(stat, nulls.vec)
  } else {
    pi0 <- 1
  }
  locfdr <- makeLocFDR(stat, nulls, pi0)
  qvalue <- makeQvalue(stat, nulls, pi0)
  df <- data.frame(stat, locfdr, qvalue)
  y <- postprocess(y, df)
  y
}

getSamStat <- function(infRepsArray, condition) {
  dims <- dim(infRepsArray)
  ranks <- array(dim=dims)
  for (k in seq_len(dims[3])) {
    # modified from samr:::resample
    ranks[,,k] <- matrixStats::rowRanks(infRepsArray[,,k] + 0.1 * runif(dims[1]*dims[2]))
  }
  fit <- samr:::wilcoxon.unpaired.seq.func(xresamp=ranks, y=condition)
  fit$tt
}

estimatePi0 <- function(stat, nulls.vec) {
  # modified from samr::samr
  qq <- quantile(nulls.vec, c(0.25, 0.75))
  2 * sum(stat > qq[1] & stat < qq[2])/length(stat)
}

makeLocFDR <- function(stat, nulls, pi0) {
  samr.obj <- makeSamrObj(stat, nulls, pi0)
  locfdr.obj <- samr:::localfdr(samr.obj, min.foldchange=0)
  locfdr.up <- samr:::predictlocalfdr(locfdr.obj$smooth.object, samr.obj$tt[samr.obj$tt >= 0])
  locfdr.lo <- samr:::predictlocalfdr(locfdr.obj$smooth.object, samr.obj$tt[samr.obj$tt < 0])
  locfdr <- numeric(length(stat))
  locfdr[stat >= 0] <- locfdr.up/100
  locfdr[stat < 0] <- locfdr.lo/100
  locfdr
}

makeQvalue <- function(stat, nulls, pi0) {
  samr.obj <- makeSamrObj(stat, nulls, pi0)
  delta.table <- samr:::samr.compute.delta.table.seq(samr.obj)
  sig <- list(pup=which(stat >= 0), plo=which(stat < 0))
  qlist <- samr:::qvalue.func(samr.obj, sig, delta.table)
  qvalue <- numeric(length(stat))
  qvalue[stat >= 0] <- qlist$qvalue.up/100
  qvalue[stat < 0] <- qlist$qvalue.lo/100
  qvalue <- pmin(qvalue, 1)
  qvalue
}

makeSamrObj <- function(stat, nulls, pi0) {
  samr.const.twoclass.unpaired.response <- "Two class unpaired"
  samr.obj <- list(
    resp.type=samr.const.twoclass.unpaired.response,
    tt=stat,
    ttstar0=nulls,
    foldchange.star=2^sign(stat),
    evo=rowMeans(apply(nulls, 2, sort)),
    pi0=pi0,
    assay.type="seq")
}
