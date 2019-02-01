#' swish: SAMseq With Inferential Samples Helps
#'
#' @param y a SummarizedExperiment containing the inferential replicate
#' matrices of median-ratio-scaled TPM as assays 'infRep1', 'infRep2', etc.
#' @param x the name of the condition variable
#' @param cov the name of the covariate for adjustment.
#' if provided a stratified Wilcoxon in performed
#' @param nperms the number of permutations
#' @param wilcoxP the quantile of the Wilcoxon statistics across
#' inferential replicates to use as the test statistic.
#' If set to NULL this will use the mean over inferential replicates
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
#' @examples
#'
#' library(SummarizedExperiment)
#' set.seed(1)
#' y <- makeSimSwishData()
#' y <- scaleInfReps(y)
#' y <- labelKeep(y)
#' y <- swish(y, "condition")
#' stat <- mcols(y)$stat
#' hist(stat,breaks=40,col="grey")
#' cols = rep(c("blue","purple","red"),each=2)
#' for (i in 1:6) {
#'   arrows(stat[i], 20, stat[i], 10, col=cols[i], length=.1, lwd=2)
#' }
#' plotInfReps(y, "condition", 1)
#' plotInfReps(y, "condition", 3)
#' plotInfReps(y, "condition", 5)
#' 
#' @export
swish <- function(y, x, cov=NULL, nperms=30, wilcoxP=0.25, estPi0=FALSE) {
  if (is.null(metadata(y)$preprocessed) || !metadata(y)$preprocessed) {
    y <- preprocess(y)
  }
  ys <- y[mcols(y)$keep,]
  infReps <- assays(ys)[grep("infRep",assayNames(ys))]
  infRepsArray <- abind::abind(as.list(infReps), along=3)
  condition <- colData(y)[[x]] 
  if (is.null(cov)) {
    stat <- getSamStat(infRepsArray, condition, wilcoxP)
    perms <- samr:::getperms(condition, nperms)
    nulls <- matrix(nrow=nrow(ys), ncol=nperms)
    for (p in seq_len(nperms)) {
      cat(p, "")
      nulls[,p] <- getSamStat(infRepsArray, condition[perms$perms[p,]], wilcoxP)
    }
    cat("\n")
  } else {
    covariate <- colData(y)[[cov]]
    out <- swish.strat(infRepsArray, condition, covariate, nperms=nperms, wilcoxP)
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

getSamStat <- function(infRepsArray, condition, p=NULL) {
  dims <- dim(infRepsArray)
  ranks <- array(dim=dims)
  for (k in seq_len(dims[3])) {
    # modified from samr:::resample
    ranks[,,k] <- matrixStats::rowRanks(infRepsArray[,,k] + 0.1 * runif(dims[1]*dims[2]))
  }
  cond2 <- condition == levels(condition)[2]
  rankSums <- sapply(seq_len(dims[3]), function(i) rowSums(ranks[,cond2,i]))
  # Wilcoxon, centered on 0:
  W <- rankSums - sum(cond2) * (dims[2] + 1)/2
  if (is.null(p)) {
    stat <- rowMeans(W)
  } else {
    stopifnot(p >= 0 & p <= 1)
    medW <- matrixStats::rowMedians(W)
    stat <- numeric(dims[1])
    stat[medW >= 0] <- matrixStats::rowQuantiles(W[medW >= 0,,drop=FALSE], probs=1-p)
    stat[medW < 0] <- matrixStats::rowQuantiles(W[medW < 0,,drop=FALSE], probs=1-p)
    idx <- abs(stat) > abs(medW)
    if (sum(idx) > 0) {
      stat[idx] <- medW[idx]
    }
  }
  stat
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
