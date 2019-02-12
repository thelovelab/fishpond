#' swish: SAMseq With Inferential Samples Helps
#'
#' @param y a SummarizedExperiment containing the inferential replicate
#' matrices of median-ratio-scaled TPM as assays 'infRep1', 'infRep2', etc.
#' @param x the name of the condition variable. A factor with two
#' levels for a two group analysis (possible to adjust for covariate
#' or matched samples, see next two arguments)
#' @param cov the name of the covariate for adjustment.
#' If provided a stratified Wilcoxon in performed.
#' Cannot be used with \code{pair}
#' @param pair the name of the pair variable, which should be the
#' number of the pair. Can be an integer or factor.
#' If specified, a signed rank test is used
#' to build the statistic. All samples across \code{x} must be
#' pairs if this is specified. Cannot be used with \code{cov}.
#' @param nperms the number of permutations
#' @param wilcoxP the quantile of the Wilcoxon statistics across
#' inferential replicates to use as the test statistic.
#' If set to NULL this will use the mean over inferential replicates
#' @param estPi0 logical, whether to estimate pi0
#'
#' @return a SummarizedExperiment with metadata columns added:
#' the statistic (either a centered Wilcoxon Mann-Whitney
#' or a signed rank statistic, aggregated over inferential replicates),
#' a log2 fold change (the median over inferential replicates,
#' and averaged over pairs or groups (if groups, weighted by sample size),
#' the local FDR and q-value, as estimated by the \code{samr} package.
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
#'
#' # histogram of the swish statistics
#' hist(mcols(y)$stat, breaks=40, col="grey")
#' cols = rep(c("blue","purple","red"),each=2)
#' for (i in 1:6) {
#'   arrows(mcols(y)$stat[i], 20,
#'          mcols(y)$stat[i], 10,
#'          col=cols[i], length=.1, lwd=2)
#' }
#'
#' # plot inferential replicates
#' plotInfReps(y, 1, "condition")
#' plotInfReps(y, 3, "condition")
#' plotInfReps(y, 5, "condition")
#'
#' @importFrom graphics axis boxplot segments
#' @importFrom stats median quantile rpois runif
#' @importFrom utils head tail
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames
#' assays assays<- colData colData<- mcols mcols<-
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' 
#' @export
swish <- function(y, x, cov=NULL, pair=NULL,
                  nperms=30, wilcoxP=0.25,
                  estPi0=FALSE) {
  # 'cov' or 'pair' or neither, but not both
  stopifnot(is.null(cov) | is.null(pair))
  if (is.null(metadata(y)$preprocessed) || !metadata(y)$preprocessed) {
    y <- labelKeep(y)
  }
  ys <- y[mcols(y)$keep,]
  infRepIdx <- grep("infRep",assayNames(y))
  infRepError(infRepIdx)
  infReps <- assays(ys)[infRepIdx]
  infRepsArray <- abind::abind(as.list(infReps), along=3)
  condition <- colData(y)[[x]]
  stopifnot(is.factor(condition))
  stopifnot(nlevels(condition) == 2)
  if (is.null(cov) & is.null(pair)) {
    #######################
    ## simple two groups ##
    #######################
    stat <- getSamStat(infRepsArray, condition, wilcoxP)
    log2FC <- getLog2FC(infRepsArray, condition)
    perms <- samr:::getperms(condition, nperms)
    nperms <- permsNote(perms, nperms)
    nulls <- matrix(nrow=nrow(ys), ncol=nperms)
    for (p in seq_len(nperms)) {
      cat(p, "")
      nulls[,p] <- getSamStat(infRepsArray,
                              condition[perms$perms[p,]], wilcoxP)
    }
    cat("\n")
  } else if (is.null(pair)) {
    #########################
    ## stratified analysis ##
    #########################
    covariate <- colData(y)[[cov]]
    out <- swish.strat(infRepsArray, condition, covariate,
                       nperms=nperms, wilcoxP)
    stat <- out$stat
    log2FC <- out$log2FC
    nulls <- out$nulls
  } else {
    #####################
    ## paired analysis ##
    #####################
    pair <- colData(y)[[pair]]
    stopifnot(is.numeric(pair) | is.character(pair) | is.factor(pair)) 
    pair <- as.integer(factor(pair))
    if (!all(table(pair, condition) == 1)) {
      stop("'pair' should have a single sample for both levels of condition")
    }
    stat <- getSignedRank(infRepsArray, condition, pair, wilcoxP)
    log2FC <- getLog2FCPair(infRepsArray, condition, pair)
    cond.sign <- ifelse(condition == levels(condition)[1], 1, -1)
    perms <- samr:::compute.block.perms(cond.sign * pair, pair, nperms)
    nperms <- permsNote(perms, nperms)
    perms <- fixPerms(perms, condition, pair)
    nulls <- matrix(nrow=nrow(ys), ncol=nperms)
    for (p in seq_len(nperms)) {
      cat(p, "")
      nulls[,p] <- getSignedRank(infRepsArray, condition[perms[p,]],
                                 pair[perms[p,]], wilcoxP)
    }
    cat("\n")
  }
  nulls.vec <- as.vector(nulls)
  if (estPi0) {
    pi0 <- estimatePi0(stat, nulls.vec)
  } else {
    pi0 <- 1
  }
  locfdr <- makeLocFDR(stat, nulls, pi0)
  qvalue <- makeQvalue(stat, nulls, pi0)
  df <- data.frame(stat, log2FC, locfdr, qvalue)
  y <- postprocess(y, df)
  y
}

getSamStat <- function(infRepsArray, condition, p=NULL) {
  dims <- dim(infRepsArray)
  ranks <- array(dim=dims)
  for (k in seq_len(dims[3])) {
    # modified from samr:::resample
    ranks[,,k] <- matrixStats::rowRanks(infRepsArray[,,k] +
                    0.1 * runif(dims[1]*dims[2]))
  }
  cond2 <- condition == levels(condition)[2]
  rankSums <- sapply(seq_len(dims[3]), function(i) rowSums(ranks[,cond2,i]))
  # Wilcoxon, centered on 0:
  W <- rankSums - sum(cond2) * (dims[2] + 1)/2
  if (is.null(p)) {
    stat <- rowMeans(W)
  } else {
    stat <- rowQuantilesTowardZero(W, p)
  }
  stat
}

getLog2FC <- function(infRepsArray, condition, pc=5) {
  dims <- dim(infRepsArray)
  cond1 <- condition == levels(condition)[1]
  cond2 <- condition == levels(condition)[2]
  log2Cond1 <- log2(apply(infRepsArray[,cond1,],c(1,3),mean) + pc)
  log2Cond2 <- log2(apply(infRepsArray[,cond2,],c(1,3),mean) + pc)
  # median over inferential replicates
  matrixStats::rowMedians(log2Cond2 - log2Cond1)
}

getLog2FCPair <- function(infRepsArray, condition, pair, pc=5) {
  dims <- dim(infRepsArray)
  o <- order(condition, pair)
  if (!all(o == seq_along(condition))) {
    infRepsArray <- infRepsArray[,o,]
  }
  n <- dims[2]
  cond1 <- (1):(n/2)
  cond2 <- (n/2 + 1):(n)
  lfc.mat <- log2(infRepsArray[,cond2,] + pc) -
               log2(infRepsArray[,cond1,] + pc)
  # median over inferential replicates
  apply(lfc.mat, 1, median)
}

rowQuantilesTowardZero <- function(W, p) {
  stopifnot(p >= 0 & p <= 1)
  medW <- matrixStats::rowMedians(W)
  stat <- numeric(nrow(W))
  stat[medW >= 0] <- matrixStats::rowQuantiles(
                       W[medW >= 0,,drop=FALSE], probs=1-p
                     )
  stat[medW < 0] <- matrixStats::rowQuantiles(
                      W[medW < 0,,drop=FALSE], probs=1-p
                    )
  # prefer the median, if this is closer to zero
  idx <- abs(stat) > abs(medW)
  if (sum(idx) > 0) {
    stat[idx] <- medW[idx]
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
  locfdr.up <- samr:::predictlocalfdr(locfdr.obj$smooth.object,
                 samr.obj$tt[samr.obj$tt >= 0])
  locfdr.lo <- samr:::predictlocalfdr(locfdr.obj$smooth.object,
                 samr.obj$tt[samr.obj$tt < 0])
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

permsNote <- function(perms, nperms) {
  if (perms$nperms.act < nperms) {
    message("note: less permuatations are available than requested")
    perms$nperms.act
  } else {
    nperms
  }
}

fixPerms <- function(perms, condition, pair) {
  perms.in <- perms$perms
  perms <- matrix(rep(seq_along(condition), nrow(perms.in)),
                  nrow=nrow(perms.in), byrow=TRUE)
  for (i in seq_len(max(pair))) {
    idx1 <- perms.in[,2*i - 1] < 0
    idx2 <- pair == i
    perms[idx1,idx2] <- perms[idx1,idx2,drop=FALSE][,2:1]
  }
  perms
}
