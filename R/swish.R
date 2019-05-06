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
#' @param interaction logical, whether to perform a test of an interaction
#' between \code{x} and \code{cov}. These are different than the other
#' tests produced by the software, in that they focus on a difference
#' in the log2 fold change across levels of \code{x} when comparing
#' the two levels in \code{cov}. If \code{pair} is specified, this
#' will perform a Wilcoxon rank sum test on the two groups
#' of matched sample LFCs. If \code{pair} is not included, multiple
#' random pairs of samples within the two groups are chosen,
#' and again a Wilcoxon rank sum test compared the LFCs across groups.
#' @param nperms the number of permutations
#' @param estPi0 logical, whether to estimate pi0
#' @param qvaluePkg character, which package to use for q-value estimation,
#' \code{samr} or \code{qvalue}
#' @param pc pseudocount for finite estimation of \code{log2FC}, not used
#' in calculation of test statistics, \code{locfdr} or \code{qvalue}
#' @param nRandomPairs the number of random pseudo-pairs (only used with
#' \code{interaction=TRUE} and un-matched samples) to use to calculate
#' the test statistic
#' @param quiet display no messages
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
#' \code{swish} is described in the following reference:
#'
#' Anqi Zhu, Avi Srivastava, Joseph G Ibrahim, Rob Patro, Michael I Love
#' "Nonparametric expression analysis using inferential replicate counts"
#' bioRxiv (2019).
#' 
#' The \code{swish} method builds upon the \code{SAMseq} method,
#' and extends it by incorporating inferential uncertainty.
#' \code{swish} internally calls functions from the \code{samr}
#' package, for example, the calculation of local FDR and q-value.
#'
#' The citation for \code{SAMseq} is:
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
#' y <- swish(y, x="condition")
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
#' @importFrom graphics axis segments plot rect abline
#' @importFrom stats median quantile rpois runif
#' @importFrom utils head tail capture.output
#' @importFrom methods is
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames
#' assays assays<- colData colData<- mcols mcols<-
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' @importFrom gtools permutations
#' 
#' @export
swish <- function(y, x, cov=NULL, pair=NULL,
                  interaction=FALSE, nperms=30, 
                  estPi0=FALSE, qvaluePkg="qvalue",
                  pc=5, nRandomPairs=30,
                  quiet=FALSE) {

  stopifnot(is(y, "SummarizedExperiment"))
  # 'cov' or 'pair' or neither, but not both
  if (!interaction) stopifnot(is.null(cov) | is.null(pair))
  # interactions require a two level covariate
  if (interaction) stopifnot(!is.null(cov))
  if (!interactive()) { quiet <- TRUE }
  if (is.null(metadata(y)$preprocessed) || !metadata(y)$preprocessed) {
    y <- labelKeep(y)
  }
  
  if (!qvaluePkg %in% c("qvalue","samr")) {
    stop("'qvaluePkg' should be 'qvalue' or 'samr'")
  }
  if (qvaluePkg == "samr") {
    if (!requireNamespace("samr", quietly=TRUE)) {
      stop("first install the 'samr' package")
    }
  }
  
  ys <- y[mcols(y)$keep,]
  infRepsArray <- getInfReps(ys)
  stopifnot(x %in% names(colData(y)))
  condition <- colData(y)[[x]]
  stopifnot(is.factor(condition))
  stopifnot(nlevels(condition) == 2)

  if (!interaction & is.null(cov) & is.null(pair)) {
    # basic two group
    out <- swishTwoGroup(infRepsArray, condition,
                         nperms, pc, quiet)
    
  } else if (!interaction & !is.null(cov)) {
    # two group with covariate stratification
    stopifnot(cov %in% names(colData(y)))
    covariate <- colData(y)[[cov]] # covariate, e.g. batch effects
    out <- swishStrat(infRepsArray, condition, covariate,
                      nperms, pc, quiet)
    
  } else if (!interaction & !is.null(pair)) {
    # two group with matched samples
    stopifnot(pair %in% names(colData(y)))
    pair <- colData(y)[[pair]] # sample pairing
    out <- swishPair(infRepsArray, condition, pair,
                     nperms, pc, quiet)
    
  } else if (interaction & !is.null(pair)) {
    # two group 'x', two group 'cov', with matched samples
    stopifnot(cov %in% names(colData(y)))
    stopifnot(pair %in% names(colData(y)))
    covariate <- colData(y)[[cov]]
    pair <- colData(y)[[pair]]
    out <- swishInterxPair(infRepsArray, condition,
                           covariate, pair, nperms,
                           pc, quiet)
    
  } else if (interaction & is.null(pair)) {
    # two group 'x', two group 'cov', samples not matched
    stopifnot(cov %in% names(colData(y)))
    covariate <- colData(y)[[cov]]
    out <- swishInterx(infRepsArray, condition,
                       covariate, nperms,
                       pc, nRandomPairs,
                       quiet)
    
  }

  # gather results from functions above
  stat <- out$stat
  log2FC <- out$log2FC
  nulls <- out$nulls
  nulls.vec <- as.vector(nulls)
  if (qvaluePkg == "qvalue") {
    pvalue <- qvalue::empPvals(abs(stat), abs(nulls))
    pi0 <- if (estPi0) NULL else 1
    q.res <- qvalue::qvalue(pvalue, pi0=pi0)
    locfdr <- q.res$lfdr
    qvalue <- q.res$qvalues
    df <- data.frame(stat, log2FC, pvalue, locfdr, qvalue)
  } else if (qvaluePkg == "samr") {
    pi0 <- if (estPi0) estimatePi0(stat, nulls.vec) else 1
    locfdr <- makeLocFDR(stat, nulls, pi0)
    qvalue <- makeQvalue(stat, nulls, pi0, quiet)
    df <- data.frame(stat, log2FC, locfdr, qvalue)
  }
  postprocess(y, df)
}

getInfReps <- function(ys) {
  infRepIdx <- grep("infRep",assayNames(ys))
  infRepError(infRepIdx)
  infReps <- assays(ys)[infRepIdx]
  abind::abind(as.list(infReps), along=3)
}

swishTwoGroup <- function(infRepsArray, condition,
                          nperms=30, pc=5, quiet=FALSE) {
  dims <- dim(infRepsArray)
  stat <- getSamStat(infRepsArray, condition)
  log2FC <- getLog2FC(infRepsArray, condition, pc)
  perms <- getPerms(condition, nperms)
  nperms <- permsNote(perms, nperms)
  nulls <- matrix(nrow=dims[1], ncol=nperms)
  if (!quiet) message("Generating test statistics over permutations")
  for (p in seq_len(nperms)) {
    if (!quiet) svMisc::progress(p, max.value=nperms, init=(p==1), gui=FALSE)
    nulls[,p] <- getSamStat(infRepsArray,
                            condition[perms$perms[p,]])
  }
  if (!quiet) message("")
  list(stat=stat, log2FC=log2FC, nulls=nulls)
}

getSamStat <- function(infRepsArray, condition) {
  dims <- dim(infRepsArray)
  ranks <- array(dim=dims)
  for (k in seq_len(dims[3])) {
    # modified from samr:::resample
    ranks[,,k] <- matrixStats::rowRanks(infRepsArray[,,k] +
                    0.1 * runif(dims[1]*dims[2]))
  }
  cond2 <- condition == levels(condition)[2]
  rankSums <- vapply(seq_len(dims[3]), function(i)
    rowSums(ranks[,cond2,i]), numeric(dims[1]))
  # Wilcoxon, centered on 0:
  W <- rankSums - sum(cond2) * (dims[2] + 1)/2
  rowMeans(W)
}

getLog2FC <- function(infRepsArray, condition, pc=5, array=FALSE) {
  dims <- dim(infRepsArray)
  cond1 <- condition == levels(condition)[1]
  cond2 <- condition == levels(condition)[2]
  diffs <- matrix(nrow=dims[1],ncol=dims[3])
  for (k in seq_len(dims[3])) {
    diffs[,k] <- log2(rowMeans(infRepsArray[,cond2,k]) + pc) -
                 log2(rowMeans(infRepsArray[,cond1,k]) + pc)
  }
  if (array) {
    return(diffs)
  }
  # median over inferential replicates
  matrixStats::rowMedians(diffs)
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

makeQvalue <- function(stat, nulls, pi0, quiet=FALSE) {
  samr.obj <- makeSamrObj(stat, nulls, pi0)
  out <- capture.output({
    delta.table <- samr:::samr.compute.delta.table.seq(samr.obj)
  })
  if (!quiet) message(paste(out,collapse="\n"))
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
    message("note: less permutations are available than requested")
    perms$nperms.act
  } else {
    nperms
  }
}

getPerms <- function(condition, nperms) {
  if (length(condition) <= 6) {
    #perms <- samr:::getperms(condition, nperms)
    # get all possible permutations from gtools
    # for n=6 this is 720, and nperms is likely < 100
    out0 <- gtools::permutations(n=length(condition), r=length(condition))
    if (nrow(out0) > nperms) {
      idx <- sample(nrow(out0), nperms)
      out <- out0[idx,]
    } else {
      out <- out0
    }
    perms <- list(perms = out,
                  all.perms.flag = as.integer(nrow(out0) <= nperms),
                  nperms.act = nrow(out))
  } else {
    # with > 6 samples we have a good number of permutations (> 5,000),
    # just do random sampling
    x <- t(replicate(nperms, sample(length(condition))))
    perms <- list(perms = x,
                  all.perms.flag = 0,
                  nperms.act = nperms)
  }
  perms
}
