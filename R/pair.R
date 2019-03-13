swishPair <- function(infRepsArray, condition, pair,
                        nperms=30, pc=5, wilcoxP, quiet=FALSE) {
  stopifnot(is.numeric(pair) | is.character(pair) | is.factor(pair)) 
  pair <- as.integer(factor(pair))
  if (!all(table(pair, condition) == 1))
    stop("'pair' should have a single sample for both levels of condition")
  stat <- getSignedRank(infRepsArray, condition, pair, wilcoxP)
  log2FC <- getLog2FCPair(infRepsArray, condition, pair, pc)
  cond.sign <- ifelse(condition == levels(condition)[1], 1, -1)
  perms <- samr:::compute.block.perms(cond.sign * pair, pair, nperms)
  nperms <- permsNote(perms, nperms)
  perms <- fixPerms(perms, condition, pair)
  nulls <- matrix(nrow=nrow(infRepsArray), ncol=nperms)
  if (!quiet) message("Generating test statistics over permutations")
  for (p in seq_len(nperms)) {
    if (!quiet) progress(p, max.value=nperms, init=(p==1), gui=FALSE)
    nulls[,p] <- getSignedRank(infRepsArray, condition[perms[p,]],
                               pair[perms[p,]], wilcoxP)
  }
  if (!quiet) message("")
  list(stat=stat, log2FC=log2FC, nulls=nulls)
}

getSignedRank <- function(infRepsArray, condition, pair, p=NULL) {
  dims <- dim(infRepsArray)
  sgn.ranks <- array(dim=c(dims[1],dims[2]/2,dims[3]))
  o <- order(condition, pair)
  grp1 <- head(o, length(condition)/2)
  grp2 <- tail(o, length(condition)/2)
  for (k in seq_len(dims[3])) {
    diff <- infRepsArray[,grp2,k] - infRepsArray[,grp1,k]
    sgn.ranks[,,k] <- sign(diff) * matrixStats::rowRanks(abs(diff) +
                        0.1 * runif(dims[1]*dims[2]/2))
  }
  # sums of signed rank, expectation is 0
  W <- apply(sgn.ranks, c(1,3), sum)
  if (is.null(p)) {
    stat <- rowMeans(W)
  } else {
    stat <- rowQuantilesTowardZero(W, p)
  }
  stat
}

getLog2FCPair <- function(infRepsArray, condition, pair, pc=5, array=FALSE) {
  dims <- dim(infRepsArray)
  o <- order(condition, pair)
  if (!all(o == seq_along(condition))) {
    infRepsArray <- infRepsArray[,o,]
  }
  n <- dims[2]
  cond1 <- (1):(n/2)
  cond2 <- (n/2 + 1):(n)
  lfcArray <- log2(infRepsArray[,cond2,] + pc) -
              log2(infRepsArray[,cond1,] + pc)
  if (array) {
    return(lfcArray)
  }
  # mean over samples
  lfcMat <- apply(lfcArray, c(1,3), mean)
  # median over inferential replicates
  matrixStats::rowMedians(lfcMat)
}
