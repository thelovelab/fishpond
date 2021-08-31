swishCor <- function(infRepsArray, condition,
                     nperms=100, pc=5, quiet=FALSE) {
  dims <- dim(infRepsArray)
  out <- getCorStat(infRepsArray, condition)
  stat <- out$stat
  ranks <- out$ranks
  # log2FC for a numeric x is calculated differently
  log2FC <- getLog2FCNumX(infRepsArray, condition, pc)
  perms <- getPerms(condition, nperms)
  nperms <- permsNote(perms, nperms)
  if (!quiet) message("Generating test statistics over permutations")
  nulls <- matrix(nrow=dims[1], ncol=nperms)
  for (p in seq_len(nperms)) {
    if (!quiet) svMisc::progress(p, max.value=nperms, init=(p==1), gui=FALSE)
    nulls[,p] <- getCorStat(infRepsArray,
                            condition[perms$perms[p,]],
                            ranks=ranks)
  }
  if (!quiet) message("")
  list(stat=stat, log2FC=log2FC, nulls=nulls)
}

getCorStat <- function(infRepsArray, condition, ranks=NULL) {
  dims <- dim(infRepsArray)
  rcondition <- rank(condition)
  if (dims[2] == 2) stop("too few samples to compute the correlation statistic")
  # calculate ranks if they are not provided...
  ranksMissing <- is.null(ranks)
  if (ranksMissing) {
    ranks <- array(dim=dims)
    for (k in seq_len(dims[3])) {
      # modified from samr:::resample
      ranks[,,k] <- matrixStats::rowRanks(infRepsArray[,,k] +
                                          0.1 * runif(dims[1]*dims[2]),
                                          ties.method = "average")
    }
  }
  corrs <- rowMeans(sapply(seq_len(dims[3]), function(k) cor(t(ranks[,,k]), rcondition)))
  if (ranksMissing) {
    return(list(stat=corrs, ranks=ranks))
  } else {
    return(corrs)
  }
}

getLog2FCNumX <- function(infRepsArray, condition, pc=5) {
  dims <- dim(infRepsArray)
  coefs <- matrix(nrow=dims[1],ncol=dims[3])
  # center condition
  condition <- condition - mean(condition)
  xtxi <- solve(t(condition) %*% condition)
  for (k in seq_len(dims[3])) {
    logCounts <- log2(infRepsArray[,,k] + pc)
    logCounts <- t(t(logCounts) - rowMeans(logCounts))
    coefs[,k] <- t(xtxi %*% t(condition) %*% t(logCounts))
  }
  # median over inferential replicates
  matrixStats::rowMedians(coefs)
}
