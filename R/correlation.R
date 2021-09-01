swishCor <- function(infRepsArray, condition, cor,
                     nperms=100, pc=5, quiet=FALSE) {
  dims <- dim(infRepsArray)
  out <- getCorStat(infRepsArray, condition, cor, pc)
  stat <- out$stat
  ranks <- out$ranks
  # log2FC for a numeric x is calculated differently
  log2FC <- getLog2FCNumX(infRepsArray, condition, pc)
  perms <- getPerms(condition, nperms)
  nperms <- permsNote(perms, nperms)
  if (!quiet) message("generating test statistics over permutations")
  nulls <- matrix(nrow=dims[1], ncol=nperms)
  for (p in seq_len(nperms)) {
    if (!quiet) svMisc::progress(p, max.value=nperms, init=(p==1), gui=FALSE)
    nulls[,p] <- getCorStat(infRepsArray, condition[perms$perms[p,]],
                            cor, pc, ranks=ranks)
  }
  if (!quiet) message("")
  list(stat=stat, log2FC=log2FC, nulls=nulls)
}

# correlation of covariate 'cov' with LFC of pairs
swishCorPair <- function(infRepsArray, condition,
                         pair, covariate, cor,
                         nperms=100, pc=5, quiet=FALSE) {
  checkPair(pair, condition)
  # check that covariate has same value across the pairs
  o <- order(condition, pair)
  n <- dim(infRepsArray)[2]
  cond1 <- (1):(n/2)
  cond2 <- (n/2 + 1):(n)
  stopifnot(all(covariate[o][cond1] == covariate[o][cond2]))
  # values for the first of the pair, re-ordered by condition and pair,
  # as occurs within getLog2FCPair()
  covariate1 <- covariate[o][cond1] 
  # get the LFC array, both for computing the test statistic...
  lfcArray <- getLog2FCPair(infRepsArray, condition, pair, pc, array=TRUE)
  out <- getCorStat(lfcArray, covariate1, cor, pc, noise=FALSE)
  stat <- out$stat
  ranks <- out$ranks
  # ...and for computing the mean LFC over samples
  lfcMat <- apply(lfcArray, c(1,3), mean)
  log2FC <- matrixStats::rowMedians(lfcMat)
  perms <- getPerms(covariate1, nperms)
  nperms <- permsNote(perms, nperms)
  # now just work with the permutations matrix
  perms <- perms$perms
  nulls <- matrix(nrow=nrow(lfcArray), ncol=nperms)
  if (!quiet) message("generating test statistics over permutations")
  for (p in seq_len(nperms)) {
    if (!quiet) svMisc::progress(p, max.value=nperms, init=(p==1), gui=FALSE)
    nulls[,p] <- getCorStat(lfcArray, covariate1[perms[p,]],
                            cor, pc, ranks=ranks, noise=FALSE)
  }
  if (!quiet) message("")
  list(stat=stat, log2FC=log2FC, nulls=nulls)
}

# get correlation test statistics
getCorStat <- function(infRepsArray, condition, cor, pc, ranks=NULL, noise=TRUE) {
  dims <- dim(infRepsArray)
  if (dims[2] == 2) stop("too few samples to compute the correlation statistic")
  ranksMissing <- is.null(ranks)
  if (cor == "spearman") {
    rcondition <- rank(condition)  
    # calculate ranks if they are not provided...
    if (ranksMissing) {
      ranks <- array(dim=dims)
      for (k in seq_len(dims[3])) {
        # modified from samr:::resample
        if (noise) {
          ranks[,,k] <- matrixStats::rowRanks(infRepsArray[,,k] +
                                              0.1 * runif(dims[1]*dims[2]),
                                              ties.method = "average")
        } else {
          ranks[,,k] <- matrixStats::rowRanks(infRepsArray[,,k],
                                              ties.method = "average")
        }
      }
    }
  } else {
    # we just pass along original data and covariate for cor="pearson"
    rcondition <- condition
    ranks <- log2(infRepsArray + pc)
  }
  corrs <- rowMeans(sapply(seq_len(dims[3]), function(k) cor(t(ranks[,,k]), rcondition)))
  if (ranksMissing) {
    return(list(stat=corrs, ranks=ranks))
  } else {
    return(corrs)
  }
}

# define LFC for numeric 'x' using regression of log counts
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
