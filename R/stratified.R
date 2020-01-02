swishStrat <- function(infRepsArray, condition, covariate,
                        nperms=100, pc=5, fast, quiet=FALSE) {
  stopifnot(is.factor(covariate))
  ngroups <- nlevels(covariate)
  groups <- levels(covariate)
  nr <- nrow(infRepsArray)
  stats <- matrix(nrow=nr, ncol=ngroups)
  lfc.mat <- matrix(nrow=nr, ncol=ngroups)
  nulls.big <- array(dim=list(nr, nperms, ngroups))
  for (i in seq_len(ngroups)) {
    g <- groups[i]
    infRepsArray.sub <- infRepsArray[,covariate == g,]
    cond.sub <- condition[covariate == g]
    perms <- getPerms(cond.sub, nperms)
    nperms <- permsNote(perms, nperms)

    # if fast==1, avoid re-computing the ranks for the permutation distribution
    if (fast == 1) {
      out <- getSamStat(infRepsArray.sub, cond.sub, returnRanks=TRUE)
      stats[,i] <- out$stat
      ranks <- out$ranks
    } else {
      stats[,i] <- getSamStat(infRepsArray.sub, cond.sub)
      ranks <- NULL
    }
    lfc.mat[,i] <- getLog2FC(infRepsArray.sub, cond.sub, pc)
    
    if (!quiet) message(paste0(
                  "Generating test statistics over permutations: ",
                  i,"/",ngroups," groups"))
    for (p in seq_len(nperms)) {
      if (!quiet) svMisc::progress(p, max.value=nperms, init=(p==1), gui=FALSE)
      nulls.big[,p,i] <- getSamStat(infRepsArray.sub,
                                    cond.sub[perms$perms[p,]],
                                    ranks=ranks)
    }
    if (!quiet) message("")
  }
  ns <- unname(table(covariate))
  wts <- 1/(ns + 1)
  stat <- as.vector(stats %*% wts)
  lfc.wts <- ns/sum(ns)
  log2FC <- as.vector(lfc.mat %*% lfc.wts)
  nulls <- matrix(nrow=nr, ncol=nperms)
  for (p in seq_len(nperms)) {
    nulls[,p] <- nulls.big[,p,] %*% wts
  }
  list(stat=stat, log2FC=log2FC, nulls=nulls)
}
