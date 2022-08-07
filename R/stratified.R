swishStrat <- function(infRepsArray, condition, covariate,
                        nperms=100, pc=5, quiet=FALSE) {
  stopifnot(is.factor(covariate))
  stopifnot(!anyNA(covariate))
  ngroups <- nlevels(covariate)
  groups <- levels(covariate)
  nr <- nrow(infRepsArray)
  stats <- matrix(nrow=nr, ncol=ngroups)
  lfc_mat <- matrix(nrow=nr, ncol=ngroups)
  nulls_big <- array(dim=list(nr, nperms, ngroups))
  for (i in seq_len(ngroups)) {
    g <- groups[i]
    infRepsArray.sub <- infRepsArray[,covariate == g,]
    cond.sub <- condition[covariate == g]
    perms <- getPerms(cond.sub, nperms)
    nperms <- permsNote(perms, nperms)

    # avoid re-computing the ranks for the permutation distribution
    out <- getSamStat(infRepsArray.sub, cond.sub, returnRanks=TRUE)
    stats[,i] <- out$stat
    ranks <- out$ranks

    # old code would re-compute ranks for every inf rep
    ## stats[,i] <- getSamStat(infRepsArray.sub, cond.sub)
    ## ranks <- NULL

    lfc_mat[,i] <- getLog2FC(infRepsArray.sub, cond.sub, pc)
    
    if (!quiet) message(paste0(
                  "Generating test statistics over permutations: ",
                  i,"/",ngroups," groups"))
    for (p in seq_len(nperms)) {
      if (!quiet) svMisc::progress(p, max.value=nperms, init=(p==1), gui=FALSE)
      nulls_big[,p,i] <- getSamStat(infRepsArray.sub,
                                    cond.sub[perms$perms[p,]],
                                    ranks=ranks)
    }
    if (!quiet) message("")
  }
  ns <- unname(table(covariate))
  wts <- 1/(ns + 1)
  stat <- as.vector(stats %*% wts)
  lfc_weights <- ns/sum(ns)
  log2FC <- as.vector(lfc_mat %*% lfc_weights)
  nulls <- matrix(nrow=nr, ncol=nperms)
  for (p in seq_len(nperms)) {
    nulls[,p] <- nulls_big[,p,] %*% wts
  }
  list(stat=stat, log2FC=log2FC, nulls=nulls)
}
