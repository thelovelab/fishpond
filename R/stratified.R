swish.strat <- function(infRepsArray, condition, covariate, nperms=30) {
  ngroups <- nlevels(covariate)
  groups <- levels(covariate)
  nr <- nrow(infRepsArray)
  stats <- matrix(nrow=nr, ncol=ngroups)
  nulls.big <- array(dim=list(nr, nperms, ngroups))
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
  nulls <- matrix(nrow=nr, ncol=nperms)
  for (p in seq_len(nperms)) {
    nulls[,p] <- nulls.big[,p,] %*% wts
  }
  browser()
  list(stat=stat, nulls=nulls)
}
