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
    perms <- samr:::getperms(cond.sub, nperms)
    stats[,i] <- getSamStat(infRepsArray.sub, cond.sub)
    for (p in seq_len(nperms)) {
      cat(p, "")
      nulls.big[,p,i] <- getSamStat(infRepsArray.sub, cond.sub[perms$perms[p,]])
    }
    cat("\n")
  }
  ns <- unname(table(covariate))
  wts <- 1/(ns + 1)
  stat <- as.vector(stats %*% wts)
  nulls <- matrix(nrow=nr, ncol=2*nperms)
  nulls <- t(apply(nulls.big, 1, cbind)) * sum(wts)
  list(stat=stat, nulls=nulls)
}
