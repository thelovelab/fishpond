swish.pair <- function(infRepsArray, condition, pair,
                        nperms=30, wilcoxP, pc=5, quiet=FALSE) {
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
