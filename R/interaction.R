swish.interx.pair <- function(infRepsArray, condition, covariate, pair,
                              nperms=30, wilcoxP, pc=5, quiet=FALSE) {
  stopifnot(is.numeric(pair) | is.character(pair) | is.factor(pair)) 
  pair <- as.integer(factor(pair))
  if (!all(table(pair, condition) == 1))
    stop("'pair' should have a single sample for both levels of 'condition'")
  stopifnot(nlevels(covariate) == 2)
  if (!all(table(pair, covariate) %in% c(0,2)))
    stop("'pair' should be nested within 'covariate'")
  # 'lfc.array' is an array of genes x pair x inf rep
  # it is in the order of the pair (1,2,3,...)
  dims <- dim(infRepsArray)
  lfcArray <- getLog2FCPair(infRepsArray, condition, pair, pc, array=TRUE)
  dat <- data.frame(pair, covariate, stringsAsFactors=FALSE)
  dat <- dat[!duplicated(dat$pair),]
  dat <- dat[order(dat$pair),]
  group <- dat$covariate # this is now along 'lfc.array'
  stopifnot(length(group) == dim(lfcArray)[2])
  stat <- getSamStat(lfcArray, group, p=wilcoxP)
  grp1 <- group == levels(group)[1]
  grp2 <- group == levels(group)[2]
  lfcMat <- apply(lfcArray[,grp2,], c(1,3), mean) -
            apply(lfcArray[,grp1,], c(1,3), mean)
  log2FC <- matrixStats::rowMedians(lfcMat)
  perms <- getPerms(group, nperms)
  nperms <- permsNote(perms, nperms)
  nulls <- matrix(nrow=dims[1], ncol=nperms)
  if (!quiet) message("Generating test statistics over permutations")
  for (p in seq_len(nperms)) {
    if (!quiet) progress(p, max.value=nperms, init=(p==1), gui=FALSE)
    nulls[,p] <- getSamStat(lfcArray,
                            group[perms$perms[p,]], wilcoxP)
  }
  if (!quiet) message("")
  list(stat=stat, log2FC=log2FC, nulls=nulls)
}
