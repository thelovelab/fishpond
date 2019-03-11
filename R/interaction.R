swishInterxPair <- function(infRepsArray, condition, covariate, pair,
                              nperms=30, wilcoxP, pc=5, quiet=FALSE) {
  stopifnot(is.numeric(pair) | is.character(pair) | is.factor(pair)) 
  pair <- as.integer(factor(pair))
  if (!all(table(pair, condition) == 1))
    stop("'pair' should have a single sample for both levels of 'condition'")
  stopifnot(nlevels(covariate) == 2)
  if (!all(table(pair, covariate) %in% c(0,2)))
    stop("'pair' should be nested within 'covariate'")
  dims <- dim(infRepsArray)
  # 'lfcArray' is an array of genes x pair x inf rep
  # it is in the order of the pair (1,2,3,...)
  lfcArray <- getLog2FCPair(infRepsArray, condition, pair, pc, array=TRUE)
  dat <- data.frame(pair, covariate, stringsAsFactors=FALSE)
  dat <- dat[!duplicated(dat$pair),]
  dat <- dat[order(dat$pair),]
  group <- dat$covariate # this is now along 'lfcArray'
  stopifnot(length(group) == dim(lfcArray)[2])
  # here we perform Wilcoxon rank sum testing of the condition LFCs across group
  stat <- getSamStat(lfcArray, group, p=wilcoxP)
  grp1 <- group == levels(group)[1]
  grp2 <- group == levels(group)[2]
  lfcMat <- apply(lfcArray[,grp2,], c(1,3), mean) -
            apply(lfcArray[,grp1,], c(1,3), mean)
  # the reported log2FC is the difference in the mean LFC between the two groups
  # the median here is taken over inferential replicates
  log2FC <- matrixStats::rowMedians(lfcMat)
  # the permutation framework is to permute which pairs are in which group
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

swishInterx <- function(infRepsArray, condition, covariate,
                        nperms=30, wilcoxP, pc=5, quiet=FALSE) {
  stopifnot(nlevels(covariate) == 2)
  if (!all(table(condition, covariate) > 0))
    stop("swish with interaction across two variables requires samples for each combination")
  dims <- dim(infRepsArray)
  return(NULL)
  #if (!quiet) message("")
  #list(stat=stat, log2FC=log2FC, nulls=nulls)
}
