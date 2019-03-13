swishInterxPair <- function(infRepsArray, condition, covariate, pair,
                              nperms=30, pc=5, wilcoxP, quiet=FALSE) {
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
  stat <- getSamStat(lfcArray, group, wilcoxP=wilcoxP)
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
                            group[perms$perms[p,]], wilcoxP=wilcoxP)
  }
  if (!quiet) message("")
  list(stat=stat, log2FC=log2FC, nulls=nulls)
}

swishInterx <- function(infRepsArray, condition, covariate,
                        nperms=30, pc=5, wilcoxP, quiet=FALSE) {
  stopifnot(nlevels(covariate) == 2)
  if (!all(table(condition, covariate) > 0))
    stop("swish with interaction across two variables requires samples for each combination")
  dims <- dim(infRepsArray)
  # here the statistic is the difference between condition xLFC across covariate groups
  stat <- getDeltaLFC(infRepsArray, condition, covariate, pc, wilcoxP)
  # in this case, the stat and the reported log2FC are the same
  log2FC <- stat
  perms <- getPerms(covariate, nperms)
  nperms <- permsNote(perms, nperms)
  nulls <- matrix(nrow=dims[1], ncol=nperms)
  if (!quiet) message("Generating test statistics over permutations")
  for (p in seq_len(nperms)) {
    if (!quiet) progress(p, max.value=nperms, init=(p==1), gui=FALSE)
    nulls[,p] <- getDeltaLFC(infRepsArray, condition,
                             covariate[perms$perms[p,]],
                             pc, wilcoxP)
  }
  if (!quiet) message("")
  list(stat=stat, log2FC=log2FC, nulls=nulls)
}

getDeltaLFC <- function(infRepsArray, condition, covariate, pc, wilcoxP) {
  grp1 <- covariate == levels(covariate)[1]
  grp2 <- covariate == levels(covariate)[2]
  lfc1 <- getLog2FC(infRepsArray[,grp1,], condition[grp1], pc=pc, array=TRUE)
  lfc2 <- getLog2FC(infRepsArray[,grp2,], condition[grp2], pc=pc, array=TRUE)
  # here our statistic is the difference in LFC between the two groups
  if (is.null(wilcoxP)) {
    stat <- rowMeans(lfc2 - lfc1)
  } else {
    stat <- rowQuantilesTowardZero(lfc2 - lfc1, wilcoxP)
  }
  stat
}
