swishInterxPair <- function(infRepsArray, condition, covariate, pair,
                              nperms=100, pc=5, quiet=FALSE) {
  stopifnot(is.numeric(pair) | is.character(pair) | is.factor(pair))
  stopifnot(!anyNA(covariate))
  stopifnot(!anyNA(pair))
  pair <- as.integer(factor(pair))
  if (!all(table(pair, condition) == 1))
    stop("'pair' should have a single sample for both levels of 'condition'")
  stopifnot(nlevels(covariate) == 2)
  if (!all(table(pair, covariate) %in% c(0,2)))
    stop("'pair' should be nested within 'cov'")
  dims <- dim(infRepsArray)

  out <- getInterxPairStat(infRepsArray, condition, covariate, pair, pc)
  stat <- out$stat
  group <- out$group
  lfcArray <- out$lfcArray

  ranks <- out$ranks
  # old code re-computed ranks per permutation
  # ranks <- NULL
  
  grp1 <- group == levels(group)[1]
  grp2 <- group == levels(group)[2]
  lfc_mat <- apply(lfcArray[,grp2,], c(1,3), mean) -
            apply(lfcArray[,grp1,], c(1,3), mean)

  # the reported log2FC is the difference in the mean LFC between the two groups
  # the median here is taken over inferential replicates
  log2FC <- rowMedians(lfc_mat)

  # the permutation framework is to permute which pairs are in which group
  perms <- getPerms(group, nperms)
  nperms <- permsNote(perms, nperms)
  nulls <- matrix(nrow=dims[1], ncol=nperms)
  if (!quiet) message("generating test statistics over permutations")
  for (p in seq_len(nperms)) {
    if (!quiet) svMisc::progress(p, max.value=nperms, init=(p==1), gui=FALSE)
    nulls[,p] <- getSamStat(lfcArray,
                            group[perms$perms[p,]],
                            ranks=ranks)
  }
  if (!quiet) message("")
  list(stat=stat, log2FC=log2FC, nulls=nulls)
}

swishInterx <- function(infRepsArray, condition, covariate,
                        nperms=100, pc=5, nRandomPairs=30,
                        quiet=FALSE) {
  stopifnot(nlevels(covariate) == 2)
  stopifnot(!anyNA(covariate))
  if (!all(table(condition, covariate) > 0))
    stop("swish with interaction across two variables requires samples for each combination")
  dims <- dim(infRepsArray)
  
  tab <- table(condition, covariate)
  # if sizes are equal, don't need to remove samples for pairing
  all_equal <- all(tab[1,] == tab[2,])

  # don't have pairs, but instead we use pseudo-pairs multiple times
  stats <- matrix(nrow=dims[1], ncol=nRandomPairs)
  for (r in seq_len(nRandomPairs)) {
    # the easy case: balanced numbers of samples across condition
    if (all_equal) {
      pair <- getPseudoPair(condition, covariate)
      out <- getInterxPairStat(infRepsArray, condition, covariate,
                               pair, pc)
    } else {
      # random subsetting to balance groups, then use pseudo pairs
      idx <- randomSamplesToRemove(tab, condition, covariate)
      pair <- getPseudoPair(condition[-idx], covariate[-idx])
      out <- getInterxPairStat(infRepsArray[,-idx,],
                               condition[-idx], covariate[-idx],
                               pair, pc)
    }
    stats[,r] <- out$stat
  }
  stat <- rowMeans(stats)
  
  # any of the iterations works to define group
  # (this is an ordered integer vector of pairs
  # across covariate, with condition collapsed as LFCs)
  group <- out$group

  # log2FC is the difference between condition LFC across covariate groups
  log2FC <- getDeltaLFC(infRepsArray, condition, covariate, pc)

  # this permutation scheme is different than others in swish
  # (and slower) because we reform lfcArray inside the
  # permutation loop - necessary because of the random
  # pseudo pairs
  perms <- getPerms(group, nperms)
  nperms <- permsNote(perms, nperms)
  nulls <- matrix(nrow=dims[1], ncol=nperms)
  if (!quiet) message("generating test statistics over permutations")
  for (p in seq_len(nperms)) {
    if (!quiet) svMisc::progress(p, max.value=nperms, init=(p==1), gui=FALSE)
    if (all_equal) {
      # first draw a pseudo-pairing
      pair <- getPseudoPair(condition, covariate)
      out <- getInterxPairStat(infRepsArray, condition, covariate,
                               pair, pc)
    } else {
      idx <- randomSamplesToRemove(tab, condition, covariate)
      pair <- getPseudoPair(condition[-idx], covariate[-idx])
      out <- getInterxPairStat(infRepsArray[,-idx,],
                               condition[-idx], covariate[-idx],
                               pair, pc)
    }
    lfcArray <- out$lfcArray
    # then permute the pseudo-pairs across covariate
    nulls[,p] <- getSamStat(lfcArray,
                            group[perms$perms[p,]])
  }
  if (!quiet) message("")
  list(stat=stat, log2FC=log2FC, nulls=nulls)
}

getInterxPairStat <- function(infRepsArray, condition, covariate, pair, pc) {
  # 'lfcArray' is an array of genes x pair x inf rep
  # it is in the order of the pair (1,2,3,...)
  lfcArray <- getLog2FCPair(infRepsArray, condition, pair, pc, array=TRUE)
  dat <- data.frame(pair, covariate, stringsAsFactors=FALSE)
  dat <- dat[!duplicated(dat$pair),]
  dat <- dat[order(dat$pair),]
  group <- dat$covariate # this is now along 'lfcArray'
  stopifnot(length(group) == dim(lfcArray)[2])
  # here we perform Wilcoxon rank sum testing of the condition LFCs across group
  out <- getSamStat(lfcArray, group, returnRanks=TRUE)
  list(stat=out$stat, ranks=out$ranks, group=group, lfcArray=lfcArray)
}

getDeltaLFC <- function(infRepsArray, condition, covariate, pc) {
  grp1 <- covariate == levels(covariate)[1]
  grp2 <- covariate == levels(covariate)[2]
  lfc1 <- getLog2FC(infRepsArray[,grp1,], condition[grp1], pc=pc, array=TRUE)
  lfc2 <- getLog2FC(infRepsArray[,grp2,], condition[grp2], pc=pc, array=TRUE)
  # the difference in LFC between the two groups, median over inf reps
  rowMedians(lfc2 - lfc1)
}

getPseudoPair <- function(condition, covariate) {
  pair <- integer(length(condition))
  cond1 <- condition == levels(condition)[1]
  cond2 <- condition == levels(condition)[2]
  pair[cond2] <- seq_len(sum(cond2))
  for (i in 1:2) {
    grp <- covariate == levels(covariate)[i]
    pair[cond1 & grp] <- sample(x=pair[cond2 & grp],
                                size=sum(cond1 & grp),
                                replace=FALSE)
  }
  pair
}

randomSamplesToRemove <- function(tab, condition, covariate) {
  cond1 <- condition == levels(condition)[1]
  cond2 <- condition == levels(condition)[2]
  cov_lvls <- levels(covariate)
  idx <- numeric()
  # find covariate levels which are imbalanced
  for (j in which(tab[1,] != tab[2,])) {
    cond1small <- tab[1,j] < tab[2,j]
    if (cond1small) {
      idx <- c(idx, sample(which(cond2 & covariate == cov_lvls[j]),
                           tab[2,j] - tab[1,j],
                           replace=FALSE))
    } else {
      idx <- c(idx, sample(which(cond1 & covariate == cov_lvls[j]),
                           tab[1,j] - tab[2,j],
                           replace=FALSE))
    }
  }
  idx
}
