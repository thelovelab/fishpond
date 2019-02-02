getSignedRank <- function(infRepsArray, condition, pair, p=NULL) {
  dims <- dim(infRepsArray)
  sgn.ranks <- array(dim=c(dims[1],dims[2]/2,dims[3]))
  o <- order(condition, pair)
  grp1 <- head(o, length(condition)/2)
  grp2 <- tail(o, length(condition)/2)
  for (k in seq_len(dims[3])) {
    diff <- infRepsArray[,grp2,k] - infRepsArray[,grp1,k]
    sgn.ranks[,,k] <- sign(diff) * matrixStats::rowRanks(abs(diff) + 0.1 * runif(dims[1]*dims[2]/2))
  }
  # sums of signed rank, expectation is 0
  W <- apply(sgn.ranks, c(1,3), sum)
  if (is.null(p)) {
    stat <- rowMeans(W)
  } else {
    stat <- rowQuantilesTowardZero(W, p)
  }
  stat
}
