pond <- function(y, x, c) {
  rep.idx <- seq_along(assayNames(y))
  keep.mat <- sapply(rep.idx, function(i) rowSums(assays(y)[[i]] >= 10) >= 3)
  keep <- apply(keep.mat, 1, all)
  ys <- y[keep,] # y subset
  tt <- sapply(rep.idx, function(i) genefilter::rowttests(assays(ys)[[i]], ys[[x]])$statistic)
  mcols(y)$stat <- numeric(nrow(y))
  mcols(y)$stat[keep] <- rowMeans(tt)
  mcols(y)$stat[!keep] <- NA
  y
}
