pond <- function(x, y, c) {
  rep.idx <- seq_along(assayNames(x))
  keep.mat <- sapply(rep.idx, function(i) rowSums(assays(x)[[i]] >= 10) >= 3)
  keep <- apply(keep.mat, 1, all)
  xs <- x[keep,] # x subset
  tt <- sapply(rep.idx, function(i) genefilter::rowttests(assays(xs)[[i]], xs[[y]])$statistic)
  mcols(x)$stat <- numeric(nrow(x))
  mcols(x)$stat[keep] <- rowMeans(tt)
  mcols(x)$stat[!keep] <- NA
  x
}
