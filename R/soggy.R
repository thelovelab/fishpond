soggy <- function(y, x, c) {
  y <- preprocess(y)
  ys <- y[mcols(y)$keep,]
  tt <- sapply(rep.idx, function(i) genefilter::rowttests(assays(ys)[[i]], ys[[x]])$statistic)
  stat <- rowMeans(tt)
  y <- postprocess(y, stat)
  y
}
