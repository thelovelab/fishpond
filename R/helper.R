preprocess <- function(y) {
  rep.idx <- seq_along(assayNames(y))
  keep.mat <- sapply(rep.idx, function(i) rowSums(assays(y)[[i]] >= 10) >= 3)
  keep <- apply(keep.mat, 1, all)
  mcols(y)$keep <- keep
  y
}

postprocess <- function(y, stat) {
  mcols(y)$stat <- numeric(nrow(y))
  mcols(y)$stat[mcols(y)$keep] <- stat
  mcols(y)$stat[!mcols(y)$keep] <- NA
}
