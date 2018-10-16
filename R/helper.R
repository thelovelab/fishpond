medianRatioScaledTPM <- function(infReps, length, meanDepth) {
  nreps <- length(infReps)
  for (k in seq_len(nreps)) {
    cat(k,"")
    tpm <- infReps[[k]] / length
    tpm <- t(t(tpm) / colSums(tpm)) * meanDepth
    # filtering here mostly for speed
    use <- rowSums(infReps[[k]] >= 10) >= 3
    loggeomeans <- rowMeans(log(tpm[use,]))
    sf <- apply(tpm[use,], 2, function(s) {
      exp(median((log(s) - loggeomeans)[is.finite(loggeomeans)]))
    })
    infReps[[k]] <- t( t(tpm)/sf )
  }
  infReps
}

preprocess <- function(y) {
  rep.idx <- seq_along(assayNames(y))
  # filtering: remains to see what we want to put here
  nreps <- length(assayNames(y))
  keep.mat <- matrix(nrow=nrow(y), ncol=nreps)
  for (k in seq_len(nreps)) {
    cat(k,"")
    keep.mat[,k] <- rowSums(assays(y)[[k]] >= 10) >= 3
  }
  keep <- apply(keep.mat, 1, all)
  mcols(y)$keep <- keep
  y
}

postprocess <- function(y, stat) {
  mcols(y)$stat <- numeric(nrow(y))
  mcols(y)$stat[mcols(y)$keep] <- stat
  mcols(y)$stat[!mcols(y)$keep] <- NA
}
