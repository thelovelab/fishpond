#' pool...
#'
#' @param y a SummarizedExperiment containing the inferential replicate
#' matrices of median-ratio-scaled TPM as assays
#' @param x the condition
#' @param cov a covariate for adjustment
#' @param nperms the number of permutations of x
#'
#' @return a SummarizedExperiment with metadata columns added
#'
#' @export
pool <- function(y, x, cov, nperms=5) {
  y <- preprocess(y)
  ys <- y[mcols(y)$keep,]
  getOrmStat(ys, x, cov)
  
  stat <- rowMeans(tt)
  y <- postprocess(y, stat)
  y
}

getOrmStat <- function(ys, x, cov) {
  nreps <- length(assayNames(ys))
  mm <- model.matrix(~cov + x)
  wald <- matrix(nrow=nrow(ys),ncol=nreps)
  for (k in seq_len(nreps)) {
    rowmaxs1 <- matrixStats::rowMaxs(assays(ys)[[k]][,x=="1"])
    rowmaxs2 <- matrixStats::rowMaxs(assays(ys)[[k]][,x=="2"])
    rowmins1 <- matrixStats::rowMins(assays(ys)[[k]][,x=="1"])
    rowmins2 <- matrixStats::rowMins(assays(ys)[[k]][,x=="2"])
    separation <- rowmaxs1 < rowmins2 | rowmaxs2 < rowmins1
    for (i in 1:20) {#seq_len(nrow(dat))) {
      cat(i, "")
      # singularity problem when complete separation 
      if (separation[i]) {
        # TODO fix... regularize coefs
        wald[i,k] <- 100
      } else {
        z <- assays(ys)[[k]][i,] + 0.1 * runif(ncol(ys))
        # ...this is too slow and often has singularity problems
        suppressWarnings({
          fit <- rms::orm.fit(x=mm, y=z)
        })
        coef <- unname(fit$coefficients)
        V <- unname(fit$var)
        wald[i,k] <- coef[length(coef)] / sqrt(V[nrow(V),nrow(V)])
      }
    }
  }
  rowMeans(wald)
}
