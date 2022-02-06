#' Make pseudo-inferential replicates from mean and variance
#'
#' Makes pseudo-inferential replicate counts from
#' \code{mean} and \code{variance} assays. The simulated
#' counts are drawn from a negative binomial distribution,
#' with \code{mu=mean} and \code{size} set using a method
#' of moments estimator for dispersion.
#'
#' Note that these simulated counts only reflect marginal
#' variance (one transcript or gene at a time),
#' and do not capture the covariance of counts across
#' transcripts or genes, unlike imported inferential
#' replicate data. Therefore, \code{makeInfReps} should
#' not be used with \code{summarizeToGene} to create
#' gene-level inferential replicates if inferential
#' replicates were originally created on the transcript
#' level. Instead, import the original inferential
#' replicates.
#'
#' @param y a SummarizedExperiment
#' @param numReps how many inferential replicates
#' @param minDisp the minimal dispersion value,
#' set after method of moments estimation from
#' inferential mean and variance
#'
#' @return a SummarizedExperiment
#'
#' @references
#'
#' Van Buren, S., Sarkar, H., Srivastava, A., Rashid, N.U.,
#' Patro, R., Love, M.I. (2020)
#' Compression of quantification uncertainty for scRNA-seq counts.
#' bioRxiv.
#' \url{https://doi.org/10.1101/2020.07.06.189639}
#' 
#' @examples
#'
#' library(SummarizedExperiment)
#' mean <- matrix(1:4,ncol=2)
#' variance <- mean
#' se <- SummarizedExperiment(list(mean=mean, variance=variance))
#' se <- makeInfReps(se, numReps=50)
#' 
#' @export
makeInfReps <- function(y, numReps, minDisp=1e-3) {
  stopifnot(is(y, "SummarizedExperiment"))
  stopifnot("mean" %in% assayNames(y))
  stopifnot("variance" %in% assayNames(y))
  stopifnot(numReps > 1)
  stopifnot(round(numReps) == numReps)
  if (any(grepl("infRep", assayNames(y)))) {
    stop("infReps already exist, remove these first")
  }
  # drop sparsity here
  m <- as.matrix(assays(y)[["mean"]])
  v <- as.matrix(assays(y)[["variance"]])
  disp <- ifelse(m > 0, pmax(minDisp, (v - m)/m^2), minDisp)
  infReps <- list()
  for (k in seq_len(numReps)) {
    infReps[[k]] <- matrix(rnbinom(n=nrow(y)*ncol(y), mu=m, size=1/disp),
                           ncol=ncol(y), dimnames=dimnames(m))
  }
  names(infReps) <- paste0("infRep", 1:numReps)
  assays(y) <- c(assays(y), infReps)
  metadata(y)$infRepsScaled <- FALSE
  y
}

#' Function for splitting SummarizedExperiment into separate RDS files
#'
#' The \code{splitSwish} function splits up the \code{y} object
#' along genes and writes a \code{Snakefile} that can be used with
#' Snakemake to distribute running \code{swish} across genes.
#' This workflow is primarily designed for large single cell datasets,
#' and so the default is to not perform length correction
#' within the distributed jobs.
#' See the alevin section of the vignette for an example. See
#' the Snakemake documention for details on how to run and customize
#' a \code{Snakefile}: \url{https://snakemake.readthedocs.io}
#' 
#' @param y a SummarizedExperiment
#' @param nsplits integer, how many pieces to break \code{y} into
#' @param prefix character, the path of the RDS files to write out,
#' e.g. \code{prefix="/path/to/swish"} will generate \code{swish.rds}
#' files at this path
#' @param snakefile character, the path of a Snakemake file, e.g.
#' \code{Snakefile}, that should be written out. If \code{NULL},
#' then no \code{Snakefile} is written out
#' @param overwrite logical, whether the \code{snakefile} and
#' RDS files (\code{swish1.rds}, ...) should overwrite existing files
#'
#' @references
#'
#' Compression and splitting across jobs:
#' 
#' Van Buren, S., Sarkar, H., Srivastava, A., Rashid, N.U.,
#' Patro, R., Love, M.I. (2020)
#' Compression of quantification uncertainty for scRNA-seq counts.
#' bioRxiv.
#' \url{https://doi.org/10.1101/2020.07.06.189639}
#'
#' Snakemake:
#' 
#' Koster, J., Rahmann, S. (2012)
#' Snakemake - a scalable bioinformatics workflow engine.
#' Bioinformatics.
#' \url{https://doi.org/10.1093/bioinformatics/bts480}
#' 
#' @return nothing, files are written out
#'
#' @export
splitSwish <- function(y, nsplits, prefix="swish",
                       snakefile=NULL, overwrite=FALSE) {
  stopifnot(nsplits > 1)
  stopifnot(nsplits == round(nsplits))
  stopifnot(nsplits < nrow(y))
  stopifnot(!is.null(rownames(y)))
  stopifnot(all(mcols(y)$keep))
  stopifnot(basename(prefix) == "swish")
  if (!is.null(snakefile)) {
    stopifnot(is(snakefile, "character"))
    if (file.exists(snakefile) & !overwrite)
      stop("snakefile already exists at specified location, see 'overwrite'")
    snake <- scan(system.file("extdata/Snakefile", package="fishpond"),
                  what="character", sep="\n", blank.lines.skip=FALSE,
                  quiet=TRUE)
    write(snake, file=snakefile)
  }
  # how many leading 0's
  width <- floor(log10(nsplits)) + 1
  nums <- formatC(seq_len(nsplits), width=max(2, width),
                  format="d", flag="0")
  files <- paste0(prefix, nums, ".rds")
  if (any(file.exists(files)) & !overwrite)
    stop("swish RDS files exist at specified locations, see 'overwrite'")
  idx <- sort(rep(seq_len(nsplits), length.out=nrow(y)))
  for (i in seq_len(nsplits)) {
    saveRDS(y[idx == i,], file=files[i])
  }
}

#' Helper function for distributing Swish on a subset of data
#'
#' This function is called by the \code{Snakefile} that is generated
#' by \code{\link{splitSwish}}. See alevin example in the vignette.
#' As such, it doesn't need to be run by users in an interactive
#' R session.
#' 
#' Note that the default for length correction is FALSE, as
#' opposed to the default in \code{\link{scaleInfReps}} which
#' is TRUE. The default for \code{numReps} here is 20.
#' 
#' @param infile path to an RDS file of a SummarizedExperiment
#' @param outfile a CSV file to write out
#' @param numReps how many inferential replicates to generate
#' @param lengthCorrect logical, see \code{\link{scaleInfReps}},
#' and Swish vignette. As this function is primarily for alevin,
#' the default is \code{FALSE}
#' @param overwrite logical, whether \code{outfile}
#' should overwrite an existing file
#' @param ... arguments passed to \code{\link{swish}}
#'
#' @return nothing, files are written out
#'
#' @export
miniSwish <- function(infile, outfile, numReps=20,
                      lengthCorrect=FALSE, overwrite=FALSE, ...) {
  stopifnot(all(is(c(infile, outfile), "character")))
  stopifnot(file.exists(infile))
  if (file.exists(outfile) & !overwrite)
    stop("outfile already exists at specified location, see 'overwrite'")
  y <- readRDS(infile)
  stopifnot(!is.null(rownames(y)))
  stopifnot(all(mcols(y)$keep))
  y <- makeInfReps(y, numReps=numReps)
  if (is.null(colData(y)$sizeFactor))
    stop("miniSwish requires pre-estimated sizeFactors stored in colData(...)")
  y <- scaleInfReps(y, lengthCorrect=lengthCorrect, sfFun=colData(y)$sizeFactor)
  out <- swish(y=y, returnNulls=TRUE, ...)
  mat <- cbind(out$stat, out$log2FC, out$nulls)
  rownames(mat) <- rownames(y)
  write.table(mat, file=outfile, col.names=FALSE, sep=",")
}

#' Read statistics and nulls from CSV file
#'
#' After running \code{\link{splitSwish}} and the associated
#' \code{Snakefile}, this function can be used to gather and
#' add the results to the original object. See the alevin
#' section of the vignette for an example.
#'
#' @param y a SummarizedExperiment (if NULL, function will
#' output a data.frame)
#' @param infile character, path to the \code{summary.csv} file
#' @param estPi0 logical, see \code{\link{swish}}
#'
#' @return the SummarizedExperiment with metadata columns added,
#' or if \code{y} is NULL, a data.frame of compiled results
#'
#' @export
addStatsFromCSV <- function(y=NULL, infile, estPi0=FALSE) {
  res <- read.table(infile, header=FALSE, row.names=1, sep=",")
  stat <- res[,1]
  log2FC <- res[,2]
  nulls <- as.matrix(res[,-c(1:2)])
  pvalue <- qvalue::empPvals(abs(stat), abs(nulls))
  pi0 <- if (estPi0) NULL else 1
  q.res <- qvalue::qvalue(pvalue, pi0=pi0)
  locfdr <- q.res$lfdr
  qvalue <- q.res$qvalues
  df <- data.frame(stat, log2FC, pvalue, locfdr, qvalue, row.names=rownames(res))
  if (is.null(y)) {
    return(df)
  } else {
    stopifnot(all(rownames(y) == rownames(df)))
    return(postprocess(y, df))
  }
}
