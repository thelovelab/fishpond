#' deswish: DESeq2-apeglm With Inferential Samples Helps
#'
#' The DESeq2-apeglm With Inferential Samples implementation supposes
#' a hierarchical distribution of log2 fold changes.
#' The final posterior standard deviation is calculated by
#' adding the posterior variance from modeling biological replicates
#' computed by \code{apeglm}, and the observed variance on the posterior mode
#' over inferential replicates. This function requires the DESeq2 and
#' apeglm packages to be installed and will print an error if they are
#' not found.
#' 
#' @param y a SummarizedExperiment containing the inferential
#' replicate matrices, as output by \code{tximeta}, and
#' then with \code{labelKeep} applied. One does not need to
#' run \code{scaleInfReps} as scaling is done internally via
#' DESeq2.
#' @param x the design matrix
#' @param coef the coefficient to test (see \code{lfcShrink})
#'
#' @return a SummarizedExperiment with metadata columns added:
#' the log2 fold change and posterior SD using inferential replicates,
#' and the original log2 fold change (apeglm) and its posterior SD
#'
#' @references
#'
#' The \code{DESeq} and \code{lfcShrink} function in the \code{DESeq2} package:
#'
#' Zhu, Ibrahim, Love "Heavy-tailed prior distributions for sequence count data:
#' removing the noise and preserving large differences" Bioinformatics (2018).
#' 
#' Love, Huber, Anders "Moderated estimation of fold change and dispersion
#' for RNA-seq data with DESeq2" Genome Biology (2014).
#'
#' @examples
#'
#' y <- makeSimSwishData()
#' y <- labelKeep(y)
#' y <- deswish(y, ~condition, "condition_2_vs_1")
#' 
#' @export
deswish <- function(y, x, coef) {
  if (is.null(mcols(y)$keep)) stop("first run labelKeep before deswish")
  ys <- y[mcols(y)$keep,]

  ys.cts <- ys
  assays(ys.cts) <- assays(ys.cts)[c("counts","length")]

  if (!requireNamespace("DESeq2", quietly=TRUE)) {
    stop("deswish requires the DESeq2 package")
  }
  if (!requireNamespace("apeglm", quietly=TRUE)) {
    stop("deswish requires the apeglm package")
  }
  dds <- DESeq2::DESeqDataSet(ys.cts, design=x)
  dds <- DESeq2::DESeq(dds, minReplicatesForReplace=Inf, quiet=TRUE)
  res <- DESeq2::lfcShrink(dds, coef=coef, type="apeglm", quiet=TRUE)

  nreps <- length(grep("infRep", assayNames(ys)))
  lfcs <- matrix(nrow=nrow(ys), ncol=nreps)
  for (i in seq_len(nreps)) {
    cat(i,"")
    rep.cts <- round(assays(ys)[[paste0("infRep",i)]])
    mode(rep.cts) <- "integer"
    DESeq2::counts(dds) <- rep.cts
    rep.res <- DESeq2::lfcShrink(dds, coef=coef, type="apeglm",
                                 apeMethod="nbinomC", quiet=TRUE)
    lfcs[,i] <- rep.res$log2FoldChange
  }
  cat("\n")

  infVarLFC <- matrixStats::rowVars(lfcs)
  log2FC <- rowMeans(lfcs)
  # suppose a hierarchical dist'n of LFC:
  # add the posterior variance from modeling biological replicates
  # and the variance on posterior mode over inferential replicates
  lfcSD <- sqrt(res$lfcSE^2 + infVarLFC)
  df <- data.frame(log2FC, lfcSD,
                   origLFC=res$log2FoldChange,
                   origSD=res$lfcSE)

  y <- postprocess(y, df)
  y
}
