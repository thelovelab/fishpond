#' deswish: DESeq2 With Inferential Samples Helps
#'
#' The DESeq2 With Inferential Samples implementation supposes
#' a hierarchical distribution of log2 fold changes.
#' The final posterior standard deviation is calculated by
#' adding the posterior variance from modeling biological replicates
#' (apeglm) and the observed variance on the posterior mode
#' over inferential replicates.
#' 
#' @param y a SummarizedExperiment containing the inferential replicate matrices,
#' as output by \code{tximeta}
#' @param x the design matrix
#'
#' @return a SummarizedExperiment with metadata columns added
#'
#' @references
#'
#' The \code{DESeq} and \code{lfcShrink} function in the \code{DESeq2} package.
#'
#' Zhu, Ibrahim, Love "Heavy-tailed prior distributions for sequence count data:
#' removing the noise and preserving large differences" Bioinformatics (2018).
#' 
#' Love, Huber, Anders "Moderated estimation of fold change and dispersion
#' for RNA-seq data with DESeq2" Genome Biology (2014).
#' 
#' @export
deswish <- function(y, x, coef) {
  ys <- y[mcols(y)$keep,]

  ys.cts <- ys
  assays(ys.cts) <- assays(ys.cts)[c("counts","length")]
  dds <- DESeq2::DESeqDataSet(ys.cts, design=x)
  dds <- DESeq2::DESeq(dds, minReplicatesForReplace=Inf, quiet=TRUE)
  res <- DESeq2::lfcShrink(dds, coef=coef, type="apeglm", quiet=TRUE)

  nreps <- length(grep("infRep", assayNames(ys)))
  lfcs <- matrix(nrow=nrow(ys), ncol=nreps)
  for (i in seq_len(nreps)) {
    cat(i,"")
    rep.cts <- round(assays(ys)[[paste0("infRep",i)]])
    mode(rep.cts) <- "integer"
    counts(dds) <- rep.cts
    rep.res <- DESeq2::lfcShrink(dds, coef=coef, type="apeglm",
                                 apeMethod="nbinomC", quiet=TRUE)
    lfcs[,i] <- rep.res$log2FoldChange
  }

  infVarLFC <- rowVars(lfcs)
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
