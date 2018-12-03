#' deswish: DESeq2 With Inferential Samples Helps
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
  mcols(y)$keep <- rowSums(assays(y)[["counts"]] >= 10) >= 6
  ys <- y[mcols(y)$keep,]

  ys.cts <- ys
  assays(ys.cts) <- assays(ys.cts)[c("counts","length")]
  dds <- DESeq2::DESeqDataSet(ys.cts, design=x)
  dds <- DESeq2::DESeq(dds, minReplicatesForReplace=Inf)
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
  stat <- rowMeans(lfcs)
  df <- data.frame(stat)

  y <- postprocess(y, df)
  y
}
