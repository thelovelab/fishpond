#' Create isoform proportions from scaled data
#'
#' Takes output of scaled (and optionally filtered) counts
#' and returns isoform proportions by dividing out the
#' total scaled count for the gene for each sample.
#' The operation is performed on the \code{counts} assay,
#' then creating a new assay called \code{isoProp},
#' and on all of the inferential replicates, turning them
#' from counts into isoform proportions. Any transcripts
#' (rows) from single isoform genes are removed, and the
#' transcripts will be re-ordered by gene ID.
#'
#' @param y a SummarizedExperiment
#' @param geneCol the name of the gene ID column in the
#' metadata columns for the rows of \code{y}
#' @param quiet display no messages
#'
#' @return a SummarizedExperiment, with single-isoform
#' transcripts removed, and transcripts now ordered by
#' gene
#' 
#' @export
isoformProportions <- function(y, geneCol="gene_id", quiet=FALSE) {
  if (!interactive()) {
    quiet <- TRUE
  }
  if (is.null(metadata(y)$infRepsScaled)) {
    stop("first run scaleInfReps()")
  }
  if (!is.null(metadata(y)$infRepsProportions)) {
    if (metadata(y)$infRepsProportions) stop("inferential replicates already proportions")
  }
  stopifnot(geneCol %in% names(mcols(y)))

  gene <- mcols(y)[[geneCol]]
  stopifnot(all(lengths(gene) == 1))
  mcols(y)$gene <- unlist(gene)
  gene.tbl <- table(mcols(y)$gene)
  # remove single isoform genes
  keep <- mcols(y)$gene %in% names(gene.tbl)[gene.tbl > 1]
  stopifnot(sum(keep) > 0)
  y <- y[keep,]
  y <- y[order(mcols(y)$gene),]
  assays(y, withDimnames=FALSE)[["isoProp"]] <- makeIsoProp(
    assays(y)[["counts"]],
    mcols(y)$gene
  )
  
  infRepIdx <- grep("infRep",assayNames(y))
  infRepError(infRepIdx)
  infReps <- assays(y)[infRepIdx]
  nreps <- length(infReps)
  for (k in seq_len(nreps)) {
    if (!quiet) svMisc::progress(k, max.value=nreps, init=(k==1), gui=FALSE)
    infReps[[k]] <- makeIsoProp(infReps[[k]],
                                mcols(y)$gene)
  }
  if (!quiet) message("")
  
  assays(y)[grep("infRep",assayNames(y))] <- infReps
  metadata(y)$infRepsProportions <- TRUE
  y
}

# not exported
makeIsoProp <- function(counts, gene) {
  totals <- rowsum(counts, gene)
  # if a sample has a gene total of 0, replace here with 1 to avoid division by 0
  totals[totals == 0] <- 1
  ngene <- length(unique(gene))
  gene.tbl <- table(gene)
  idx <- rep(seq_len(ngene), gene.tbl)
  big.totals <- totals[idx,]
  counts / big.totals
}
