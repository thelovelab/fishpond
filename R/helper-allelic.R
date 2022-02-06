#' Import allelic counts as a SummarizedExperiment
#'
#' Read in Salmon quantification of allelic counts from a
#' diploid transcriptome. Assumes that diploid transcripts
#' are marked with the following suffix: an underscore and
#' a consistent symbol for each of the two alleles,
#' e.g. \code{ENST123_M} and \code{ENST123_P},
#' or \code{ENST123_alt} and \code{ENST123_ref}.
#' There must be exactly two alleles for each transcript, 
#' and the \code{--keep-duplicates} option should be used in
#' Salmon indexing to avoid removing transcripts with identical sequence.
#' The output object has half the number of transcripts,
#' with the two alleles either stored in a \code{"wide"} object,
#' or as re-named \code{"assays"}. Note carefully that the symbol
#' provided to \code{a1} is used as the effect allele,
#' and \code{a2} is used as the non-effect allele
#' (see the \code{format} argument description and Value
#' description below).
#'
#' Requires the tximeta package.
#' \code{skipMeta=TRUE} is used, as it is assumed
#' the diploid transcriptome does not match any reference
#' transcript collection. This may change in future iterations
#' of the function, depending on developments in upstream
#' software.
#'
#' @param coldata a data.frame as used in \code{tximeta}
#' @param a1 the symbol for the effect allele
#' @param a2 the symbol for the non-effect allele
#' @param format either \code{"wide"} or \code{"assays"} for whether
#' to combine the allelic counts as columns (wide) or put the allelic
#' count information in different assay slots (assays).
#' For wide output, the data for the non-effect allele (a2) comes first,
#' then the effect allele (a1), e.g. \code{[a2 | a1]}. The \code{ref} level
#' of the factor variable \code{se$allele} will be \code{"a2"}
#' (so by default comparisons will be: a1 vs a2).
#' For assays output, all of the original matrices are renamed with a prefix,
#' either \code{a1-} or \code{a2-}.
#' @param tx2gene optional, a data.frame with first column indicating
#' transcripts, second column indicating genes (or any other transcript
#' grouping). Either this should include the \code{a1} and \code{a2}
#' suffix for the transcripts and genes, or those will be added internally,
#' if it is detected that the first transcript does not have these suffices.
#' For example if \code{_alt} or \code{_ref}, or \code{_M} or \code{_P}
#' (as indicated by the \code{a1} and \code{a2} arguments) are not present
#' in the table, the table rows will be duplicated with those suffices
#' added on behalf of the user.
#' If not provided, the output object will be transcript-level.
#' Note: do not attempt to set the \code{txOut} argument, it will
#' conflict with internal calls to downstream functions.
#' Note: if the a1/a2 suffices are not at the end of the transcript name
#' in the quantification files, e.g. \code{ENST123_M|<metadata>},
#' then \code{ignoreAfterBar=TRUE} can be used to match regardless of
#' the string following \code{|} in the quantification files.
#' @param ... any arguments to pass to tximeta
#'
#' @return a SummarizedExperiment, with allele counts (and other data)
#' combined into a wide matrix \code{[a2 | a1]}, or as assays (a1, then a2).
#' The original strings associated with a1 and a2 are stored in the
#' metadata of the object, in the \code{alleles} list element.
#' Note the \code{ref} level of \code{se$allele} will be \code{"a2"}, 
#' such that comparisons by default will be a1 vs a2 (effect vs non-effect).
#' 
#' @export
importAllelicCounts <- function(coldata, a1, a2,
                                format=c("wide","assays"),
                                tx2gene=NULL, ...) {
  format <- match.arg(format)
  if (!requireNamespace("tximeta", quietly=TRUE)) {
    stop("this function requires installing the Bioconductor package 'tximeta'")
  }

  a1match <- paste0("_",a1,"$")
  a2match <- paste0("_",a2,"$")

  txOut <- is.null(tx2gene) # output transcripts if no tx2gene provided  
  if (!txOut) {
    stopifnot(ncol(tx2gene) == 2)
    # see if tx2gene already has tagged the txps and genes
    if (grepl(a1match, tx2gene[1,1]) | grepl(a2match, tx2gene[1,1])) {
      # ensure same number of a1 and a2 alleles in txps and genes
      sum1txp <- sum(grepl(a1match, tx2gene[,1]))
      sum2txp <- sum(grepl(a2match, tx2gene[,1]))
      sum1gene <- sum(grepl(a1match, tx2gene[,2]))
      sum2gene <- sum(grepl(a2match, tx2gene[,2]))
      stopifnot(sum1txp > 0 & sum2txp > 0 & sum1gene > 0 & sum2gene > 0)
      stopifnot(sum1txp == sum2txp)
      stopifnot(sum1gene == sum2gene)
    } else {
      a2a1_vec <- rep(c(a2, a1), each=nrow(tx2gene))
      tx2gene <- data.frame(
        tx=paste0(rep(tx2gene[,1], 2), "_", a2a1_vec),
        gene=paste0(rep(tx2gene[,2], 2), "_", a2a1_vec)
      )
    }
  }
  
  se <- tximeta::tximeta(coldata, skipMeta=TRUE, tx2gene=tx2gene, txOut=txOut, ...)

  # remove any characters after "|"
  rownames(se) <- sub("\\|.*", "", rownames(se))

  ntxp <- nrow(se)/2
  n <- ncol(se)

  # gather transcript names for a1 and a2 alleles
  txp_nms_a1 <- grep(a1match, rownames(se), value=TRUE)
  stopifnot(length(txp_nms_a1) == ntxp)
  txp_nms_a2 <- sub(paste0("_",a1),paste0("_",a2),txp_nms_a1)
  stopifnot(all(txp_nms_a2 %in% rownames(se)))
  stopifnot(length(txp_nms_a1) == length(txp_nms_a2))
  txp_nms <- sub(paste0("_",a1),"",txp_nms_a1)
  
  if (format == "wide") {
    coldata_wide <- data.frame(
      allele=factor(rep(c("a2","a1"), each=n), levels=c("a2","a1"))
    )
    # add any other covariate data
    for (v in setdiff(names(coldata), c("files","names"))) {
      coldata_wide[[v]] <- coldata[[v]]
    }
    rownames(coldata_wide) <- paste0(colnames(se), "-", coldata_wide$allele)

    assays_wide <- lapply(assays(se), function(a) {
      a_wide <- cbind(a[txp_nms_a2,], a[txp_nms_a1,])
      rownames(a_wide) <- txp_nms
      colnames(a_wide) <- rownames(coldata_wide)
      a_wide
    })
    # make a new SE
    wide <- SummarizedExperiment(assays=assays_wide,
                                 colData=coldata_wide)
    metadata(wide) <- c(metadata(se), list(alleles=c(a1=a1, a2=a2)))
    return(wide)
  } else if (format == "assays") {
    se_a1 <- se[txp_nms_a1,]
    se_a2 <- se[txp_nms_a2,]
    rownames(se_a1) <- txp_nms
    rownames(se_a2) <- txp_nms
    # rename the assays
    assayNames(se_a1) <- paste0("a1-", assayNames(se_a1))
    assayNames(se_a2) <- paste0("a2-", assayNames(se_a2))
    # add the a2 matrices to the a1 SE object
    assays(se_a1) <- c(assays(se_a1), assays(se_a2))
    metadata(se_a1) <- c(metadata(se_a1), list(alleles=c(a1=a1, a2=a2)))
    return(se_a1)
  }
}

#' Make a GRanges linking transcripts to TSS within gene
#'
#' This helper function takes either a TxDb/EnsDb or
#' GRanges object as input and outputs a GRanges object
#' where transcripts are aggregated to the gene + TSS
#' (transcription start site).
#'
#' @param x either TxDb/EnsDb or GRanges object. The GRanges
#' object should have metadata columns \code{tx_id} and
#' \code{gene_id}
#'
#' @return GRanges with columns \code{tx_id}, \code{tss}, and
#' \code{group_id}
#'
#' @examples
#' \dontrun{
#' library(EnsDb.Hsapiens.v86)
#' edb <- EnsDb.Hsapiens.v86
#' t2t <- makeTx2Tss(edb)
#' }
#' 
#' @export
makeTx2Tss <- function(x) {
  if (is(x, "GRanges")) {
    txps <- x
    if (!all(c("tx_id","gene_id") %in% names(mcols(txps)))) {
      stop("'tx_id' and 'gene_id' must be mcols of 'x'")
    }
  } else if (is(x, "TxDb")) {
    if (!requireNamespace("GenomicFeatures", quietly=TRUE)) {
      stop("'x' as TxDb requires GenomicFeatures")
    }
    if (!requireNamespace("AnnotationDbi", quietly=TRUE)) {
      stop("'x' as TxDb requires AnnotationDbi")
    }
    txps <- GenomicFeatures::transcripts(x)
    mcols(txps)$tx_id <- mcols(txps)$tx_name
    tx_id_stripped <- sub("\\..*", "", mcols(txps)$tx_id)
    mcols(txps)$gene_id <- AnnotationDbi::mapIds(x, tx_id_stripped, "GENEID", "TXNAME")
    missing.idx <- is.na(mcols(txps)$gene_id)
    mcols(txps)$gene_id[ missing.idx ] <- mcols(txps)$tx_id[ missing.idx ]
  } else if (is(x, "EnsDb")) {
    if (!requireNamespace("ensembldb", quietly=TRUE)) {
      stop("'x' as EnsDb requires ensembldb")
    }
    txps <- ensembldb::transcripts(x)
  } else {
    stop("'x' should be a GRanges or TxDb/EnsDb")
  }
  tss <- GenomicRanges::start(GenomicRanges::resize(txps, width=1))
  mcols(txps)$tss <- tss
  mcols(txps)$group_id <- paste0(mcols(txps)$gene_id, "-", mcols(txps)$tss)
  txps
}
