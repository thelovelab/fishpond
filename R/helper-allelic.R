#' Import allelic counts as a SummarizedExperiment
#'
#' Read in Salmon quantification of allelic counts from a
#' diploid transcriptome. Assumes that diploid transcripts
#' are marked with the following suffix: an underscore and
#' a consistent symbol for each of the two alleles,
#' e.g. \code{ENST123_M} and \code{ENST123_P},
#' or \code{ENST123_alt} and \code{ENST123_ref}, etc.
#' \code{importAllelicCounts} requires the tximeta package.
#' Further information in Details below.
#'
#' \strong{Requirements} - There must be exactly two alleles for each
#' transcript, and the \code{--keep-duplicates} option should be used
#' in Salmon indexing to avoid removal of transcripts with identical
#' sequence. The output object has half the number of transcripts,
#' with the two alleles either stored in a \code{"wide"} object, or as
#' re-named \code{"assays"}. Note carefully that the symbol provided
#' to \code{a1} is used as the effect allele, and \code{a2} is used as
#' the non-effect allele (see the \code{format} argument description
#' and Value description below).
#'
#' \strong{tx2gene} - The two columns should include the \code{a1} and
#' \code{a2} suffix for the transcripts and genes/groups, or those
#' will be added internally, if it is detected that the first
#' transcript does not have these suffices. For example if
#' \code{_alt} or \code{_ref}, or \code{_M} or \code{_P} (as indicated
#' by the \code{a1} and \code{a2} arguments) are not present in the
#' table, the table rows will be duplicated with those suffices added
#' on behalf of the user. If \code{tx2gene} is not provided, the
#' output object will be transcript-level. Do not attempt to set the
#' \code{txOut} argument, it will conflict with internal calls to
#' downstream functions. If the a1/a2 suffices are not at the end of
#' the transcript name in the quantification files,
#' e.g. \code{ENST123_M|<metadata>}, then \code{ignoreAfterBar=TRUE}
#' can be used to match regardless of the string following \code{|} in
#' the quantification files.
#'
#' \code{skipMeta=TRUE} is used, as it is assumed the diploid
#' transcriptome does not match any reference transcript
#' collection. This may change in future iterations of the function,
#' depending on developments in upstream annotations and software.
#'
#' If \code{tx2gene} is a GRanges object, the rowRanges of the output
#' will be the reduced ranges of the grouped input ranges, with
#' \code{tx_id} collapsed into a CharacterList. Other metadata columns
#' are not manipulated, just the metadata for the first range is
#' returned.
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
#' grouping). Alternatively, this can be a GRanges object with
#' required columns \code{tx_id}, and \code{group_id}
#' (see \code{makeTx2Tss}). For more information on this argument,
#' see Details.
#' @param ... any arguments to pass to tximeta
#' 
#' @return a SummarizedExperiment, with allele counts (and other data)
#' combined into a wide matrix \code{[a2 | a1]}, or as assays (a1, then a2).
#' The original strings associated with a1 and a2 are stored in the
#' metadata of the object, in the \code{alleles} list element.
#' Note the reference level of \code{se$allele} will be \code{"a2"}, 
#' such that comparisons by default will be a1 vs a2
#' (effect vs non-effect). 
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

  # output transcripts if no tx2gene provided
  txOut <- is.null(tx2gene) 
  # if tx2gene is ranges, some extra steps to modify output rowRanges
  t2gRanges <- is(tx2gene, "GRanges")
  if (!txOut) {
    # if tx2gene is ranges, then pull out the table for collapsing
    if (t2gRanges) {
      cols <- c("tx_id","group_id")
      stopifnot(all(cols %in% names(mcols(tx2gene))))
      # swap around variable names to run the data.frame-based code
      txps <- tx2gene
      tx2gene <- mcols(txps)[,cols]
    }
    stopifnot(ncol(tx2gene) == 2)
    # see if tx2gene already has tagged the txps and genes
    diploid <- FALSE # save whether tx2gene is diploid
    if (grepl(a1match, tx2gene[1,1]) | grepl(a2match, tx2gene[1,1])) {
      # ensure same number of a1 and a2 alleles in txps and genes
      sum1txp <- sum(grepl(a1match, tx2gene[,1]))
      sum2txp <- sum(grepl(a2match, tx2gene[,1]))
      sum1gene <- sum(grepl(a1match, tx2gene[,2]))
      sum2gene <- sum(grepl(a2match, tx2gene[,2]))
      stopifnot(sum1txp > 0 & sum2txp > 0 & sum1gene > 0 & sum2gene > 0)
      stopifnot(sum1txp == sum2txp)
      stopifnot(sum1gene == sum2gene)
      diploid <- TRUE
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
    out <- SummarizedExperiment(assays=assays_wide,
                                colData=coldata_wide)
    metadata(out) <- c(metadata(se), list(alleles=c(a1=a1, a2=a2)))
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
    out <- se_a1
  }
  # extra steps to assign rowRanges to output
  if (t2gRanges) {
    if (diploid) {
      txps <- txps[grepl(a1match, names(txps))]
      names(txps) <- sub(a1match,"",names(txps))
      mcols(txps)$tx_id <- sub(a1match,"",mcols(txps)$tx_id)
      mcols(txps)$group_id <- sub(a1match,"",mcols(txps)$group_id)
    }
    stopifnot(all(rownames(out) %in% mcols(txps)$group_id))
    tx_list <- CharacterList(split(mcols(txps)$tx_id, mcols(txps)$group_id))
    # detect if the TSS have been grouped by `maxgap`
    # (the group_id will not be equal to "gene-tss")
    tss_grouped <- mcols(txps)$group_id[1] != paste0(mcols(txps)$gene_id[1],
                                                     "-", mcols(txps)$tss[1])
    # want to save individual TSS information in a list
    if (tss_grouped) {
      tss_list <- IntegerList(split(mcols(txps)$tss, mcols(txps)$group_id))
    }
    tx_starts <- sapply(split(start(txps), mcols(txps)$group_id), min)
    tx_ends <- sapply(split(end(txps), mcols(txps)$group_id), max)
    txps <- txps[!duplicated(mcols(txps)$group_id)]
    names(txps) <- mcols(txps)$group_id
    start(txps) <- tx_starts[names(txps)]
    end(txps) <- tx_ends[names(txps)]
    txps <- txps[rownames(out)]
    mcols(txps)$tx_id <- tx_list[rownames(out)]
    if (tss_grouped) {
      mcols(txps)$tss <- tss_list[rownames(out)]
    }
    SummarizedExperiment::rowRanges(out) <- txps
  }
  out
}

#' Make a GRanges linking transcripts to TSS within gene
#'
#' This helper function takes either a TxDb/EnsDb or
#' GRanges object as input and outputs a GRanges object
#' where transcripts are aggregated to the gene + TSS
#' (transcription start site). For nearby TSS that should
#' be grouped together, see \code{maxgap}.
#'
#' @param x either TxDb/EnsDb or GRanges object. The GRanges
#' object should have metadata columns \code{tx_id} and
#' \code{gene_id}
#' @param maxgap integer, number of basepairs to use determining
#' whether to combine nearby TSS
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
makeTx2Tss <- function(x, maxgap=0) {
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
    # if any missing gene IDs, just propagate the txp IDs
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
  tss_ranges <- GenomicRanges::resize(txps, width=1)
  tss <- GenomicRanges::start(tss_ranges)
  mcols(txps)$tss <- tss
  if (maxgap == 0) {
    # group ID is just gene ID + TSS
    mcols(txps)$group_id <- paste0(mcols(txps)$gene_id, "-", mcols(txps)$tss)
  } else {
    # in this case, we will number the TSS groups with 1,2,...
    # nearby TSS are grouped together according to maxgap
    tss_ranges <- sort(tss_ranges)
    txps <- txps[ names(tss_ranges) ]
    txps <- txps[ order(mcols(txps)$gene_id) ]
    sp_tss <- split(mcols(txps)$tss, mcols(txps)$gene_id)
    tss_groups <- lapply(sp_tss, joinNearbyTss, maxgap)
    mcols(txps)$group_id <- paste0(mcols(txps)$gene_id, "-", unlist(tss_groups))
  }
  txps
}

joinNearbyTss <- function(tss, maxgap) {
  delta <- diff(tss, lag=1)
  new_tss_loc <- as.integer(delta > maxgap)
  tss_loc <- c(1,cumsum(new_tss_loc) + 1)
  tss_loc
}

#' Plot allelic counts in a gene context using Gviz
#'
#' Plot allelic data (allelic proportions, isoform propostions)
#' in a gene context leveraging the Gviz package. See the allelic
#' vignette for example usage. TPM and count filters are used by
#' default to clean up the plot of features with minimal signal;
#' note that the isoform proportion displayed at the bottom of the
#' plot is among the features that pass the filtering steps.
#' If the function is not responding, it is likely due to issues
#' connecting to UCSC servers to draw the ideogram, in this case
#' set \code{ideogram=FALSE}.
#' 
#' @param y a SummarizedExperiment (see \code{swish})
#' @param gene the name of the gene of interest, requires
#' a column \code{gene_id} in the metadata columns of the
#' rowRanges of y
#' @param db either a TxDb or EnsDb object to use for the gene model
#' @param region GRanges, the region to be displayed in the Gviz plot.
#' if not specified, will be set according to the gene plus 20%
#' of the total gene extent on either side
#' @param symbol alternative to \code{gene}, to specify
#' the gene of interest according to a column \code{symbol}
#' in the metadata columns of the rowRanges of y
#' @param genome UCSC genome code (e.g. \code{"hg38"},
#' if not specified it will use the \code{GenomeInfoDb::genome()}
#' of the rowRanges of \code{y}
#' @param tpmFilter minimum TPM value (mean over samples) to keep a feature
#' @param countFilter minimum count value (mean over samples) to keep a feature
#' @param pc pseudocount to avoid dividing by zero in allelic proportion calculation
#' @param transcriptAnnotation argument passed to Gviz::GeneRegionTrack
#' (\code{"symbol"}, \code{"gene"}, \code{"transcript"}, etc.)
#' @param labels list, labels for a2 (non-effect) and a1 (effect) alleles
#' @param qvalue logical, whether to inclue qvalue track
#' @param log2FC logical, whether to include log2FC track
#' @param ideogram logical, whether to include ideogram track
#' @param cov character specifying a factor or integer variable to use
#' to facet the allelic proportion plots, should be a column in
#' \code{colData(y)}
#' @param covFacetIsoform logical, if \code{cov} is provided,
#' should it also be used to facet the isoform proportion track,
#' in addition to the allelic proportion track
#' @param allelicCol the colors of the points and lines for allelic proportion
#' @param isoformCol the colors of the points and lines for isoform proportion
#' @param statCol the color of the lollipops for q-value and log2FC
#' @param gridCol the color of the grid in the data tracks
#' @param baselineCol the color of the horizontal baseline for q-value and lo2gFC
#' @param titleCol font color of the side titles (track labels)
#' @param titleAxisCol axis color of the side titles (track labels)
#' @param titleBgCol background color of the side titles (track labels)
#' @param geneBorderCol the color of the borders and font in gene region track
#' @param geneFillCol the color of the fill in the gene region track
#' @param genomeAxisCol line color of the genome axis track
#' @param innerFontCol font color of genome axis track, ideogram, and
#' allelic proportion legend
#' @param ... additional arguments passed to \code{Gviz::plotTracks()}
#'
#' @return nothing, a plot is displayed
#'
#' @export
plotAllelicGene <- function(y, gene, db,
                            region=NULL, symbol=NULL, genome=NULL,
                            tpmFilter=1, countFilter=10, pc=1,
                            transcriptAnnotation="symbol",
                            labels=list(a2="a2",a1="a1"),
                            qvalue=TRUE,
                            log2FC=TRUE,
                            ideogram=FALSE,
                            cov=NULL,
                            covFacetIsoform=FALSE,
                            allelicCol=c("dodgerblue","goldenrod1"),
                            isoformCol="firebrick",
                            statCol="black",
                            gridCol="grey80",
                            baselineCol="black",
                            titleCol="black",
                            titleAxisCol="black",
                            titleBgCol="white",
                            geneBorderCol="darkblue",
                            geneFillCol="darkblue",
                            genomeAxisCol="black",
                            innerFontCol="black", ...) {
  if (!requireNamespace("Gviz", quietly=TRUE)) {
    stop("plotAllelicGene() requires 'Gviz' Bioconductor package")
  }
  if (!requireNamespace("GenomeInfoDb", quietly=TRUE)) {
    stop("plotAllelicGene() requires 'GenomeInfoDb' Bioconductor package")
  }
  stopifnot("gene_id" %in% names(mcols(y)))
  if (!is.null(symbol)) {
    stopifnot("symbol" %in% names(mcols(y)))
    gene <- mcols(y)$gene_id[match(symbol, mcols(y)$symbol)]
  }
  stopifnot(is.character(gene))
  stopifnot(is(db, "TxDb") | is(db, "EnsDb"))
  stopifnot(length(allelicCol) == 2)
  stopifnot(all(c("a2","a1") %in% names(labels)))
  if (!is.null(cov)) {
    stopifnot(cov %in% names(colData(y)))
  } else {
    stopifnot(!covFacetIsoform)
  }
  # pull out the ranges of y to build the Gviz plot
  gr <- rowRanges(y)
  # standard chr's and UCSC for compatibility w Gviz
  gr <- GenomeInfoDb::keepStandardChromosomes(
                        gr, pruning.mode="coarse")
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  # gene must be in the metadata
  stopifnot(gene %in% gr$gene_id)
  # pull out the ranges for the gene of interest
  gr <- gr[gr$gene_id == gene]
  stopifnot(!any(duplicated(names(gr))))
  # use TPM and count filtering to remove lowly expressed features
  if (!is.null(tpmFilter) | !is.null(countFilter)) {
    keep <- rep(TRUE, length(gr))
    if (!is.null(tpmFilter)) {
      keep <- keep & rowMeans(assay(y[names(gr),,drop=FALSE],
                                    "abundance")) >= tpmFilter
    }
    if (!is.null(countFilter)) {
      keep <- keep & rowMeans(assay(y[names(gr),,drop=FALSE],
                                    "counts")) >= countFilter
    }
    stopifnot(sum(keep) > 0)
    gr <- gr[keep,]
  }
  if (!"qvalue" %in% names(mcols(gr))) {
    stop("expecting qvalue and log2FC, first run swish()")
  }
  regionProvided <- TRUE
  if (is.null(region)) {
    regionProvided <- FALSE
    # the region to be displayed
    region <- range(gr)
    total_width <- width(region)
  }  
  # plot data at the TSS
  gr <- flank(gr, width=1)
  strand(gr) <- "*"
  # so data plots work, need to be sorted
  gr <- sort(gr)
  gr$minusLogQ <- -log10(gr$qvalue)
  # define upper bounds for the q-value and LFC
  qUpper <- 1.2 * max(gr$minusLogQ)
  lfcUpper <- 1.2 * max(abs(gr$log2FC))
  chr <- as.character(seqnames(region))
  ##########################################
  ## store allelic data in GRanges object ##
  ##########################################
  gr_allelic <- gr
  mcols(gr_allelic) <- NULL
  stopifnot("counts" %in% assayNames(y))
  allelic_counts <- assay(y[names(gr),,drop=FALSE], "counts")
  stopifnot(rownames(allelic_counts) == names(gr_allelic))
  n <- ncol(allelic_counts)/2
  stopifnot(all(y$allele == rep(c("a2","a1"),each=n)))
  total_counts <- (
    allelic_counts[,1:n,drop=FALSE] +
    allelic_counts[,(n+1):(2*n),drop=FALSE])
  allelic_prop <- (allelic_counts+pc) /
    cbind(total_counts+2*pc, total_counts+2*pc)
  rowMeans(allelic_prop, na.rm=TRUE)
  mcols(gr_allelic) <- allelic_prop
  ##########################################
  ## store isoform data in GRanges object ##
  ##########################################
  gr_isoform <- gr
  mcols(gr_isoform) <- NULL
  stopifnot("abundance" %in% assayNames(y))
  allelic_tpm <- assay(y[names(gr),,drop=FALSE], "abundance")
  total_tpm <- (
    allelic_tpm[,1:n,drop=FALSE] +
    allelic_tpm[,(n+1):(2*n),drop=FALSE])
  gene_tpm <- colSums(total_tpm)
  isoform_prop <- t(t(total_tpm) / gene_tpm)
  isoUpper <- 1.1 * max(isoform_prop)
  mcols(gr_isoform) <- isoform_prop
  ####################
  ## ideogram track ##
  ####################
  if (ideogram) {
    if (is.null(genome)) {
      genome <- GenomeInfoDb::genome(gr)[1]
    }
    ideo_track <- Gviz::IdeogramTrack(genome=genome, chromosome=chr,
                                      fontcolor=innerFontCol)
  } else {
    ideo_track <- NULL
  }
  ##################
  ## genome track ##
  ##################
  genome_track <- Gviz::GenomeAxisTrack(col=genomeAxisCol,
                                        fontcolor=innerFontCol)
  # for ensembldb EnsDb
  if (is(db, "EnsDb")) {
    if (!requireNamespace("ensembldb", quietly=TRUE)) {
      stop("db=EnsDb requires 'ensembldb' Bioconductor package")
    }
    GenomeInfoDb::seqlevelsStyle(db) <- "UCSC"
    # warning is also showing up on ensembldb vignette...
    suppressWarnings({
      gene_region <- ensembldb::getGeneRegionTrackForGviz(
        db,
        chromosome=chr,
        start=start(region)+10,
        end=end(region)-10)
    })
    gene_track <- Gviz::GeneRegionTrack(
      range=gene_region,
      name="gene model",
      transcriptAnnotation=transcriptAnnotation,
      col=geneBorderCol, col.line=geneBorderCol,
      fontcolor.group=geneBorderCol, fill=geneFillCol)
  } else {
    message("using TxDb, assuming UCSC seqlevelsStyle")
    gene_track <- Gviz::GeneRegionTrack(
      range=db,
      chromosome=chr,
      start=start(region)+10,
      end=end(region)-10,
      name="gene model",
      transcriptAnnotation=transcriptAnnotation,
      col=geneBorderCol, col.line=geneBorderCol,
      fontcolor.group=geneBorderCol, fill=geneFillCol)
  }
  ##################################
  ## put together the data tracks ##
  ##################################
  if (qvalue) {
    qvalue_track <- Gviz::DataTrack(
      grSelect(gr, "minusLogQ"),
      type=c("p","h","g"), name="-log10 qvalue",
      col=statCol, col.grid=gridCol, cex=1.5, lwd=2,
      ylim=c(0,qUpper), baseline=0)
  } else {
    qvalue_track <- NULL
  }
  if (log2FC) {
    lfc_track <- Gviz::DataTrack(
      grSelect(gr, "log2FC"),
      type=c("p","h","g"), name="log2FC",
      col=statCol, col.grid=gridCol, cex=1.5, lwd=2,
      baseline=0,
      ylim=c(-lfcUpper,lfcUpper))
  } else {
    lfc_track <- NULL
  }
  allele <- y$allele
  # optionally relabelling of alleles
  if (!(labels$a2 == "a2" & labels$a1 == "a1")) {
    levels(allele) <- c(labels$a2, labels$a1)
  }
  # case where we are not faceting across a covariate
  if (is.null(cov)) {
    allele_track <- list(Gviz::DataTrack(
      gr_allelic, type=c("p","a","g"), name="allelic prop.",
      groups=allele, col=allelicCol, col.grid=gridCol, lwd=2,
      baseline=0.5, fontcolor.legend=innerFontCol))
  } else {
    # faceting case
    covariate <- colData(y)[[cov]]
    if (!is(covariate, "factor")) {
      covariate <- factor(covariate)
    }
    allele_track <- lapply(levels(covariate), function(i) {
      gr_allelic_sub <- gr_allelic
      idx <- covariate == i
      mcols(gr_allelic_sub) <- mcols(gr_allelic_sub)[,idx]
      Gviz::DataTrack(
        gr_allelic_sub, type=c("p","a","g"), name=i,
        groups=allele[idx], col=allelicCol, col.grid=gridCol, lwd=2,
        baseline=0.5, fontcolor.legend=innerFontCol)
    })
  }
  # case where we are not faceting across a covariate
  if (!covFacetIsoform) {
    isoform_track <- list(Gviz::DataTrack(
      gr_isoform, type=c("p","a","g"), name="isoform prop.",
      col=isoformCol, col.grid=gridCol, lwd=2,
      baseline=0, ylim=c(0, isoUpper)))
  } else {
    # faceting case
    isoform_track <- lapply(levels(covariate), function(i) {
      gr_isoform_sub <- gr_isoform
      # subset to one half of the covariate vector
      idx <- covariate[1:n] == i
      mcols(gr_isoform_sub) <- mcols(gr_isoform_sub)[,idx]
      Gviz::DataTrack(
        gr_isoform_sub, type=c("p","a","g"), name=i,
        col=isoformCol, col.grid=gridCol, lwd=2,
        baseline=0, ylim=c(0, isoUpper))
    })
  }
  if (!regionProvided) {
    eps <- round(.2 * total_width)
    gvizFrom <- start(region) - eps
    gvizTo <- end(region) + eps
  } else {
    gvizFrom <- start(region)
    gvizTo <- end(region)
  }
  ##############################
  ## finally, plot the tracks ##
  ##############################
  tracks <- c(list(ideo_track, genome_track, gene_track,
                   qvalue_track, lfc_track),
              allele_track, # already list
              isoform_track) # already list
  tracks <- tracks[!sapply(tracks, is.null)]
  Gviz::plotTracks(
    tracks, from=gvizFrom, to=gvizTo,
    col.title=titleCol,
    col.axis=titleAxisCol,
    background.title=titleBgCol,
    col.baseline=baselineCol, ...)
}

grSelect <- function(gr, col) {
  out <- gr
  mcols(out) <- NULL
  mcols(out)[col] <- mcols(gr)[col]
  out
}

#' Plot allelic ratio heatmap 
#'
#' Plot allelic ratio heatmap over features and samples
#' using the pheatmap package. The a1/(a2 + a1) ratio
#' is displayed.
#' 
#' @param y a SummarizedExperiment (see \code{swish})
#' @param idx a numeric or logical vector of which features
#' to plot
#' @param breaks breaks passed along to pheatmap
#' @param cluster_cols logical, passed to pheatmap
#' @param main title of the plot
#' @param stripAfterChar for the column names, if specified
#' will strip allelic identifiers after this character,
#' default is hyphen. set to NULL to avoid this action
#' @param ... other arguments passed to pheatmap
#'
#' @return nothing, a plot is displayed
#'
#' @export
plotAllelicHeatmap <- function(y, idx,
                               breaks=NULL,
                               cluster_cols=FALSE,
                               main="Allelic ratio",
                               stripAfterChar="-",
                               ...) {
  if (!requireNamespace("GenomeInfoDb", quietly=TRUE)) {
    stop("plotAllelicHeatmap() requires 'pheatmap' CRAN package")
  }
  stopifnot(ncol(y) > 1)
  if ("mean" %in% assayNames(y)) {
    cts <- assay(y, "mean")[idx,,drop=FALSE]
    message("using posterior mean for calculating ratio")
  } else {
    cts <- assay(y, "counts")[idx,,drop=FALSE]
    message("using counts, for posterior mean, run computeInfRV")
  }
  stopifnot("allele" %in% names(colData(y)))
  cts_a2 <- cts[,y$allele == "a2",drop=FALSE]
  cts_a1 <- cts[,y$allele == "a1",drop=FALSE]
  tot <- cts_a2 + cts_a1
  ratio <- cts_a1 / tot
  if (!is.null(breaks)) {
    delta <- max(abs(ratio - 0.5))
    breaks <- seq(from=0.5 - delta, 0.5 + delta, length.out=101)
  }
  if (!is.null(stripAfterChar)) {
    colnames(ratio) <- sub(paste0(stripAfterChar,".*"),
                           "",colnames(ratio))
  }
  pheatmap::pheatmap(ratio, breaks=breaks,
                     cluster_cols=cluster_cols,
                     main=main, ...)
}
