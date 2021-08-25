#' Load in data from alevin-fry USA mode
#'
#' Enables easy loading of sparse data matrices provided by alevin-fry USA mode.
#' Alevin-fry - https://www.biorxiv.org/content/10.1101/2021.06.29.450377v1
#'
#'
#' @param fry.dir A path to the output directory returned by alevin-fry quant command. This directory
#' should contain a metainfo.json, and an `alevin` folder which contains quants_mat.mtx,
#' quants_mat_cols.txt and  quants_mat_rows.txt.
#' @param which_counts  A vector specifying which kinds of counts will be considered to build the final count matrix.
#  There are three options: U (Unspliced), S (Spliced) and A (Ambiguous).
#' For single-cell RNA-seq, `c('S', 'A')` is recommended (the default value); for single nucleus RNA-seq, `c('U', 'S', 'A')`
#' is recommended.
#' @param velocity A boolean (default: `FALSE`) specifying whether you want to return a `SingleCellExperiment with
#' that includes the values from the `"U"`-nspliced counts stored in an additional `"unspliced"` matrix in the assay
#' slot. When `TRUE`, the kinds specification provided in the `which_counts` parameter determines the types of counts
#' used for the normal, gene-level, qauntitation returned in the default `"counts"` assay.
#' @param verbose A boolean specifying if showing messages when running the function
#'
#' @details
#' This function consumes the result folder returned by running alevin-fry quant in unspliced, spliced, ambiguous (USA) quantification mode, and returns a `SingleCellExperiement` object
#' that contains a final count for each gene within each cell. In USA mode, alevin-fry quant returns a count matrix contains three types of count for each feature (gene)
#' within each sample (cell or nucleus), which represent the spliced mRNA count of the gene, the unspliced mRNA count of the gene, and the count of mRNA whose splicing status is ambiguous.
#' In this function, these three counts of a gene within a cell will either be summed or discarded to get the final count of the gene by specifying `which_counts`. The returned object will
#' contains a gene by cell count matrix, with rownames as the barcode of samples and colnames as the feature names.
#'
#' @return A `SingleCellExperiment` object contains a gene by cell count matrix. The row names are feature names, and the column names are cell barcodes.
#' When `velocity = TRUE`, a second assay matrix named `"unspliced"` will be included. This allows the returned object to house both quantitation matrices
#' required for velocity analysis using tools like scvelo/velociraptor.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom Matrix readMM t
#' @importFrom utils read.table
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' # Get path for minimal example avelin-fry output dir
#' testdat <- fishpond:::readExampleFryData("fry-usa-basic")
#'
#' # Load avelein-fry gene quantification by summing spliced and ambiguous events
#' sce <- loadFry(testdat$parent_dir, which_counts = c('S', 'A'))
#' SummarizedExperiment::assayNames(sce)
#'
#' # Load the same data and include an "unspliced" assay matrix to pass down
#' # to velociraptor, perhaps
#' scev <- loadFry(testdat$parent_dir, which_counts = c('S', 'A'), velocity = TRUE)
#' SummarizedExperiment::assayNames(scev)
loadFry <- function(fry.dir, which_counts = c('S', 'A'), velocity = FALSE,
                    verbose = FALSE) {
    # Check `fry.dir` is legit
    quant_file <- file.path(fry.dir, "alevin", "quants_mat.mtx")
    if (!file.exists(quant_file)) {
      stop("The `fry.dir` directory provided does not look like a directory generated from alevin-fry:\n",
           sprintf("Missing quant file: %s", quant_file))
    }
    # in alevin-fry 0.4.1, meta_info.json is changed to quant.json, so check both
    # read in metadata
    qfile <- file.path(fry.dir, "quant.json")
    if (!file.exists(qfile)) {
      qfile <- file.path(fry.dir, "meta_info.json")
    }

    # read in metadata
    meta_info <- fromJSON(qfile)
    ng <- meta_info$num_genes
    usa_mode <- meta_info$usa_mode

    # figure out how we're going to summarize expression counts
    wc_opts <- c('U', 'S', 'A')
    if (usa_mode) {
      stopifnot(
        "`which_counts` must be a character vector with length() > 1" = {
          is.character(which_counts) && length(which_counts) > 1
        },
        "`which_counts` can only include elements from c('U', 'S', 'A')" = {
          all(which_counts %in% wc_opts)
        },
        "`velocity` must be a boolean flag" = {
          is.logical(velocity) && length(velocity) == 1L
        })
      if (velocity) {
        if ("U" %in% which_counts) {
          if (verbose) {
            message("'U' removed from `which_counts` when `velocity` set to `TRUE`")
          }
          which_counts <- setdiff(which_counts, "U")
        }
      }
      if (length(which_counts) == 0) {
        if (velocity) {
          stop("Pleae provide at least one of c('S', 'A') in `which_counts` when `velocity = TRUE`")
        } else {
          stop("Pleae provide at least one of c('U', 'S', 'A') in `which_counts`")
        }
      }
      if (verbose) {
        message("processing input in USA mode, will return ", paste(which_counts, collapse = '+'))
        if (velocity) {
          message("unspliced counts will be included in the `'unspliced'` assay matrix")
        }
      }
    } else {
      if (velocity) {
        if (verbose) {
          message("`velocity` is ignored when processing data in standard mode, ",
                  "will return spliced count")
        }
        velocity <- FALSE
      } else if (verbose) {
        message("processing input in standard mode, will return spliced count")
      }
    }

    # if usa mode, each gene gets 3 rows, so the actual number of genes is ng/3
    if (usa_mode) {
      if (ng %% 3 != 0) {
        stop("The number of quantified targets is not a multiple of 3")
      }
      ng <- as.integer(ng/3)
    }

    # read in count matrix, gene names, and barcodes
    af_raw <-  readMM(file = file.path(fry.dir, "alevin", "quants_mat.mtx"))
    afg <-  read.table(file.path(fry.dir, "alevin", "quants_mat_cols.txt"),
                       strip.white = TRUE, header = FALSE, nrows = ng,
                       col.names = c("gene_ids"), row.names = 1)
    afc <-  read.table(file.path(fry.dir, "alevin", "quants_mat_rows.txt"),
                       strip.white = TRUE, header = FALSE,
                       col.names = c("barcodes"), row.names = 1)

    # if in usa_mode, sum up counts in different status according to which_counts
    if (usa_mode) {
      rd <- list("S" = seq(1, ng), "U" =  seq(ng + 1, 2 * ng),
                 "A" =  seq(2 * ng + 1, 3 * ng))
      o <- af_raw[, rd[[which_counts[1]]], drop = FALSE]
      for (wc in which_counts[-1]) {
        o <- o + af_raw[, rd[[wc]], drop = FALSE]
      }
    } else {
      o <- af_raw
    }

    alist <- list(counts = t(o))
    if (velocity) {
      alist[["unspliced"]] <- t(af_raw[, rd[["U"]], drop = FALSE])
    }

    # create SingleCellExperiment object
    sce <- SingleCellExperiment(alist, colData = afc, rowData = afg)
    sce
}
