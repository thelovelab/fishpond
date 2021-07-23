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
#' is recommanded.  
#' @param verbose A boolean specifying if showing messages when running the function 
#'
#' @details
#' This function consumes the result folder returned by running alevin-fry quant in unspliced, spliced, ambiguous (USA) quantification mode, and returns a `SingleCellExperiement` object
#' that contains a final count for each gene within each cell. In USA mode, alevin-fry quant returns a count matrix contains three types of count for each feature (gene) 
#' within each sample (cell or nucleus), which represent the spliced mRNA count of the gene, the unspliced mRNA count of the gene, and the count of mRNA whose splicing status is ambiguous.
#' In this function, these three counts of a gene within a cell will either be summed or discarded to get the final count of the gene by specifying `which_counts`. The returned object will 
#' contains a gene by cell count matrix, with rownames as the barcode of samples and colnames as the feature names.
#'
#'
#' @return A `SingleCellExperiment` object contains a gene by cell count matrix. The row names are feature names, and the column names are cell barcodes.
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
#' \dontrun{
#' # For output from alevin-fry USA mode 
#' fry.dir <- system.file("inst", "extdata", "alevin", "test_loadFry", package="fishpond")
#' is.element('meta_info.json', list.files('fry.dir')) # Should reture TRUE
#' list.files(file.path(fry.dir, 'alevin')) 
#' # Should show quants_mat.mtx, quants_mat_cols.txt and quants_mat_rows.txt
#' 
#' sce <- loadFry(fry.dir = fry.dir, which_counts = c('S', 'A'))
#' }
#'

loadFry <- function(fry.dir, which_counts = c('S', 'A'), verbose = FALSE) {
    # in alevin-fry 0.4.1, meta_info.json is changed to quant.json, so check both
    # read in metadata
    qfile <- file.path(fry.dir, "quant.json")
    if (!file.exists(qfile)) {
      qfile <- file.path(fry.dir, "meta_info.json")
    }

    # read in metadata
    meta_info <- fromJSON(file = qfile)
    ng <- meta_info$num_genes
    usa_mode <- meta_info$usa_mode
    
    if (usa_mode) {
      if (length(which_counts) == 0) {
        stop("Please at least provide one status in 'U' 'S' 'A' ")
      }
      if (verbose) {
        message("processing input in USA mode, will return ", paste(which_counts, collapse = '+'))
      }
    } else if (verbose) {
      message("processing input in standard mode, will return spliced count")
    }
    
    # read in count matrix
    af_raw <-  readMM(file = file.path(fry.dir, "alevin", "quants_mat.mtx"))
    # if usa mode, each gene gets 3 rows, so the actual number of genes is ng/3
    if (usa_mode) {
      if (ng %% 3 != 0) {
        stop("The number of quantified targets is not a multiple of 3")
      }
      ng <- as.integer(ng/3)
    }
    
    # read in gene name file and cell barcode file
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

    # create SingleCellExperiment object
    sce <-  SingleCellExperiment(list(counts = t(o)),
                                colData = afc,
                                rowData = afg
    )
    sce
}