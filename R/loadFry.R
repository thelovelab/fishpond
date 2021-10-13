#' Load in data from alevin-fry USA mode
#'
#' Enables easy loading of sparse data matrices provided by alevin-fry USA mode.
#' Alevin-fry - \url{https://www.biorxiv.org/content/10.1101/2021.06.29.450377v1}
#'
#' @param fryDir A path to the output directory returned by
#' alevin-fry quant command. This directory should contain a
#' \code{metainfo.json}, and an alevin folder which contains
#' \code{quants_mat.mtx}, \code{quants_mat_cols.txt} and
#' \code{quants_mat_rows.txt}
#' @param outputFormat This argument can \emph{either} be a list that defines the desired format of the output 
#' \code{SingleCellExperiment} object \emph{or} a string that represents one of 
#' the pre-defined output formats, which are "scRNA", "snRNA", "scVelo" and "velocity". 
#' See details for the explainations of the pre-defined formats and how to define custom format.
#' @param nonzero A boolean specifying if filtering
#' genes with non-zero expression value across all cells in the output \code{SingleCellExperiment} object.
#' @param verbose A boolean specifying if showing
#' messages when running the function
#'
#' @section Details about \code{loadFry}:
#' This function consumes the result folder returned by running
#' alevin-fry quant in unspliced, spliced, ambiguous (USA) 
#' quantification mode, and returns a \code{SingleCellExperiement} object
#' that contains a final count for each gene within each cell. In
#' USA mode, alevin-fry quant returns a count matrix contains three
#' types of count for each feature (gene) within each sample (cell
#' or nucleus), which represent the spliced mRNA count of the gene (U),
#' the unspliced mRNA count of the gene (S), and the count of mRNA whose
#' splicing status is ambiguous (A). For each slot defined by \code{outputFormat},
#' these three counts of a gene within a cell will be summed 
#' to get the final count of the gene according to the rule defined in the \code{outputFormat}. 
#' The returned object will contains the desired slots defined by \code{outputFormat}, 
#' with rownames as the barcode of samples and colnames as the feature names.
#' 
#' 
#' @section Details about the output format:
#' The \code{outputFormat} argument takes \emph{either} be a list that defines 
#' the desired format of the output 
#' \code{SingleCellExperiment} object \emph{or} a string that represents one of 
#' the pre-defined output format. 
#' 
#' Currently the pre-defined formats 
#' of the output \code{SingleCellExperiment} object are: 
#' \describe{
#' \item{"scRNA":}{This format is recommended for single cell experiments. 
#' It returns a \code{counts} assay slot that contains the S+A count of each gene in each cell.}
#' \item{"snRNA":}{This format is recommended for single nucleus experiments. 
#' It returns a \code{counts} assay slot that contains the U+S+A count of each gene in each cell.}
#' \item{"raw":}{This format put the three kinds of counts into three separate assay slots, 
#' which are \code{unspliced}, \code{spliced} and \code{ambiguous}.}
#' \item{"velocity":}{This format contains two assay slots. 
#' The \code{spliced} slot contains the S+A count of each gene in each cell.
#' The \code{unspliced} slot contains the U counts of each gene in each cell.}
#' \item{"scVelo":}{This format is for direct entry into velociraptor R package or 
#' other scVelo downstream analysis pipeline for velocity
#' analysis in R with Bioconductor. It adds the expected 
#' "S"-pliced assay and removes errors for size factors being
#' non-positive.}
#' }
#' 
#' A custom output format can be defined using a list. Each element in the list 
#' defines an assay slot in the output \code{SingleCellExperiment} object. 
#' The name of an element in the list will be the name of the corresponding 
#' assay slot in the output object. Each element in the list should be defined as 
#' a vector that takes at least one of the three kinds of count, which are U, S and A.
#' See the provided toy example for defining a custom output format.
#' 
#' @section Details about \code{load_fry_raw}:
#' This function processes alevin-fry's quantification result contained within the input folder.
#' This function returns a list that consists of the gene count matrix, the gene names list, the barcode list, 
#' and some metadata, such as the number of genes in the experiment and whether alevin-fry was executed 
#' in USA mode. In the returned list, the all-in-one count matrix, \code{count_mat}, 
#' returned from the USA mode of alevin-fry consists of the spliced count of genes defined in \code{gene.names}
#' for all barcodes defined in \code{barcodes}, followed by the unspliced count of genes in the same order 
#' for all cells, then followed by the ambiguous count of genes in the same order for all cells.  
#'
#' @return A \code{SingleCellExperiment} object that contains one or more assay slots.
#' Each assay slot consists of a gene by cell
#' count matrix. The row names are feature names, and the column
#' names are cell barcodes.
#'
#' @concept preprocessing
#'
#' @examples
#' 
#' # Get path for minimal example avelin-fry output dir
#' testdat <- fishpond:::readExampleFryData("fry-usa-basic")
#' 
#' # This is exactly how the velocity format defined internally.
#' custom_velocity_format = list("spliced" = c("S", "A"), "unspliced" = c("U"))
#'
#' # Load alevin-fry gene quantification in velocity format
#' sce <- loadFry(fryDir = testdat$parent_dir, outputFormat = custom_velocity_format, verbose = TRUE)
#' SummarizedExperiment::assayNames(sce)
#'
#' # Load the same data but use pre-defined, velociraptor R pckage desired format
#' scvelo_format = "scVelo"
#' 
#' scev <- loadFry(fryDir = testdat$parent_dir, outputFormat = scvelo_format, nonzero = TRUE)
#' SummarizedExperiment::assayNames(scev)
#' 
#' @name loadFry
NULL

#' @export
#' @rdname loadFry
#' @importFrom jsonlite fromJSON
#' @importFrom Matrix readMM
#' @importFrom utils read.table
load_fry_raw <- function(fryDir, verbose = FALSE) {
  # Check `fryDir` is legit
  quant.file <- file.path(fryDir, "alevin", "quants_mat.mtx")
  if (!file.exists(quant.file)) {
    stop("The `fryDir` directory provided does not look like a directory generated from alevin-fry:\n",
      sprintf("Missing quant file: %s", quant.file)
    )
  }
  if (verbose) {
    message("Quant file found.")
  }
  
  # Since alevin-fry 0.4.1, meta_info.json is changed to quant.json, we check both
  # read in metadata
  qfile <- file.path(fryDir, "quant.json")
  if (!file.exists(qfile)) {
    qfile <- file.path(fryDir, "meta_info.json")
  }
  
  # read in metadata
  meta.info <- fromJSON(qfile)
  ng <- meta.info$num_genes
  usa.mode <- meta.info$usa_mode
  
  if (verbose) {
    message("Meta data read.")
    message(paste0("USA mode: ", usa.mode))
  }
  
  # if usa mode, each gene gets 3 rows, so the actual number of genes is ng/3
  if (usa.mode) {
    if (ng %% 3 != 0) {
      stop("The number of quantified targets is not a multiple of 3")
    }
    ng <- as.integer(ng/3)
  }
  
  # read in count matrix, gene names, and barcodes
  count.mat <- as(readMM(file = quant.file), "dgCMatrix")
  afg <-  read.table(file.path(fryDir, "alevin", "quants_mat_cols.txt"),
                     strip.white = TRUE, header = FALSE, nrows = ng,
                     col.names = c("gene_ids"))
  rownames(afg) = afg$gene_ids
  afc <-  read.table(file.path(fryDir, "alevin", "quants_mat_rows.txt"),
                     strip.white = TRUE, header = FALSE,
                     col.names = c("barcodes"))
  rownames(afc) = afc$barcodes
  
  if (verbose) {
    message(paste("Processing", ng, "genes", "and", nrow(count.mat), "barcodes."))
  }
  
  list(count.mat = count.mat, gene.names = afg, barcodes = afc, meta.data = list(num.genes = ng, usa.mode = usa.mode))
}

#' @export
#' @rdname loadFry
#' @importFrom Matrix t colSums
#' @importFrom SingleCellExperiment SingleCellExperiment
loadFry <- function(fryDir, 
                    outputFormat = "scRNA", 
                    nonzero = FALSE,
                    verbose = FALSE) {

  # scvelo doesn't actually need this. but sce /w velociraptor /w defaults
  # requires size factors should be positive and this was the easiest solution I could think of,
  # with the assumption, that BioC users generally build their pipelines with other BioC
  # packages with up or downstream analysis in mind.
  if(outputFormat == "scVelo" && nonzero == FALSE) {
    message("velociraptor R package filters genes with zero expression by default. To mimic this behavior, please set nonzero = TRUE.")
  }
  
  # load in fry result
  fry.raw = load_fry_raw(fryDir, verbose)
  meta.data = fry.raw$meta.data
  
  
  # if in usa.mode, sum up counts in different status according to which.counts
  if (meta.data$usa.mode) {
    # preparation
    predefined.format = list("scRNA" = list("counts" = c("S", "A")),
                             "snRNA" = list("counts" = c("S", "A")),
                             "velocity" = list("spliced" = c("S", "A"), "unspliced" = c("U")),
                             "scVelo" = list("counts" = c("S", "A"), "spliced" = c("S", "A"), "unspliced" = c("U")),
                             "raw" = list("spliced" = c("S"), "unspliced" = c("U"), "ambiguous" = c("A"))
    )
    valid.counts = c("U", "S", "A")
    
    # check outputFormat
    if (is.character(outputFormat)) {
      # Check whether outputFormat is a predefined format
      if (! (outputFormat %in% names(predefined.format))) {
        stop("Provided outputFormat string is invalid. Please check the function description for the list of predifined format.")
      }
      
      if (verbose) {
        message("Using pre-defined output format: ", outputFormat)
      }
      
      # get the slots
      output.slots = predefined.format[[outputFormat]]
      
    } else if (is.list(outputFormat)) {
      # check whether user-defined slots are all 
      for (slot.name in names(outputFormat)) {
        if (is.null(outputFormat[[slot.name]])) {              
          stop(paste0("The provided slot \"", slot.name, "\" is empty. Please remove it."))
        }
        else if (!all(outputFormat[[slot.name]] %in% valid.counts)) {
          stop("Please use U, S and A ONLY to indicate which counts will be considered to build a slot")
        }
      }
      if (!all(unlist(outputFormat) %in% valid.counts)) {
        stop("Please use U, S and A ONLY to indicate which counts will be considered to build a slot")
      }
      
      output.slots = outputFormat
      
      if (verbose) {
        message("Using user-defined output slots.")
      }
    }
    # If we are here, the output.slots is valid.

    alist <- vector(mode = "list", length = length(output.slots))
    names(alist) = names(output.slots)
    ng = meta.data$num.genes
    rd <- list("S" = seq(1, ng), "U" =  seq(ng + 1, 2 * ng),
               "A" =  seq(2 * ng + 1, 3 * ng))
    # fill in each slot
    for (slot.name in names(output.slots)) {
      which.counts = output.slots[[slot.name]]
      alist[[slot.name]] <- fry.raw$count.mat[, rd[[which.counts[1]]], drop = FALSE]
      if (length(which.counts) > 1) {
        # build slot
        for (wc in which.counts[-1]) {
          alist[[slot.name]] <- alist[[slot.name]] + fry.raw$count.mat[, rd[[wc]], drop = FALSE]
        }
      }
      alist[[slot.name]] = t(alist[[slot.name]])
      if (verbose) {
        message(paste(c(paste0("Building the \"", slot.name, "\" slot, which contains"), which.counts), collapse = " "))
      }
    }
  } else {
    if(verbose) {
      message("Not in USA mode, ignore argument outputFormat.")
    }

    # define output matrix
    alist = list(counts = t(fry.raw$count.mat))
  }
  
  if (verbose) {
    message("Constructing output SingleCellExperiment object.")
  }
  
  # create SingleCellExperiment object
  sce <- SingleCellExperiment(alist, colData = fry.raw$barcodes, rowData = fry.raw$gene.names)
  
  # filter all zero genes
  if (nonzero) {
    for(assay.name in names(sce@assays)) {
      sce <- sce[, colSums(sce@assays@data[[assay.name]]) > 0]
    }
  }
                                                    
  if (verbose) {
    message("Done.")
  }
  
  sce
}
