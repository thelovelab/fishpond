#' Load in data from alevin-fry USA mode
#'
#' Enables easy loading of sparse data matrices provided by
#' alevin-fry USA mode.
#'
#' @param fryDir path to the output directory returned by
#' alevin-fry quant command. This directory should contain a
#' \code{metainfo.json}, and an alevin folder which contains
#' \code{quants_mat.mtx}, \code{quants_mat_cols.txt} and
#' \code{quants_mat_rows.txt}
#' @param outputFormat can be \emph{either} be a list that defines the
#' desired format of the output \code{SingleCellExperiment} object
#' \emph{or} a string that represents one of the pre-defined output
#' formats, which are "scRNA", "snRNA", "all", "scVelo", "velocity", "U+S+A" and "S+A".
#' See details for the explanations of the pre-defined formats and
#' how to define custom format.
#' @param nonzero whether to filter cells with non-zero expression
#' value across all genes (default \code{FALSE}).
#' If \code{TRUE}, this will filter based on all assays.
#' If a string vector of assay names, it will filter based
#' on the matching assays in the vector.
#' If not in USA mode, it must be TRUE/FALSE/counts.
#' @param quiet logical whether to display no messages
#'
#' @section Details about \code{loadFry}:
#' This function consumes the result folder returned by running
#' alevin-fry quant in unspliced, spliced, ambiguous (USA) 
#' quantification mode, and returns a \code{SingleCellExperiment} object
#' that contains a final count for each gene within each cell. In
#' USA mode, alevin-fry quant returns a count matrix contains three
#' types of count for each feature (gene) within each sample (cell
#' or nucleus), which represent the spliced mRNA count of the gene (S),
#' the unspliced mRNA count of the gene (U), and the count of UMIs whose
#' splicing status is ambiguous for the gene (A). For each assay
#' defined by \code{outputFormat}, these three counts of a gene
#' within a cell will be summed to get the final count of the gene
#' according to the rule defined in the \code{outputFormat}.  The
#' returned object will contains the desired assays defined by
#' \code{outputFormat}, with rownames as the barcode of samples and
#' colnames as the feature names.
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
#' It returns a \code{counts} assay that contains the S+A count of each gene in each cell,
#' and a \code{unspliced} assay that contains the U count of each gene in each cell.}
#' \item{"snRNA", "all" and "U+S+A":}{These three formats are the same.
#' They return a \code{counts} assay that contains the U+S+A count of each gene in 
#' each cell without any extra layers. "snRNA" is recommended for single-nucleus 
#' RNA-sequencing experiments. CellRanger 7 returns this format for both single-cell 
#' and single-nucleus experiments.}
#' \item{"S+A":}{returns a \code{counts} assay that contains the S+A 
#' count of each gene in each cell.}
#' \item{"raw":}{This format put the three kinds of counts into three separate assays, 
#' which are \code{unspliced}, \code{spliced} and \code{ambiguous}.}
#' \item{"velocity":}{This format contains two assays. 
#' The \code{spliced} assay contains the S+A count of each gene in each cell.
#' The \code{unspliced} assay contains the U counts of each gene in each cell.}
#' \item{"scVelo":}{This format is for direct entry into velociraptor R package or 
#' other scVelo downstream analysis pipeline for velocity
#' analysis in R with Bioconductor. It adds the expected 
#' "S"-pliced assay and removes errors for size factors being
#' non-positive.}
#' }
#' 
#' A custom output format can be defined using a list. Each element in the list 
#' defines an assay in the output \code{SingleCellExperiment} object. 
#' The name of an element in the list will be the name of the corresponding 
#' assay in the output object. Each element in the list should be defined as 
#' a vector that takes at least one of the three kinds of count, which are U, S and A.
#' See the provided toy example for defining a custom output format.
#' 
#' @return A \code{SingleCellExperiment} object that contains one
#' or more assays. Each assay consists of a gene by cell count matrix.
#' The row names are feature names, and the column names are cell
#' barcodes
#'
#' @references
#'
#' alevin-fry publication:
#'
#' He, D., Zakeri, M., Sarkar, H. et al. "Alevin-fry unlocks rapid, accurate and 
#' memory-frugal quantification of single-cell RNA-seq data."
#' Nature Methods 19, 316–322 (2022).
#' \url{https://doi.org/10.1038/s41592-022-01408-3}
#' 
#' @examples
#' 
#' # Get path for minimal example avelin-fry output dir
#' testdat <- fishpond:::readExampleFryData("fry-usa-basic")
#' 
#' # This is exactly how the velocity format defined internally.
#' custom_velocity_format <- list("spliced"=c("S","A"), "unspliced"=c("U"))
#'
#' # Load alevin-fry gene quantification in velocity format
#' sce <- loadFry(fryDir=testdat$parent_dir, outputFormat=custom_velocity_format)
#' SummarizedExperiment::assayNames(sce)
#'
#' # Load the same data but use pre-defined, velociraptor R pckage desired format
#' scvelo_format <- "scVelo"
#' 
#' scev <- loadFry(fryDir=testdat$parent_dir, outputFormat=scvelo_format, nonzero=TRUE)
#' SummarizedExperiment::assayNames(scev)
#'
#' @author Dongze He, with contributions from Steve Lianoglou, Wes Wilson
#' 
#' @export
loadFry <- function(fryDir, 
                    outputFormat = "scRNA", 
                    nonzero = FALSE,
                    quiet = FALSE) {
  
  # load in fry result
  fry.raw <- load_fry_raw(fryDir, quiet)
  meta.data <- fry.raw$meta.data
  
  
  # if in usa.mode, sum up counts in different status according to which.counts
  if (meta.data$usa.mode) {
    # preparation
    predefined.format <- list("scrna" = list("counts" = c("S", "A"), "unspliced" = c("U")),
                             "snrna" = list("counts" = c("U", "S", "A")),
                             "all" = list("counts" = c("U", "S", "A")),
                             "U+S+A" = list("counts" = c("U", "S", "A")),
                             "S+A" = list("counts" = c("S", "A")),
                             "velocity" = list("spliced" = c("S", "A"), "unspliced" = c("U")),
                             "scvelo" = list("counts" = c("S", "A"), "spliced" = c("S", "A"), "unspliced" = c("U")),
                             "raw" = list("spliced" = c("S"), "unspliced" = c("U"), "ambiguous" = c("A"))
    )
    valid.counts <- c("U", "S", "A")
    
    # check outputFormat
    if (is.character(outputFormat)) {
      outputFormat = tolower(outputFormat)
      # Check whether outputFormat is a predefined format
      if (! (outputFormat %in% names(predefined.format))) {
        stop("Provided outputFormat string is invalid. Please check the function description
for the list of predifined format")
      }
      
      if (!quiet) {
        message("Using pre-defined output format: ", outputFormat)
      }
      
      # get the assays
      output.assays <- predefined.format[[outputFormat]]
      
    } else if (is.list(outputFormat)) {
      # check whether user-defined assays are all 
      for (assay.name in names(outputFormat)) {
        if (is.null(outputFormat[[assay.name]])) {              
          stop(paste0("The provided assay '", assay.name, "' is empty. Please remove it"))
        }
        else if (!all(outputFormat[[assay.name]] %in% valid.counts)) {
          stop("Please use U, S and A ONLY to indicate which counts will be considered to build a assay")
        }
      }
      if (!all(unlist(outputFormat) %in% valid.counts)) {
        stop("Please use U, S and A ONLY to indicate which counts will be considered to build a assay")
      }
      
      output.assays <- outputFormat
      
      if (!quiet) {
        message("Using user-defined output assays")
      }
    }
    
    # If we are here, the output.assays is valid.
    # then we check the assay names in nonzero
    if (is.logical(nonzero)) {
      if (nonzero) {
        nonzero <- names(output.assays)
      } else {
        if (is.character(outputFormat) && outputFormat == "scvelo") {
          nonzero <- c("counts")
        } else {
          nonzero <- c()
        }
      }
    } else if (is.character(nonzero)) {
      if (length(nonzero) > 0) {
        for (idx in seq_along(nonzero)) {
          if (!nonzero[idx] %in% names(output.assays)) {
            warning(paste0("In the provided nonzero vector, '",
                           nonzero[idx],
                           "' is not one of the output assays, ignored"))
            nonzero <- nonzero[-idx]
          }
        }
      }
    } else {
      warning("Invalid nonzero, ignored")
      nonzero <- c()
    }
    
    # assembly
    alist <- vector(mode = "list", length = length(output.assays))
    names(alist) <- names(output.assays)
    ng <- meta.data$num.genes
    rd <- list("S" = seq(1, ng), "U" =  seq(ng + 1, 2 * ng),
               "A" =  seq(2 * ng + 1, 3 * ng))
    # fill in each assay
    for (assay.name in names(output.assays)) {
      which.counts <- output.assays[[assay.name]]
      alist[[assay.name]] <- fry.raw$count.mat[, rd[[which.counts[1]]], drop = FALSE]
      if (length(which.counts) > 1) {
        # build assay
        for (wc in which.counts[-1]) {
          alist[[assay.name]] <- alist[[assay.name]] + fry.raw$count.mat[, rd[[wc]], drop = FALSE]
        }
      }
      alist[[assay.name]] <- t(alist[[assay.name]])
      if (!quiet) {
        message(paste(c(paste0("Building the '",
                               assay.name,
                               "' assay, which contains"),
                        which.counts),
                      collapse = " "))
      }
    }
  } else {
    # not in USA mode
    if (!quiet) {
      message("Not in USA mode, set assay name as \"counts\"")
    }
    
    if (is.logical(nonzero)) {
      if (nonzero) {
        nonzero <- c("counts")
      } else {
        nonzero <- c()
      }
    } else {
      if (nonzero != "counts") {
        message("Not in USA mode, nonzero must be TRUE/FALSE/counts")
        message("Set nonzero as FALSE")
        nonzero <- c()
      } else {
        nonzero <- c("counts")
      }
    }

    # define output matrix
    alist <- list(counts = t(fry.raw$count.mat))
  }
  
  if (!quiet) {
    message("Constructing output SingleCellExperiment object")
  }
  
  # create SingleCellExperiment object
  sce <- SingleCellExperiment(alist,
                              colData = fry.raw$barcodes,
                              rowData = fry.raw$gene.names)
  
  # filter all zero cells in zero, one or multiple assays
  for (assay.name in nonzero) {
    sce <- sce[, colSums(assay(sce, assay.name)) > 0]
  }
  
  if (!quiet) {
    message("Done")
  }
  
  sce
}

load_fry_raw <- function(fryDir, quiet = FALSE) {
  # Check `fryDir` is legit
  if (!quiet) {
    message("locating quant file")
  }
  quant.file <- file.path(fryDir, "alevin", "quants_mat.mtx")
  if (!file.exists(quant.file)) {
    stop("The `fryDir` directory provided does not look like a directory generated from alevin-fry:\n",
         sprintf("Missing quant file: %s", quant.file)
    )
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
  
  if (!quiet) {
    message("Reading meta data")
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
  count.mat <- my_as_dgcmatrix(readMM(file = quant.file))
  afg <-  read.table(file.path(fryDir, "alevin", "quants_mat_cols.txt"),
                     strip.white = TRUE, header = FALSE, nrows = ng,
                     col.names = c("gene_ids"))
  rownames(afg) <- afg$gene_ids
  afc <-  read.table(file.path(fryDir, "alevin", "quants_mat_rows.txt"),
                     strip.white = TRUE, header = FALSE,
                     col.names = c("barcodes"))
  rownames(afc) <- afc$barcodes
  
  if (!quiet) {
    message(paste("Processing", ng, "genes", "and", nrow(count.mat), "barcodes"))
  }
  
  list(count.mat = count.mat, gene.names = afg, barcodes = afc, meta.data = list(num.genes = ng, usa.mode = usa.mode))
}

# Functions to help with running or creating test data
# none of these functions are exported on purpose

# Example alevin-fry quant dataset ---------------------------------------------
#
# These methods provide an orthologous way to create and read in test examples
# from alevin outputs than what is implemented in loadFry().
#
# The functions provide a mechanism that enables you to define the sample output
# in a plain text matrix format (data.frame) named `example-dat.csv` and
# place that in a new directory under extdata/alevin/example-quants.
#
# Look at the `extdata/alevin/example-quants/fry-usa-basic/example-dat.csv` for
# an example data file that created a sample `salmon alevin-fry quant` directory
# we can use to test `loadFry` against.
#
# Once we created the extdata/alevin/example-quants/fry-usa-basic directory
# and put the the example-dat.csv file in it, we then run the following command
# create all of the *.mtx and other files to make this look like an alevin
# output directoy.
#
# ```{r}
# devtools::load_all(".")
# writeExampleFryDat("fry-usa-basic")
# ```

#' Loads an example data matrix from a csv data from a top-level example
#' fry output directory.
#'
#' @noRd
#' @param example_name the name of the folder that holds the example data.
#' @return a list of primitive data types you can use to serialize a mock
#' output dir from a `salmon alevin-fry quant` run.
readExampleFryData <- function(example_name = "fry-usa-basic", usa = TRUE, ...) {
  example_dir <- system.file("extdata", "alevin", "example-quants",
                             example_name, package = "fishpond",
                             mustWork = TRUE)
  dat <- read.csv(file.path(example_dir, "example-dat.csv"),
                  strip.white = TRUE, comment.char = "#")
  m <- as.matrix(dat)
  dimnames(m) <- list(NULL, NULL)

  if (usa) {
    genes <- unique(sub("_.*$", "", colnames(dat)))
  } else {
    genes <- colnames(dat)
  }

  out <- list(
    parent_dir = example_dir,
    matrix = m,
    barcodes = rownames(dat),
    genes = colnames(dat))

  if (usa) {
    out$genes <- unique(sub("_.*$", "", out$genes))
    out$usa <- list(
      U = grep("_U$", colnames(dat)),
      S = grep("_S$", colnames(dat)),
      A = grep("_A$", colnames(dat)))
  }

  out
}

#' Serializes example fry output data from an `example-fry-dat.csv` as if
#' it were produced from an alevin-fry quant run
#'
#' @noRd
#' @param x the name of the top level directory in
#' `extdata/alevin/example-fry-USA-quants` or a list result from calling
#' `readExampleFryData`
writeExampleFryDat <- function(x = "fry-usa-basic", ...) {
  if (is.character(x)) {
    x <- readExampleFryData(x)
  }
  stopifnot(
    is.list(x),
    all(c("matrix", "genes", "barcodes", "parent_dir") %in% names(x)))
  if (is.null(x$usa)) {
    stop("Haven't kicked the tires on a non USA like output")
  }

  out.dir <- file.path(x$parent_dir, "alevin")

  if (!dir.exists(out.dir)) {
    stopifnot(dir.create(out.dir))
  }

  m <- my_as_dgcmatrix(Matrix::Matrix(x$matrix, sparse = TRUE))
  
  Matrix::writeMM(m, file.path(out.dir, "quants_mat.mtx"))

  write.table(
    x$genes,
    file = file.path(out.dir, "quants_mat_cols.txt"),
    quote = FALSE,col.names = FALSE, row.names = FALSE, sep = "\t")

  write.table(
    x$barcodes,
    file = file.path(out.dir, "quants_mat_rows.txt"),
    quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

  # create metadata
  meta_info = list()
  meta_info[["alt_resolved_cell_numbers"]] = list()
  meta_info[["cmd"]] = ""
  meta_info[["dump_eq"]] = FALSE
  meta_info[["num_genes"]] = length(x$genes) * ifelse(is.null(x$usa), 1, 3)
  meta_info[["num_quantified_cells"]] = length(x$barcodes)
  meta_info[["resolution_strategy"]] = "CellRangerLike"
  meta_info[["usa_mode"]] = TRUE

  write(
    jsonlite::toJSON(meta_info, pretty=TRUE),
    file = file.path(x$parent_dir, "meta_info.json"))
}

# added August 2022
# see https://cran.r-project.org/web/packages/Matrix/vignettes/Design-issues.pdf
my_as_dgcmatrix <- function(matrix) {
  as(as(as(matrix, "dMatrix"), "generalMatrix"), "CsparseMatrix")
}
