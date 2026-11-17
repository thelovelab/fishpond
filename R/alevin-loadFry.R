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
#' formats, which are "scRNA" (or "S+A"), "snRNA" (or "all", "U+S+A"), "scVelo", and "velocity".
#' See details for the explanations of the pre-defined formats and
#' how to define custom format.
#' @param aux_columns vector of characters containing the column names of the
#' auxiliary information in the barcodes file starting from the second column.
#' The first column is assumed to be the barcodes and is named as "barcodes".
#' Extra auxiliary columns in the barcodes file without a specified name will be ignored.
#' @param gene_id_to_name path to a file that contains the mapping from gene names to gene ids. 
#' It is only needed if\cr
#'   1. you are not using the simpleaf pipeline (`simpleaf index` + `simpleaf quant`),\cr
#'   2. you have such a file, and,\cr
#'   3. you want to add this information to the colData of your output object.\cr
#' If you do, please ensure it is a tab-separated, two-column file without a header, and the first column is the gene ids and the second column is the gene names.
#' @param nonzero whether to filter cells with non-zero expression
#' value across all genes (default \code{FALSE}).
#' If \code{TRUE}, this will filter based on all assays.
#' If a string vector of assay names, it will filter based
#' on the matching assays in the vector.
#' If not in USA mode, it must be TRUE/FALSE/counts.
#' @param quiet logical whether to display no messages
#' @param featureLevel \code{"gene"} (the default) reads the cell by gene
#' matrix \code{quants_mat.mtx}. \code{"probe"} reads the cell by probe
#' matrix \code{probe_quants_mat.mtx}, which \code{alevin-fry quant-probe
#' --probe-mtx} (or \code{simpleaf multiplex-quant --probe-mtx}) writes.
#' The function stops with an error if \code{quant.json} does not have
#' \code{"probe_mtx": true}.
#' @param probeRef only for \code{featureLevel = "probe"}: the path to the
#' 10x probe-set CSV of the run, which gives the gene of each probe. If it
#' is \code{NULL} (the default), the function uses \code{probe_set.csv} in
#' \code{fryDir}. The function stops with an error if the file does not
#' exist. See the Flex section.
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
#' RNA-sequencing experiments. "raw" is recommended for mimicing CellRanger 7's behavior, 
#' which returns this format for both single-cell and single-nucleus experiments.}
#' \item{"S+A":}{This format returns a \code{counts} assay that contains the S+A 
#' count of each gene in each cell.}
#' \item{"raw":}{This format puts the three kinds of counts into three separate assays, 
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
#' @section Details about Flex (multiplex) output:
#' \code{loadFry} also reads the output of \code{simpleaf multiplex-quant}
#' and \code{alevin-fry quant-probe} (10x Flex). In this output, the
#' row label of a cell in \code{quants_mat_rows.txt} is
#' \code{<sample_name>_<cell_barcode>}, because the same barcode can occur
#' in two samples. The function adds the two parts to the \code{colData}
#' as \code{sample_name} and \code{cell_barcode}. They come from the
#' \code{sample_name} and \code{CB} columns of \code{featureDump.txt},
#' whose line i describes cell i; the function stops with an error if a
#' line does not match its cell. If there is no \code{featureDump.txt}
#' file, but the directory has a \code{collation_manifest.bin} or a
#' \code{sample_info.json} file (a multiplex run), the function gives a
#' warning and adds no columns. Other outputs get no new columns.
#'
#' The \code{gene_id_to_name.tsv} file of a Flex run can have more gene ids
#' than the matrix. The function adds the gene name of each gene in the
#' matrix (a left join); a gene that is not in the file gets \code{NA}.
#'
#' With \code{featureLevel = "probe"}, the rows are the probes in
#' \code{probe_quants_mat_cols.txt}, and the columns are the same cells as
#' at gene level. The \code{rowData} is the probe-set CSV as it is: all its
#' columns, with their names (for example \code{gene_id}, \code{probe_seq},
#' \code{probe_id}, \code{included}, \code{region} and \code{gene_name}),
#' read with the \code{read.csv} defaults (for example, \code{included} is
#' logical). The function skips the lines at the top of
#' the CSV that start with \code{#}, and uses the CSV rows of the matrix
#' probes, in the order of the matrix (by \code{probe_id}, because the probe
#' order of the matrix is not the order of the CSV). It stops with an error
#' if a probe of the matrix is not in the CSV, or if the CSV has a
#' \code{probe_id} more than one time. A 10x \code{"unspliced"} probe does
#' not span a splice junction; it does not mean pre-mRNA. The sum of the
#' probe counts of a gene (\code{gene_id}) is the gene count of the gene
#' matrix.
#'
#' @return A \code{SingleCellExperiment} object that contains one
#' or more assays. Each assay consists of a gene by cell count matrix.
#' The row names are feature names, and the column names are cell
#' barcodes (for a multiplex run, \code{<sample_name>_<cell_barcode>})
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
#' # Load a minimal Flex (multiplex) output, at gene level and at probe level
#' flexdir <- system.file("extdata", "alevin", "example-quants", "flex",
#'                        "flex_probe", package = "fishpond")
#' sce <- loadFry(flexdir)
#' table(sce$sample_name)
#' scep <- loadFry(flexdir, featureLevel = "probe")
#' SummarizedExperiment::rowData(scep)
#'
#' @author Dongze He, with contributions from Steve Lianoglou, Wes Wilson
#' 
#' @export
loadFry <- function(fryDir, 
                    outputFormat = "scRNA", 
                    aux_columns = c("X", "Y"),
                    gene_id_to_name = NULL,
                    nonzero = FALSE,
                    quiet = FALSE,
                    featureLevel = c("gene", "probe"),
                    probeRef = NULL) {
  featureLevel <- match.arg(featureLevel)

  # load in fry result
  fry.raw <- load_fry_raw(fryDir, aux_columns, gene_id_to_name, quiet,
                          featureLevel, probeRef)
  meta.data <- fry.raw$meta.data
  
  
  # if in usa.mode, sum up counts in different status according to which.counts
  if (meta.data$usa.mode) {
    # preparation
    predefined.format <- list("scrna" = list("counts" = c("S", "A"), "unspliced" = c("U")),
                             "snrna" = list("counts" = c("U", "S", "A")),
                             "all" = list("counts" = c("U", "S", "A")),
                             "u+s+a" = list("counts" = c("U", "S", "A")),
                             "s+a" = list("counts" = c("S", "A")),
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
        stop("Provided outputFormat string is invalid. Please check the function description for the list of predifined format")
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

load_fry_raw <- function(fryDir, aux_columns, gene_id_to_name, quiet = FALSE,
                         featureLevel = "gene", probeRef = NULL) {
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
    if (!file.exists(qfile)) {
      stop("The `fryDir` directory provided does not look like a directory generated from alevin-fry:\n",
           sprintf("Missing meta_info.json or quant.json file: %s", qfile)
      )
    }
  }

  # read in metadata
  meta.info <- fromJSON(qfile)
  ng <- meta.info$num_genes
  usa.mode <- meta.info$usa_mode

  # the probe matrix of alevin-fry quant-probe --probe-mtx
  probe.level <- featureLevel == "probe"
  if (probe.level) {
    if (!isTRUE(meta.info$probe_mtx)) {
      stop(sprintf(paste0("featureLevel = \"probe\" needs a probe matrix, but %s does not ",
                          "have \"probe_mtx\": true. Run alevin-fry quant-probe (or simpleaf ",
                          "multiplex-quant) with --probe-mtx."), qfile))
    }
    quant.file <- file.path(fryDir, "alevin", "probe_quants_mat.mtx")
  }
  mat.prefix <- if (probe.level) "probe_quants_mat" else "quants_mat"
  
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
  if (probe.level) {
    afg <- .readProbeVar(fryDir, .findProbeRef(fryDir, probeRef), quiet)
    ng <- nrow(afg)
  } else {
    afg <-  read.table(file.path(fryDir, "alevin", "quants_mat_cols.txt"),
                       strip.white = TRUE, header = FALSE, nrows = ng,
                       col.names = c("gene_ids"))
    rownames(afg) <- afg$gene_ids
  }
  afc <- read.table(file.path(fryDir, "alevin", paste0(mat.prefix, "_rows.txt")),
                     strip.white = TRUE, header = FALSE)
  
  # assign colnames to afc
  if (ncol(afc) == 1) {
    colnames(afc) <- "barcodes"
  } else {
    afc_colnames <- c("barcodes", aux_columns)
    num_cols = min(ncol(afc), length(afc_colnames))

    if (ncol(afc) != length(afc_colnames)) {
      warning(paste0("Number of auxiliary columns in barcodes file does not match the provided column names. Using the first", num_cols, "columns in the barcode file with names ", paste(afc_colnames[1:num_cols], collapse = ", ")))
    }

    afc = afc[, 1:num_cols]
    colnames(afc) <- afc_colnames[1:num_cols]
  }

  # assign rownames to afc
  rownames(afc) <- afc$barcodes

  # for a multiplex (Flex) run, add the cell barcode and the sample name
  afc <- .addFrySampleInfo(fryDir, afc, quiet)

  # read in gene id to gene name mapping
  gene_id_to_name_path = NULL
  if (!is.null(gene_id_to_name)) {
    if (file.exists(gene_id_to_name)) {
      gene_id_to_name_path <- gene_id_to_name
    } else {
      warning("The provided gene_id_to_name file does not exist; ignored")
    }
  } else {
    default_gene_id_to_name_path <- file.path(fryDir, "gene_id_to_name.tsv")
    if (file.exists(default_gene_id_to_name_path)) {
      gene_id_to_name_path <- default_gene_id_to_name_path
      .say(quiet, "Using simpleaf default gene id to name file gene_id_to_name.tsv in the given fryDir")
    }
  }

  # see if we can add gene names to the coldata
  # (not at probe level: the rowData is the probe-set CSV as it is)
  if (!probe.level && !is.null(gene_id_to_name_path)) {
    gene_id_to_name <- read.table(gene_id_to_name_path, sep = "\t", header = FALSE)
    names(gene_id_to_name) <- c("gene_ids", "gene_names")
    rownames(gene_id_to_name) <- gene_id_to_name$gene_ids
    afg$gene_names <- gene_id_to_name[afg$gene_ids, "gene_names"]
  }

  if (!quiet) {
    message(paste("Processing", ng, if (probe.level) "probes" else "genes",
                  "and", nrow(count.mat), "barcodes"))
  }
  
  list(count.mat = count.mat, gene.names = afg, barcodes = afc, meta.data = list(num.genes = ng, usa.mode = usa.mode))
}

# Return the path of the 10x probe-set CSV: `probeRef`, or else
# <fryDir>/probe_set.csv. Stop if the file does not exist.
.findProbeRef <- function(fryDir, probeRef) {
  path <- if (is.null(probeRef)) file.path(fryDir, "probe_set.csv") else probeRef
  if (!file.exists(path)) {
    stop(sprintf(paste0("Cannot find the probe-set CSV: %s does not exist. Give the 10x ",
                        "probe-set CSV of the run with the probeRef argument."), path))
  }
  path
}

# Return the rowData of the probe matrix: the rows of the 10x probe-set CSV,
# as they are (all columns, read with the read.csv defaults), in the column
# order of probe_quants_mat.mtx.
.readProbeVar <- function(fryDir, path, quiet) {
  .say(quiet, "Using the probe-set CSV ", path)
  # skip the "#key=value" lines at the top of the CSV
  lines <- readLines(path)
  n.meta <- match(FALSE, startsWith(lines, "#")) - 1
  csv <- read.csv(path, skip = n.meta)
  dup <- unique(csv$probe_id[duplicated(csv$probe_id)])
  if (length(dup) > 0) {
    stop(sprintf("The probe-set CSV %s has %d probe ids more than one time, for example '%s'.",
                 path, length(dup), dup[1]))
  }
  probes <- readLines(file.path(fryDir, "alevin", "probe_quants_mat_cols.txt"))
  # join by name: the probe order of the matrix is not the order of the CSV
  idx <- match(probes, csv$probe_id)
  if (anyNA(idx)) {
    stop(sprintf("%d probes of the probe matrix are not in %s, for example '%s'.",
                 sum(is.na(idx)), path, probes[is.na(idx)][1]))
  }
  afg <- csv[idx, , drop = FALSE]
  rownames(afg) <- afg$probe_id
  afg
}

# For a multiplex (Flex) run, add the columns cell_barcode and sample_name to
# the barcode table `afc`. The row label of a multiplex run is
# <sample_name>_<cell_barcode>.
# If featureDump.txt has a sample_name column, use its CB and sample_name
# columns: line i describes row i. If featureDump.txt does not exist and the
# directory is from a multiplex run, give a warning and add no columns.
# Other outputs do not change.
.addFrySampleInfo <- function(fryDir, afc, quiet) {
  labels <- as.character(afc$barcodes)
  fd.file <- file.path(fryDir, "featureDump.txt")

  if (file.exists(fd.file)) {
    header <- strsplit(readLines(fd.file, n = 1), "\t", fixed = TRUE)[[1]]
    if (!("sample_name" %in% header)) {
      return(afc)
    }
    .say(quiet, "Reading the sample of each cell from ", fd.file)
    # read only CB and sample_name, as strings (a sample name such as "1e5" stays a string)
    col.classes <- ifelse(header %in% c("CB", "sample_name"), "character", "NULL")
    fd <- read.table(fd.file, sep = "\t", header = TRUE, quote = "",
                     comment.char = "", na.strings = character(0),
                     check.names = FALSE, colClasses = col.classes)
    if (nrow(fd) != length(labels)) {
      stop(sprintf("%s has %d cells, but the count matrix has %d rows.",
                   fd.file, nrow(fd), length(labels)))
    }
    expected <- paste0(fd$sample_name, "_", fd$CB)
    bad <- which(expected != labels)
    if (length(bad) > 0) {
      stop(sprintf(paste0("%s does not have the same row order as the rows file of ",
                          "the count matrix: %d rows are different. The first is row %d: ",
                          "'%s' in featureDump.txt, '%s' in the rows file."),
                   fd.file, length(bad), bad[1], expected[bad[1]], labels[bad[1]]))
    }
    afc$cell_barcode <- fd$CB
    afc$sample_name <- fd$sample_name
    return(afc)
  }

  multiplex.files <- file.path(fryDir, c("collation_manifest.bin", "sample_info.json"))
  if (any(file.exists(multiplex.files))) {
    warning(sprintf(paste0("%s is missing in a multiplex run, so sample_name and ",
                           "cell_barcode are not added to the colData."), fd.file))
  }

  afc
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

.say <- function(quiet, ...) {
  if (!quiet) {
    message(...)
  }
}