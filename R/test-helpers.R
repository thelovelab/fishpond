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

  m <- Matrix::Matrix(x$matrix, sparse = TRUE)
  m <- as(m, "dgCMatrix")
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
