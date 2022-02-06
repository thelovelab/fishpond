#' readEDS - a utility function for quickly reading in Alevin's EDS format
#'
#' @param numOfGenes number of genes
#' @param numOfOriginalCells number of cells
#' @param countMatFilename pointer to the EDS file, \code{quants_mat.gz}
#' @param tierImport whether the \code{countMatFilename} refers to a quants tier file
#'
#' @return a genes x cells sparse matrix, of the class \code{dgCMatrix}
#'
#' @useDynLib fishpond
#' @export
readEDS <- function(numOfGenes, numOfOriginalCells, countMatFilename, tierImport=FALSE) {
  getSparseMatrix(numOfGenes, numOfOriginalCells, path.expand(countMatFilename), tierImport)
}
