#' readEDS - a utility function for quickly reading in Alevin's EDS format
#'
#' @param numOfGenes number of genes
#' @param numOfOriginalCells number of cells
#' @param countMatFilename pointer to the EDS file, \code{quants_mat.gz}
#'
#' @return a genes x cells sparse matrix, of the class \code{dgCMatrix}
#'
#' @useDynLib fishpond
#' @export
readEDS <- function(numOfGenes, numOfOriginalCells, countMatFilename) {
  mat <- getSparseMatrix(numOfGenes, numOfOriginalCells, path.expand(countMatFilename))
}
