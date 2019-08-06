#' readEDS - a utility function for quickly reading in Alevin's EDS format
#'
#' @param files pointer to the EDS file
#' @param num.genes number of genes
#' @param num.cells number of cells
#'
#' @return a list, the positions and the expression values of
#' the non-zero expressions in the Alevin EDS file.
#' \strong{NOTE:} the positions are 0-based, i.e. ranging from
#' [0, num.cells-1]
#'
#' @useDynLib fishpond
#' @export
readEDS <- function(files, num.genes, num.cells) {
  files <- path.expand(files)
  pos <- getPositions(num.genes, num.cells, files)
  exp <- getExpression(num.genes, num.cells, pos, files)
  list(pos=pos, exp=exp)
}
