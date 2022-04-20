#' Fishpond: downstream methods and tools for expression data
#'
#' This package provides statistical methods and other tools for
#' working with Salmon and Alevin quantification of RNA-seq data.
#' Fishpond contains the Swish non-parametric method for
#' detecting differential transcript expression (DTE). Swish can
#' also be used to detect differential gene expresion (DGE),
#' to perform allelic analysis, or to assess changes in isoform
#' proportions.
#'
#' The main Swish functions are:
#' \itemize{
#' \item \code{\link{scaleInfReps}} - scaling transcript or gene expression data
#' \item \code{\link{labelKeep}} - labelling which features have sufficient counts
#' \item \code{\link{swish}} - perform non-parametric differential analysis
#' \item Plots, e.g., \code{\link{plotMASwish}}, \code{\link{plotInfReps}}
#' }
#' 
#' All software-related questions should be posted to the Bioconductor Support Site:
#' 
#' \url{https://support.bioconductor.org}
#'
#' The code can be viewed at the GitHub repository,
#' which also lists the contributor code of conduct:
#'
#' \url{https://github.com/mikelove/fishpond}
#' 
#' @references
#'
#' Swish method:
#' 
#' Zhu, A., Srivastava, A., Ibrahim, J.G., Patro, R., Love, M.I. (2019)
#' Nonparametric expression analysis using inferential replicate counts.
#' Nucleic Acids Research.
#' \url{https://doi.org/10.1093/nar/gkz622}
#'
#' Compression, makeInfReps and splitSwish:
#'
#' Van Buren, S., Sarkar, H., Srivastava, A., Rashid, N.U.,
#' Patro, R., Love, M.I. (2020)
#' Compression of quantification uncertainty for scRNA-seq counts.
#' bioRxiv.
#' \url{https://doi.org/10.1101/2020.07.06.189639}
#'
#' @importFrom graphics plot points segments rect abline title axis legend
#' @importFrom stats median quantile qnorm rpois runif rnbinom var cor
#' @importFrom utils head tail capture.output read.table write.table read.csv
#' @importFrom methods is as slot
#' @importFrom SummarizedExperiment SummarizedExperiment
#' assayNames assayNames<- assay assay<- assays assays<-
#' colData colData<- mcols mcols<- rowRanges rowRanges<-
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom IRanges IntegerList
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' @importFrom GenomicRanges start end strand width
#' start<- end<- strand<- resize flank sort seqnames
#' @importFrom gtools permutations
#' @importFrom Matrix rowSums readMM t colSums
#' @importFrom matrixStats rowRanks rowMedians rowVars rowQuantiles
#' @importFrom svMisc progress
#' @importFrom jsonlite fromJSON
#' 
#' @docType package
#' @name fishpond-package
#' @aliases fishpond-package
#' @keywords package
NULL
