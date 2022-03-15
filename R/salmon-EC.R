#' Construct a matrix of transcript compatibility counts from salmon output
#'
#' Constructs a count matrix with equivalence class identifiers
#' in the rows and barcode identifiers in the columns. The count matrix is 
#' generated from one or multiple `eq_classes.txt` files that have been created
#' by running salmon with the --dumpEq flag. Salmon - 
#' \url{https://doi.org/10.1038/nmeth.4197}
#'
#' @rdname salmonEC
#' 
#' @param paths `Charachter` or `character vector`, full path specifying the 
#' location of the `eq_classes.txt` files generated with salmon.
#' @param tx2gene A `dataframe` linking transcript identifiers to their 
#' corresponding gene identifiers. Transcript identifiers must be in a column 
#' `isoform_id`. Corresponding gene identifiers must be in a column `gene_id`.
#' @param multigene `Logical`, should equivalence classes that are compatible 
#' with multiple genes be retained? Default is `FALSE`, removing such ambiguous
#' equivalence classes.
#' @param quiet `Logical`, set `TRUE` to avoid displaying messages.
#'
#' @author Jeroen Gilis
#' 
#' @return a sparse count matrix 
#' 
#' @importFrom Matrix sparseMatrix sparseVector
#' @export
salmonEC <- function(paths, tx2gene, multigene = FALSE, quiet = FALSE){
  
  # get line number where the TCCs start (same for each file, get from first file)
  startread <- data.table::fread(paths[1], nrows=1, sep = " ", 
                     quote = "", header = FALSE)$V1 + 2
  
  # order tx2gene like transcripts in ECC files
  tx_lookup <- data.table::fread(paths[1], skip = 2, nrows = startread-2, 
                     sep = " ", quote = "", header = FALSE)$V1
  tx2gene <- tx2gene[match(tx_lookup, tx2gene$isoform_id),]
  
  # For each sample/cell generate a list of 2:
  # (1) names of the equivalence classes and (2) corresponding counts
  unique_ecs <- c()
  sample_list <- vector("list", length = length(paths))
  
  for (i in seq_along(paths)) {
    if (!quiet) svMisc::progress(i, max.value=length(paths),
                                 init=(i==1), gui=FALSE)
    sample_list[[i]] <- readEq(file = paths[i],
                               geneSet = tx2gene[,"gene_id"],
                               startread = startread,
                               multigene = multigene)
    unique_ecs <- c(unique_ecs, setdiff(sample_list[[i]][[1]],unique_ecs))
  }
  
  # Constructing sparseMatrix output
  # message("Constructing sparseMatrix output \n")
  vec_list <- vector("list", length = length(sample_list))
  for (i in seq_along(sample_list)) {
    # if (!quiet) svMisc::progress(i, max.value=length(paths),
    #                               init=(i==1), gui=FALSE)
    # wrangle elements sample_list to use as input for sparseMatrix
    ## wrangle rows
    EC_names <- sample_list[[i]][[1]]
    RowIdx <- match(EC_names, unique_ecs)
    
    ## wrangle entries
    X <- sample_list[[i]][[2]]
    
    vec_list[[i]] <- sparseVector(
      i = RowIdx,
      x = X,
      length = length(unique_ecs)
    )
  }
  ECC_mat <- do.call(sv.cbind, vec_list)
  return(ECC_mat)
}

# Helper function that reads eq_classes.txt files and converts them to lists
# with two elements:
# (i) the names of the equivalence classes
# (ii) the count for each equivalence class
# Additionally allows for remove equivalence classes that are compatible with
# multiple genes.
readEq <- function(file, geneSet, startread, multigene){
  
  ec_df <- data.table::fread(file, skip=startread, sep = " ", quote = "", header = FALSE)
  eccs <- strsplit(ec_df$V1,"\t",fixed=TRUE)
  
  # remove first and last -> keep ECCs and not their corresponding TXs and counts
  eccs_hlp <- lapply(eccs, function(x) x[-c(1,length(x))])
  
  if(!multigene){ # remove ECs corresponding to multiple genes
    hlp <- unlist(lapply(eccs_hlp, function(element) {
      length(unique(geneSet[as.integer(element)+1]))!=1
    }))
    eccs_hlp[hlp] <- NULL
    eccs[hlp] <- NULL
  }
  
  #ec_txs <- unlist(lapply(eccs_hlp, stringi::stri_c, collapse = "|"))
  ec_txs <- unlist(lapply(eccs_hlp, paste, collapse = "|"))
  ec_counts <- as.integer(vapply(eccs, data.table::last, FUN.VALUE = character(1)))
  
  return(list(ec_txs, ec_counts))
}

# Helper function to convert a list of sparseVectors to a sparseMatrix
# obtained from https://stackoverflow.com/questions/8843700, 
# response by petrelharp
sv.cbind <- function(...) {
  input <- lapply(list(...), as, "dsparseVector")
  thelength <- unique(sapply(input,length))
  stopifnot(length(thelength)==1)
  return(sparseMatrix(
    x = unlist(lapply(input,slot,"x")), 
    i = unlist(lapply(input,slot,"i")), 
    p = c(0, cumsum(sapply(input,function(x){length(x@x)}))),
    dims = c(thelength, length(input))
  ))
}


