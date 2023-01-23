#' Construct a sparse matrix of transcript compatibility counts from salmon 
#' output
#'
#' Constructs a count matrix with equivalence class identifiers
#' in the rows. The count matrix is generated from one or multiple 
#' `eq_classes.txt` files that have been created by running salmon with the 
#' --dumpEq flag. Salmon - \url{https://doi.org/10.1038/nmeth.4197}
#'
#' @rdname salmonEC
#' 
#' @param paths `Charachter` or `character vector`, path specifying the 
#' location of the `eq_classes.txt` files generated with salmon.
#' @param tx2gene A `dataframe` linking transcript identifiers to their 
#' corresponding gene identifiers. Transcript identifiers must be in a column
#' `isoform_id`. Corresponding gene identifiers must be in a column `gene_id`.
#' @param multigene `Logical`, should equivalence classes that are compatible 
#' with multiple genes be retained? Default is `FALSE`, removing such ambiguous
#' equivalence classes.
#' @param ignoreTxVersion logical, whether to split the isoform id on the '.' 
#' character to remove version information to facilitate matching with the
#' isoform id in `tx2gene` (default FALSE).
#' @param ignoreAfterBar logical, whether to split the isoform id on the '|' 
#' character to facilitate matching with the isoform id in `tx2gene` 
#' (default FALSE).
#' @param quiet `Logical`, set `TRUE` to avoid displaying messages.
#'
#' @author Jeroen Gilis
#' 
#' @return A list with two elements. The first element `counts` is a sparse 
#' count matrix with equivalence class identifiers in the rows. If multiple 
#' paths are specified, the columns are in the same order as the paths. The 
#' second element `tx2gene_matched` allows for linking those identifiers to
#' their respective transcripts and genes.
#' 
#' @section Details:
#' The resulting count matrix uses equivalence class identifiers as rownames.
#' These can be linked to respective transcripts and genes using the 
#' `tx2gene_matched` element of the output. Specifically, if the equivalence 
#' class identifier reads 1|2|8, then the equivalence class is compatible with
#' the transcripts and their respective genes in rows 1, 2 and 8 of 
#' `tx2gene_matched`.
#' 
#' @importFrom Matrix sparseMatrix sparseVector
#' @export
salmonEC <- function(paths, tx2gene, multigene = FALSE, ignoreTxVersion = FALSE, 
                      ignoreAfterBar = FALSE, quiet = FALSE){
  
  if (!requireNamespace("data.table", quietly=TRUE)) {
    stop("salmonEC() requires CRAN package data.table")
  }
  
  if(!all(file.exists(paths))){
    stop("The following paths do not exist: ",
         paste(paths[which(!file.exists(paths))], "\n")
    )
  }
  
  if(!isTRUE(all.equal(colnames(tx2gene),
                       c("isoform_id","gene_id")))){
    stop("tx2gene does not contain columns gene_id and isoform_id")
  }
  
  # get line number where the TCCs start (same for each file, get from first file)
  startread <- data.table::fread(paths[1], nrows=1, sep = " ",
                                 quote = "", header = FALSE)$V1 + 2
  
  # order tx2gene like transcripts in ECC files
  tx_lookup <- data.table::fread(paths[1], skip = 2, nrows = startread-2, 
                                 sep = " ", quote = "", header = FALSE)$V1
  
  if (ignoreTxVersion) {
    tx2gene$isoform_id <- sub("\\..*", "", tx2gene$isoform_id)
    tx_lookup <- sub("\\..*", "", tx_lookup)
  } else if (ignoreAfterBar) {
    tx2gene$isoform_id <- sub("\\|.*", "", tx2gene$isoform_id)
    tx_lookup <- sub("\\|.*", "", tx_lookup)
  }
  
  if (!any(tx_lookup %in% tx2gene$isoform_id)) {
    txFromFile <- paste0("Example IDs (file): [", 
                         paste(head(tx_lookup,3),collapse=", "),
                         ", ...]")
    txFromTable <- paste0("Example IDs (tx2gene): [", 
                          paste(head(tx2gene$isoform_id,3),
                                collapse=", "),", ...]")
    stop(paste0("
  None of the transcripts in the quantification files are present
  in the first column of tx2gene. Check to see that you are using
  the same annotation for both.\n\n",txFromFile,"\n\n",txFromTable,
                "\n\n  This can sometimes (not always) be fixed using 'ignoreTxVersion' or 'ignoreAfterBar'.\n\n"))
  }

  # matching of tx_lookup and tx2gene
  tx2gene <- tx2gene[match(tx_lookup, tx2gene$isoform_id),]
  
  ntxmissing <- sum(!tx_lookup %in% tx2gene$isoform_id)
  if (ntxmissing > 0) message("transcripts missing from tx2gene: ", ntxmissing)
  
  # For each sample/cell generate a list of 2:
  # (1) names of the equivalence classes and (2) corresponding counts
  unique_ecs <- c()
  sample_list <- vector("list", length = length(paths))
  
  for (i in seq_along(paths)) {
    if (!quiet) svMisc::progress(i, max.value=length(paths),
                                 init=(i==1), gui=FALSE)
    sample_list[[i]] <- readEq(file = paths[i],
                               geneSet = as.character(tx2gene[,"gene_id"]),
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
  rownames(ECC_mat) <- unique_ecs
  
  rownames(tx2gene) <- seq_len(nrow(tx2gene))
  
  return(list(counts = ECC_mat, tx2gene_matched = tx2gene))
}

# Helper function that reads eq_classes.txt files and converts them to lists
# with two elements:
# (i) the names of the equivalence classes
# (ii) the count for each equivalence class
# Additionally allows for remove equivalence classes that are compatible with
# multiple genes.
readEq <- function(file, geneSet, startread, multigene){
  
  ec_df <- data.table::fread(file, skip=startread, sep = " ", quote = "", 
                             header = FALSE)
  eccs <- strsplit(ec_df$V1,"\t",fixed=TRUE)
  
  # remove first and last -> keep ECCs and not their corresponding TXs and counts
  eccs_hlp <- lapply(eccs, function(x) as.integer(x[-c(1, length(x))])+1)
  # +1 to account for going from 0-indexing to 1-indexing
  
  if(!multigene){ # remove ECs corresponding to multiple genes
    hlp <- unlist(lapply(eccs_hlp, function(element) {
      if (any(is.na(geneSet[element]))) {
        return(TRUE)
      }
      else {
        return(length(unique(geneSet[element])) != 1)
      }
    }))
    eccs_hlp[hlp] <- NULL
    eccs[hlp] <- NULL
  }
  
  #ec_txs <- unlist(lapply(eccs_hlp, stringi::stri_c, collapse = "|"))
  ec_txs <- unlist(lapply(eccs_hlp, paste, collapse = "|"))
  ec_counts <- as.integer(vapply(eccs, data.table::last, 
                                 FUN.VALUE = character(1)))
  
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


