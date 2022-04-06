#' Construct a sparse matrix of transcript compatibility counts from alevin 
#' output
#'
#' Constructs a UMI count matrix with equivalence class identifiers
#' in the rows and barcode identifiers in the columns. The count matrix is 
#' generated from one or multiple `bfh.txt` files that have been created by 
#' running alevin-fry with the --dumpBFH flag. Alevin-fry - 
#' \url{https://doi.org/10.1186/s13059-019-1670-y}
#'
#' @rdname alevinEC
#' 
#' @param paths `Charachter` or `character vector`, path specifying the 
#' location of the `bfh.txt` files generated with alevin-fry.
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
#' count matrix with equivalence class identifiers in the rows and barcode 
#' identifiers in the columns. The second element `tx2gene_matched` allows for 
#' linking the equivalence class identifiers to their respective transcripts 
#' and genes.
#' 
#' @section Details:
#' The resulting count matrix uses equivalence class identifiers as rownames.
#' These can be linked to respective transcripts and genes using the 
#' `tx2gene_matched` element of the output. Specifically, if the equivalence 
#' class identifier reads 1|2|8, then the equivalence class is compatible with
#' the transcripts and their respective genes in rows 1, 2 and 8 of 
#' `tx2gene_matched`.
#' 
#' @importFrom Matrix sparseMatrix
#' @export
alevinEC <- function(paths, tx2gene, multigene = FALSE, ignoreTxVersion = FALSE,
                     ignoreAfterBar = FALSE, quiet = FALSE){

  if (!requireNamespace("data.table", quietly=TRUE)) {
    stop("alevinEC() requires CRAN package data.table")
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
  
  # Read and wrangle link between genes and transcripts
  n_tx <- data.table::fread(paths[1], nrows=1,
                sep = " ",quote = "", header = FALSE)$V1
  tname <- data.table::fread(paths[1], skip = 3, nrows = n_tx, 
                 sep = " ", quote = "", header = FALSE)$V1
  if (ignoreTxVersion) {
    tx2gene$isoform_id <- sub("\\..*", "", tx2gene$isoform_id)
    tname <- sub("\\..*", "", tname)
  } else if (ignoreAfterBar) {
    tx2gene$isoform_id <- sub("\\|.*", "", tx2gene$isoform_id)
    tname <- sub("\\|.*", "", tname)
  }
  
  if (!any(tname %in% tx2gene$isoform_id)) {
    txFromFile <- paste0("Example IDs (file): [", paste(head(tname,3),collapse=", "),", ...]")
    txFromTable <- paste0("Example IDs (tx2gene): [", paste(head(tx2gene$isoform_id,3),collapse=", "),", ...]")
    stop(paste0("
None of the transcripts in the quantification files are present
in the first column of tx2gene. Check to see that you are using
the same annotation for both.\n\n",txFromFile,"\n\n",txFromTable,
                "\n\n  This can sometimes (not always) be fixed using 'ignoreTxVersion' or 'ignoreAfterBar'.\n\n"))
  }
  
  # matching of tname and tx2gene
  tx2gene <- tx2gene[match(tname, tx2gene$isoform_id),]
  
  ntxmissing <- sum(!tname %in% tx2gene$isoform_id)
  if (ntxmissing > 0) message("transcripts missing from tx2gene: ", ntxmissing)
  
  # Reading and wrangling alevin output
  # message("Reading and wrangling alevin output")
  sample_list <- vector("list", length = length(paths))
  for (i in seq_along(paths)) {
    if (!quiet) svMisc::progress(i, max.value=length(paths),
                                 init=(i==1), gui=FALSE)
    sample_list[[i]] <- readBFH(file = paths[i],
                                tx2gene = tx2gene,
                                multigene = multigene)
  }
  
  # Get unique equivalence classes to later construct sparse matrix
  unique_ecs <- unique(unlist(lapply(sample_list, "[[", 1), use.names = FALSE))
  
  # Constructing sparseMatrix output
  # message("Constructing sparseMatrix output \n")
  mat_list <- vector("list", length = length(sample_list))
  for (i in seq_along(sample_list)) {
    # if (!quiet) svMisc::progress(i, max.value=length(paths),
    #                               init=(i==1), gui=FALSE)
    # wrangle elements sample_list to use as input for sparseMatrix
    ## wrangle rows
    EC_names <- sample_list[[i]][[1]]
    times <- sample_list[[i]][[4]]
    RowIdx <- match(EC_names, unique_ecs)
    RowIdx <- rep(RowIdx, times=times)
    
    ## wrangle columns
    bc_id_to_name <- sample_list[[i]][[5]]
    ColIdx <- sample_list[[i]][[2]]
    
    ## wrangle entries
    X <- sample_list[[i]][[3]]
    
    # construct sparse matrix
    mat_list[[i]] <- sparseMatrix(
      i = RowIdx,
      j = ColIdx,
      x = X,
      dims = c(length(unique_ecs), length(bc_id_to_name)),
      dimnames = list(unique_ecs, bc_id_to_name)
    )
  }
  matrix_umis <- do.call(cbind, mat_list)
  
  rownames(tx2gene) <- seq_len(nrow(tx2gene))
  
  return(list(counts = matrix_umis, tx2gene_matched = tx2gene))
}


# readBFH: helper to read bfh.txt files, creating lists with 5 elements: 
# (i) the names of the ECs
# (ii) the indices of the barcodes in which each of the ECs are expressed
# (iii) the UMI expression level for each EC/barcode pair
# (iv) number of barcodes in which each EC is expressed (to construct sparseMatrix)
# (v) the names of the barcodes
readBFH <- function(file, tx2gene, multigene){
  
  numbers <- data.table::fread(file, nrows=3, sep = " ", quote = "", header = FALSE)
  n_tx <- numbers$V1[1]
  n_bc <- numbers$V1[2]
  
  # read barcode names
  bc_id_to_name <- data.table::fread(file, skip = 3+n_tx, nrows = n_bc, 
                         sep = " ", quote = "", header = FALSE)$V1
  
  # read ECs
  startread <- sum(n_tx, n_bc, 3) # first line with EC info
  ec_df <- data.table::fread(file, skip=startread, sep = " ", quote = "", header = FALSE)
  eccs <- strsplit(ec_df$V1,"\t",fixed=TRUE)
  
  # code strongly based on read_bfh() from 
  # https://github.com/k3yavi/vpolo/blob/master/vpolo/alevin/parser.py#L374-L444
  hlpList <- lapply(eccs, function(toks) {
    
    num_labels <- as.integer(toks[1])
    txps <- as.integer(toks[2:(num_labels+1)])+1
    genes <-  tx2gene$gene_id[txps]
    if(!multigene){
      if(all(is.na(genes))){
        return(NULL)
      }
      if(length(unique(genes))>1){
        return(NULL)
      }
    }
    
    idx <- num_labels+3
    num_bcs <- as.integer(toks[idx]) # number of cells in which this ECs 
    # is present in this sample
    
    bc_names <- vector(mode = "numeric", length=num_bcs)
    num_umis <- vector(mode = "numeric", length=num_bcs)
    
    for (bc in seq_len(num_bcs)) {
      idx <- idx+1
      bc_names[bc] <- as.integer(toks[idx])+1
      idx <- idx+1
      num_umi <- as.integer(toks[idx])
      num_umis[bc] <- num_umi
      idx <- idx+(2*num_umi)
    }
    
    times <- length(bc_names)
    return(list(txps, bc_names, num_umis, times))
  })
  hlpList <- hlpList[!sapply(hlpList, is.null)]
  
  # wrangle
  bc_names <- unlist(lapply(hlpList, `[[`, 2))
  num_umis <- unlist(lapply(hlpList, `[[`, 3))
  times <- unlist(lapply(hlpList, `[[`, 4))
  
  EC_names <- unlist(lapply(lapply(hlpList, `[[`, 1), paste, collapse="|"))
  
  return(list(EC_names, bc_names, num_umis, times, bc_id_to_name))
}
