context("salmonEC")
library(data.table)
library(fishpond)

test_that("Importing transcript compatibility counts from salmon output works",{


  if (packageVersion("tximportData") >= "1.23.4") {
  
    # import test data
    dir <- system.file("extdata", package="tximportData")
    files <- c(file.path(dir,"salmon_ec/SRR7311351/aux_info/eq_classes.txt"),
               file.path(dir,"salmon_ec/SRR7311386/aux_info/eq_classes.txt"))
    file.exists(files)
    
    tx2gene <- read.csv2(file.path(dir, "salmon_ec/tx2gene_tasic.csv"),
                         header = TRUE, row.names = NULL)
    tx2gene <- tx2gene[,c(2,3)]
    colnames(tx2gene) <- c("isoform_id", "gene_id")

    slow <- FALSE
    
    # use salmonEC ~9 seconds
    EC_mat <- salmonEC(paths = files,
                       tx2gene = tx2gene,
                       ignoreTxVersion = TRUE,
                       ignoreAfterBar = FALSE,
                       multigene = FALSE,
                       quiet = FALSE)
    
    # test output type
    expect_true(validObject(EC_mat))
    expect_true(length(EC_mat) == 2) # list of two
    expect_true(is(EC_mat$counts, "dgCMatrix"))
    expect_true(is(EC_mat$counts, "Matrix"))
    expect_true(is(EC_mat$tx2gene_matched, "data.frame"))
    
    # test output dimensions
    expect_true(ncol(EC_mat$counts) == 2) # 2 cells
    
    # test colnames
    expect_true(is.null(colnames(EC_mat$counts)))
    
    # test rownames type
    expect_true(is(rownames(EC_mat$counts), "character"))
    
    # test multigene = FALSE and equivalence class to gene matching
    eccs <- strsplit(rownames(EC_mat$counts),"\\|",fixed=FALSE)
    nrGeneForEachEcc <- unlist(lapply(eccs, function(ecc){
      length(unique(EC_mat$tx2gene_matched$gene_id[as.integer(ecc)]))
    }))
    # each ecc must only be compatible with a single gene if multigene = FALSE
    expect_true(all(nrGeneForEachEcc == 1))
    
    if (slow) {
      # test quiet argument
      expect_silent(salmonEC(paths = files,
                             tx2gene = tx2gene,
                             ignoreTxVersion = TRUE,
                             ignoreAfterBar = FALSE,
                             multigene = FALSE,
                             quiet = TRUE))

      # test multigene argument
      EC_mat_multi <- salmonEC(paths = files,
                               tx2gene = tx2gene,
                               ignoreTxVersion = TRUE,
                               ignoreAfterBar = FALSE,
                               multigene = TRUE,
                               quiet = TRUE)
      
      # test output type
      expect_true(validObject(EC_mat_multi))
      expect_true(length(EC_mat_multi) == 2) # list of two
      expect_true(is(EC_mat_multi$counts, "dgCMatrix"))
      expect_true(is(EC_mat_multi$tx2gene_matched, "data.frame"))
      
      expect_true(nrow(EC_mat_multi$counts) >= nrow(EC_mat$counts)) # multigene=FALSE removes ECs
      
      # test patial but correct tx2gene list
      expect_message(EC_mat_2 <- salmonEC(paths = files,
                                          tx2gene = tx2gene[1:50000,],
                                          ignoreTxVersion = TRUE,
                                          ignoreAfterBar = FALSE,
                                          multigene = FALSE,
                                          quiet = TRUE))
      
      expect_true(validObject(EC_mat_2))
      expect_true(nrow(EC_mat$counts) >= nrow(EC_mat_2$counts)) # missing transcripts
    }
    
    # test faulty tx2gene list (e.g. when not ignoring tx version)
    expect_error(EC_mat_3 <- salmonEC(paths = files,
                                      tx2gene = tx2gene,
                                      ignoreTxVersion = FALSE,
                                      ignoreAfterBar = FALSE,
                                      multigene = FALSE,
                                      quiet = TRUE))

  }
})
