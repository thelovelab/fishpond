context("alevinEC")
library(data.table)
library(fishpond)

test_that("Importing transcript compatibility counts from alevin output works",{

  if (packageVersion("tximportData") >= "1.23.4") {
  
    # import test data
    dir <- system.file("extdata", package="tximportData")
    files <- c(file.path(dir,"alevin/mouse1_unst_50/alevin/bfh.txt"),
               file.path(dir,"alevin/mouse1_LPS2_50/alevin/bfh.txt"))
    file.exists(files)
    
    tx2gene <- data.table::fread(file.path(dir, "tx2gene_alevin.tsv"),
                                 header = FALSE)
    colnames(tx2gene) <- c("isoform_id", "gene_id")

    slow <- FALSE
    
    # use alevinEC ~10 seconds
    EC_mat <- alevinEC(paths = files, 
                       tx2gene = tx2gene, 
                       multigene = FALSE,
                       ignoreTxVersion = FALSE,
                       ignoreAfterBar = FALSE,
                       quiet = TRUE)
    
    # test output type
    expect_true(validObject(EC_mat))
    expect_true(is(EC_mat$counts, "dgCMatrix"))
    expect_true(is(EC_mat$counts, "Matrix"))
    expect_true(is(EC_mat$tx2gene_matched, "data.frame"))
    
    # test output dimensions
    expect_true(ncol(EC_mat$counts) == 100) # 2 files with 50 cells
    
    # test barcode names / colnames
    expect_true(is(colnames(EC_mat$counts), "character"))
    barcodes1 <- data.table::fread(file.path(dir,"alevin/mouse1_unst_50/alevin/quants_mat_rows.txt"),
                                   header = FALSE)[[1]]
    barcodes1 <- paste0(barcodes1, "_1")
    barcodes2 <- data.table::fread(file.path(dir,"alevin/mouse1_LPS2_50/alevin/quants_mat_rows.txt"),
                                   header = FALSE)[[1]]
    barcodes2 <- paste0(barcodes2, "_2")

    ### TODO: this test fails on my end (Mike, Aug 12 2022)
    expect_true(all(colnames(EC_mat$counts) == c(barcodes1, barcodes2)))
    
    # test rownames type
    expect_true(is(rownames(EC_mat$counts), "character"))

    if (slow) {
      # test quiet argument
      expect_silent(alevinEC(paths = files, 
                             tx2gene = tx2gene,
                             ignoreTxVersion = FALSE,
                             ignoreAfterBar = FALSE,
                             multigene = FALSE, 
                             quiet = TRUE))

      # test multigene argument
      EC_mat_multi <- alevinEC(paths = files, 
                               tx2gene = tx2gene, 
                               multigene = TRUE,
                               ignoreTxVersion = FALSE,
                               ignoreAfterBar = FALSE,
                               quiet = TRUE)
      
      ## test output type
      expect_true(validObject(EC_mat_multi))
      expect_true(is(EC_mat_multi$counts, "dgCMatrix"))
      expect_true(is(EC_mat_multi$counts, "Matrix"))
      expect_true(is(EC_mat_multi$tx2gene_matched, "data.frame"))
      
      ## test output dimensions
      expect_true(ncol(EC_mat_multi$counts) == 100) # 2 files with 50 cells
      expect_true(nrow(EC_mat_multi$counts) >= nrow(EC_mat$counts)) # multigene=FALSE removes ECs
      
      ## test barcode names
      expect_true(all(colnames(EC_mat_multi$counts) == c(barcodes1, barcodes2)))
      
      ## test rownames type
      expect_true(is(rownames(EC_mat_multi$counts), "character"))
    }
    
    # faulty path specification
    files2 <- c(file.path(dir,"alevin/mouse1_unst_50/alevin/wrong.txt"),
                file.path(dir,"alevin/mouse1_LPS2_50/alevin/wrangAgain.txt"))
    
    expect_error(alevinEC(paths = files2, 
                          tx2gene = tx2gene,
                          ignoreTxVersion = FALSE,
                          ignoreAfterBar = FALSE,
                          multigene = FALSE, 
                          quiet = TRUE))
    
    # faulty tx2gene colnames
    colnames(tx2gene) <- c("transcripts", "genes")
    expect_error(alevinEC(paths = files, 
                          tx2gene = tx2gene, 
                          multigene = FALSE,
                          ignoreTxVersion = FALSE,
                          ignoreAfterBar = FALSE,
                          quiet = TRUE))
    
    tx2gene <- as.data.frame.list(tx2gene)
    colnames(tx2gene) <- NULL
    expect_error(alevinEC(paths = files, 
                          tx2gene = tx2gene,
                          ignoreTxVersion = FALSE,
                          ignoreAfterBar = FALSE,
                          multigene = FALSE, 
                          quiet = TRUE))
    
    colnames(tx2gene) <- c("isoform_id", "gene_id")

    if (slow) {
      # test incomplete tx2gene annotation
      expect_message(EC_mat2 <- alevinEC(paths = files, 
                                         tx2gene = tx2gene[1:50000,], 
                                         multigene = FALSE,
                                         ignoreTxVersion = FALSE,
                                         ignoreAfterBar = FALSE,
                                         quiet = TRUE))
      
      expect_true(nrow(EC_mat$counts) >= nrow(EC_mat2$counts)) # incomplete annotation
    }
    
    # test faulty tx2gene list (e.g. when not ignoring tx version)
    tx2gene_faulty <- tx2gene
    tx2gene_faulty$isoform_id <- sub("\\..*", "", tx2gene_faulty$isoform_id)
    expect_error(EC_mat_3 <- alevinEC(paths = files,
                                      tx2gene = tx2gene_faulty,
                                      ignoreTxVersion = FALSE,
                                      ignoreAfterBar = FALSE,
                                      multigene = FALSE,
                                      quiet = TRUE))

  }
  
})
