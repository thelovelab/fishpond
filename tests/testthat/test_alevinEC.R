context("alevinEC")
library(data.table)
library(fishpond)

test_that("Importing transcript compatibility counts from alevin output works",{
  
  # import test data
  dir <- system.file("extdata", package="tximportData")
  files <- c(file.path(dir,"alevin/mouse1_unst_50/alevin/bfh.txt"),
             file.path(dir,"alevin/mouse1_LPS2_50/alevin/bfh.txt"))
  file.exists(files)
  
  tx2gene <- data.table::fread(file.path(dir, "tx2gene_alevin.tsv"),
                               header = FALSE)
  colnames(tx2gene) <- c("isoform_id", "gene_id")
  
  # use alevinEC
  EC_mat <- alevinEC(paths = files, 
                     tx2gene = tx2gene, 
                     multigene = FALSE, 
                     quiet = TRUE)
  
  # test output type
  expect_true(validObject(EC_mat))
  expect_true(is(EC_mat, "dgCMatrix"))
  expect_true(is(EC_mat, "Matrix"))
  
  # test output dimensions
  expect_true(ncol(EC_mat) == 100) # 2 files with 50 cells
  
  # test barcode names / colnames
  expect_true(is(colnames(EC_mat), "character"))
  barcodes1 <- data.table::fread(file.path(dir,"alevin/mouse1_unst_50/alevin/quants_mat_rows.txt"),
                                 header = FALSE)[[1]]
  barcodes2 <- data.table::fread(file.path(dir,"alevin/mouse1_LPS2_50/alevin/quants_mat_rows.txt"),
                                 header = FALSE)[[1]]
  expect_true(all(colnames(EC_mat) == c(barcodes1, barcodes2)))
  
  # test rownames type
  expect_true(is(rownames(EC_mat), "character"))
  
  # test quiet argument
  expect_silent(alevinEC(paths = files, 
                         tx2gene = tx2gene, 
                         multigene = FALSE, 
                         quiet = TRUE))
  
  expect_error(expect_silent(alevinEC(paths = files, 
                                      tx2gene = tx2gene, 
                                      multigene = FALSE, 
                                      quiet = FALSE)))
  # test multigene argument
  EC_mat_multi <- alevinEC(paths = files, 
                           tx2gene = tx2gene, 
                           multigene = TRUE, 
                           quiet = TRUE)
  
  ## test output type
  expect_true(validObject(EC_mat_multi))
  expect_true(is(EC_mat_multi, "dgCMatrix"))
  expect_true(is(EC_mat_multi, "Matrix"))
  
  ## test output dimensions
  expect_true(ncol(EC_mat_multi) == 100) # 2 files with 50 cells
  expect_true(nrow(EC_mat_multi) >= nrow(EC_mat)) # multigene=FALSE removes ECs
  
  ## test barcode names
  expect_true(all(colnames(EC_mat_multi) == c(barcodes1, barcodes2)))
  
  ## test rownames type
  expect_true(is(rownames(EC_mat_multi), "character"))
  
  # faulty path specification
  files2 <- c(file.path(dir,"alevin/mouse1_unst_50/alevin/wrong.txt"),
              file.path(dir,"alevin/mouse1_LPS2_50/alevin/wrangAgain.txt"))
  
  expect_error(alevinEC(paths = files2, 
                        tx2gene = tx2gene, 
                        multigene = FALSE, 
                        quiet = TRUE))
  
  # faulty tx2gene colnames
  colnames(tx2gene) <- c("transcripts", "genes")
  expect_error(alevinEC(paths = files, 
                        tx2gene = tx2gene, 
                        multigene = FALSE, 
                        quiet = TRUE))
  
  tx2gene <- as.data.frame.list(tx2gene)
  colnames(tx2gene) <- NULL
  expect_error(alevinEC(paths = files, 
                        tx2gene = tx2gene, 
                        multigene = FALSE, 
                        quiet = TRUE))
})