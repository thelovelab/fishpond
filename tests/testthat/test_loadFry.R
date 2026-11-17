context("loadFry")

# Notes on how to generate test data sets for use here can be round in the
# R/test-helpers.R file

test_that("Reading in Alevin-fry USA count matrix works", {
  dat <- fishpond:::readExampleFryData("fry-usa-basic")
  expect_true(dir.exists(dat$parent_dir))

  # read default quantification (S + A)
  sce <- loadFry(dat$parent_dir)
  expect_equal(nrow(sce), length(dat$genes))
  expect_equal(ncol(sce), length(dat$barcodes))
  expect_equal(sort(SummarizedExperiment::assayNames(sce)), sort(c("counts", "unspliced")))

  cts <- SummarizedExperiment::assay(sce, "counts")
  expect_s4_class(cts, "dgCMatrix")

  # Add the spliced and ambiguous reads manualy and convert into a matrix
  M <- local({
    m <- dat$matrix[, dat$usa$S, drop = FALSE] + dat$matrix[, dat$usa$A, drop = FALSE]
    dimnames(m) <- list(dat$barcodes, dat$genes)
    Matrix::Matrix(t(m), sparse = TRUE)
  })

  expect_equal(cts, M)
  
})

# The paths of the test directories in inst/extdata/alevin/example-quants
fryExample <- function(...) {
  system.file("extdata", "alevin", "example-quants", ..., package = "fishpond",
              mustWork = TRUE)
}

# The base R content of a SingleCellExperiment, for a comparison with saved values
snapSCE <- function(sce) {
  list(assays = lapply(SummarizedExperiment::assays(sce), as.matrix),
       rowData = as.data.frame(SummarizedExperiment::rowData(sce)),
       colData = as.data.frame(SummarizedExperiment::colData(sce)),
       dimnames = dimnames(sce), metadata = S4Vectors::metadata(sce))
}

test_that("loadFry gives the same output as before for directories that are not multiplex", {
  # loadFry_regression.rds was made with the loadFry of fishpond 2.19.2
  # (commit 111f6c0), with the same directories and arguments as below.
  expected <- readRDS(test_path("loadFry_regression.rds"))
  dirs <- c(fry_usa_basic = fryExample("fry-usa-basic"),
            legacy_usa = fryExample("flex", "legacy_usa"))
  cases <- list(
    scRNA = list(outputFormat = "scRNA"), snRNA = list(outputFormat = "snRNA"),
    all = list(outputFormat = "all"), `S+A` = list(outputFormat = "S+A"),
    `U+S+A` = list(outputFormat = "U+S+A"), velocity = list(outputFormat = "velocity"),
    scVelo = list(outputFormat = "scVelo"), raw = list(outputFormat = "raw"),
    custom = list(outputFormat = list(spliced = c("S", "A"), unspliced = "U")),
    nonzero = list(outputFormat = "scRNA", nonzero = TRUE))
  for (d in names(dirs)) {
    for (k in names(cases)) {
      sce <- do.call(loadFry, c(list(fryDir = dirs[[d]], quiet = TRUE), cases[[k]]))
      expect_identical(snapSCE(sce), expected[[d]][[k]], info = paste(d, k))
    }
  }
  # a single-sample run: featureDump.txt has no sample_name column
  single <- fryExample("flex", "flex_single")
  expect_identical(snapSCE(loadFry(single, quiet = TRUE)), expected$flex_single$default)
  expect_identical(snapSCE(loadFry(single, nonzero = TRUE, quiet = TRUE)),
                   expected$flex_single$nonzero)
})

# Copy a test directory into a new temporary directory, so that a test can change it
copyFryExample <- function(name) {
  tmp <- tempfile("fry")
  dir.create(tmp)
  file.copy(fryExample("flex", name), tmp, recursive = TRUE)
  file.path(tmp, name)
}

readLinesTab <- function(path) {
  do.call(rbind, strsplit(readLines(path), "\t", fixed = TRUE))
}

test_that("loadFry adds the cell barcode and the sample of each cell of a Flex run", {
  d <- fryExample("flex", "flex_gene")
  sce <- loadFry(d, quiet = TRUE)
  rows <- readLines(file.path(d, "alevin", "quants_mat_rows.txt"))
  cd <- SummarizedExperiment::colData(sce)

  # the column names are the row labels <sample_name>_<cell_barcode>
  expect_identical(colnames(sce), rows)
  expect_identical(colnames(cd), c("barcodes", "cell_barcode", "sample_name"))
  expect_identical(cd$barcodes, rows)
  expect_identical(paste0(cd$sample_name, "_", cd$cell_barcode), rows)

  # sample names keep "_", and a name such as "1e5" stays a string
  expect_identical(sort(unique(cd$sample_name)), c("1e5", "S_1", "S_1_A"))
  info <- jsonlite::fromJSON(file.path(d, "sample_info.json"))$samples
  expect_identical(as.vector(table(cd$sample_name)[info$name]), as.integer(info$num_cells))

  # the same bare barcode in two samples gives two cells
  dup <- duplicated(cd$cell_barcode) | duplicated(cd$cell_barcode, fromLast = TRUE)
  expect_equal(sum(dup), 2)
  expect_false(anyDuplicated(colnames(sce)) > 0)

  # the counts are the transposed matrix of quants_mat.mtx
  m <- Matrix::readMM(file.path(d, "alevin", "quants_mat.mtx"))
  expect_equal(unname(as.matrix(SummarizedExperiment::assay(sce, "counts"))),
               unname(as.matrix(Matrix::t(m))))
})

test_that("loadFry left-joins gene_id_to_name.tsv onto the genes of a Flex run", {
  d <- fryExample("flex", "flex_gene")
  sce <- loadFry(d, quiet = TRUE)
  cols <- readLines(file.path(d, "alevin", "quants_mat_cols.txt"))
  g2n <- readLinesTab(file.path(d, "gene_id_to_name.tsv"))
  expect_gt(nrow(g2n), length(cols))
  rd <- SummarizedExperiment::rowData(sce)
  expect_identical(rownames(sce), cols)
  expect_identical(colnames(rd), c("gene_ids", "gene_names"))
  expect_identical(rd$gene_ids, cols)
  expect_identical(rd$gene_names, g2n[match(cols, g2n[, 1]), 2])
})

test_that("loadFry gives NA as the name of a gene that is not in gene_id_to_name.tsv", {
  d <- copyFryExample("flex_gene")
  path <- file.path(d, "gene_id_to_name.tsv")
  lines <- readLines(path)
  writeLines(lines[!startsWith(lines, "GENE000030\t")], path)
  rd <- SummarizedExperiment::rowData(loadFry(d, quiet = TRUE))
  expect_identical(is.na(rd$gene_names), rd$gene_ids == "GENE000030")
})

test_that("loadFry stops for a bad line or a gene id two times in gene_id_to_name.tsv", {
  d <- copyFryExample("flex_gene")
  path <- file.path(d, "gene_id_to_name.tsv")
  lines <- readLines(path)
  writeLines(c(lines, "GENE000099"), path)
  expect_error(loadFry(d, quiet = TRUE), "did not have 2 elements")
  writeLines(c(lines, "GENE000099\tA\tB"), path)
  expect_error(loadFry(d, quiet = TRUE), "did not have 2 elements")
  writeLines(c(lines, "GENE000010\tOther"), path)
  # rownames<- gives a warning before the error
  expect_error(suppressWarnings(loadFry(d, quiet = TRUE)), "duplicate 'row.names'")
})

test_that("loadFry stops if featureDump.txt and the rows file do not align", {
  expect_error(loadFry(fryExample("flex", "flex_bad_featuredump"), quiet = TRUE),
               "same row order")
  d <- copyFryExample("flex_gene")
  path <- file.path(d, "featureDump.txt")
  lines <- readLines(path)
  writeLines(lines[-length(lines)], path)
  expect_error(loadFry(d, quiet = TRUE), "has 9 cells, but the count matrix has 10 rows")
  lines[1] <- sub("^CB", "Barcode", lines[1])
  writeLines(lines, path)
  expect_error(loadFry(d, quiet = TRUE), "same row order")
})

test_that("loadFry splits the row labels on the last '_' without featureDump.txt", {
  full <- loadFry(fryExample("flex", "flex_gene"), quiet = TRUE)
  sce <- loadFry(fryExample("flex", "flex_no_featuredump"), quiet = TRUE)
  cd <- SummarizedExperiment::colData(sce)
  expect_identical(colnames(cd), c("barcodes", "cell_barcode", "sample_name"))
  expect_identical(cd$cell_barcode, full$cell_barcode)
  expect_identical(cd$sample_name, full$sample_name)
})

test_that("loadFry stops if a multiplex row label has no '_' and barcode", {
  d <- copyFryExample("flex_no_featuredump")
  path <- file.path(d, "alevin", "quants_mat_rows.txt")
  lines <- readLines(path)
  lines[1] <- gsub("_", "-", lines[1])
  writeLines(lines, path)
  expect_error(loadFry(d, quiet = TRUE), "Cannot split")
})

test_that("loadFry does not split the row labels of a run that is not multiplex", {
  d <- copyFryExample("flex_no_featuredump")
  file.remove(file.path(d, "sample_info.json"))
  sce <- loadFry(d, quiet = TRUE)
  expect_identical(colnames(SummarizedExperiment::colData(sce)), "barcodes")
})

test_that("loadFry with nonzero = TRUE removes the Flex cell with no counts", {
  sce <- loadFry(fryExample("flex", "flex_gene"), nonzero = TRUE, quiet = TRUE)
  expect_equal(ncol(sce), 9)
  expect_true(all(Matrix::colSums(SummarizedExperiment::assay(sce, "counts")) > 0))
})

# The probe-set CSV rows (without the "#" lines) as a data.frame of strings
readProbeSet <- function(path) {
  lines <- readLines(path)
  read.csv(text = lines[!startsWith(lines, "#")], colClasses = "character")
}

test_that("loadFry reads the probe matrix with the probe_set.csv file in fryDir", {
  d <- fryExample("flex", "flex_probe")
  sce <- loadFry(d, featureLevel = "probe", quiet = TRUE)
  cols <- readLines(file.path(d, "alevin", "probe_quants_mat_cols.txt"))
  csv <- readProbeSet(file.path(d, "probe_set.csv"))
  rd <- SummarizedExperiment::rowData(sce)

  # the rows are the probes, in the column order of the matrix
  expect_identical(rownames(sce), cols)
  expect_identical(colnames(rd), c("probe_id", "gene_ids", "gene_names", "region"))
  expect_identical(rd$probe_id, cols)
  i <- match(cols, csv$probe_id)
  expect_identical(rd$gene_ids, csv$gene_id[i])
  expect_identical(rd$gene_names, csv$gene_name[i])
  expect_identical(rd$region, csv$region[i])
  # the CSV also has probes that are not in the matrix (included = FALSE)
  expect_gt(nrow(csv), length(cols))

  # the counts are the transposed matrix of probe_quants_mat.mtx
  m <- Matrix::readMM(file.path(d, "alevin", "probe_quants_mat.mtx"))
  expect_equal(unname(as.matrix(SummarizedExperiment::assay(sce, "counts"))),
               unname(as.matrix(Matrix::t(m))))
})

test_that("loadFry gives the same cells at probe level and at gene level", {
  d <- fryExample("flex", "flex_probe")
  probe <- loadFry(d, featureLevel = "probe", quiet = TRUE)
  gene <- loadFry(d, quiet = TRUE)
  expect_identical(SummarizedExperiment::colData(probe), SummarizedExperiment::colData(gene))
  expect_identical(colnames(SummarizedExperiment::rowData(gene)), c("gene_ids", "gene_names"))
})

test_that("loadFry probe counts summed by gene_ids are the gene counts", {
  d <- fryExample("flex", "flex_probe")
  probe <- loadFry(d, featureLevel = "probe", quiet = TRUE)
  gene <- loadFry(d, quiet = TRUE)
  g <- SummarizedExperiment::rowData(probe)$gene_ids
  p2g <- Matrix::sparseMatrix(i = match(g, rownames(gene)), j = seq_along(g), x = 1,
                              dims = c(nrow(gene), length(g)))
  summed <- as.matrix(p2g %*% SummarizedExperiment::assay(probe, "counts"))
  expect_equal(unname(summed), unname(as.matrix(SummarizedExperiment::assay(gene, "counts"))))
})

test_that("loadFry uses gene_id_to_name.tsv for a probe-set CSV without gene_name", {
  d <- fryExample("flex", "flex_probe")
  ref <- loadFry(d, featureLevel = "probe", quiet = TRUE)
  sce <- loadFry(d, featureLevel = "probe", probeRef = fryExample("flex", "probe_set_v1.csv"),
                 quiet = TRUE)
  rd <- SummarizedExperiment::rowData(sce)
  expect_identical(colnames(rd), c("probe_id", "gene_ids", "gene_names"))
  expect_identical(rd$gene_ids, SummarizedExperiment::rowData(ref)$gene_ids)
  expect_identical(rd$gene_names, SummarizedExperiment::rowData(ref)$gene_names)
})

test_that("loadFry gives NA for an empty gene_name or region in the probe-set CSV", {
  d <- fryExample("flex", "flex_probe")
  path <- tempfile(fileext = ".csv")
  lines <- readLines(file.path(d, "probe_set.csv"))
  i <- grep("^GENE000050,", lines)
  lines[i] <- sub(",spliced,Epsilon$", ",,", lines[i])
  writeLines(lines, path)
  rd <- SummarizedExperiment::rowData(loadFry(d, featureLevel = "probe", probeRef = path,
                                              quiet = TRUE))
  k <- rd$gene_ids == "GENE000050"
  expect_true(all(is.na(rd$gene_names[k])) && all(is.na(rd$region[k])))
  expect_false(anyNA(rd$gene_names[!k]))
})

test_that("loadFry stops if it finds no probe-set CSV", {
  d <- copyFryExample("flex_probe")
  file.remove(file.path(d, "probe_set.csv"))
  expect_error(loadFry(d, featureLevel = "probe", quiet = TRUE),
               file.path(d, "probe_set.csv"), fixed = TRUE)
  expect_error(loadFry(fryExample("flex", "flex_probe"), featureLevel = "probe",
                       probeRef = file.path(tempdir(), "missing.csv"), quiet = TRUE),
               "missing.csv does not exist")
  # the gene level does not need the CSV
  expect_identical(dim(loadFry(d, quiet = TRUE)), c(6L, 10L))
})

test_that("loadFry stops for a probe-set CSV that does not fit the probe matrix", {
  d <- fryExample("flex", "flex_probe")
  lines <- readLines(file.path(d, "probe_set.csv"))
  path <- tempfile(fileext = ".csv")
  writeLines(lines[!grepl("^GENE000010,", lines)], path)
  expect_error(loadFry(d, featureLevel = "probe", probeRef = path, quiet = TRUE),
               "1 probes of the probe matrix are not in")
  writeLines(c(lines, lines[grep("^GENE000010,", lines)]), path)
  expect_error(loadFry(d, featureLevel = "probe", probeRef = path, quiet = TRUE),
               "more than one time")
  writeLines(sub("probe_id", "probe", lines), path)
  expect_error(loadFry(d, featureLevel = "probe", probeRef = path, quiet = TRUE),
               "does not have the columns probe_id")
  writeLines(sub(",spliced,Alpha$", ",S,Alpha", lines), path)
  expect_error(loadFry(d, featureLevel = "probe", probeRef = path, quiet = TRUE),
               "must have only 'spliced'")
})

test_that("loadFry with nonzero = TRUE at probe level removes the cell with no counts", {
  sce <- loadFry(fryExample("flex", "flex_probe"), featureLevel = "probe",
                 nonzero = TRUE, quiet = TRUE)
  expect_equal(dim(sce), c(7, 9))
})

test_that("loadFry stops at probe level without a probe matrix or for a USA mode run", {
  for (name in c("flex_gene", "legacy_usa")) {
    expect_error(loadFry(fryExample("flex", name), featureLevel = "probe",
                         probeRef = fryExample("flex", "flex_probe", "probe_set.csv"),
                         quiet = TRUE),
                 "probe_mtx", info = name)
  }
  expect_error(loadFry(fryExample("flex", "flex_probe"), featureLevel = "transcript",
                       quiet = TRUE), "should be one of")
  d <- copyFryExample("flex_probe")
  qpath <- file.path(d, "quant.json")
  writeLines(sub('"usa_mode": false', '"usa_mode": true', readLines(qpath)), qpath)
  expect_error(loadFry(d, featureLevel = "probe", quiet = TRUE), "USA mode")
})
