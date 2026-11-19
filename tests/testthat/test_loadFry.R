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

test_that("loadFry warns and adds no columns for a multiplex run without featureDump.txt", {
  d <- fryExample("flex", "flex_no_featuredump")
  expect_warning(sce <- loadFry(d, quiet = TRUE),
                 "featureDump.txt is missing in a multiplex run")
  expect_identical(colnames(SummarizedExperiment::colData(sce)), "barcodes")
})

test_that("loadFry is silent without featureDump.txt for a run that is not multiplex", {
  d <- copyFryExample("flex_no_featuredump")
  file.remove(file.path(d, "sample_info.json"))
  expect_silent(sce <- loadFry(d, quiet = TRUE))
  expect_identical(colnames(SummarizedExperiment::colData(sce)), "barcodes")
})

test_that("loadFry with nonzero = TRUE removes the Flex cell with no counts", {
  sce <- loadFry(fryExample("flex", "flex_gene"), nonzero = TRUE, quiet = TRUE)
  expect_equal(ncol(sce), 9)
  expect_true(all(Matrix::colSums(SummarizedExperiment::assay(sce, "counts")) > 0))
})

# The probe-set CSV rows (without the "#" lines), read as loadFry reads them
readProbeSet <- function(path) {
  lines <- readLines(path)
  read.csv(text = lines[!startsWith(lines, "#")], na.strings = c("NA", ""))
}

# The CSV rows of the matrix probes, in matrix order, with the probe ids as row names
probeRows <- function(csv, probes) {
  x <- csv[match(probes, csv$probe_id), , drop = FALSE]
  rownames(x) <- probes
  x
}

test_that("loadFry reads the probe matrix with the probe_set.csv file in fryDir", {
  d <- fryExample("flex", "flex_probe")
  sce <- loadFry(d, featureLevel = "probe", quiet = TRUE)
  cols <- readLines(file.path(d, "alevin", "probe_quants_mat_cols.txt"))
  csv <- readProbeSet(file.path(d, "probe_set.csv"))
  rd <- SummarizedExperiment::rowData(sce)

  # the rowData is the CSV as it is: all its columns and values, in matrix order
  expect_identical(rownames(sce), cols)
  expect_identical(colnames(rd), c("gene_id", "probe_seq", "probe_id", "included",
                                   "region", "gene_name"))
  expect_identical(as.data.frame(rd), probeRows(csv, cols))
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

test_that("loadFry probe counts summed by gene_id are the gene counts", {
  d <- fryExample("flex", "flex_probe")
  probe <- loadFry(d, featureLevel = "probe", quiet = TRUE)
  gene <- loadFry(d, quiet = TRUE)
  g <- SummarizedExperiment::rowData(probe)$gene_id
  p2g <- Matrix::sparseMatrix(i = match(g, rownames(gene)), j = seq_along(g), x = 1,
                              dims = c(nrow(gene), length(g)))
  summed <- as.matrix(p2g %*% SummarizedExperiment::assay(probe, "counts"))
  expect_equal(unname(summed), unname(as.matrix(SummarizedExperiment::assay(gene, "counts"))))
})

test_that("loadFry adds no gene names to a probe-set CSV without gene_name", {
  d <- fryExample("flex", "flex_probe")
  path <- fryExample("flex", "probe_set_v1.csv")
  sce <- loadFry(d, featureLevel = "probe", probeRef = path, quiet = TRUE)
  cols <- readLines(file.path(d, "alevin", "probe_quants_mat_cols.txt"))
  csv <- readProbeSet(path)
  expect_false("gene_name" %in% colnames(csv))
  expect_identical(as.data.frame(SummarizedExperiment::rowData(sce)), probeRows(csv, cols))
})

test_that("loadFry gives NA for empty CSV values and copies other values and columns as they are", {
  d <- fryExample("flex", "flex_probe")
  path <- tempfile(fileext = ".csv")
  lines <- readLines(file.path(d, "probe_set.csv"))
  lines <- sub(",spliced,Epsilon$", ",,", lines)
  lines <- sub(",spliced,Alpha$", ",S,Alpha", lines)
  body <- !startsWith(lines, "#")
  lines[body] <- paste0(lines[body], c(",extra", rep(",x", sum(body) - 1)))
  writeLines(lines, path)
  rd <- SummarizedExperiment::rowData(loadFry(d, featureLevel = "probe", probeRef = path,
                                              quiet = TRUE))
  expect_identical(colnames(rd)[7], "extra")
  expect_identical(rd$region[rd$gene_id == "GENE000050"], NA_character_)
  expect_identical(rd$gene_name[rd$gene_id == "GENE000050"], NA_character_)
  expect_identical(rd$region[rd$gene_id == "GENE000010"], "S")
  # TRUE/FALSE is logical, and an empty value is NA
  expect_identical(unique(rd$included), TRUE)
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
               "7 probes of the probe matrix are not in")
})

test_that("loadFry with nonzero = TRUE at probe level removes the cell with no counts", {
  sce <- loadFry(fryExample("flex", "flex_probe"), featureLevel = "probe",
                 nonzero = TRUE, quiet = TRUE)
  expect_equal(dim(sce), c(7, 9))
})

test_that("loadFry stops at probe level without a probe matrix", {
  for (name in c("flex_gene", "legacy_usa")) {
    expect_error(loadFry(fryExample("flex", name), featureLevel = "probe",
                         probeRef = fryExample("flex", "flex_probe", "probe_set.csv"),
                         quiet = TRUE),
                 "probe_mtx", info = name)
  }
  expect_error(loadFry(fryExample("flex", "flex_probe"), featureLevel = "transcript",
                       quiet = TRUE), "should be one of")
})

# The C reader of the count matrix (src/readFryMtx.c)

# Write a MatrixMarket file. 'entries' is a character vector of entry lines.
writeMtx <- function(entries, nr, nc, nnz = length(entries),
                     banner = "%%MatrixMarket matrix coordinate real general",
                     comments = "%", eol = "\n") {
  path <- tempfile(fileext = ".mtx")
  lines <- c(banner, comments, paste(nr, nc, nnz), entries)
  con <- file(path, "wb")
  writeLines(lines, con, sep = eol)
  close(con)
  path
}

# The reader that loadFry used before: a features x cells dgCMatrix
oldReader <- function(path) {
  Matrix::t(fishpond:::my_as_dgcmatrix(Matrix::readMM(path)))
}

# A random cells x features matrix as entry lines, in the order 'ord'
randomEntries <- function(nr, nc, density, values) {
  set.seed(nr * 1000 + nc)
  k <- sample(nr * nc, ceiling(density * nr * nc))
  data.frame(i = (k - 1) %% nr + 1, j = (k - 1) %/% nr + 1,
             x = sample(values, length(k), replace = TRUE))
}

entryLines <- function(e) paste(e$i, e$j, e$x)

test_that("the C reader gives the same matrix as readMM for entries in any order", {
  e <- randomEntries(30, 20, 0.3, c("1", "2", "7", "0.5", "1.25", "1234"))
  orders <- list(cellFeature = order(e$i, e$j), random = sample(nrow(e)),
                 reversed = rev(order(e$i, e$j)), featureCell = order(e$j, e$i))
  for (o in names(orders)) {
    f <- writeMtx(entryLines(e[orders[[o]], ]), 30, 20)
    m <- fishpond:::.readFryMtx(f)
    expect_s4_class(m, "dgCMatrix")
    expect_true(validObject(m))
    expect_identical(m, oldReader(f), info = o)
  }
})

test_that("the C reader reads a cell in two runs far apart, an empty cell and a full cell", {
  # cell 2 has two runs; the cell index goes down; cell 3 has no entries;
  # cell 4 has all features
  entries <- c("2 1 1", "2 2 2", "4 1 1", "4 2 1", "4 3 1", "4 4 1", "4 5 1",
               "1 3 5", "2 4 3", "2 5 4")
  f <- writeMtx(entries, 4, 5)
  m <- fishpond:::.readFryMtx(f)
  expect_identical(m, oldReader(f))
  expect_identical(diff(m@p), c(1L, 4L, 0L, 5L))
})

test_that("the C reader sorts a cell with features that do not increase", {
  entries <- c("1 5 1", "1 3 2", "1 1 3", "2 2 4", "2 1 5")
  f <- writeMtx(entries, 2, 5)
  expect_identical(fishpond:::.readFryMtx(f), oldReader(f))
})

test_that("the C reader reads small and special files as readMM does", {
  files <- list(
    oneByOne = writeMtx("1 1 3", 1, 1),
    noEntries = writeMtx(character(0), 3, 4),
    explicitZero = writeMtx(c("1 1 0", "2 2 1"), 2, 2),
    integerBanner = writeMtx(c("1 2 3", "2 1 4"), 2, 2,
                             banner = "%%MatrixMarket matrix coordinate integer general"),
    upperCaseBanner = writeMtx(c("1 2 3", "2 1 4"), 2, 2,
                               banner = "%%MatrixMarket MATRIX Coordinate Real General"),
    comments = writeMtx(c("1 2 3", "", "2 1 4  "), 2, 2, nnz = 2,
                        comments = c("% a comment", "%", "", paste0("%", strrep("x", 3000)))),
    trailingSpaces = writeMtx(c("1 2 3 ", "2\t1\t4"), 2, 2),
    crlf = writeMtx(c("1 2 3", "2 1 4.5"), 2, 2, eol = "\r\n"))
  sizeLineSpaces <- tempfile(fileext = ".mtx")
  writeLines(c("%%MatrixMarket matrix coordinate real general", "2 2 2   ", "1 2 3", "2 1 4"),
             sizeLineSpaces)
  files$sizeLineSpaces <- sizeLineSpaces
  for (k in names(files)) {
    expect_identical(fishpond:::.readFryMtx(files[[k]]), oldReader(files[[k]]), info = k)
  }
})

test_that("the C reader gives the same values as readMM for numbers in many formats", {
  set.seed(1)
  u <- c(runif(300), runif(300) * 10^sample(-8:8, 300, replace = TRUE),
         sample(1e6, 100), 2^52 + 1, 2^53 + 1)
  values <- c(sprintf("%.17g", u), sprintf("%.15g", u), sprintf("%.6g", u),
              sprintf("%e", u), sprintf("%E", u), sprintf("%.3f", u),
              ".25", "-0.5", "+2", "1e3", "1E-3", "0.1000000000000000055511151231257827",
              "12345678901234567", "1234567890123456789", "007", "-0")
  nc <- 50
  nr <- ceiling(length(values) / nc)
  k <- seq_along(values)
  f <- writeMtx(paste((k - 1) %/% nc + 1, (k - 1) %% nc + 1, values), nr, nc)
  expect_identical(fishpond:::.readFryMtx(f), oldReader(f))
})

test_that("the C reader gives the same matrix as readMM for the example files", {
  files <- list.files(system.file("extdata", "alevin", "example-quants", package = "fishpond"),
                      pattern = "\\.mtx$", recursive = TRUE, full.names = TRUE)
  expect_gt(length(files), 5)
  for (f in files) {
    expect_identical(fishpond:::.readFryMtx(f), oldReader(f), info = f)
  }
})

test_that("the C reader gives the same matrix for any block size", {
  e <- randomEntries(30, 20, 0.3, c("1", "2", "7", "0.5", "1.25", "1234"))
  f <- writeMtx(entryLines(e[sample(nrow(e)), ]), 30, 20,
                comments = c("% a comment", paste0("%", strrep("y", 2500))))
  for (block in c(1L, 2L, 7L, 1023L, 1024L, 1025L, 4194304L)) {
    expect_identical(fishpond:::.readFryMtx(f, block = block), oldReader(f), info = block)
  }
})

test_that("the C reader reads a last line without a line end", {
  f <- tempfile(fileext = ".mtx")
  writeBin(charToRaw("%%MatrixMarket matrix coordinate real general\n2 2 2\n1 2 3\n2 1 4"), f)
  expect_identical(fishpond:::.readFryMtx(f), oldReader(f))
})

test_that("the C reader ignores a long comment line, also with a NUL character", {
  f <- tempfile(fileext = ".mtx")
  long <- c(charToRaw(paste0("%", strrep("x", 1022))), as.raw(0), charToRaw(strrep("z", 2000)))
  writeBin(c(charToRaw("%%MatrixMarket matrix coordinate real general\n"), long,
             charToRaw("\n2 2 2\n1 2 3\n2 1 4\n")), f)
  ref <- writeMtx(c("1 2 3", "2 1 4"), 2, 2)
  for (block in c(1L, 1023L, 4194304L)) {
    expect_identical(fishpond:::.readFryMtx(f, block = block), oldReader(ref), info = block)
  }
})

test_that("the C reader reads gzip, bzip2 and xz compressed files as readMM does", {
  e <- randomEntries(30, 20, 0.3, c("1", "2", "7", "0.5"))
  lines <- c("%%MatrixMarket matrix coordinate real general", "%",
             paste(30, 20, nrow(e)), entryLines(e[sample(nrow(e)), ]))
  for (type in c("gzip", "bzip2", "xz")) {
    f <- tempfile(fileext = ".mtx")
    con <- switch(type, gzip = gzfile(f, "wb"), bzip2 = bzfile(f, "wb"), xz = xzfile(f, "wb"))
    writeLines(lines, con)
    close(con)
    expect_identical(fishpond:::.readFryMtx(f), oldReader(f), info = type)
  }
})

test_that("the C reader stops with a clear message for a file that is not correct", {
  rd <- function(...) fishpond:::.readFryMtx(writeMtx(...))
  expect_error(rd(c("1 1 1", "2 3 1", "2 1 1", "2 3 2"), 2, 3),
               "more than one entry for row 2, column 3")
  expect_error(rd(c("1 1 1", "3 1 1"), 2, 3), "row index that is not in 1..2")
  expect_error(rd(c("1 1 1", "0 1 1"), 2, 3), "row index that is not in 1..2")
  expect_error(rd(c("1 1 1", "1 4 1"), 2, 3), "column index that is not in 1..3")
  expect_error(rd(c("1 1 1", "1 -1 1"), 2, 3), "column index that is not in 1..3")
  expect_error(rd(c("1 1 1", "2 2 1"), 2, 3, nnz = 3),
               "has 2 entries, but its size line gives 3 entries")
  expect_error(rd(c("1 1 1", "2 2 1"), 2, 3, nnz = 1),
               "has more entries than its size line gives \\(1\\)")
  expect_error(rd("1 1 1", 2, 3, nnz = "2147483648"),
               "has 2147483648 entries. A dgCMatrix can have at most 2147483647 entries.")
  for (b in c("%%MatrixMarket matrix array real general",
              "%%MatrixMarket matrix coordinate pattern general",
              "%%MatrixMarket matrix coordinate complex general",
              "%%MatrixMarket matrix coordinate real symmetric",
              "%%matrixmarket matrix coordinate real general",
              "%%MatrixMarket matrix coordinate real general extra")) {
    expect_error(rd("1 1 1", 2, 3, banner = b),
                 "is not a MatrixMarket 'coordinate real general' file. The first line is: %%",
                 info = b)
  }
  expect_error(rd(c("1 1", "2 2 1"), 2, 3), "line 4 does not have 3 numbers")
  expect_error(rd(c("1 1 1 1", "2 2 1"), 2, 3), "line 4 does not have 3 numbers")
  expect_error(rd(c("1 1 x", "2 2 1"), 2, 3), "line 4 does not have 3 numbers")
  expect_error(rd(c("1.0 1 1", "2 2 1"), 2, 3), "line 4 does not have 3 numbers")
  expect_error(rd(paste("1 1", strrep("1", 1100)), 2, 3), "line 4 is too long")
  expect_error(rd(character(0), 2, 3, nnz = ""), "does not have a correct MatrixMarket size line: 2 3")
  expect_error(rd(character(0), "2147483648", 3, 0),
               "does not have a correct MatrixMarket size line")
  empty <- tempfile(fileext = ".mtx")
  file.create(empty)
  expect_error(fishpond:::.readFryMtx(empty), "is empty. It is not a MatrixMarket file.")
  noSize <- tempfile(fileext = ".mtx")
  writeLines(c("%%MatrixMarket matrix coordinate real general", "% comment"), noSize)
  expect_error(fishpond:::.readFryMtx(noSize), "does not have a MatrixMarket size line")
  nul <- tempfile(fileext = ".mtx")
  writeBin(c(charToRaw("%%MatrixMarket matrix coordinate real general\n2 2 2\n1 2 3"),
             as.raw(0), charToRaw(" 9\n2 1 4\n")), nul)
  expect_error(fishpond:::.readFryMtx(nul), "line 3 has a NUL character")
  expect_error(rd("1 1 NA", 2, 3), "line 4 does not have 3 numbers")
  expect_error(suppressWarnings(fishpond:::.readFryMtx(file.path(tempdir(), "no-such-file.mtx"))),
               "cannot open")
})

test_that("loadFry gives the number of barcodes in its message", {
  dat <- fishpond:::readExampleFryData("fry-usa-basic")
  expect_message(loadFry(dat$parent_dir),
                 paste("and", length(dat$barcodes), "barcodes"))
})

test_that("loadFry reads a directory with a name that is not ASCII", {
  d <- file.path(tempfile("fry"), "été_中")
  ok <- tryCatch(dir.create(d, recursive = TRUE), warning = function(w) FALSE,
                 error = function(e) FALSE)
  if (!isTRUE(ok) || !dir.exists(d)) {
    skip("this platform cannot make a directory with a name that is not ASCII")
  }
  src <- fryExample("flex", "flex_probe")
  file.copy(list.files(src, full.names = TRUE), d, recursive = TRUE)
  if (!file.exists(file.path(d, "alevin", "probe_quants_mat.mtx"))) {
    skip("this platform cannot copy files into a directory with a name that is not ASCII")
  }
  for (level in c("gene", "probe")) {
    expect_identical(loadFry(d, featureLevel = level, quiet = TRUE),
                     loadFry(src, featureLevel = level, quiet = TRUE), info = level)
  }
})
