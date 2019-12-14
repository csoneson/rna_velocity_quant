suppressPackageStartupMessages({
  library(tximeta)
  library(dplyr)
  library(ggplot2)
  library(zeallot)
  library(BUSpaRse)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(velocyto.R)
  library(scater)
  library(readr)
})

read_velocyto <- function(loomfile, sampleid) {
  loom <- velocyto.R::read.loom.matrices(loomfile)
  velocyto <- SingleCellExperiment(
    assays = list(counts = loom$spliced,
                  spliced = loom$spliced,
                  unspliced = loom$unspliced)
  )
  colnames(velocyto) <- gsub(paste0(sampleid, ":"), "", gsub("x$", "", colnames(velocyto)))
  velocyto
}

read_starsolo <- function(solodir, sampleid) {
  genes <- readr::read_delim(file.path(solodir, "features.tsv"), col_names = FALSE, delim = "\t")
  cells <- readr::read_delim(file.path(solodir, "barcodes.tsv"), col_names = FALSE, delim = "\t")
  n_genes <- nrow(genes)
  n_cells <- nrow(cells)
  ## Read matrix. col 1 = gene, col 2 = cell, col 3 = spliced, col 4 = unspliced, col 5 = ambiguous
  matr <- readr::read_delim(file.path(solodir, "matrix.mtx"), skip = 3, delim = " ", col_names = FALSE)
  scounts <- Matrix::sparseMatrix(
    i = matr$X1,
    j = matr$X2,
    x = matr$X3,
    dims = c(n_genes, n_cells)
  )
  ucounts <- Matrix::sparseMatrix(
    i = matr$X1,
    j = matr$X2,
    x = matr$X4,
    dims = c(n_genes, n_cells)
  )
  rownames(scounts) <- rownames(ucounts) <- genes$X2
  colnames(scounts) <- colnames(ucounts) <- cells$X1
  SingleCellExperiment(
    assays = list(counts = scounts,
                  spliced = scounts,
                  unspliced = ucounts)
  )
}

read_starsolo_subtract <- function(solodir, sampleid) {
  ## Read gene body and exon counts, infer intronic counts as the difference
  ## Gene body
  fgenes <- readr::read_delim(file.path(solodir, "GeneFull/raw/features.tsv"), col_names = FALSE, delim = "\t")
  fcells <- readr::read_delim(file.path(solodir, "GeneFull/raw/barcodes.tsv"), col_names = FALSE, delim = "\t")
  fn_genes <- nrow(fgenes)
  fn_cells <- nrow(fcells)
  ## Read matrix. col 1 = gene, col 2 = cell, col 3 = spliced, col 4 = unspliced, col 5 = ambiguous
  fmatr <- readr::read_delim(file.path(solodir, "GeneFull/raw/matrix.mtx"), skip = 3, delim = " ", col_names = FALSE)
  fcounts <- Matrix::sparseMatrix(
    i = fmatr$X1,
    j = fmatr$X2,
    x = fmatr$X3,
    dims = c(fn_genes, fn_cells)
  )
  rownames(fcounts) <- fgenes$X2
  colnames(fcounts) <- fcells$X1
  
  ## Exons
  sgenes <- readr::read_delim(file.path(solodir, "Gene/raw/features.tsv"), col_names = FALSE, delim = "\t")
  scells <- readr::read_delim(file.path(solodir, "Gene/raw/barcodes.tsv"), col_names = FALSE, delim = "\t")
  sn_genes <- nrow(sgenes)
  sn_cells <- nrow(scells)
  ## Read matrix. col 1 = gene, col 2 = cell, col 3 = spliced, col 4 = unspliced, col 5 = ambiguous
  smatr <- readr::read_delim(file.path(solodir, "Gene/raw/matrix.mtx"), skip = 3, delim = " ", col_names = FALSE)
  scounts <- Matrix::sparseMatrix(
    i = smatr$X1,
    j = smatr$X2,
    x = smatr$X3,
    dims = c(sn_genes, sn_cells)
  )
  rownames(scounts) <- sgenes$X2
  colnames(scounts) <- scells$X1
  
  stopifnot(all(fgenes$X1 == sgenes$X1))
  stopifnot(all(fgenes$X2 == sgenes$X2))
  stopifnot(all(fcells$X1 == scells$X1))
  
  umatr <- fmatr %>% dplyr::rename(f = X3) %>%
    dplyr::full_join(smatr %>% dplyr::rename(s = X3)) %>%
    dplyr::mutate(f = replace(f, is.na(f), 0),
                  s = replace(s, is.na(s), 0)) %>%
    dplyr::mutate(u = f - s) %>%
    dplyr::mutate(u = replace(u, u < 0, 0))
  ucounts <- Matrix::sparseMatrix(
    i = umatr$X1,
    j = umatr$X2,
    x = umatr$u,
    dims = c(sn_genes, sn_cells)
  )
  rownames(ucounts) <- sgenes$X2
  colnames(ucounts) <- scells$X1
  
  stopifnot(all(rownames(scounts) == rownames(fcounts)))
  stopifnot(all(rownames(scounts) == rownames(ucounts)))
  stopifnot(all(colnames(scounts) == colnames(fcounts)))
  stopifnot(all(colnames(scounts) == colnames(ucounts)))
  
  SingleCellExperiment(
    assays = list(counts = scounts,
                  spliced = scounts,
                  unspliced = ucounts)
  )
}

sce_from_scounts_ucounts <- function(scounts, ucounts) {
  ucounts <- ucounts[rownames(ucounts) %in% rownames(scounts), colnames(ucounts) %in% colnames(scounts)]
  ucounts2 <- scounts
  ucounts2[] <- 0
  tm <- as(ucounts, "dgTMatrix")
  ucounts2[cbind(match(rownames(ucounts)[tm@i + 1], rownames(ucounts2)), 
                 match(colnames(ucounts)[tm@j + 1], colnames(ucounts2)))] <- tm@x
  stopifnot(all(rownames(ucounts2) == rownames(scounts)))
  stopifnot(all(colnames(ucounts2) == colnames(scounts)))
  stopifnot(sum(ucounts2) == sum(ucounts))
  SingleCellExperiment(
    assays = list(counts = as(scounts, "dgCMatrix"),
                  spliced = as(scounts, "dgCMatrix"),
                  unspliced = as(ucounts2, "dgCMatrix"))
  )
}

read_alevin_with_decoys <- function(spliceddir, unspliceddir, sampleid, tx2gene) {
  cdna <- tximeta(coldata = data.frame(
    names = sampleid,
    files = file.path(spliceddir, "quants_mat.gz"),
    stringsAsFactors = FALSE
  ), type = "alevin")
  cdna <- cdna[grep("I\\.*$", rownames(cdna), invert = TRUE), ]
  introns <- tximeta(coldata = data.frame(
    names = sampleid,
    files = file.path(unspliceddir, "quants_mat.gz"),
    stringsAsFactors = FALSE
  ), type = "alevin")
  introns <- introns[grep("I\\.*$", rownames(introns)), ]
  ucounts <- round(assay(introns, "counts"))
  scounts <- round(assay(cdna, "counts"))
  rownames(ucounts) <- gsub("\\.*I\\.*$", "", rownames(ucounts))
  rownames(scounts) <- gsub("\\.*$", "", rownames(scounts))
  cdna_introns_decoy <- sce_from_scounts_ucounts(scounts, ucounts)
  rownames(cdna_introns_decoy) <- scater::uniquifyFeatureNames(
    ID = rownames(cdna_introns_decoy),
    names = tx2gene$gene_name[match(rownames(cdna_introns_decoy), tx2gene$gene_id)]
  )
  cdna_introns_decoy
}

read_alevin_cdna_introns <- function(alevindir, sampleid, tx2gene) {
  cdna_introns <- tximeta(coldata = data.frame(
    names = sampleid,
    files = file.path(alevindir, "quants_mat.gz"),
    stringsAsFactors = FALSE
  ), type = "alevin")
  uidx <- grep("\\.*I\\.*$", rownames(cdna_introns))
  sidx <- grep("\\.*I\\.*$", rownames(cdna_introns), invert = TRUE)
  ucounts <- round(assay(cdna_introns, "counts")[uidx, ])
  scounts <- round(assay(cdna_introns, "counts")[sidx, ])
  rownames(ucounts) <- gsub("\\.*I\\.*$", "", rownames(ucounts))
  rownames(scounts) <- gsub("\\.*$", "", rownames(scounts))
  cdna_introns <- sce_from_scounts_ucounts(scounts, ucounts)
  rownames(cdna_introns) <- scater::uniquifyFeatureNames(
    ID = rownames(cdna_introns),
    names = tx2gene$gene_name[match(rownames(cdna_introns), tx2gene$gene_id)]
  )
  cdna_introns
}

read_alevin_spliced_unspliced <- function(alevindir, sampleid, tx2gene) {
  ## spliced + unspliced
  spliced_unspliced <- tximeta(coldata = data.frame(
    names = sampleid,
    files = file.path(alevindir, "quants_mat.gz"),
    stringsAsFactors = FALSE
  ), type = "alevin")
  uidx <- grep("_unspliced$", rownames(spliced_unspliced))
  sidx <- grep("_unspliced$", rownames(spliced_unspliced), invert = TRUE)
  ucounts <- round(assay(spliced_unspliced, "counts")[uidx, ])
  scounts <- round(assay(spliced_unspliced, "counts")[sidx, ])
  rownames(ucounts) <- gsub("_unspliced$", "", rownames(ucounts))
  spliced_unspliced <- sce_from_scounts_ucounts(scounts, ucounts)
  rownames(spliced_unspliced) <- scater::uniquifyFeatureNames(
    ID = rownames(spliced_unspliced),
    names = tx2gene$gene_name[match(rownames(spliced_unspliced), tx2gene$gene_id)]
  )
  spliced_unspliced
}

read_alevin_spliced <- function(alevindir, sampleid, tx2gene) {
  ## spliced
  spliced <- tximeta(coldata = data.frame(
    names = sampleid,
    files = file.path(alevindir, "quants_mat.gz"),
    stringsAsFactors = FALSE
  ), type = "alevin")
  rownames(spliced) <- scater::uniquifyFeatureNames(
    ID = rownames(spliced),
    names = tx2gene$gene_name[match(rownames(spliced), tx2gene$gene_id)]
  )
  spliced <- as(spliced, "SingleCellExperiment")
  spliced
}

read_kallisto_bustools <- function(kallistodir, splicedname, unsplicedname) {
  c(spliced_bus, unspliced_bus) %<-% 
    read_velocity_output(spliced_dir = kallistodir,
                         spliced_name = splicedname,
                         unspliced_dir = kallistodir,
                         unspliced_name = unsplicedname)
  rownames(spliced_bus) <- gsub("\\.*$", "", rownames(spliced_bus))
  rownames(unspliced_bus) <- gsub("\\.*$", "", rownames(unspliced_bus))
  stopifnot(all(rownames(spliced_bus) == rownames(unspliced_bus)))
  kallistobus_cells <- intersect(colnames(spliced_bus), colnames(unspliced_bus))
  kallistobus <- SingleCellExperiment(
    assays = list(counts = spliced_bus[, kallistobus_cells],
                  spliced = spliced_bus[, kallistobus_cells],
                  unspliced = unspliced_bus[, kallistobus_cells])
  )
  rownames(kallistobus) <- scater::uniquifyFeatureNames(
    ID = rownames(kallistobus),
    names = tx2gene$gene_name[match(rownames(kallistobus), tx2gene$gene_id)]
  )
  kallistobus
}

