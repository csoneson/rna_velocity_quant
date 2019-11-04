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
  rownames(ucounts) <- gsub("\\.I\\.*$", "", rownames(ucounts))
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
  uidx <- grep("\\.I\\.*$", rownames(cdna_introns))
  sidx <- grep("\\.I\\.*$", rownames(cdna_introns), invert = TRUE)
  ucounts <- round(assay(cdna_introns, "counts")[uidx, ])
  scounts <- round(assay(cdna_introns, "counts")[sidx, ])
  rownames(ucounts) <- gsub("\\.I\\.*$", "", rownames(ucounts))
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
                         unspliced_name = splicedname)
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

