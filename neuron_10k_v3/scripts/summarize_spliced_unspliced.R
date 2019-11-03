args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

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

print(topdir)
print(tx2gene)
print(outrds)

tx2gene <- readRDS(tx2gene)

## ========================================================================= ##
## Read quantifications and create SingleCellExperiment objects
## ========================================================================= ##
## CellRanger + velocyto
loom <- velocyto.R::read.loom.matrices(file.path(topdir, "quants/cellranger/neuron_10k_v3/velocyto/neuron_10k_v3.loom"))
velocyto <- SingleCellExperiment(
  assays = list(counts = loom$spliced,
                spliced = loom$spliced,
                unspliced = loom$unspliced)
)
colnames(velocyto) <- gsub("neuron_10k_v3:", "", gsub("x$", "", colnames(velocyto)))

## cDNA/introns separately (with decoys)
cdna <- tximeta(coldata = data.frame(
  names = "neuron_10k_v3",
  files = file.path(topdir, "quants/salmon_gencodevM21_cDNA_intronsasdecoy/alevin/quants_mat.gz"),
  stringsAsFactors = FALSE
), type = "alevin")
cdna <- cdna[grep("I\\.$", rownames(cdna), invert = TRUE), ]
introns <- tximeta(coldata = data.frame(
  names = "neuron_10k_v3",
  files = file.path(topdir, "quants/salmon_gencodevM21_introns_cDNAasdecoy/alevin/quants_mat.gz"),
  stringsAsFactors = FALSE
), type = "alevin")
introns <- introns[grep("I\\.$", rownames(introns)), ]
ucounts <- round(assay(introns, "counts"))
scounts <- round(assay(cdna, "counts"))
rownames(ucounts) <- gsub("\\.I\\.$", "", rownames(ucounts))
rownames(scounts) <- gsub("\\.$", "", rownames(scounts))
ucounts <- ucounts[rownames(ucounts) %in% rownames(scounts), colnames(ucounts) %in% colnames(scounts)]
ucounts2 <- scounts
ucounts2[] <- 0
tm <- as(ucounts, "dgTMatrix")
ucounts2[cbind(match(rownames(ucounts)[tm@i + 1], rownames(ucounts2)), 
               match(colnames(ucounts)[tm@j + 1], colnames(ucounts2)))] <- tm@x
stopifnot(all(rownames(ucounts2) == rownames(scounts)))
stopifnot(all(colnames(ucounts2) == colnames(scounts)))
stopifnot(sum(ucounts2) == sum(ucounts))
cdna_introns_decoy <- SingleCellExperiment(
  assays = list(counts = as(scounts, "dgCMatrix"),
                spliced = as(scounts, "dgCMatrix"),
                unspliced = as(ucounts2, "dgCMatrix"))
)
rownames(cdna_introns_decoy) <- scater::uniquifyFeatureNames(
  ID = rownames(cdna_introns_decoy),
  names = tx2gene$gene_name[match(rownames(cdna_introns_decoy), tx2gene$gene_id)]
)

## cDNA + introns
cdna_introns <- tximeta(coldata = data.frame(
  names = "neuron_10k_v3",
  files = file.path(topdir, "quants/salmon_gencodevM21_cDNA_introns/alevin/quants_mat.gz"),
  stringsAsFactors = FALSE
), type = "alevin")
uidx <- grep("\\.I\\.$", rownames(cdna_introns))
sidx <- grep("\\.I\\.$", rownames(cdna_introns), invert = TRUE)
ucounts <- round(assay(cdna_introns, "counts")[uidx, ])
scounts <- round(assay(cdna_introns, "counts")[sidx, ])
rownames(ucounts) <- gsub("\\.I\\.$", "", rownames(ucounts))
rownames(scounts) <- gsub("\\.$", "", rownames(scounts))
ucounts <- ucounts[rownames(ucounts) %in% rownames(scounts), ]
ucounts2 <- scounts
ucounts2[] <- 0
tm <- as(ucounts, "dgTMatrix")
ucounts2[cbind(match(rownames(ucounts)[tm@i + 1], rownames(ucounts2)), 
               match(colnames(ucounts)[tm@j + 1], colnames(ucounts2)))] <- tm@x
stopifnot(all(rownames(ucounts2) == rownames(scounts)))
stopifnot(all(colnames(ucounts2) == colnames(scounts)))
stopifnot(sum(ucounts2) == sum(ucounts))
cdna_introns <- SingleCellExperiment(
  assays = list(counts = as(scounts, "dgCMatrix"),
                spliced = as(scounts, "dgCMatrix"),
                unspliced = as(ucounts2, "dgCMatrix"))
)
rownames(cdna_introns) <- scater::uniquifyFeatureNames(
  ID = rownames(cdna_introns),
  names = tx2gene$gene_name[match(rownames(cdna_introns), tx2gene$gene_id)]
)

## cDNA + collapsed introns
cdna_intronscollapsed <- tximeta(coldata = data.frame(
  names = "neuron_10k_v3",
  files = file.path(topdir, "quants/salmon_gencodevM21_cDNA_intronscollapsed/alevin/quants_mat.gz"),
  stringsAsFactors = FALSE
), type = "alevin")
uidx <- grep("\\.I\\.$", rownames(cdna_intronscollapsed))
sidx <- grep("\\.I\\.$", rownames(cdna_intronscollapsed), invert = TRUE)
ucounts <- round(assay(cdna_intronscollapsed, "counts")[uidx, ])
scounts <- round(assay(cdna_intronscollapsed, "counts")[sidx, ])
rownames(ucounts) <- gsub("\\.I\\.$", "", rownames(ucounts))
rownames(scounts) <- gsub("\\.$", "", rownames(scounts))
ucounts <- ucounts[rownames(ucounts) %in% rownames(scounts), ]
ucounts2 <- scounts
ucounts2[] <- 0
tm <- as(ucounts, "dgTMatrix")
ucounts2[cbind(match(rownames(ucounts)[tm@i + 1], rownames(ucounts2)), 
               match(colnames(ucounts)[tm@j + 1], colnames(ucounts2)))] <- tm@x
stopifnot(all(rownames(ucounts2) == rownames(scounts)))
stopifnot(all(colnames(ucounts2) == colnames(scounts)))
stopifnot(sum(ucounts2) == sum(ucounts))
cdna_intronscollapsed <- SingleCellExperiment(
  assays = list(counts = as(scounts, "dgCMatrix"),
                spliced = as(scounts, "dgCMatrix"),
                unspliced = as(ucounts2, "dgCMatrix"))
)
rownames(cdna_intronscollapsed) <- scater::uniquifyFeatureNames(
  ID = rownames(cdna_intronscollapsed),
  names = tx2gene$gene_name[match(rownames(cdna_intronscollapsed), tx2gene$gene_id)]
)

## spliced + unspliced
spliced_unspliced <- tximeta(coldata = data.frame(
  names = "neuron_10k_v3",
  files = file.path(topdir, "quants/salmon_gencodevM21_spliced_unspliced/alevin/quants_mat.gz"),
  stringsAsFactors = FALSE
), type = "alevin")
uidx <- grep("_unspliced$", rownames(spliced_unspliced))
sidx <- grep("_unspliced$", rownames(spliced_unspliced), invert = TRUE)
ucounts <- round(assay(spliced_unspliced, "counts")[uidx, ])
scounts <- round(assay(spliced_unspliced, "counts")[sidx, ])
rownames(ucounts) <- gsub("_unspliced$", "", rownames(ucounts))
ucounts <- ucounts[rownames(ucounts) %in% rownames(scounts), ]
ucounts2 <- scounts
ucounts2[] <- 0
tm <- as(ucounts, "dgTMatrix")
ucounts2[cbind(match(rownames(ucounts)[tm@i + 1], rownames(ucounts2)), 
               match(colnames(ucounts)[tm@j + 1], colnames(ucounts2)))] <- tm@x
stopifnot(all(rownames(ucounts2) == rownames(scounts)))
stopifnot(all(colnames(ucounts2) == colnames(scounts)))
stopifnot(sum(ucounts2) == sum(ucounts))
spliced_unspliced <- SingleCellExperiment(
  assays = list(counts = as(scounts, "dgCMatrix"),
                spliced = as(scounts, "dgCMatrix"),
                unspliced = as(ucounts2, "dgCMatrix"))
)
rownames(spliced_unspliced) <- scater::uniquifyFeatureNames(
  ID = rownames(spliced_unspliced),
  names = tx2gene$gene_name[match(rownames(spliced_unspliced), tx2gene$gene_id)]
)


## spliced
spliced <- tximeta(coldata = data.frame(
  names = "neuron_10k_v3",
  files = file.path(topdir, "quants/salmon_gencodevM21_spliced/alevin/quants_mat.gz"),
  stringsAsFactors = FALSE
), type = "alevin")
rownames(spliced) <- scater::uniquifyFeatureNames(
  ID = rownames(spliced),
  names = tx2gene$gene_name[match(rownames(spliced), tx2gene$gene_id)]
)
spliced <- as(spliced, "SingleCellExperiment")

## kallisto/bustools
## separate
kallistodir <- file.path(topdir, "quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns")
c(spliced_bus, unspliced_bus) %<-% 
  read_velocity_output(spliced_dir = kallistodir,
                       spliced_name = "spliced",
                       unspliced_dir = kallistodir,
                       unspliced_name = "unspliced")
rownames(spliced_bus) <- gsub("\\.$", "", rownames(spliced_bus))
rownames(unspliced_bus) <- gsub("\\.$", "", rownames(unspliced_bus))
stopifnot(all(rownames(spliced_bus) == rownames(unspliced_bus)))
kallistobus_cells <- intersect(colnames(spliced_bus), colnames(unspliced_bus))
kallistobus_separate <- SingleCellExperiment(
  assays = list(counts = spliced_bus[, kallistobus_cells],
                spliced = spliced_bus[, kallistobus_cells],
                unspliced = unspliced_bus[, kallistobus_cells])
)
rownames(kallistobus_separate) <- scater::uniquifyFeatureNames(
  ID = rownames(kallistobus_separate),
  names = tx2gene$gene_name[match(rownames(kallistobus_separate), tx2gene$gene_id)]
)

## separate, permissive
kallistodir <- file.path(topdir, "quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns")
c(spliced_bus, unspliced_bus) %<-% 
  read_velocity_output(spliced_dir = kallistodir,
                       spliced_name = "spliced.permissive",
                       unspliced_dir = kallistodir,
                       unspliced_name = "unspliced.permissive")
rownames(spliced_bus) <- gsub("\\.$", "", rownames(spliced_bus))
rownames(unspliced_bus) <- gsub("\\.$", "", rownames(unspliced_bus))
stopifnot(all(rownames(spliced_bus) == rownames(unspliced_bus)))
kallistobus_cells <- intersect(colnames(spliced_bus), colnames(unspliced_bus))
kallistobus_sep_permissive <- SingleCellExperiment(
  assays = list(counts = spliced_bus[, kallistobus_cells], 
                spliced = spliced_bus[, kallistobus_cells],
                unspliced = unspliced_bus[, kallistobus_cells])
)
rownames(kallistobus_sep_permissive) <- scater::uniquifyFeatureNames(
  ID = rownames(kallistobus_sep_permissive),
  names = tx2gene$gene_name[match(rownames(kallistobus_sep_permissive), tx2gene$gene_id)]
)

## collapse
kallistodir <- file.path(topdir, "quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns")
c(spliced_bus, unspliced_bus) %<-% 
  read_velocity_output(spliced_dir = kallistodir,
                       spliced_name = "spliced",
                       unspliced_dir = kallistodir,
                       unspliced_name = "unspliced")
rownames(spliced_bus) <- gsub("\\.$", "", rownames(spliced_bus))
rownames(unspliced_bus) <- gsub("\\.$", "", rownames(unspliced_bus))
stopifnot(all(rownames(spliced_bus) == rownames(unspliced_bus)))
kallistobus_cells <- intersect(colnames(spliced_bus), colnames(unspliced_bus))
kallistobus_collapse <- SingleCellExperiment(
  assays = list(counts = spliced_bus[, kallistobus_cells], 
                spliced = spliced_bus[, kallistobus_cells],
                unspliced = unspliced_bus[, kallistobus_cells])
)
rownames(kallistobus_collapse) <- scater::uniquifyFeatureNames(
  ID = rownames(kallistobus_collapse),
  names = tx2gene$gene_name[match(rownames(kallistobus_collapse), tx2gene$gene_id)]
)

## collapse, permissive
kallistodir <- file.path(topdir, "quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns")
c(spliced_bus, unspliced_bus) %<-% 
  read_velocity_output(spliced_dir = kallistodir,
                       spliced_name = "spliced.permissive",
                       unspliced_dir = kallistodir,
                       unspliced_name = "unspliced.permissive")
rownames(spliced_bus) <- gsub("\\.$", "", rownames(spliced_bus))
rownames(unspliced_bus) <- gsub("\\.$", "", rownames(unspliced_bus))
stopifnot(all(rownames(spliced_bus) == rownames(unspliced_bus)))
kallistobus_cells <- intersect(colnames(spliced_bus), colnames(unspliced_bus))
kallistobus_coll_permissive <- SingleCellExperiment(
  assays = list(counts = spliced_bus[, kallistobus_cells], 
                spliced = spliced_bus[, kallistobus_cells],
                unspliced = unspliced_bus[, kallistobus_cells])
)
rownames(kallistobus_coll_permissive) <- scater::uniquifyFeatureNames(
  ID = rownames(kallistobus_coll_permissive),
  names = tx2gene$gene_name[match(rownames(kallistobus_coll_permissive), tx2gene$gene_id)]
)

## ========================================================================= ##
## subset to shared cells/genes
## ========================================================================= ##
shared_cells <- Reduce(intersect, list(colnames(velocyto),
                                       colnames(cdna_introns_decoy),
                                       colnames(cdna_introns),
                                       colnames(cdna_intronscollapsed),
                                       colnames(spliced_unspliced),
                                       colnames(spliced),
                                       colnames(kallistobus_separate),
                                       colnames(kallistobus_collapse),
                                       colnames(kallistobus_sep_permissive),
                                       colnames(kallistobus_coll_permissive)))
shared_genes <- Reduce(intersect, list(rownames(velocyto),
                                       rownames(cdna_introns_decoy),
                                       rownames(cdna_introns),
                                       rownames(cdna_intronscollapsed),
                                       rownames(spliced_unspliced),
                                       rownames(spliced),
                                       rownames(kallistobus_separate),
                                       rownames(kallistobus_collapse),
                                       rownames(kallistobus_sep_permissive),
                                       rownames(kallistobus_coll_permissive)))

velocyto <- velocyto[shared_genes, shared_cells]
cdna_introns_decoy <- cdna_introns_decoy[shared_genes, shared_cells]
cdna_introns <- cdna_introns[shared_genes, shared_cells]
cdna_intronscollapsed <- cdna_intronscollapsed[shared_genes, shared_cells]
spliced_unspliced <- spliced_unspliced[shared_genes, shared_cells]
spliced <- spliced[shared_genes, shared_cells]
kallistobus_separate <- kallistobus_separate[shared_genes, shared_cells]
kallistobus_collapse <- kallistobus_collapse[shared_genes, shared_cells]
kallistobus_sep_permissive <- kallistobus_sep_permissive[shared_genes, shared_cells]
kallistobus_coll_permissive <- kallistobus_coll_permissive[shared_genes, shared_cells]

## ========================================================================= ##
## Save
## ========================================================================= ##
saveRDS(velocyto, file = file.path(dirname(outrds), "sce_velocyto.rds"))
saveRDS(cdna_introns_decoy, file = file.path(dirname(outrds), "sce_cdna_introns_decoy.rds"))
saveRDS(cdna_introns, file = file.path(dirname(outrds), "sce_cdna_introns.rds"))
saveRDS(cdna_intronscollapsed, file = file.path(dirname(outrds), "sce_cdna_intronscollapsed.rds"))
saveRDS(spliced_unspliced, file = file.path(dirname(outrds), "sce_spliced_unspliced.rds"))
saveRDS(spliced, file = file.path(dirname(outrds), "sce_spliced.rds"))
saveRDS(kallistobus_separate, file = file.path(dirname(outrds), "sce_kallistobus_separate.rds"))
saveRDS(kallistobus_collapse, file = file.path(dirname(outrds), "sce_kallistobus_collapse.rds"))
saveRDS(kallistobus_sep_permissive, file = file.path(dirname(outrds), "sce_kallistobus_sep_permissive.rds"))
saveRDS(kallistobus_coll_permissive, file = file.path(dirname(outrds), "sce_kallistobus_coll_permissive.rds"))
saveRDS(NULL, file = outrds)

date()
sessionInfo()
