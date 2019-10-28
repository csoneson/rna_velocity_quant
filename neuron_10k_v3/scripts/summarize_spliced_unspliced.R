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
  library(velocyto.R)
  library(scater)
})

print(topdir)
print(tx2gene)
print(outrds)

tx2gene <- readRDS(tx2gene)

## ========================================================================= ##
## Read quantifications and create SummarizedExperiment objects
## ========================================================================= ##
## CellRanger + velocyto
loom <- velocyto.R::read.loom.matrices(file.path(topdir, "quants/cellranger/neuron_10k_v3/velocyto/neuron_10k_v3.loom"))
velocyto <- SummarizedExperiment(
  assays = list(spliced = loom$spliced,
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
cdna_introns_decoy <- SummarizedExperiment(
  assays = list(spliced = scounts,
                unspliced = ucounts2)
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
cdna_introns <- SummarizedExperiment(
  assays = list(spliced = scounts,
                unspliced = ucounts2)
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
cdna_intronscollapsed <- SummarizedExperiment(
  assays = list(spliced = scounts,
                unspliced = ucounts2)
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
spliced_unspliced <- SummarizedExperiment(
  assays = list(spliced = scounts,
                unspliced = ucounts2)
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
kallistobus_separate <- SummarizedExperiment(
  assays = list(spliced = spliced_bus[, kallistobus_cells],
                unspliced = unspliced_bus[, kallistobus_cells])
)
rownames(kallistobus_separate) <- scater::uniquifyFeatureNames(
  ID = rownames(kallistobus_separate),
  names = tx2gene$gene_name[match(rownames(kallistobus_separate), tx2gene$gene_id)]
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
kallistobus_collapse <- SummarizedExperiment(
  assays = list(spliced = spliced_bus[, kallistobus_cells],
                unspliced = unspliced_bus[, kallistobus_cells])
)
rownames(kallistobus_collapse) <- scater::uniquifyFeatureNames(
  ID = rownames(kallistobus_collapse),
  names = tx2gene$gene_name[match(rownames(kallistobus_collapse), tx2gene$gene_id)]
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
                                       colnames(kallistobus_collapse)))
shared_genes <- Reduce(intersect, list(rownames(velocyto),
                                       rownames(cdna_introns_decoy),
                                       rownames(cdna_introns),
                                       rownames(cdna_intronscollapsed),
                                       rownames(spliced_unspliced),
                                       rownames(spliced),
                                       rownames(kallistobus_separate),
                                       rownames(kallistobus_collapse)))

velocyto <- velocyto[shared_genes, shared_cells]
cdna_introns_decoy <- cdna_introns_decoy[shared_genes, shared_cells]
cdna_introns <- cdna_introns[shared_genes, shared_cells]
cdna_intronscollapsed <- cdna_intronscollapsed[shared_genes, shared_cells]
spliced_unspliced <- spliced_unspliced[shared_genes, shared_cells]
spliced <- spliced[shared_genes, shared_cells]
kallistobus_separate <- kallistobus_separate[shared_genes, shared_cells]
kallistobus_collapse <- kallistobus_collapse[shared_genes, shared_cells]

## ========================================================================= ##
## Plot
## ========================================================================= ##
df0 <- dplyr::bind_rows(
  data.frame(method = "velocyto",
             spliced_count = rowSums(assay(velocyto, "spliced")),
             unspliced_count = rowSums(assay(velocyto, "unspliced")),
             total_count = rowSums(assay(velocyto, "spliced")) + 
               rowSums(assay(velocyto, "unspliced")),
             stringsAsFactors = FALSE),
  data.frame(method = "cdna_introns_decoy",
             spliced_count = rowSums(assay(cdna_introns_decoy, "spliced")),
             unspliced_count = rowSums(assay(cdna_introns_decoy, "unspliced")),
             total_count = rowSums(assay(cdna_introns_decoy, "spliced")) + 
               rowSums(assay(cdna_introns_decoy, "unspliced")),
             stringsAsFactors = FALSE),
  data.frame(method = "cdna_introns",
             spliced_count = rowSums(assay(cdna_introns, "spliced")),
             unspliced_count = rowSums(assay(cdna_introns, "unspliced")),
             total_count = rowSums(assay(cdna_introns, "spliced")) + 
               rowSums(assay(cdna_introns, "unspliced")),
             stringsAsFactors = FALSE),
  data.frame(method = "cdna_intronscollapsed",
             spliced_count = rowSums(assay(cdna_intronscollapsed, "spliced")),
             unspliced_count = rowSums(assay(cdna_intronscollapsed, "unspliced")),
             total_count = rowSums(assay(cdna_intronscollapsed, "spliced")) + 
               rowSums(assay(cdna_intronscollapsed, "unspliced")),
             stringsAsFactors = FALSE),
  data.frame(method = "spliced_unspliced",
             spliced_count = rowSums(assay(spliced_unspliced, "spliced")),
             unspliced_count = rowSums(assay(spliced_unspliced, "unspliced")),
             total_count = rowSums(assay(spliced_unspliced, "spliced")) + 
               rowSums(assay(spliced_unspliced, "unspliced")),
             stringsAsFactors = FALSE),
  data.frame(method = "spliced",
             spliced_count = rowSums(assay(spliced, "counts")),
             total_count = rowSums(assay(spliced, "counts")),
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_separate",
             spliced_count = rowSums(assay(kallistobus_separate, "spliced")),
             unspliced_count = rowSums(assay(kallistobus_separate, "unspliced")),
             total_count = rowSums(assay(kallistobus_separate, "spliced")) + 
               rowSums(assay(kallistobus_separate, "unspliced")),
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_collapse",
             spliced_count = rowSums(assay(kallistobus_collapse, "spliced")),
             unspliced_count = rowSums(assay(kallistobus_collapse, "unspliced")),
             total_count = rowSums(assay(kallistobus_collapse, "spliced")) + 
               rowSums(assay(kallistobus_collapse, "unspliced")),
             stringsAsFactors = FALSE)
) %>% 
  tidyr::gather(key = "type", value = "numi", -method)

pdf(gsub("rds", "pdf", outrds), width = 9, height = 9)
ggplot(df0 %>% dplyr::mutate(
  type = factor(type, levels = c("spliced_count", "unspliced_count", "total_count"))), 
       aes(x = method, y = numi + 1)) + 
  geom_violin(aes(fill = method), alpha = 0.25, scale = "width") + 
  geom_boxplot(width = 0.025) + 
  facet_wrap(~ type, ncol = 1) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_log10() + 
  ggtitle("Total count per gene, across all cells")

df1 <- dplyr::bind_rows(
  data.frame(method = "velocyto",
             fraction_unspliced = rowSums(assay(velocyto, "unspliced"))/
               (rowSums(assay(velocyto, "unspliced")) + 
                  rowSums(assay(velocyto, "spliced"))),
             stringsAsFactors = FALSE),
  data.frame(method = "cdna_introns_decoy", 
             fraction_unspliced = rowSums(assay(cdna_introns_decoy, "unspliced"))/
               (rowSums(assay(cdna_introns_decoy, "unspliced")) + 
                  rowSums(assay(cdna_introns_decoy, "spliced"))),
             stringsAsFactors = FALSE),
  data.frame(method = "cdna_introns", 
             fraction_unspliced = rowSums(assay(cdna_introns, "unspliced"))/
               (rowSums(assay(cdna_introns, "unspliced")) + 
                  rowSums(assay(cdna_introns, "spliced"))),
             stringsAsFactors = FALSE),
  data.frame(method = "cdna_intronscollapsed", 
             fraction_unspliced = rowSums(assay(cdna_intronscollapsed, "unspliced"))/
               (rowSums(assay(cdna_intronscollapsed, "unspliced")) + 
                  rowSums(assay(cdna_intronscollapsed, "spliced"))),
             stringsAsFactors = FALSE),
  data.frame(method = "spliced_unspliced", 
             fraction_unspliced = rowSums(assay(spliced_unspliced, "unspliced"))/
               (rowSums(assay(spliced_unspliced, "unspliced")) + 
                  rowSums(assay(spliced_unspliced, "spliced"))),
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_separate", 
             fraction_unspliced = rowSums(assay(kallistobus_separate, "unspliced"))/
               (rowSums(assay(kallistobus_separate, "unspliced")) + 
                  rowSums(assay(kallistobus_separate, "spliced"))),
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_collapse", 
             fraction_unspliced = rowSums(assay(kallistobus_collapse, "unspliced"))/
               (rowSums(assay(kallistobus_collapse, "unspliced")) + 
                  rowSums(assay(kallistobus_collapse, "spliced"))),
             stringsAsFactors = FALSE)
)

ggplot(df1, aes(x = method, y = fraction_unspliced)) + 
  geom_violin(aes(fill = method), alpha = 0.25, scale = "width") + 
  geom_boxplot(width = 0.025) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ggtitle("Fraction of total gene counts that are 'unspliced'")

## Correlation between "spliced" counts and counts from 
## quantification using only spliced transcripts
spliced_counts <- assay(spliced, "counts")
cdna_introns_counts <- assay(cdna_introns, "spliced")
cdna_intronscollapsed_counts <- assay(cdna_intronscollapsed, "spliced")
cdna_introns_decoy_counts <- assay(cdna_introns_decoy, "spliced")
spliced_unspliced_counts <- assay(spliced_unspliced, "spliced")
kallistobus_separate_counts <- assay(kallistobus_separate, "spliced")
kallistobus_collapse_counts <- assay(kallistobus_collapse, "spliced")
velocyto_counts <- assay(velocyto, "spliced")
ncells <- 100
df2 <- dplyr::bind_rows(
  data.frame(method = "velocyto", 
             correlation = sapply(seq_len(ncol(spliced))[1:ncells], function(i) {
               cor(spliced_counts[, i],
                   velocyto_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE),
  data.frame(method = "cdna_introns_decoy", 
             correlation = sapply(seq_len(ncol(spliced))[1:ncells], function(i) {
               cor(spliced_counts[, i],
                   cdna_introns_decoy_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE),
  data.frame(method = "cdna_introns", 
             correlation = sapply(seq_len(ncol(spliced))[1:ncells], function(i) {
               cor(spliced_counts[, i],
                   cdna_introns_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE),
  data.frame(method = "cdna_intronscollapsed", 
             correlation = sapply(seq_len(ncol(spliced))[1:ncells], function(i) {
               cor(spliced_counts[, i],
                   cdna_intronscollapsed_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE),
  data.frame(method = "spliced_unspliced", 
             correlation = sapply(seq_len(ncol(spliced))[1:ncells], function(i) {
               cor(spliced_counts[, i],
                   spliced_unspliced_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_separate", 
             correlation = sapply(seq_len(ncol(spliced))[1:ncells], function(i) {
               cor(spliced_counts[, i],
                   kallistobus_separate_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_collapse", 
             correlation = sapply(seq_len(ncol(spliced))[1:ncells], function(i) {
               cor(spliced_counts[, i],
                   kallistobus_collapse_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE)
)

ggplot(df2, aes(x = method, y = correlation)) + 
  geom_violin(aes(fill = method), alpha = 0.5, scale = "width") + 
  geom_boxplot(width = 0.025) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  labs(title = "Cell-wise correlation of 'spliced' counts with\ncounts obtained when only considering spliced transcripts",
       subtitle = paste0(ncells, " cells"))
dev.off()

saveRDS(list(velocyto = velocyto, cdna_introns_decoy = cdna_introns_decoy, 
             cdna_introns = cdna_introns, cdna_intronscollapsed = cdna_intronscollapsed,
             spliced_unspliced = spliced_unspliced, spliced = spliced, 
             kallistobus_separate = kallistobus_separate, 
             kallistobus_collapse = kallistobus_collapse), 
        file = outrds)

date()
sessionInfo()
