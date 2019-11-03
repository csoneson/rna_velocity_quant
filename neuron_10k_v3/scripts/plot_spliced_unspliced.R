args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(tximeta)
  library(dplyr)
  library(ggplot2)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(scater)
})

print(topdir)
print(outrds)

velocyto <- readRDS(file.path(topdir, "output/sce_velocyto.rds"))
cdna_introns_decoy <- readRDS(file.path(topdir, "output/sce_cdna_introns_decoy.rds"))
cdna_introns <- readRDS(file.path(topdir, "output/sce_cdna_introns.rds"))
cdna_intronscollapsed <- readRDS(file.path(topdir, "output/sce_cdna_intronscollapsed.rds"))
spliced_unspliced <- readRDS(file.path(topdir, "output/sce_spliced_unspliced.rds"))
spliced <- readRDS(file.path(topdir, "output/sce_spliced.rds"))
kallistobus_separate <- readRDS(file.path(topdir, "output/sce_kallistobus_separate.rds"))
kallistobus_collapse <- readRDS(file.path(topdir, "output/sce_kallistobus_collapse.rds"))
kallistobus_sep_permissive <- readRDS(file.path(topdir, "output/sce_kallistobus_sep_permissive.rds"))
kallistobus_coll_permissive <- readRDS(file.path(topdir, "output/sce_kallistobus_coll_permissive.rds"))

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
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_sep_permissive",
             spliced_count = rowSums(assay(kallistobus_sep_permissive, "spliced")),
             unspliced_count = rowSums(assay(kallistobus_sep_permissive, "unspliced")),
             total_count = rowSums(assay(kallistobus_sep_permissive, "spliced")) + 
               rowSums(assay(kallistobus_sep_permissive, "unspliced")),
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_coll_permissive",
             spliced_count = rowSums(assay(kallistobus_coll_permissive, "spliced")),
             unspliced_count = rowSums(assay(kallistobus_coll_permissive, "unspliced")),
             total_count = rowSums(assay(kallistobus_coll_permissive, "spliced")) + 
               rowSums(assay(kallistobus_coll_permissive, "unspliced")),
             stringsAsFactors = FALSE)
) %>% 
  tidyr::gather(key = "type", value = "numi", -method)

pdf(gsub("rds", "pdf", outrds), width = 11, height = 9)
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
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_sep_permissive", 
             fraction_unspliced = rowSums(assay(kallistobus_sep_permissive, "unspliced"))/
               (rowSums(assay(kallistobus_sep_permissive, "unspliced")) + 
                  rowSums(assay(kallistobus_sep_permissive, "spliced"))),
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_coll_permissive", 
             fraction_unspliced = rowSums(assay(kallistobus_coll_permissive, "unspliced"))/
               (rowSums(assay(kallistobus_coll_permissive, "unspliced")) + 
                  rowSums(assay(kallistobus_coll_permissive, "spliced"))),
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
kallistobus_sep_permissive_counts <- assay(kallistobus_sep_permissive, "spliced")
kallistobus_coll_permissive_counts <- assay(kallistobus_coll_permissive, "spliced")
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
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_sep_permissive", 
             correlation = sapply(seq_len(ncol(spliced))[1:ncells], function(i) {
               cor(spliced_counts[, i],
                   kallistobus_sep_permissive_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_coll_permissive", 
             correlation = sapply(seq_len(ncol(spliced))[1:ncells], function(i) {
               cor(spliced_counts[, i],
                   kallistobus_coll_permissive_counts[, i],
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

## Correlation between spliced and unspliced counts
cdna_introns_ucounts <- assay(cdna_introns, "unspliced")
cdna_intronscollapsed_ucounts <- assay(cdna_intronscollapsed, "unspliced")
cdna_introns_decoy_ucounts <- assay(cdna_introns_decoy, "unspliced")
spliced_unspliced_ucounts <- assay(spliced_unspliced, "unspliced")
kallistobus_separate_ucounts <- assay(kallistobus_separate, "unspliced")
kallistobus_collapse_ucounts <- assay(kallistobus_collapse, "unspliced")
kallistobus_sep_permissive_ucounts <- assay(kallistobus_sep_permissive, "unspliced")
kallistobus_coll_permissive_ucounts <- assay(kallistobus_coll_permissive, "unspliced")
velocyto_ucounts <- assay(velocyto, "spliced")
ncells <- 100
df3 <- dplyr::bind_rows(
  data.frame(method = "velocyto", 
             correlation = sapply(seq_len(ncol(velocyto))[1:ncells], function(i) {
               cor(velocyto_ucounts[, i],
                   velocyto_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE),
  data.frame(method = "cdna_introns_decoy", 
             correlation = sapply(seq_len(ncol(cdna_introns_decoy))[1:ncells], function(i) {
               cor(cdna_introns_decoy_ucounts[, i],
                   cdna_introns_decoy_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE),
  data.frame(method = "cdna_introns", 
             correlation = sapply(seq_len(ncol(cdna_introns))[1:ncells], function(i) {
               cor(cdna_introns_ucounts[, i],
                   cdna_introns_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE),
  data.frame(method = "cdna_intronscollapsed", 
             correlation = sapply(seq_len(ncol(cdna_intronscollapsed))[1:ncells], function(i) {
               cor(cdna_intronscollapsed_ucounts[, i],
                   cdna_intronscollapsed_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE),
  data.frame(method = "spliced_unspliced", 
             correlation = sapply(seq_len(ncol(spliced_unspliced))[1:ncells], function(i) {
               cor(spliced_unspliced_ucounts[, i],
                   spliced_unspliced_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_separate", 
             correlation = sapply(seq_len(ncol(kallistobus_separate))[1:ncells], function(i) {
               cor(kallistobus_separate_ucounts[, i],
                   kallistobus_separate_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_collapse", 
             correlation = sapply(seq_len(ncol(kallistobus_collapse))[1:ncells], function(i) {
               cor(kallistobus_collapse_ucounts[, i],
                   kallistobus_collapse_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_sep_permissive", 
             correlation = sapply(seq_len(ncol(kallistobus_sep_permissive))[1:ncells], function(i) {
               cor(kallistobus_sep_permissive_ucounts[, i],
                   kallistobus_sep_permissive_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE),
  data.frame(method = "kallistobus_coll_permissive", 
             correlation = sapply(seq_len(ncol(kallistobus_coll_permissive))[1:ncells], function(i) {
               cor(kallistobus_coll_permissive_ucounts[, i],
                   kallistobus_coll_permissive_counts[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE)
)

ggplot(df3, aes(x = method, y = correlation)) + 
  geom_violin(aes(fill = method), alpha = 0.5, scale = "width") + 
  geom_boxplot(width = 0.025) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  labs(title = "Cell-wise correlation of 'spliced' and 'unspliced' counts",
       subtitle = paste0(ncells, " cells"))
dev.off()


saveRDS(NULL, file = outrds)

date()
sessionInfo()