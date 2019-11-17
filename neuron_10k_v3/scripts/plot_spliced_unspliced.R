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

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(topdir)
print(gtf)
print(methods)
print(outrds)

sces <- lapply(methods, function(nm) {
  readRDS(file.path(topdir, paste0("output/sce_", nm, ".rds")))
})
for (nm in names(sces)) {
  metadata(sces[[nm]]) <- list(dataset = nm)
}
alevin_spliced <- readRDS(file.path(topdir, "output/sce_alevin_spliced.rds"))
metadata(alevin_spliced) <- list(dataset = "alevin_spliced")

gtf <- rtracklayer::import(gtf)

## ========================================================================= ##
## Plot
## ========================================================================= ##
df0 <- do.call(dplyr::bind_rows, lapply(sces, function(w) {
  data.frame(method = metadata(w)$dataset,
             spliced_count = rowSums(assay(w, "spliced")),
             unspliced_count = rowSums(assay(w, "unspliced")),
             total_count = rowSums(assay(w, "spliced")) + 
               rowSums(assay(w, "unspliced")),
             stringsAsFactors = FALSE)
})) %>% 
  tidyr::gather(key = "type", value = "numi", -method)

pdf(gsub("rds", "pdf", outrds), width = 11, height = 11)
ggplot(df0 %>% dplyr::mutate(
  type = factor(type, levels = c("spliced_count", "unspliced_count", "total_count"))), 
  aes(x = method, y = numi + 1)) + 
  geom_violin(aes(fill = method), alpha = 0.25, scale = "width") + 
  geom_boxplot(width = 0.025) + 
  facet_wrap(~ type, ncol = 1) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_log10() + 
  ggtitle("Total count per gene, across all cells")

df1 <- do.call(dplyr::bind_rows, lapply(sces, function(w) {
  data.frame(method = metadata(w)$dataset,
             fraction_unspliced = rowSums(assay(w, "unspliced"))/
               (rowSums(assay(w, "unspliced")) + 
                  rowSums(assay(w, "spliced"))),
             stringsAsFactors = FALSE)
}))

ggplot(df1, aes(x = method, y = fraction_unspliced)) + 
  geom_violin(aes(fill = method), alpha = 0.25, scale = "width") + 
  geom_boxplot(width = 0.025) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  ggtitle("Fraction of total gene counts that are 'unspliced'")

## Correlation between "spliced" counts and counts from 
## quantification using only spliced transcripts
alevin_spliced_counts <- assay(alevin_spliced, "counts")
ncells <- 100
df2 <- do.call(dplyr::bind_rows, lapply(sces, function(w) {
  tmp <- assay(w, "spliced")
  data.frame(method = metadata(w)$dataset, 
             correlation = sapply(seq_len(ncol(alevin_spliced_counts))[1:ncells], function(i) {
               cor(alevin_spliced_counts[, i],
                   tmp[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE)
}))

ggplot(df2, aes(x = method, y = correlation)) + 
  geom_violin(aes(fill = method), alpha = 0.5, scale = "width") + 
  geom_boxplot(width = 0.025) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(title = "Cell-wise correlation of 'spliced' counts with\ncounts obtained when only considering spliced transcripts",
       subtitle = paste0(ncells, " cells"))


## Correlation between spliced and unspliced counts
ncells <- 100
df3 <- do.call(dplyr::bind_rows, lapply(sces, function(w) {
  sct <- assay(w, "spliced")
  uct <- assay(w, "unspliced")
  data.frame(method = metadata(w)$dataset,
             correlation = sapply(seq_len(ncol(w))[1:ncells], function(i) {
               cor(sct[, i],
                   uct[, i],
                   method = "spearman")
             }),
             stringsAsFactors = FALSE)
}))

ggplot(df3, aes(x = method, y = correlation)) + 
  geom_violin(aes(fill = method), alpha = 0.5, scale = "width") + 
  geom_boxplot(width = 0.025) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(title = "Cell-wise correlation of 'spliced' and 'unspliced' counts",
       subtitle = paste0(ncells, " cells"))


dev.off()
saveRDS(NULL, file = outrds)

date()
sessionInfo()