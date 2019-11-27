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
  library(cowplot)
})

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(topdir)
print(gtf)
print(methods)
print(outrds)

sces <- lapply(methods, function(nm) {
  readRDS(file.path(topdir, paste0("output/sce/sce_", nm, ".rds")))
})
for (nm in names(sces)) {
  metadata(sces[[nm]]) <- list(dataset = nm)
}
alevin_spliced <- readRDS(file.path(topdir, "output/sce/sce_alevin_spliced.rds"))
metadata(alevin_spliced) <- list(dataset = "alevin_spliced")

gtf <- rtracklayer::import(gtf)

## ========================================================================= ##
## Plot
## ========================================================================= ##
pdf(gsub("rds", "pdf", outrds), width = 11, height = 11)

## Total number of assigned reads
df0 <- do.call(dplyr::bind_rows, lapply(sces, function(w) {
  data.frame(method = metadata(w)$dataset, 
             spliced_count = sum(assay(w, "spliced")),
             unspliced_count = sum(assay(w, "unspliced")),
             stringsAsFactors = FALSE) %>%
    dplyr::mutate(total_count = spliced_count + unspliced_count)
})) %>%
  tidyr::gather(key = "type", value = "numi", -method)

ggplot(df0 %>% dplyr::mutate(
  type = factor(type, levels = c("spliced_count", "unspliced_count", "total_count"))), 
  aes(x = method, y = numi, fill = type)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  facet_wrap(~ type, nrow = 1) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Total count across all genes and cells")

## Total number of reads, by gene, across all cells
df0 <- do.call(dplyr::bind_rows, lapply(sces, function(w) {
  data.frame(method = metadata(w)$dataset,
             spliced_count = rowSums(assay(w, "spliced")),
             unspliced_count = rowSums(assay(w, "unspliced")),
             total_count = rowSums(assay(w, "spliced")) + 
               rowSums(assay(w, "unspliced")),
             stringsAsFactors = FALSE)
})) %>% 
  tidyr::gather(key = "type", value = "numi", -method)

ggplot(df0 %>% dplyr::mutate(
  type = factor(type, levels = c("spliced_count", "unspliced_count", "total_count"))), 
  aes(x = method, y = numi + 1)) + 
  geom_violin(aes(fill = method), alpha = 0.25) + 
  facet_wrap(~ type, ncol = 1) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_log10() + 
  ggtitle("Total count per gene, across all cells")

## Total number of reads, by cell, across all genes
df0 <- do.call(dplyr::bind_rows, lapply(sces, function(w) {
  data.frame(method = metadata(w)$dataset,
             spliced_count = colSums(assay(w, "spliced")),
             unspliced_count = colSums(assay(w, "unspliced")),
             total_count = colSums(assay(w, "spliced")) + 
               colSums(assay(w, "unspliced")),
             stringsAsFactors = FALSE)
})) %>% 
  tidyr::gather(key = "type", value = "numi", -method)

ggplot(df0 %>% dplyr::mutate(
  type = factor(type, levels = c("spliced_count", "unspliced_count", "total_count"))), 
  aes(x = method, y = numi)) + 
  geom_violin(aes(fill = method), alpha = 0.25) + 
  facet_wrap(~ type, ncol = 1) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Total count per cell, across all genes")

## Fraction spliced, by gene, across all cells
df1 <- do.call(dplyr::bind_rows, lapply(sces, function(w) {
  data.frame(method = metadata(w)$dataset,
             fraction_unspliced = rowSums(assay(w, "unspliced"))/
               (rowSums(assay(w, "unspliced")) + 
                  rowSums(assay(w, "spliced"))),
             stringsAsFactors = FALSE)
}))

ggplot(df1, aes(x = method, y = fraction_unspliced)) + 
  geom_violin(aes(fill = method), alpha = 0.25) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  ggtitle("Fraction of total gene counts that are 'unspliced'")

## Fraction spliced, by cell, across all gene
df1 <- do.call(dplyr::bind_rows, lapply(sces, function(w) {
  data.frame(method = metadata(w)$dataset,
             fraction_unspliced = colSums(assay(w, "unspliced"))/
               (colSums(assay(w, "unspliced")) + 
                  colSums(assay(w, "spliced"))),
             stringsAsFactors = FALSE)
}))

ggplot(df1, aes(x = method, y = fraction_unspliced)) + 
  geom_violin(aes(fill = method), alpha = 0.25) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  ggtitle("Fraction of total cell counts that are 'unspliced'")


## Correlation between spliced and unspliced counts
df3 <- do.call(dplyr::bind_rows, lapply(sces, function(w) {
  sct <- as.matrix(assay(w, "spliced"))
  uct <- as.matrix(assay(w, "unspliced"))
  a <- scale(sct, center = TRUE, scale = TRUE)
  b <- scale(uct, center = TRUE, scale = TRUE)
  corrs <- colSums(a * b)/(nrow(a) - 1)
  data.frame(method = metadata(w)$dataset,
             correlation = corrs,
             stringsAsFactors = FALSE)
}))

ggplot(df3, aes(x = method, y = correlation)) + 
  geom_violin(aes(fill = method), alpha = 0.5) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(title = "Cell-wise correlation of 'spliced' and 'unspliced' counts")


dev.off()
saveRDS(NULL, file = outrds)

date()
sessionInfo()