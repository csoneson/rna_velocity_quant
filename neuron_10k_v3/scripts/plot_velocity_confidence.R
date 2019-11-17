args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(scater)
  library(readr)
  library(UpSetR)
  library(pheatmap)
})

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(topdir)
print(methods)
print(outrds)

cellinfo <- lapply(methods, function(nm) {
  readr::read_csv(paste0(topdir, "/plots/velocity/anndata_", nm, "/anndata_", nm, 
                         "_cell_info.csv")) %>%
    dplyr::mutate(method = nm)
})
geneinfo <- lapply(methods, function(nm) {
  readr::read_csv(paste0(topdir, "/plots/velocity/anndata_", nm, "/anndata_", nm, 
                         "_gene_info.csv")) %>%
    dplyr::mutate(method = nm)
})

velocity_confidence <- do.call(dplyr::bind_rows, lapply(cellinfo, function(w) {
  w %>% dplyr::select(index, method, velocity_confidence, 
                      velocity_confidence_transition, velocity_self_transition, 
                      velocity_length)
}))

velocity_genes <- do.call(dplyr::bind_rows, lapply(geneinfo, function(w) {
  w %>% dplyr::select(index, method, velocity_r2, 
                      velocity_score)
}))

pdf(gsub("rds$", "pdf", outrds), width = 10, height = 10)
pheatmap(cor(velocity_confidence %>% dplyr::select(index, method, velocity_length) %>%
               tidyr::spread(key = method, value = velocity_length) %>%
               tibble::column_to_rownames("index")), 
         cluster_rows = TRUE, cluster_cols = TRUE, main = "Velocity length")

ggplot(velocity_confidence, aes(x = method, y = velocity_confidence)) + 
  geom_violin(aes(fill = method), alpha = 0.25, scale = "width") + 
  geom_boxplot(width = 0.025) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Velocity confidence, per cell",
       subtitle = paste0("The velocity confidence checks whether the velocities ", 
                         "agree with the neighboring velocities, \ni.e. ", 
                         "a kind of smoothness score"))

ggplot(velocity_confidence, aes(x = method, y = velocity_confidence_transition)) + 
  geom_violin(aes(fill = method), alpha = 0.25, scale = "width") + 
  geom_boxplot(width = 0.025) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Velocity confidence (transition), per cell",
       subtitle = paste0("The velocity transition confidence checks whether ", 
                         "the embedded velocities obtained from the \ntransition ", 
                         "probabilities truly reflect the velocities ", 
                         "in high dimensional space."))

ggplot(velocity_confidence, aes(x = method, y = velocity_self_transition)) + 
  geom_violin(aes(fill = method), alpha = 0.25, scale = "width") + 
  geom_boxplot(width = 0.025) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Velocity self-transition, per cell",
       subtitle = paste0(""))

ggplot(velocity_genes, aes(x = method, y = velocity_r2)) + 
  geom_violin(aes(fill = method), alpha = 0.25, scale = "width") + 
  geom_boxplot(width = 0.025) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Velocity r2, per gene",
       subtitle = paste0(""))

ggplot(velocity_genes, aes(x = method, y = velocity_score)) + 
  geom_violin(aes(fill = method), alpha = 0.25, scale = "width") + 
  geom_boxplot(width = 0.025) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Velocity score, per gene",
       subtitle = paste0(""))

allgenes <- as.character(Reduce(union, lapply(geneinfo, function(w) w$index)))
selgenes <- do.call(cbind, lapply(geneinfo, function(w) as.numeric(allgenes %in% w$index)))
rownames(selgenes) <- allgenes
selgenes <- data.frame(selgenes)
upset(selgenes, nsets = length(geneinfo), keep.order = TRUE, 
      order.by = "freq", decreasing = TRUE)

upset(selgenes %>% dplyr::select(-contains("busparse")), 
      nsets = length(geneinfo), keep.order = TRUE, 
      order.by = "freq", decreasing = TRUE)

dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
