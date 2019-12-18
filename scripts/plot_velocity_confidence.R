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
print(plothelperscript)
print(dataset)
print(methods)
print(outrds)

source(plothelperscript)

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
  w %>% dplyr::select(index, method, contains("clusters"), velocity_confidence, 
                      velocity_confidence_transition, velocity_self_transition, 
                      velocity_length, velocity_pseudotime, latent_time)
}))
if ("clusters" %in% colnames(velocity_confidence) && !is.null(cluster_levels[[dataset]])) {
  velocity_confidence$clusters <- factor(
    as.character(velocity_confidence$clusters), 
    levels = cluster_levels[[dataset]]
  )
}

velocity_genes <- do.call(dplyr::bind_rows, lapply(geneinfo, function(w) {
  w %>% dplyr::select(index, method, 
                      velocity_score, fit_alpha, fit_beta, fit_gamma, 
                      fit_likelihood)
}))

pdf(gsub("rds$", "pdf", outrds), width = 10, height = 10)

if ("clusters" %in% colnames(velocity_confidence)) {
  ggplot(velocity_confidence, aes(x = clusters, y = latent_time)) + 
    geom_violin(aes(fill = clusters), alpha = 0.25) + theme_bw() + 
    facet_wrap(~ method) + 
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    labs(title = "Latent time vs cluster")
}

ggplot(velocity_confidence, aes(x = method, y = velocity_confidence)) + 
  geom_boxplot(aes(fill = method), alpha = 0.25) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Velocity confidence, per cell",
       subtitle = paste0("The velocity confidence checks whether the velocities ", 
                         "agree with the neighboring velocities, \ni.e. ", 
                         "a kind of smoothness score"))

ggplot(velocity_confidence, aes(x = method, y = velocity_length)) + 
  geom_boxplot(aes(fill = method), alpha = 0.25) + 
  theme_bw() + 
  scale_y_log10() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Velocity length, per cell",
       subtitle = "")

ggplot(velocity_confidence, aes(x = method, y = velocity_confidence_transition)) + 
  geom_boxplot(aes(fill = method), alpha = 0.25) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Velocity confidence (transition), per cell",
       subtitle = paste0("The velocity transition confidence checks whether ", 
                         "the embedded velocities obtained from the \ntransition ", 
                         "probabilities truly reflect the velocities ", 
                         "in high dimensional space."))

ggplot(velocity_confidence, aes(x = method, y = velocity_self_transition)) + 
  geom_boxplot(aes(fill = method), alpha = 0.25) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Velocity self-transition, per cell",
       subtitle = paste0(""))

ggplot(velocity_genes, aes(x = method, y = velocity_score)) + 
  geom_boxplot(aes(fill = method), alpha = 0.25) + 
  scale_y_sqrt() + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Velocity score, per gene",
       subtitle = paste0(""))

ggplot(velocity_genes, aes(x = method, y = fit_likelihood)) + 
  geom_boxplot(aes(fill = method), alpha = 0.25, scale = "width") + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Fit likelihood, per gene",
       subtitle = paste0(""))

dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()