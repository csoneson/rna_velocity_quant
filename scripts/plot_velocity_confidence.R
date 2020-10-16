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
  library(Seurat)
})

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(topdir)
print(velosuffix)
print(plothelperscript)
print(dataset)
print(methods)
print(outrds)

source(plothelperscript)

methods_short <- shorten_methods(methods)

cellinfo <- lapply(methods, function(nm) {
  readr::read_csv(paste0(topdir, "/plots/velocity", velosuffix, 
                         "/anndata_", nm, "/anndata_", nm, 
                         "_cell_info.csv")) %>%
    dplyr::mutate(method = nm)
})
geneinfo <- lapply(methods, function(nm) {
  gi <- readr::read_csv(paste0(topdir, "/plots/velocity", velosuffix, 
                         "/anndata_", nm, "/anndata_", nm, 
                         "_gene_info.csv")) %>%
    dplyr::mutate(method = nm)

  ## Estimate r2 value for line (through origin) fit to Mu vs Ms
  w <- ReadH5AD(file.path(topdir, paste0(
    "output/anndata_with_velocity", velosuffix, "/anndata_", 
    nm, "_with_velocity.h5ad")))
  Ms <- as.matrix(GetAssayData(GetAssay(w, "Ms")))
  Mu <- as.matrix(GetAssayData(GetAssay(w, "Mu")))
  calc_gamma <- vapply(seq_len(nrow(Ms)), function(i) {
    sum(Ms[i, ] * Mu[i, ])/sum(Ms[i, ] * Ms[i, ])
  }, 0)
  est_r2 <- vapply(seq_len(nrow(Ms)), function(i) {
    sstotal <- sum((Mu[i, ] - mean(Mu[i, ])) ^ 2)
    ssregr <- sum((Mu[i, ] - calc_gamma[i] * Ms[i, ]) ^ 2)
    r2 <- 1 - ssregr/sstotal
    r2
  }, 0)
  
  gi %>% dplyr::left_join(
    data.frame(index = rownames(Ms),
               est_r2 = est_r2),
    by = "index")
})


velocity_confidence <- do.call(dplyr::bind_rows, lapply(cellinfo, function(w) {
  w %>% dplyr::select(index, method, contains("clusters"), velocity_confidence, 
                      velocity_confidence_transition, max_cosine_corr, 
                      velocity_length, velocity_pseudotime, latent_time)
})) %>%
  dplyr::left_join(methods_short, by = "method")
if ("clusters" %in% colnames(velocity_confidence) && !is.null(cluster_levels[[dataset]])) {
  velocity_confidence$clusters <- factor(
    as.character(velocity_confidence$clusters), 
    levels = cluster_levels[[dataset]]
  )
}

velocity_genes <- do.call(dplyr::bind_rows, lapply(geneinfo, function(w) {
  w %>% dplyr::select(index, method, 
                      velocity_score, fit_alpha, fit_beta, fit_gamma, 
                      fit_likelihood, est_r2, contains("fit_r2"))
})) %>%
  dplyr::left_join(methods_short, by = "method")

pdf(gsub("rds$", "pdf", outrds), width = 10, height = 10)

if ("clusters" %in% colnames(velocity_confidence)) {
  ggplot(velocity_confidence, aes(x = clusters, y = latent_time)) + 
    geom_violin(aes(fill = clusters), alpha = 0.25) + theme_bw() + 
    facet_wrap(~ method) + 
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    labs(title = paste0(gsub("_", " ", dataset), ", latent time vs cluster"))
}

ggplot(velocity_confidence, aes(x = method_short, y = velocity_confidence)) + 
  geom_boxplot(aes(fill = mtype), alpha = 0.25) + 
  theme_bw() + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = paste0(gsub("_", " ", dataset), ", velocity confidence, per cell"),
       subtitle = paste0("The velocity confidence checks whether the velocities ", 
                         "agree with the neighboring velocities, \ni.e. ", 
                         "a kind of smoothness score"))

ggplot(velocity_confidence, aes(x = method_short, y = velocity_length)) + 
  geom_boxplot(aes(fill = mtype), alpha = 0.25) + 
  theme_bw() + 
  scale_y_log10() + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = paste0(gsub("_", " ", dataset), ", velocity length, per cell"),
       subtitle = "")

ggplot(velocity_confidence, aes(x = method_short, y = velocity_confidence_transition)) + 
  geom_boxplot(aes(fill = mtype), alpha = 0.25) + 
  theme_bw() + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = paste0(gsub("_", " ", dataset), ", velocity confidence (transition), per cell"),
       subtitle = paste0("The velocity transition confidence checks whether ", 
                         "the embedded velocities obtained from the \ntransition ", 
                         "probabilities truly reflect the velocities ", 
                         "in high dimensional space."))

ggplot(velocity_confidence, aes(x = method_short, y = max_cosine_corr)) + 
  geom_boxplot(aes(fill = mtype), alpha = 0.25) + 
  theme_bw() + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = paste0(gsub("_", " ", dataset), ", max cosine correlation, per cell"),
       subtitle = paste0(""))

if ("clusters" %in% colnames(velocity_confidence)) {
  print(
    ggplot(velocity_confidence, aes(x = method_short, y = max_cosine_corr)) + 
      geom_boxplot(aes(fill = mtype), alpha = 0.25) + 
      theme_bw() + 
      facet_wrap(~ clusters) + 
      scale_fill_manual(values = base_method_colors, name = "") + 
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = paste0(gsub("_", " ", dataset), ", max cosine correlation, per cell"),
           x = "",
           y = "Max cosine correlation"))
  
  print(
    ggplot(velocity_confidence, aes(x = clusters, y = max_cosine_corr)) + 
      geom_boxplot(aes(fill = mtype), alpha = 0.25) + 
      theme_bw() + 
      facet_wrap(~ method_short) + 
      scale_fill_manual(values = base_method_colors, name = "") + 
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = paste0(gsub("_", " ", dataset), ", max cosine correlation, per cell"),
           x = "",
           y = "Max cosine correlation"))
}

ggplot(velocity_genes, aes(x = method_short, y = velocity_score)) + 
  geom_boxplot(aes(fill = mtype), alpha = 0.25) + 
  scale_y_sqrt() + 
  theme_bw() + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = paste0(gsub("_", " ", dataset), ", velocity score, per gene"),
       subtitle = paste0(""))

ggplot(velocity_genes, aes(x = method_short, y = fit_likelihood)) + 
  geom_boxplot(aes(fill = mtype), alpha = 0.25, scale = "width") + 
  theme_bw() + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = paste0(gsub("_", " ", dataset), ", fit likelihood, per gene"),
       subtitle = paste0(""))

if ("fit_r2" %in% colnames(velocity_genes)) {
  print(ggplot(velocity_genes, aes(x = method_short, y = fit_r2)) + 
          geom_boxplot(aes(fill = mtype), alpha = 0.25) + 
          theme_bw() + coord_cartesian(ylim = c(-1, NA)) + 
          scale_fill_manual(values = base_method_colors, name = "") + 
          theme(legend.position = "none",
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
          labs(title = paste0(gsub("_", " ", dataset), ", fit R2, per gene"),
               subtitle = paste0("")))
}

print(ggplot(velocity_genes, aes(x = method_short, y = est_r2)) + 
        geom_boxplot(aes(fill = mtype), alpha = 0.25) + 
        theme_bw() + coord_cartesian(ylim = c(-1, NA)) + 
        scale_fill_manual(values = base_method_colors, name = "") + 
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        labs(title = paste0(gsub("_", " ", dataset), ", est R2, per gene"),
             subtitle = paste0("")))

dev.off()

saveRDS(list(cell_vel = velocity_confidence %>% dplyr::mutate(dataset = dataset),
             gene_vel = velocity_genes %>% dplyr::mutate(dataset = dataset)), 
        file = outrds)
date()
sessionInfo()
