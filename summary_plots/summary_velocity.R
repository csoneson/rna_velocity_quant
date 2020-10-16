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


topdir <- ".."
plothelperscript <- paste0(topdir, "/scripts/plot_helpers.R")
datasets <- c("spermatogenesis_mouse", "dentate_gyrus_mouse",
              "pancreas_mouse", "pfc_mouse", "oldbrain_mouse")

source(plothelperscript)

files <- lapply(datasets, function(ds) {
  readRDS(file.path(topdir, ds, "plots/velocity_confidence", 
                    paste0(gsub("_10k_v3", "", gsub("_mouse", "", ds)), 
                           "_plot_velocity_confidence.rds")))
})

cellinfo <- do.call(dplyr::bind_rows, lapply(files, function(f) f$cell_vel)) %>%
  dplyr::mutate(dataset = gsub("_", " ", dataset))
geneinfo <- do.call(dplyr::bind_rows, lapply(files, function(f) f$gene_vel)) %>%
  dplyr::mutate(dataset = gsub("_", " ", dataset))

## Average scaled dot product with velocities of nearest neighbors
vfiles <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  readRDS(file.path(topdir, ds, "plots/compare_dimred_velocities",
                    paste0(gsub("_10k_v3", "", gsub("_mouse", "", ds)),
                           "_plot_compare_dimred_velocities.rds"))) %>%
    dplyr::filter(dimred == "UMAP_alevin_spliced_gentrome") %>%
    dplyr::select(cell, dataset, method_short, mtype, rtype, ave_dot_scaled) %>%
    dplyr::rename(index = cell) %>%
    dplyr::mutate(dataset = gsub("_", "", dataset))
}))
cellinfo <- cellinfo %>% 
  dplyr::left_join(vfiles, by = c("index", "dataset", "method_short", "mtype", "rtype"))

pdf("velocity_sd_latent_time.pdf")
ggplot(cellinfo %>% dplyr::filter(dataset != "Neuron") %>%
         dplyr::group_by(dataset, method_short, mtype, clusters) %>% 
         dplyr::summarize(ltsd = sd(latent_time)), 
       aes(x = dataset, y = ltsd)) + 
  geom_boxplot(outlier.size = -1, aes(fill = mtype)) + 
  geom_jitter(width = 0.2, height = 0) + 
  facet_wrap(~ method_short) + 
  theme_bw() + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 8)) + 
  labs(x = "", 
       y = "Within-cell type standard deviation of latent time estimate")
dev.off()  

pdf("velocity_max_cosine_corr_across_datasets.pdf")
ggplot(cellinfo %>% dplyr::filter(dataset != "Neuron"), 
       aes(x = dataset, y = max_cosine_corr)) + 
  geom_boxplot(outlier.size = 0.1, aes(fill = mtype)) + 
  facet_wrap(~ method_short) + 
  theme_bw() + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 8)) + 
  labs(x = "", 
       y = "Maximal cosine correlation by cell")
dev.off()  

pdf("velocity_ave_cos_theta_velodimred_across_datasets.pdf")
ggplot(cellinfo %>% dplyr::filter(dataset != "Neuron"), 
       aes(x = dataset, y = ave_dot_scaled)) + 
  geom_boxplot(outlier.size = 0.1, aes(fill = mtype)) + 
  facet_wrap(~ method_short) + 
  theme_bw() + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 8)) + 
  labs(x = "", 
       y = "Average cosine of angle between velocity embedding\nand velocity embeddings of neighboring cells")
dev.off()  

pdf("ms_vs_mu_est_r2_across_datasets.pdf")
## Calculate K-S statistic for est_r2 between null and non-null data sets
kscomp <- geneinfo %>% dplyr::filter(dataset != "Neuron") %>% 
  dplyr::select(dataset, est_r2, method_short) %>% 
  dplyr::group_by(method_short) %>% 
  dplyr::summarize(ksstat = ks.test(est_r2[dataset %in% c("PFC", "OldBrain")], 
                                    est_r2[!(dataset %in% c("PFC", "OldBrain"))])$stat) %>%
  dplyr::mutate(ksstat = paste0("KS = ", signif(ksstat, digits = 2)))
ggplot(geneinfo %>% dplyr::filter(dataset != "Neuron"), 
       aes(x = dataset, y = est_r2)) + 
  # geom_violin(aes(fill = mtype), scale = "width") + 
  geom_boxplot(outlier.size = 0.1, aes(fill = mtype)) +
  facet_wrap(~ method_short) + 
  geom_text(data = kscomp, aes(x = 4.5, y = 1.13, label = ksstat), size = 3) + 
  theme_bw() + coord_cartesian(ylim = c(-1, 1.2)) + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 8)) + 
  labs(x = "", 
       y = "Fit R2, unspliced vs spliced abundance")
dev.off()  

pdf("summary_velocity.pdf", width = 8, height = 8)
g1 <- ggplot(cellinfo, aes(x = method_short, y = velocity_confidence)) + 
  geom_boxplot(aes(fill = mtype), alpha = 0.5, size = 0.5) + 
  theme_bw() + facet_wrap(~ dataset, ncol = 1) + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Velocity confidence, per cell", y = "Velocity confidence", x = "")

g2 <- ggplot(geneinfo, aes(x = method_short, y = fit_likelihood)) + 
  geom_boxplot(aes(fill = mtype), alpha = 0.5, size = 0.5) + 
  theme_bw() + facet_wrap(~ dataset, ncol = 1) + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Fit likelihood, per gene", y = "Fit likelihood", x = "")

cowplot::plot_grid(g1, g2, nrow = 1, labels = c("A", "B"))
dev.off()

date()
sessionInfo()

