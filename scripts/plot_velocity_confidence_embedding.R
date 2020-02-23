args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(cowplot)
  library(Seurat)
})
source(plothelperscript)

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(plothelperscript)
print(topdir)
print(methods)
print(outrds)

methods_short <- shorten_methods(methods) %>%
  dplyr::filter(method %in% methods) %>%
  dplyr::arrange(as.factor(method_short))
methods <- methods_short$method
names(methods) <- methods

## ------------------------------------------------------------------------- ##
## Read data
## ------------------------------------------------------------------------- ##
sce <- readRDS(file.path(topdir, paste0("output/sce/sce_starsolo.rds")))
umap <- reducedDim(sce, "UMAP_alevin_spliced_gentrome")

cellinfo_all <- lapply(methods, function(nm) {
  readr::read_csv(paste0(topdir, "/plots/velocity/anndata_", nm, "/anndata_", 
                         nm, "_cell_info.csv")) %>%
    dplyr::mutate(method = nm)
})
cellinfo_shared <- lapply(methods, function(nm) {
  readr::read_csv(paste0(topdir, "/plots/velocity/anndata_", nm, 
                         "_shared_genes/anndata_", nm, 
                         "_shared_genes_cell_info.csv")) %>%
    dplyr::mutate(method = nm)
})

## Shared genes
plotlist_shared <- lapply(methods, function(m) {
  sm <- methods_short$method_short[match(m, methods_short$method)]
  df <- data.frame(umap)
  df$confidence <- 
    cellinfo_shared[[m]]$velocity_confidence[match(rownames(df), cellinfo_shared[[m]]$index)]
  ggplot(df, aes(x = X1, y = X2, fill = confidence)) + 
    geom_point(shape = 21) + 
    scale_fill_gradient(low = "white", high = "red", na.value = "grey50", limits = c(0, 1),
                        name = "Velocity confidence") + 
    labs(x = "", y = "", title = sm) + 
    theme_void() + theme(legend.position = "bottom")
})

png(gsub("\\.rds", "_shared_genes.png", outrds), width = 8, height = 12, unit = "in", res = 200)
print(cowplot::plot_grid(
  cowplot::plot_grid(plotlist = lapply(plotlist_shared, function(w) w + theme(legend.position = "none")), 
                     ncol = 3),
  cowplot::get_legend(plotlist_shared[[1]] + theme(legend.position = "bottom")),
  ncol = 1, rel_heights = c(1, 0.05)
))
dev.off()

## All genes
plotlist_all <- lapply(methods, function(m) {
  sm <- methods_short$method_short[match(m, methods_short$method)]
  df <- data.frame(umap)
  df$confidence <- 
    cellinfo_all[[m]]$velocity_confidence[match(rownames(df), cellinfo_all[[m]]$index)]
  ggplot(df, aes(x = X1, y = X2, fill = confidence)) + 
    geom_point(shape = 21) + 
    scale_fill_gradient(low = "white", high = "red", na.value = "grey50", limits = c(0, 1),
                        name = "Velocity confidence") + 
    labs(x = "", y = "", title = sm) + 
    theme_void() + theme(legend.position = "bottom")
})

png(gsub("\\.rds", ".png", outrds), width = 8, height = 12, unit = "in", res = 200)
print(cowplot::plot_grid(
  cowplot::plot_grid(plotlist = lapply(plotlist_all, function(w) w + theme(legend.position = "none")), 
                     ncol = 3),
  cowplot::get_legend(plotlist_all[[1]] + theme(legend.position = "bottom")),
  ncol = 1, rel_heights = c(1, 0.05)
))
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()

