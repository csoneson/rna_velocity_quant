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
print(showgene)
print(outvelopng)

methods_short <- shorten_methods(methods)

## ------------------------------------------------------------------------- ##
## Read data
## ------------------------------------------------------------------------- ##
sces <- lapply(methods, function(nm) {
  readRDS(file.path(topdir, paste0("output/sce/sce_", nm, ".rds")))
})
fracunspliced <- lapply(sces, function(w) {
  assay(w, "unspliced")/(assay(w, "spliced") + assay(w, "unspliced"))
})
cells <- colnames(sces[[1]])
stopifnot(all(sapply(fracunspliced, function(w) all(colnames(w) == cells))))

umap <- data.frame(reducedDim(sces[[1]], "UMAP_alevin_spliced"),
                   stringsAsFactors = FALSE)
stopifnot(rownames(umap) == cells)

seurats <- lapply(methods, function(nm) {
  ReadH5AD(file.path(topdir, paste0("output/anndata_with_velocity/anndata_", 
                                    nm, "_with_velocity.h5ad")))
})
velocities <- lapply(seurats, function(w) {
  as.matrix(GetAssayData(GetAssay(w, "velocity")))
})
stopifnot(all(sapply(velocities, function(w) all(colnames(w) == cells))))

## ------------------------------------------------------------------------- ##
## Fraction unspliced
## ------------------------------------------------------------------------- ##
plotlist <- lapply(methods, function(m) {
  sm <- methods_short$method_short[match(m, methods_short$method)]
  df <- umap
  df$fracunspl <- fracunspliced[[m]][match(showgene, rownames(fracunspliced[[m]])), ]
  ggplot(df, aes(x = X1, y = X2, fill = fracunspl)) + 
    geom_point(color = "grey", shape = 21) + 
    scale_fill_gradientn(colors = viridis::viridis(21), na.value = "white", limits = c(0, 1),
                          name = paste0(showgene, " fraction unspliced")) + 
    labs(x = "", y = "", title = sm) + 
    theme_void() + theme(legend.position = "bottom")
})

png(gsub("_velocity.png", "_fracunspliced.png", outvelopng), width = 8, height = 12,
    unit = "in", res = 200)
print(cowplot::plot_grid(
  cowplot::plot_grid(plotlist = lapply(plotlist, function(w) w + theme(legend.position = "none")), ncol = 3),
  cowplot::get_legend(plotlist[[1]] + theme(legend.position = "bottom")),
  ncol = 1, rel_heights = c(1, 0.05)
))
dev.off()

## ------------------------------------------------------------------------- ##
## Velocities
## ------------------------------------------------------------------------- ##
## Get range of values for color scale
rng <- range(unlist(lapply(methods, function(m) {
  velocities[[m]][match(showgene, rownames(velocities[[m]])), ]
})), na.rm = TRUE)

plotlist <- lapply(methods, function(m) {
  sm <- methods_short$method_short[match(m, methods_short$method)]
  df <- umap
  df$velocity <- velocities[[m]][match(showgene, rownames(velocities[[m]])), ]
  ggplot(df, aes(x = X1, y = X2, fill = velocity)) + 
    geom_point(color = "grey", shape = 21) + 
    scale_fill_gradientn(colors = viridis::viridis(21), na.value = "white", limits = rng,
                         name = paste0(showgene, " velocity")) + 
    labs(x = "", y = "", title = sm) + 
    theme_void() + theme(legend.position = "bottom")
})

png(outvelopng, width = 8, height = 12,
    unit = "in", res = 200)
print(cowplot::plot_grid(
  cowplot::plot_grid(plotlist = lapply(plotlist, function(w) w + theme(legend.position = "none")), ncol = 3),
  cowplot::get_legend(plotlist[[1]] + theme(legend.position = "bottom")),
  ncol = 1, rel_heights = c(1, 0.05)
))
dev.off()

date()
sessionInfo()


