args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(cowplot)
  library(magick)
})

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(topdir)
print(methods)
print(outrds)

## Individually chosen genes
umaps <- lapply(methods, function(m) {
  file.path(topdir, "plots/velocity",
            paste0("anndata_", m), 
            paste0("anndata_", m, "_scvelo_UMAP_alevin_spliced_stream.png"))
})

pdf(gsub("\\.rds$", ".pdf", outrds), width = 8, height = 9)
cowplot::plot_grid(plotlist = lapply(umaps, function(w) {
  cowplot::ggdraw() + cowplot::draw_image(w)
}), ncol = 3)
dev.off()

## Shared genes
umapss <- lapply(methods, function(m) {
  file.path(topdir, "plots/velocity",
            paste0("anndata_", m, "_shared_genes"), 
            paste0("anndata_", m, "_shared_genes_scvelo_UMAP_alevin_spliced_stream.png"))
})

pdf(gsub("\\.rds$", "_shared_genes.pdf", outrds), width = 8, height = 9)
cowplot::plot_grid(plotlist = lapply(umapss, function(w) {
  cowplot::ggdraw() + cowplot::draw_image(w)
}), ncol = 3)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
