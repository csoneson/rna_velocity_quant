args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(cowplot)
  library(magick)
  library(ggplot2)
})
source(plothelperscript)

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(topdir)
print(plothelperscript)
print(dimred)
print(methods)
print(outrds)

methods_short <- shorten_methods(methods)

## Individually chosen genes
umaps <- lapply(methods, function(m) {
  file.path(topdir, "plots/velocity",
            paste0("anndata_", m), 
            paste0("anndata_", m, "_scvelo_", dimred,
                   "_alevin_spliced_stream.png"))
})

pdf(gsub("\\.rds$", paste0("_", dimred, ".pdf"), outrds), width = 8, height = 9)
cowplot::plot_grid(plotlist = lapply(names(umaps), function(nm) {
  title <- cowplot::ggdraw() + cowplot::draw_label(methods_short$method_short[match(nm, methods_short$method)], size = 7)
  w <- umaps[[nm]]
  ump <- cowplot::ggdraw() + cowplot::draw_image(w)
  cowplot::plot_grid(title, ump + theme(panel.background = element_rect(color = "black")), 
                     ncol = 1, rel_heights = c(0.1, 1))
}), ncol = 3)
dev.off()

## Shared genes
umapss <- lapply(methods, function(m) {
  file.path(topdir, "plots/velocity",
            paste0("anndata_", m, "_shared_genes"), 
            paste0("anndata_", m, "_shared_genes_scvelo_", dimred,
                   "_alevin_spliced_stream.png"))
})

pdf(gsub("\\.rds$", paste0("_", dimred, "_shared_genes.pdf"), outrds), width = 8, height = 9)
cowplot::plot_grid(plotlist = lapply(names(umapss), function(nm) {
  title <- cowplot::ggdraw() + cowplot::draw_label(methods_short$method_short[match(nm, methods_short$method)], size = 7)
  w <- umapss[[nm]]
  ump <- cowplot::ggdraw() + cowplot::draw_image(w)
  cowplot::plot_grid(title, ump + theme(panel.background = element_rect(color = "black")), 
                     ncol = 1, rel_heights = c(0.1, 1))
}), ncol = 3)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
