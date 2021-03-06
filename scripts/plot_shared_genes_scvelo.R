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
  library(cowplot)
})

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods
dataset <- gsub("_", " ", dataset)

print(topdir)
print(velosuffix)
print(plothelperscript)
print(dataset)
print(methods)
print(outrds)

source(plothelperscript)
methods_short <- shorten_methods(methods) %>%
  dplyr::filter(method %in% methods) %>%
  dplyr::arrange(as.factor(method_short))
methods <- methods_short$method
names(methods) <- methods

geneinfo <- lapply(methods, function(nm) {
  readr::read_csv(file.path(topdir, paste0("plots/velocity", velosuffix, "/anndata_",
                                           nm, "/anndata_", nm, 
                                           "_gene_info.csv")))
})

allgenes <- as.character(Reduce(union, lapply(geneinfo, function(w) w$index)))
selgenes <- do.call(cbind, lapply(geneinfo, 
                                  function(w) as.numeric(allgenes %in% w$index)))
rownames(selgenes) <- allgenes
selgenes <- data.frame(selgenes)
colnames(selgenes) <- methods_short$method_short[match(colnames(selgenes), methods_short$method)]

pdf(gsub("rds$", "pdf", outrds), width = 9, height = 5)
g <- upset(selgenes, nsets = length(geneinfo), keep.order = TRUE, 
           order.by = "freq", decreasing = TRUE)
print(cowplot::plot_grid(ggdraw() + draw_label(dataset, size = 13),
                         g$Main_bar, g$Matrix, ncol = 1, align = "v", 
                         rel_heights = c(0.2, 2.5, 1)))
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()

