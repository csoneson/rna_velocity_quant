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
  library(ggrepel)
})
source(plothelperscript)

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(plothelperscript)
print(topdir)
print(methods)
print(outrds)

sces <- lapply(methods, function(nm) {
  readRDS(file.path(topdir, paste0("output/sce/sce_", nm, ".rds")))
})

methods_short <- shorten_methods(methods)

## ------------------------------------------------------------------------- ##
## Overall similarity
## ------------------------------------------------------------------------- ##
## Spearman correlations of row/column sums
corrs <- do.call(dplyr::bind_rows, lapply(names(sces), function(s1) {
  do.call(dplyr::bind_rows, lapply(names(sces), function(s2) {
    print(c(s1, s2))
    stopifnot(all(rownames(sces[[s1]]) == rownames(sces[[s2]])))
    stopifnot(all(colnames(sces[[s1]]) == colnames(sces[[s2]])))
    data.frame(m1 = s1, m2 = s2, 
               corrCellSumsSpliced = cor(colSums(assay(sces[[s1]], "spliced")),
                                         colSums(assay(sces[[s2]], "spliced")), 
                                         method = "spearman"),
               corrCellSumsUnspliced = cor(colSums(assay(sces[[s1]], "unspliced")),
                                           colSums(assay(sces[[s2]], "unspliced")), 
                                           method = "spearman"),
               corrGeneSumsSpliced = cor(rowSums(assay(sces[[s1]], "spliced")),
                                         rowSums(assay(sces[[s2]], "spliced")), 
                                         method = "spearman"),
               corrGeneSumsUnspliced = cor(rowSums(assay(sces[[s1]], "unspliced")),
                                           rowSums(assay(sces[[s2]], "unspliced")), 
                                           method = "spearman"),
               stringsAsFactors = FALSE)
  }))
}))

## Create separate correlation matrices for spliced and unspliced counts
corrcells <- corrs %>% dplyr::select(m1, m2, corrCellSumsSpliced) %>% 
  tidyr::spread(key = m2, value = corrCellSumsSpliced) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("m1") %>%
  as.matrix()
corrcellu <- corrs %>% dplyr::select(m1, m2, corrCellSumsUnspliced) %>% 
  tidyr::spread(key = m2, value = corrCellSumsUnspliced) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("m1") %>%
  as.matrix()
corrgenes <- corrs %>% dplyr::select(m1, m2, corrGeneSumsSpliced) %>% 
  tidyr::spread(key = m2, value = corrGeneSumsSpliced) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("m1") %>%
  as.matrix()
corrgeneu <- corrs %>% dplyr::select(m1, m2, corrGeneSumsUnspliced) %>% 
  tidyr::spread(key = m2, value = corrGeneSumsUnspliced) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("m1") %>%
  as.matrix()

stopifnot(all(rownames(corrcells) == colnames(corrcells)))
stopifnot(all(rownames(corrcellu) == colnames(corrcellu)))
stopifnot(all(rownames(corrgenes) == colnames(corrgenes)))
stopifnot(all(rownames(corrgeneu) == colnames(corrgeneu)))

## MDS
cmdcellspliced <- as.data.frame(cmdscale(d = sqrt(2 - 2 * corrcells), k = 2)) %>%
  tibble::rownames_to_column("method") %>%
  dplyr::left_join(methods_short, by = "method")

cmdcellunspliced <- as.data.frame(cmdscale(d = sqrt(2 - 2 * corrcellu), k = 2)) %>%
  tibble::rownames_to_column("method") %>%
  dplyr::left_join(methods_short, by = "method")

cmdgenespliced <- as.data.frame(cmdscale(d = sqrt(2 - 2 * corrgenes), k = 2)) %>%
  tibble::rownames_to_column("method") %>%
  dplyr::left_join(methods_short, by = "method")

cmdgeneunspliced <- as.data.frame(cmdscale(d = sqrt(2 - 2 * corrgeneu), k = 2)) %>%
  tibble::rownames_to_column("method") %>%
  dplyr::left_join(methods_short, by = "method")

plset <- list(
  aes(x = V1, y = V2, shape = rtype, label = method_short),
  geom_point(aes(color = mtype), size = 7, alpha = 0.7),
  geom_label_repel(size = 1.5),
  theme_minimal(),
  scale_color_manual(values = base_method_colors, name = ""),
  scale_shape_discrete(name = ""),
  labs(x = "MDS1", y = "MDS2")
)
pdf(gsub("\\.rds$", "_mds_spliced_unspliced.pdf", outrds), width = 10, height = 6)
cowplot::plot_grid(
  cowplot::plot_grid(
    ggplot(cmdcellspliced) + 
      plset + ggtitle("MDS, Cells (spliced)") + 
      theme(legend.position = "none"),
    ggplot(cmdcellunspliced) + 
      plset + ggtitle("MDS, Cells (unspliced)") + 
      theme(legend.position = "none"),
    nrow = 1, labels = c("A", "B"), rel_widths = c(1, 1)
  ),
  cowplot::plot_grid(
    ggplot(cmdgenespliced) + 
      plset + ggtitle("MDS, Genes (spliced)") + 
      theme(legend.position = "none"),
    ggplot(cmdgeneunspliced) + 
      plset + ggtitle("MDS, Genes (unspliced)") + 
      theme(legend.position = "none"),
    nrow = 1, labels = c("C", "D"), rel_widths = c(1, 1)
  ),
  cowplot::get_legend(ggplot(cmdcellspliced) + 
                        plset + theme(legend.position = "bottom")),
  ncol = 1, labels = "", rel_heights = c(1, 1, 0.1)
)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
