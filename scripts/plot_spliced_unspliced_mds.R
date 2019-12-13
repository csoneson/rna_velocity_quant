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
## Overall similarity (RMSE of count matrices) 
## ------------------------------------------------------------------------- ##
rmse <- do.call(dplyr::bind_rows, lapply(names(sces), function(s1) {
  do.call(dplyr::bind_rows, lapply(names(sces), function(s2) {
    data.frame(m1 = s1, m2 = s2, 
               RMSEspliced = sqrt(mean((assay(sces[[s1]], "spliced") - 
                                          assay(sces[[s2]], "spliced"))^2)),
               RMSEunspliced = sqrt(mean((assay(sces[[s1]], "unspliced") - 
                                            assay(sces[[s2]], "unspliced"))^2)),
               stringsAsFactors = FALSE)
  }))
}))
rmses <- rmse %>% dplyr::select(m1, m2, RMSEspliced) %>% 
  tidyr::spread(key = m2, value = RMSEspliced) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("m1") %>%
  as.matrix()
rmseu <- rmse %>% dplyr::select(m1, m2, RMSEunspliced) %>% 
  tidyr::spread(key = m2, value = RMSEunspliced) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("m1") %>%
  as.matrix()

## MDS
cmdspliced <- as.data.frame(cmdscale(rmses, k = 2)) %>%
  tibble::rownames_to_column("method") %>%
  dplyr::left_join(methods_short, by = "method")

cmdunspliced <- as.data.frame(cmdscale(rmseu, k = 2)) %>%
  tibble::rownames_to_column("method") %>%
  dplyr::left_join(methods_short, by = "method")

plset <- list(
  aes(x = V1, y = V2, shape = rtype, label = method_short),
  geom_point(aes(color = mtype), size = 7, alpha = 0.7),
  geom_label_repel(size = 3),
  theme_minimal(),
  scale_color_manual(values = base_method_colors, name = ""),
  scale_shape_discrete(name = ""),
  labs(x = "MDS1", y = "MDS2")
)
pdf(gsub("\\.rds$", "_mds_spliced_unspliced.pdf", outrds), width = 10, height = 6)
cowplot::plot_grid(
  cowplot::plot_grid(
    ggplot(cmdspliced) + 
      plset + ggtitle("MDS, RMSE (spliced)") + 
      theme(legend.position = "none"),
    ggplot(cmdunspliced) + 
      plset + ggtitle("MDS, RMSE (unspliced)") + 
      theme(legend.position = "none"),
    nrow = 1, labels = c("A", "B"), rel_widths = c(1, 1)
  ),
  cowplot::get_legend(ggplot(cmdspliced) + 
                        plset + theme(legend.position = "bottom")),
  ncol = 1, labels = "", rel_heights = c(1, 0.1)
)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
