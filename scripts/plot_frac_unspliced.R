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

sumdf_bygene <- do.call(dplyr::bind_rows, lapply(sces, function(w) {
  data.frame(gene = rownames(w),
             method = metadata(w)$method, 
             spliced = rowSums(assay(w, "spliced")),
             unspliced = rowSums(assay(w, "unspliced")),
             total = rowSums(assay(w, "unspliced")) + rowSums(assay(w, "spliced")),
             frac_unspliced = rowSums(assay(w, "unspliced"))/(rowSums(assay(w, "unspliced")) + rowSums(assay(w, "spliced"))),
             stringsAsFactors = FALSE
  )
}))

## ------------------------------------------------------------------------- ##
## Velocity
## ------------------------------------------------------------------------- ##
## Read genes used for velocity calculations (2,000 per method)
geneinfo <- lapply(methods, function(nm) {
  readr::read_csv(paste0(topdir, "/plots/velocity/anndata_", nm, "/anndata_", nm, 
                         "_gene_info.csv")) %>%
    dplyr::mutate(method = nm)
})
velocitygenes_freq <- table(unlist(lapply(geneinfo, function(w) w$index)))
velocitygenes <- 
  as.character(Reduce(union, lapply(geneinfo, function(w) w$index)))
sumdf_bygene$scvelo_selected <- sumdf_bygene$gene %in% velocitygenes
sumdf_bygene$scvelo_selected_nbr <- as.numeric(velocitygenes_freq[sumdf_bygene$gene])
sumdf_bygene$scvelo_selected_nbr[is.na(sumdf_bygene$scvelo_selected_nbr)] <- 0

## Summarize across methods
sumdf_bygene_acrossmethods <- sumdf_bygene %>%
  dplyr::group_by(gene, scvelo_selected, scvelo_selected_nbr) %>%
  dplyr::summarize(sd_frac_unspliced = sd(frac_unspliced))

## Plot ----
sumdf_bygene_acrossmethods_scvelosel <- sumdf_bygene_acrossmethods %>% dplyr::filter(scvelo_selected)
qtl <- quantile(sumdf_bygene_acrossmethods_scvelosel$sd_frac_unspliced, 0.9, na.rm = TRUE)
ggplot(sumdf_bygene_acrossmethods_scvelosel, 
       aes(x = sd_frac_unspliced)) + 
  geom_histogram(bins = 100, fill = "lightgrey") + 
  geom_vline(xintercept = qtl) + 
  theme_bw() + 
  labs(x = "Standard deviation of fraction UMIs in unspliced targets",
       y = "Number of genes",
       title = "Variability of the unspliced fractions across quantifications, by gene",
       subtitle = "Only genes selected by scVelo for at least one quantification") 

## Retain top 10% ---- 
genes_to_keep <- sumdf_bygene_acrossmethods_scvelosel %>%
  dplyr::filter(sd_frac_unspliced > qtl) %>%
  dplyr::pull(gene)

## Cluster the selected genes based on their unspliced fraction pattern across methods ---- 
clstdata <- sumdf_bygene %>% dplyr::filter(gene %in% genes_to_keep) %>%
  dplyr::select(gene, method, frac_unspliced) %>%
  dplyr::left_join(methods_short) %>% 
  dplyr::select(gene, method_short, frac_unspliced) %>%
  tidyr::spread(key = "method_short", value = "frac_unspliced") %>%
  as.data.frame() %>%
  tibble::column_to_rownames("gene")
clstannot <- sumdf_bygene %>% dplyr::filter(gene %in% genes_to_keep) %>%
  dplyr::select(gene, method, total) %>%
  dplyr::left_join(methods_short) %>% 
  dplyr::select(gene, method_short, total) %>%
  tidyr::spread(key = "method_short", value = "total") %>%
  as.data.frame() %>%
  tibble::column_to_rownames("gene")

hcl <- hclust(d = as.dist(sqrt(2 - 2*cor(t(clstdata)))))
clusters <- cutree(hcl, k = 10)

pdf(gsub("\\.rds$", "_fracunspliced_clustering.pdf", outrds), 
    width = 10, height = 35)
print(pheatmap::pheatmap(
  clstdata, cluster_rows = hcl, cutree_rows = 10, 
  scale = "none", fontsize_row = 4))
dev.off()

for (i in unique(clusters)) {
  gn <- names(clusters[clusters == i])
  if (length(gn) > 2) {
    pdf(gsub("\\.rds$", paste0("_fracunspliced_clustering_cluster", i, ".pdf"), outrds),
        width = 8, height = 5 + length(gn) * 0.07)
    print(pheatmap::pheatmap(
      clstdata[gn, ], scale = "none",
      fontsize_row = 5,
      main = paste0("cluster ", i)))
    dev.off()
  }
}


saveRDS(NULL, file = outrds)

date()
sessionInfo()