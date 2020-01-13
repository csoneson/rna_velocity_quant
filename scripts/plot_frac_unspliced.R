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
  library(cdata)
})
source(plothelperscript)

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(plothelperscript)
print(topdir)
print(refdir)  ## directory where uniqueness files are
print(tx2gene)
print(methods)
print(outrds)

## ------------------------------------------------------------------------- ##
## Read data
## ------------------------------------------------------------------------- ##
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
             frac_unspliced = rowSums(assay(w, "unspliced"))/
               (rowSums(assay(w, "unspliced")) + rowSums(assay(w, "spliced"))),
             stringsAsFactors = FALSE
  )
}))

tx2gene <- readRDS(tx2gene)
uniq <- merge_uniq(refdir = refdir, tx2gene = tx2gene, 
                   keepgenes = sumdf_bygene$gene)

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
  dplyr::summarize(sd_frac_unspliced = sd(frac_unspliced),
                   mean_frac_unspliced = mean(frac_unspliced),
                   min_frac_unspliced = min(frac_unspliced),
                   max_frac_unspliced = max(frac_unspliced))

## Plot ----
sumdf_bygene_acrossmethods_scvelosel <- 
  sumdf_bygene_acrossmethods %>% dplyr::filter(scvelo_selected)
qtl <- quantile(sumdf_bygene_acrossmethods_scvelosel$sd_frac_unspliced, 0.9, na.rm = TRUE)
pdf(gsub("\\.rds$", "_fracunspliced_sd_distr.pdf", outrds))
ggplot(sumdf_bygene_acrossmethods_scvelosel, 
       aes(x = sd_frac_unspliced)) + 
  geom_histogram(bins = 100, fill = "lightgrey") + 
  geom_vline(xintercept = qtl) + 
  theme_bw() + 
  labs(x = "Standard deviation of fraction UMIs in unspliced targets",
       y = "Number of genes",
       title = "Variability of the unspliced fractions across quantifications, by gene",
       subtitle = paste0("Only genes selected by scVelo for at least one quantification (n=", 
                         nrow(sumdf_bygene_acrossmethods_scvelosel), ")")) 
dev.off()

## Retain top 10% ---- 
genes_to_keep <- sumdf_bygene_acrossmethods_scvelosel %>%
  dplyr::filter(sd_frac_unspliced > qtl) %>%
  dplyr::pull(gene)
print(length(genes_to_keep))

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
    width = 10, height = 25)
print(pheatmap::pheatmap(
  clstdata, cluster_rows = hcl, cutree_rows = 10, 
  scale = "none", fontsize_row = 4,
  annotation_row = data.frame(clusters = factor(clusters), row.names = names(clusters))))
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

## Summarize each cluster by its centroid
ms <- split(clstdata, f = clusters)
m <- do.call(dplyr::bind_rows, lapply(ms, colMeans))

## Find most similar gene in each cluster
cl_rep <- do.call(dplyr::bind_rows, lapply(unique(clusters), function(w) {
  tmp <- cor(t(ms[[w]]), t(m[w, , drop = FALSE]))
  data.frame(cluster = w, 
             gene = rownames(tmp),
             cor = tmp[, 1],
             stringsAsFactors = FALSE)
}))
print(do.call(dplyr::bind_rows, 
              lapply(split(cl_rep, cl_rep$cluster), 
                     function(w) w %>% arrange(desc(cor)) %>% head(3))))
write.table(cl_rep %>% dplyr::arrange(cluster, desc(cor)), 
            file = gsub("\\.rds$", "_fracunspliced_cluster_centroid_corrs.txt", outrds),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
      
## Plot centroids and correlations within each cluster
rn <- round(1e7 * runif(1))
tmpdir <- tempdir()
png(paste0(tmpdir, "/pheatmap", rn, ".png"), width = 6,
    height = 5, unit = "in", res = 400)
print(pheatmap::pheatmap(m, cluster_rows = FALSE, cluster_cols = TRUE, treeheight_col = 10, 
                         main = "Centroid cluster profiles"))
dev.off()
g <- ggplot(cl_rep, aes(x = factor(cluster), y = cor)) + geom_boxplot(outlier.size = -1) + 
  geom_jitter(width = 0.2, height = 0) + theme_bw() + 
  labs(x = "Cluster", y = "Correlations with centroid profile")
pdf(gsub("\\.rds$", "_fracunspliced_cluster_centroids.pdf", outrds), width = 10, height = 5)
cowplot::plot_grid(
  cowplot::ggdraw() +
    cowplot::draw_image(paste0(tmpdir, "/pheatmap", rn, ".png")),
  g,
  nrow = 1, rel_widths = c(6, 4), labels = c("A", "B")
)
dev.off()

## Uniqueness for each cluster
uniqsub <- uniq %>% 
  dplyr::left_join(
    as.data.frame(clusters) %>% tibble::rownames_to_column("gene")
  ) %>%
  dplyr::filter(!is.na(clusters)) %>%
  tidyr::unite(col = "actype", ctype, atype, sep = ".") %>%
  dplyr::mutate(clusters = factor(clusters, levels = sort(unique(as.numeric(clusters)))))
  
pdf(gsub("\\.rds$", "_uniqueness_by_cluster.pdf", outrds), width = 7, height = 5)
ggplot(uniqsub, aes(x = clusters, 
                    y = frac_unique, 
                    fill = clusters)) + 
  geom_boxplot(alpha = 0.5) + 
  facet_wrap(~ actype) + 
  theme_bw() + 
  labs(x = "Cluster", y = "Uniqueness") + 
  theme(legend.position = "none")
dev.off()

saveRDS(NULL, file = outrds)

date()
sessionInfo()