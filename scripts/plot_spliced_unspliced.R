args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(tximeta)
  library(dplyr)
  library(ggplot2)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(scater)
  library(cowplot)
  library(pheatmap)
  library(ggrepel)
})

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(topdir)
print(gtf)
print(tx2gene)
print(methods)
print(outrds)

sces <- lapply(methods, function(nm) {
  readRDS(file.path(topdir, paste0("output/sce/sce_", nm, ".rds")))
})
for (nm in names(sces)) {
  metadata(sces[[nm]]) <- list(method = nm)
}

gtf <- rtracklayer::import(gtf)
tx2gene <- readRDS(tx2gene)

methods_short <- data.frame(method = methods, 
                            stringsAsFactors = FALSE) %>%
  dplyr::mutate(method_short = 
                  gsub("ude", "", 
                       gsub("kallisto_bustools", "kallisto|bus", 
                            gsub("collapse", "coll", 
                                 gsub("separate", "sep",
                                      gsub("iso", "", 
                                           gsub("prepref_", "",
                                                gsub("_cdna_introns", "", method)))))))) %>%
  dplyr::mutate(
    mtype = stringr::str_extract(
      method_short, "alevin|kallisto\\|bus|starsolo|velocyto"
    ),
    rtype = stringr::str_extract(
      method, "separate|collapse"
    )
  ) %>%
  dplyr::mutate(rtype = replace(rtype, is.na(rtype), "N/A")) %>%
  dplyr::mutate(rtype = factor(rtype, levels = c("collapse", "separate", "N/A")))


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
  scale_color_manual(values = c(alevin = "#999999", `kallisto|bus` = "#009E73",
                                starsolo = "#0072B2", velocyto = "#CC79A7"), name = ""),
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

## ------------------------------------------------------------------------- ##
## Total number of assigned reads
## ------------------------------------------------------------------------- ##
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

## Add info about uniqueness
uniq_separate <- read.delim(
  file.path(topdir, "reference/prepref_isoseparate_uniqueness.txt"),
  header = TRUE, as.is = TRUE)
uniq_collapse <- read.delim(
  file.path(topdir, "reference/prepref_isocollapse_uniqueness.txt"),
  header = TRUE, as.is = TRUE)

uniq <- dplyr::bind_rows(
  uniq_separate %>% 
    dplyr::mutate(ctype = c("exonic", "intronic")[grepl("I\\.", gene) + 1]) %>%
    dplyr::mutate(gene = gsub("I\\.", "", gene)) %>%
    dplyr::mutate(frac_unique = unique/total) %>%
    dplyr::select(gene, ctype, frac_unique) %>%
    dplyr::mutate(atype = "separate"),
  uniq_collapse %>% 
    dplyr::mutate(ctype = c("exonic", "intronic")[grepl("I\\.", gene) + 1]) %>%
    dplyr::mutate(gene = gsub("I\\.", "", gene)) %>%
    dplyr::mutate(frac_unique = unique/total) %>%
    dplyr::select(gene, ctype, frac_unique) %>%
    dplyr::mutate(atype = "collapse")
) %>%
  dplyr::mutate(frac_unique_bin = Hmisc::cut2(frac_unique, 
                                              cuts = c(0, 0.001, 0.5, 0.999, 1))) %>%
  dplyr::mutate(gene = tx2gene$gene_name[match(gene, tx2gene$gene_id)]) %>%
  dplyr::filter(gene %in% sumdf_bygene$gene)

pdf(gsub("\\.rds$", "_uniqueness.pdf", outrds), width = 8, height = 5)
ggplot(uniq %>% dplyr::select(-frac_unique_bin) %>% dplyr::group_by(atype) %>% 
         tidyr::spread(key = "ctype", value = "frac_unique"),
       aes(x = exonic, y = intronic)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point(alpha = 0.25) + 
  facet_wrap(~ atype) + 
  theme_bw() + 
  labs(x = "Fraction unique k-mers, exonic features",
       y = "Fraction unique k-mers, intronic features",
       title = "Fraction unique k-mers")

ggplot(uniq %>% dplyr::select(-frac_unique_bin) %>% dplyr::group_by(ctype) %>% 
         tidyr::spread(key = "atype", value = "frac_unique"),
       aes(x = separate, y = collapse)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point(alpha = 0.25) + 
  facet_wrap(~ ctype) + 
  theme_bw() + 
  labs(x = "Fraction unique k-mers, separate",
       y = "Fraction unique k-mers, collapse",
       title = "Fraction unique k-mers")

ggplot(uniq %>% dplyr::select(-frac_unique_bin), aes(x = frac_unique)) + 
  geom_histogram(bins = 100) + 
  facet_grid(atype ~ ctype) + 
  theme_bw() + 
  labs(x = "Fraction unique k-mers",
       y = "count",
       title = "Fraction unique k-mers")
dev.off()

pdf(gsub("\\.rds$", "_total_count.pdf", outrds), width = 7, height = 10)
for (ct in c("exonic", "intronic")) {
  for (at in c("collapse", "separate")) {
    print(ggplot(dplyr::bind_rows(
      sumdf_bygene %>% 
        dplyr::left_join(uniq %>% dplyr::filter(ctype == ct & 
                                                  atype == at)) %>% 
        dplyr::mutate(frac_unique_bin = as.character(frac_unique_bin)) %>% 
        dplyr::group_by(method, frac_unique_bin) %>% 
        dplyr::summarize(spliced = sum(spliced), 
                         unspliced = sum(unspliced), 
                         total = sum(total)) %>% 
        tidyr::gather(key = "ctype", value = "count", spliced, unspliced, total) %>%
        dplyr::mutate(ctype = factor(ctype, 
                                     levels = c("total", "spliced", "unspliced"))) %>%
        dplyr::left_join(methods_short, by = "method"), 
      sumdf_bygene %>% group_by(method) %>% 
        dplyr::summarize(spliced = sum(spliced), 
                         unspliced = sum(unspliced), 
                         total = sum(total)) %>% 
        tidyr::gather(key = "ctype", value = "count", spliced, unspliced, total) %>%
        dplyr::mutate(frac_unique_bin = "overall") %>%
        dplyr::mutate(ctype = factor(ctype, 
                                     levels = c("total", "spliced", "unspliced"))) %>%
        dplyr::left_join(methods_short, by = "method")
    ) %>%
      dplyr::mutate(frac_unique_bin = relevel(factor(frac_unique_bin), ref = "overall")),
    aes(x = method_short, y = count, fill = mtype)) + 
      geom_bar(stat = "identity") + 
      facet_grid(frac_unique_bin ~ ctype, scale = "free_y") + 
      theme_bw() + 
      scale_fill_manual(values = c(alevin = "#999999", `kallisto|bus` = "#009E73",
                                   starsolo = "#0072B2", velocyto = "#CC79A7"), name = "") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "none") + 
      ggtitle(paste0("Total count, stratified by ", ct, " uniqueness (", at, ")")))
  }
}
dev.off()

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
  dplyr::left_join(
    uniq %>% dplyr::select(gene, ctype, atype, frac_unique) %>%
      tidyr::unite(col = "catype", ctype, atype, sep = ".") %>%
      dplyr::mutate(catype = paste0("fracunique.", catype)) %>% 
      tidyr::spread(key = "catype", value = "frac_unique") 
  ) %>% 
  as.data.frame() %>%
  tibble::column_to_rownames("gene")

hcl <- hclust(d = as.dist(sqrt(2 - 2*cor(t(clstdata)))))
clusters <- cutree(hcl, k = 10)

pdf(gsub("\\.rds$", "_fracunspliced_clustering.pdf", outrds), 
    width = 10, height = 35)
print(pheatmap::pheatmap(
  clstdata, cluster_rows = hcl, cutree_rows = 10, 
  scale = "none", fontsize_row = 4, 
  annotation_row = clstannot[match(rownames(clstdata), rownames(clstannot)), 
                             grep("fracunique", colnames(clstannot))]))
dev.off()

for (i in unique(clusters)) {
  gn <- names(clusters[clusters == i])
  if (length(gn) > 2) {
    pdf(gsub("\\.rds$", paste0("_fracunspliced_clustering_cluster", i, ".pdf"), outrds),
        width = 8, height = 5 + length(gn) * 0.07)
    print(pheatmap::pheatmap(
      clstdata[gn, ], scale = "none",
      fontsize_row = 5,
      main = paste0("cluster ", i),
      annotation_row = clstannot[match(gn, rownames(clstannot)), 
                                 grep("fracunique", colnames(clstannot))]))
    dev.off()
  }
}


saveRDS(NULL, file = outrds)

date()
sessionInfo()