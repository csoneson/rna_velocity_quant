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
  library(Seurat)
})
source(plothelperscript)

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(topdir)
print(plothelperscript)
print(methods)
print(genetxt)
print(outrds)

methods_short <- shorten_methods(methods)

sces <- lapply(methods, function(nm) {
  readRDS(file.path(topdir, paste0("output/sce/sce_", nm, ".rds")))
})

genes <- read.delim(genetxt, header = FALSE, as.is = TRUE)[, 1]
sces <- lapply(sces, function(w) w[match(genes, rownames(w)), ])

seurats <- lapply(methods, function(nm) {
  w <- ReadH5AD(file.path(topdir, paste0("output/anndata_with_velocity/anndata_", 
                                         nm, "_shared_genes_with_velocity.h5ad")))
  stopifnot(all(rownames(w) == genes))
  w
})

gene_corrs <- do.call(
  dplyr::bind_rows, 
  lapply(seq_len(length(methods) - 1), function(j) {
    jj <- methods[j]
    do.call(
      dplyr::bind_rows,  
      lapply((j + 1):(length(methods)), function(k) {
        kk <- methods[k]
        message(jj, " - ", kk)
        
        a <- as.data.frame(t(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "spliced")))))
        b <- as.data.frame(t(as.matrix(GetAssayData(GetAssay(seurats[[kk]], "spliced")))))
        corrs_spliced <- mapply(a, b, FUN = cor, 
                                MoreArgs = list(use = "na.or.complete", method = "spearman"))
        
        a <- as.data.frame(t(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "unspliced")))))
        b <- as.data.frame(t(as.matrix(GetAssayData(GetAssay(seurats[[kk]], "unspliced")))))
        corrs_unspliced <- mapply(a, b, FUN = cor, 
                                  MoreArgs = list(use = "na.or.complete", method = "spearman"))
        
        a <- as.data.frame(t(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "velocity")))))
        b <- as.data.frame(t(as.matrix(GetAssayData(GetAssay(seurats[[kk]], "velocity")))))
        corrs_velocity <- mapply(a, b, FUN = cor, 
                                 MoreArgs = list(use = "na.or.complete", method = "spearman"))
        
        a <- as.data.frame(t(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "spliced")))) + 
                             t(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "unspliced")))))
        b <- as.data.frame(t(as.matrix(GetAssayData(GetAssay(seurats[[kk]], "spliced")))) + 
                             t(as.matrix(GetAssayData(GetAssay(seurats[[kk]], "unspliced")))))
        corrs_total <- mapply(a, b, FUN = cor, 
                              MoreArgs = list(use = "na.or.complete", method = "spearman"))
        
        dplyr::bind_rows(
          data.frame(method1 = jj,
                     method2 = kk, 
                     gene = names(corrs_spliced), 
                     ctype = "spliced",
                     dtype = "gene",
                     corrs = corrs_spliced, 
                     stringsAsFactors = FALSE),
          data.frame(method1 = jj,
                     method2 = kk,
                     gene = names(corrs_unspliced), 
                     ctype = "unspliced",
                     dtype = "gene",
                     corrs = corrs_unspliced, 
                     stringsAsFactors = FALSE),
          data.frame(method1 = jj,
                     method2 = kk,
                     gene = names(corrs_velocity), 
                     ctype = "velocity",
                     dtype = "gene",
                     corrs = corrs_velocity, 
                     stringsAsFactors = FALSE),
          data.frame(method1 = jj,
                     method2 = kk,
                     gene = names(corrs_total), 
                     ctype = "total",
                     dtype = "gene",
                     corrs = corrs_total, 
                     stringsAsFactors = FALSE)
        )
      }))
  })) %>%
  dplyr::mutate(method1 = methods_short$method_short[match(method1, methods_short$method)],
                method2 = methods_short$method_short[match(method2, methods_short$method)])
gene_corrs$method1 <- factor(gene_corrs$method1, levels = methods_short$method_short[match(methods, methods_short$method)])
gene_corrs$method2 <- factor(gene_corrs$method2, levels = methods_short$method_short[match(methods, methods_short$method)])

gene_corrs_within <- do.call(
  dplyr::bind_rows, 
  lapply(seq_len(length(methods)), function(j) {
    jj <- methods[j]
    
    a <- as.data.frame(t(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "spliced")))))
    b <- as.data.frame(t(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "unspliced")))))
    corrs_spliced_unspliced <- mapply(a, b, FUN = cor, 
                                      MoreArgs = list(use = "na.or.complete", method = "spearman"))
    
    a <- as.data.frame(t(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "spliced")))) + 
                         t(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "unspliced")))))
    b <- as.data.frame(t(abs(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "velocity"))))))
    corrs_total_absvelocity <- mapply(a, b, FUN = cor, 
                                      MoreArgs = list(use = "na.or.complete", method = "spearman"))
    
    a <- as.data.frame(t(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "unspliced"))))/
                         (t(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "spliced")))) + 
                            t(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "unspliced"))))))
    b <- as.data.frame(t(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "velocity")))))
    corrs_fracunspl_velocity <- mapply(a, b, FUN = cor, 
                                       MoreArgs = list(use = "na.or.complete", method = "spearman"))
    
    dplyr::bind_rows(
      data.frame(method = jj,
                 gene = names(corrs_spliced_unspliced), 
                 ctype = "spliced-unspliced",
                 dtype = "gene",
                 corrs = corrs_spliced_unspliced, 
                 stringsAsFactors = FALSE),
      data.frame(method = jj,
                 gene = names(corrs_total_absvelocity), 
                 ctype = "total-absvelocity",
                 dtype = "gene",
                 corrs = corrs_total_absvelocity, 
                 stringsAsFactors = FALSE),
      data.frame(method = jj,
                 gene = names(corrs_fracunspl_velocity), 
                 ctype = "fracunspliced-velocity",
                 dtype = "gene",
                 corrs = corrs_fracunspl_velocity, 
                 stringsAsFactors = FALSE)
    )
  })) %>%
  dplyr::mutate(method = methods_short$method_short[match(method, methods_short$method)])
gene_corrs_within$method <- factor(gene_corrs_within$method, levels = methods_short$method_short[match(methods, methods_short$method)])

## Same but for cells
cell_corrs <- do.call(
  dplyr::bind_rows, 
  lapply(seq_len(length(methods) - 1), function(j) {
    jj <- methods[j]
    do.call(
      dplyr::bind_rows,  
      lapply((j + 1):(length(methods)), function(k) {
        kk <- methods[k]
        message(jj, " - ", kk)
        
        a <- as.data.frame(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "spliced"))))
        b <- as.data.frame(as.matrix(GetAssayData(GetAssay(seurats[[kk]], "spliced"))))
        corrs_spliced <- mapply(a, b, FUN = cor, 
                                MoreArgs = list(use = "na.or.complete", method = "spearman"))
        
        a <- as.data.frame(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "unspliced"))))
        b <- as.data.frame(as.matrix(GetAssayData(GetAssay(seurats[[kk]], "unspliced"))))
        corrs_unspliced <- mapply(a, b, FUN = cor, 
                                  MoreArgs = list(use = "na.or.complete", method = "spearman"))
        
        a <- as.data.frame(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "velocity"))))
        b <- as.data.frame(as.matrix(GetAssayData(GetAssay(seurats[[kk]], "velocity"))))
        corrs_velocity <- mapply(a, b, FUN = cor, 
                                 MoreArgs = list(use = "na.or.complete", method = "spearman"))
        
        a <- as.data.frame(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "spliced"))) + 
                             as.matrix(GetAssayData(GetAssay(seurats[[jj]], "unspliced"))))
        b <- as.data.frame(as.matrix(GetAssayData(GetAssay(seurats[[kk]], "spliced"))) + 
                             as.matrix(GetAssayData(GetAssay(seurats[[kk]], "unspliced"))))
        corrs_total <- mapply(a, b, FUN = cor, 
                              MoreArgs = list(use = "na.or.complete", method = "spearman"))
        
        if ("clusters" %in% colnames(colData(sces[[jj]]))) {
          clsce <- sces[[jj]]$clusters
          clseu <- levels(factor(sces[[jj]]$clusters))[seurats[[jj]]@meta.data$clusters + 1]
        } else {
          clsce <- rep(NA, ncol(sces[[jj]]))
          clseu <- rep(NA, ncol(sces[[jj]]))
        }
        dplyr::bind_rows(
          data.frame(method1 = jj,
                     method2 = kk, 
                     cell = names(corrs_spliced), 
                     cluster = clseu,
                     ctype = "spliced",
                     dtype = "cell",
                     corrs = corrs_spliced, 
                     stringsAsFactors = FALSE),
          data.frame(method1 = jj,
                     method2 = kk,
                     cell = names(corrs_unspliced), 
                     cluster = clseu,
                     ctype = "unspliced",
                     dtype = "cell",
                     corrs = corrs_unspliced, 
                     stringsAsFactors = FALSE),
          data.frame(method1 = jj,
                     method2 = kk,
                     cell = names(corrs_velocity), 
                     cluster = clseu,
                     ctype = "velocity",
                     dtype = "cell",
                     corrs = corrs_velocity, 
                     stringsAsFactors = FALSE),
          data.frame(method1 = jj,
                     method2 = kk,
                     cell = names(corrs_total), 
                     cluster = clseu,
                     ctype = "total",
                     dtype = "cell",
                     corrs = corrs_total, 
                     stringsAsFactors = FALSE)
        )
      }))
  })) %>%
  dplyr::mutate(method1 = methods_short$method_short[match(method1, methods_short$method)],
                method2 = methods_short$method_short[match(method2, methods_short$method)])
cell_corrs$method1 <- factor(cell_corrs$method1, levels = methods_short$method_short[match(methods, methods_short$method)])
cell_corrs$method2 <- factor(cell_corrs$method2, levels = methods_short$method_short[match(methods, methods_short$method)])

cell_corrs_within <- do.call(
  dplyr::bind_rows, 
  lapply(seq_len(length(methods)), function(j) {
    jj <- methods[j]
    
    a <- as.data.frame(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "spliced"))))
    b <- as.data.frame(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "unspliced"))))
    corrs_spliced_unspliced <- mapply(a, b, FUN = cor, 
                                      MoreArgs = list(use = "na.or.complete", method = "spearman"))
    
    a <- as.data.frame(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "spliced"))) + 
                         as.matrix(GetAssayData(GetAssay(seurats[[jj]], "unspliced"))))
    b <- as.data.frame(as.matrix(abs(GetAssayData(GetAssay(seurats[[jj]], "velocity")))))
    corrs_total_absvelocity <- mapply(a, b, FUN = cor, 
                                      MoreArgs = list(use = "na.or.complete", method = "spearman"))
    
    a <- as.data.frame(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "unspliced")))/
                         (as.matrix(GetAssayData(GetAssay(seurats[[jj]], "spliced"))) + 
                            as.matrix(GetAssayData(GetAssay(seurats[[jj]], "unspliced")))))
    b <- as.data.frame(as.matrix(abs(GetAssayData(GetAssay(seurats[[jj]], "velocity")))))
    corrs_fracunspl_velocity <- mapply(a, b, FUN = cor, 
                                       MoreArgs = list(use = "na.or.complete", method = "spearman"))
    
    if ("clusters" %in% colnames(colData(sces[[jj]]))) {
      clsce <- sces[[jj]]$clusters
      clseu <- levels(factor(sces[[jj]]$clusters))[seurats[[jj]]@meta.data$clusters + 1]
    } else {
      clsce <- rep(NA, ncol(sces[[jj]]))
      clseu <- rep(NA, ncol(sces[[jj]]))
    }
    dplyr::bind_rows(
      data.frame(method = jj,
                 cell = names(corrs_spliced_unspliced), 
                 cluster = clseu,
                 ctype = "spliced-unspliced",
                 dtype = "cell",
                 corrs = corrs_spliced_unspliced, 
                 stringsAsFactors = FALSE),
      data.frame(method = jj,
                 cell = names(corrs_total_absvelocity), 
                 cluster = clseu,
                 ctype = "total-absvelocity",
                 dtype = "cell",
                 corrs = corrs_total_absvelocity, 
                 stringsAsFactors = FALSE),
      data.frame(method = jj,
                 cell = names(corrs_fracunspl_velocity), 
                 cluster = clseu,
                 ctype = "fracunspliced-velocity",
                 dtype = "cell",
                 corrs = corrs_fracunspl_velocity, 
                 stringsAsFactors = FALSE)
    )
  })) %>%
  dplyr::mutate(method = methods_short$method_short[match(method, methods_short$method)])
cell_corrs_within$method <- factor(cell_corrs_within$method, levels = methods_short$method_short[match(methods, methods_short$method)])

## Combine gene and cell correlations
gene_cell_corrs <- dplyr::bind_rows(
  gene_corrs %>% dplyr::rename(m1 = method1) %>%
    dplyr::rename(method1 = method2,
                  method2 = m1) %>%
    tidyr::unite(dtype, ctype, col = "cdtype", sep = ": ", remove = FALSE) %>%
    dplyr::select(method1, method2, ctype, dtype, cdtype, corrs),
  cell_corrs %>% 
    tidyr::unite(dtype, ctype, col = "cdtype", sep = ": ", remove = FALSE) %>%
    dplyr::select(method1, method2, ctype, dtype, cdtype, corrs)
)
gene_cell_corrs$method1 <- factor(gene_cell_corrs$method1, levels = methods_short$method_short[match(methods, methods_short$method)])
gene_cell_corrs$method2 <- factor(gene_cell_corrs$method2, levels = methods_short$method_short[match(methods, methods_short$method)])

gene_cell_corrs_within <- dplyr::bind_rows(
  gene_corrs_within %>% 
    tidyr::unite(dtype, ctype, col = "cdtype", sep = ": ", remove = FALSE) %>%
    dplyr::select(method, ctype, dtype, cdtype, corrs),
  cell_corrs_within %>% 
    tidyr::unite(dtype, ctype, col = "cdtype", sep = ": ", remove = FALSE) %>%
    dplyr::select(method, ctype, dtype, cdtype, corrs)
)
gene_cell_corrs_within$method <- factor(gene_cell_corrs_within$method, levels = methods_short$method_short[match(methods, methods_short$method)])


pdf(gsub("rds$", "pdf", outrds), width = 15, height = 13)

ggplot(gene_corrs, aes(x = ctype, y = corrs, fill = ctype)) + 
  geom_boxplot() + 
  facet_grid(method1 ~ method2) + 
  theme_bw() + 
  theme(strip.text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        title = element_text(size = 20),
        legend.position = "none") + 
  labs(title = "Correlation, abundances and velocities, by gene",
       subtitle = "For genes selected by all methods")

ggplot(cell_corrs, aes(x = ctype, y = corrs, fill = ctype)) + 
  geom_boxplot() + 
  facet_grid(method1 ~ method2) + 
  theme_bw() + 
  theme(strip.text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        title = element_text(size = 20)) + 
  labs(title = "Correlation, abundances and velocities, by cell",
       subtitle = "Over genes selected by all methods")

if (!any(is.na(cell_corrs$cluster))) {
  for (ct in unique(cell_corrs$ctype)) {
    print(ggplot(cell_corrs %>% dplyr::filter(ctype == ct), 
                 aes(x = cluster, y = corrs, fill = cluster)) + 
            geom_boxplot() + 
            facet_grid(method1 ~ method2) + 
            theme_bw() + 
            theme(strip.text = element_text(size = 5),
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                  title = element_text(size = 20)) + 
            labs(title = paste0("Correlation, counts and velocities, by cell (", ct, ")"),
                 subtitle = "Over genes selected by all methods"))
  }
}

dev.off()

png(gsub("\\.rds$", "_genes_plus_cells.png", outrds), width = 15, height = 13,
    unit = "in", res = 300)
ggplot(gene_cell_corrs, aes(x = ctype, y = corrs, fill = dtype)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") + 
  geom_boxplot() + 
  facet_grid(method1 ~ method2) + 
  theme_bw() + 
  theme(strip.text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        title = element_text(size = 18)) + 
  labs(title = "Correlation, abundances and velocities, by gene and cell",
       subtitle = "Using genes selected by all methods",
       x = "", y = "Correlation") + 
  scale_fill_manual(values = c(cell = "red", gene = "blue"), name = "")
dev.off()

png(gsub("\\.rds$", "_genes_plus_cells_within.png", outrds), width = 9, height = 9,
    unit = "in", res = 300)
ggplot(gene_cell_corrs_within, aes(x = cdtype, y = corrs, fill = dtype)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") + 
  geom_boxplot() + 
  facet_wrap(~ method) + 
  theme_bw() + 
  theme(strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        title = element_text(size = 18)) + 
  labs(title = "Correlation, abundances and velocities, by gene and cell",
       subtitle = "Using genes selected by all methods",
       x = "", y = "Correlation") + 
  scale_fill_manual(values = c(cell = "red", gene = "blue"), name = "")
dev.off()

saveRDS(NULL, file = outrds)

date()
sessionInfo()
