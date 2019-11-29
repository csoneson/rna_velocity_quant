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

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(topdir)
print(methods)
print(genetxt)
print(outrds)

sces <- lapply(methods, function(nm) {
  readRDS(file.path(topdir, paste0("output/sce/sce_", nm, ".rds")))
})
for (nm in names(sces)) {
  metadata(sces[[nm]]) <- list(dataset = nm)
}

genes <- read.delim(genetxt, header = FALSE, as.is = TRUE)[, 1]
sces <- lapply(sces, function(w) w[match(genes, rownames(w)), ])

seurats <- lapply(methods, function(nm) {
  w <- ReadH5AD(file.path(topdir, paste0("output/anndata_with_velocity/anndata_", nm, "_shared_genes_with_velocity.h5ad")))
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
        a <- scale(t(as.matrix(assay(sces[[jj]], "spliced"))), center = TRUE, scale = TRUE)
        b <- scale(t(as.matrix(assay(sces[[kk]], "spliced"))), center = TRUE, scale = TRUE)
        corrs_spliced <- colSums(a * b)/(nrow(a) - 1)
        a <- scale(t(as.matrix(assay(sces[[jj]], "unspliced"))), center = TRUE, scale = TRUE)
        b <- scale(t(as.matrix(assay(sces[[kk]], "unspliced"))), center = TRUE, scale = TRUE)
        corrs_unspliced <- colSums(a * b)/(nrow(a) - 1)
        a <- scale(t(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "velocity")))), center = TRUE, scale = TRUE)
        b <- scale(t(as.matrix(GetAssayData(GetAssay(seurats[[kk]], "velocity")))), center = TRUE, scale = TRUE)
        corrs_velocity <- colSums(a * b)/(nrow(a) - 1)
        a <- scale(t(as.matrix(assay(sces[[jj]], "spliced"))) + 
                     t(as.matrix(assay(sces[[jj]], "unspliced"))), center = TRUE, scale = TRUE)
        b <- scale(t(as.matrix(assay(sces[[kk]], "spliced"))) + 
                     t(as.matrix(assay(sces[[kk]], "unspliced"))), center = TRUE, scale = TRUE)
        corrs_total <- colSums(a * b)/(nrow(a) - 1)
        dplyr::bind_rows(
          data.frame(method1 = jj,
                     method2 = kk, 
                     gene = names(corrs_spliced), 
                     ctype = "spliced",
                     corrs = corrs_spliced, 
                     stringsAsFactors = FALSE),
          data.frame(method1 = jj,
                     method2 = kk,
                     gene = names(corrs_unspliced), 
                     ctype = "unspliced",
                     corrs = corrs_unspliced, 
                     stringsAsFactors = FALSE),
          data.frame(method1 = jj,
                     method2 = kk,
                     gene = names(corrs_velocity), 
                     ctype = "velocity",
                     corrs = corrs_velocity, 
                     stringsAsFactors = FALSE),
          data.frame(method1 = jj,
                     method2 = kk,
                     gene = names(corrs_total), 
                     ctype = "total",
                     corrs = corrs_total, 
                     stringsAsFactors = FALSE)
        )
      }))
  }))

gene_corrs$method1 <- factor(gene_corrs$method1, levels = methods[methods %in% gene_corrs$method1])
gene_corrs$method2 <- factor(gene_corrs$method2, levels = methods[methods %in% gene_corrs$method2])

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
        a <- scale(as.matrix(assay(sces[[jj]], "spliced")), center = TRUE, scale = TRUE)
        b <- scale(as.matrix(assay(sces[[kk]], "spliced")), center = TRUE, scale = TRUE)
        tmp <- a * b
        corrs_spliced <- colSums(tmp, na.rm = TRUE)/(colSums(!is.na(tmp)) - 1)
        a <- scale(as.matrix(assay(sces[[jj]], "unspliced")), center = TRUE, scale = TRUE)
        b <- scale(as.matrix(assay(sces[[kk]], "unspliced")), center = TRUE, scale = TRUE)
        tmp <- a * b
        corrs_unspliced <- colSums(tmp, na.rm = TRUE)/(colSums(!is.na(tmp)) - 1)
        a <- scale(as.matrix(GetAssayData(GetAssay(seurats[[jj]], "velocity"))), center = TRUE, scale = TRUE)
        b <- scale(as.matrix(GetAssayData(GetAssay(seurats[[kk]], "velocity"))), center = TRUE, scale = TRUE)
        tmp <- a * b
        corrs_velocity <- colSums(tmp, na.rm = TRUE)/(colSums(!is.na(tmp)) - 1)
        a <- scale(as.matrix(assay(sces[[jj]], "spliced")) + 
                     as.matrix(assay(sces[[jj]], "unspliced")), center = TRUE, scale = TRUE)
        b <- scale(as.matrix(assay(sces[[kk]], "spliced")) + 
                     as.matrix(assay(sces[[kk]], "unspliced")), center = TRUE, scale = TRUE)
        tmp <- a * b
        corrs_total <- colSums(tmp, na.rm = TRUE)/(colSums(!is.na(tmp)) - 1)
        dplyr::bind_rows(
          data.frame(method1 = jj,
                     method2 = kk, 
                     cell = names(corrs_spliced), 
                     cluster = sces[[jj]]$clusters,
                     ctype = "spliced",
                     corrs = corrs_spliced, 
                     stringsAsFactors = FALSE),
          data.frame(method1 = jj,
                     method2 = kk,
                     cell = names(corrs_unspliced), 
                     cluster = sces[[jj]]$clusters,
                     ctype = "unspliced",
                     corrs = corrs_unspliced, 
                     stringsAsFactors = FALSE),
          data.frame(method1 = jj,
                     method2 = kk,
                     cell = names(corrs_velocity), 
                     cluster = levels(factor(sces[[jj]]$clusters))[seurats[[jj]]@meta.data$clusters + 1],
                     ctype = "velocity",
                     corrs = corrs_velocity, 
                     stringsAsFactors = FALSE),
          data.frame(method1 = jj,
                     method2 = kk,
                     cell = names(corrs_total), 
                     cluster = sces[[jj]]$clusters,
                     ctype = "total",
                     corrs = corrs_total, 
                     stringsAsFactors = FALSE)
        )
      }))
  }))

cell_corrs$method1 <- factor(cell_corrs$method1, levels = methods[methods %in% cell_corrs$method1])
cell_corrs$method2 <- factor(cell_corrs$method2, levels = methods[methods %in% cell_corrs$method2])


pdf(gsub("rds$", "pdf", outrds), width = 20, height = 20)

ggplot(gene_corrs, aes(x = ctype, y = corrs, fill = ctype)) + 
  geom_boxplot() + 
  facet_grid(method1 ~ method2) + 
  theme_bw() + 
  theme(strip.text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        title = element_text(size = 20)) + 
  labs(title = "Correlation, counts and velocities, by gene",
       subtitle = "For genes selected by all methods")

ggplot(cell_corrs, aes(x = ctype, y = corrs, fill = ctype)) + 
  geom_boxplot() + 
  facet_grid(method1 ~ method2) + 
  theme_bw() + 
  theme(strip.text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        title = element_text(size = 20)) + 
  labs(title = "Correlation, counts and velocities, by cell",
       subtitle = "Over genes selected by all methods")

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

dev.off()

saveRDS(NULL, file = outrds)

date()
sessionInfo()
