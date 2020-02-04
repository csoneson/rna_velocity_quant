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
print(outrds)

## Number of neighbors
nn <- 10

methods_short <- shorten_methods(methods)

# dimreds <- c("PCA", "TSNE", "UMAP")
# inputs <- c("alevin_spliced", "starsolo", "starsolo_unspliced",
#             "starsolo_summed", "starsolo_concatenated")
# dfiles <- expand.grid(dimreds, inputs) %>%
#   tidyr::unite(col = "dfile", Var1, Var2, sep = "_") %>% 
#   dplyr::pull(dfile)
dfiles <- c("PCA_alevin_spliced_gentrome", "TSNE_alevin_spliced_gentrome",
            "UMAP_alevin_spliced_gentrome",
            "UMAP_starsolo", "UMAP_starsolo_unspliced", "UMAP_starsolo_summed",
            "UMAP_starsolo_concatenated")

dfiles <- expand.grid(methods, dfiles) %>%
  dplyr::mutate(Var1 = as.character(Var1),
                Var2 = as.character(Var2)) %>% 
  dplyr::mutate(veloemb = paste0(topdir, "/plots/velocity/anndata_", Var1, 
                                 "/anndata_", Var1, "_velocity_", Var2, ".csv"),
                veloall = paste0(topdir, "/output/anndata_with_velocity/anndata_", 
                                 Var1, "_with_velocity.h5ad"),
                veloshared = paste0(topdir, "/output/anndata_with_velocity/anndata_",
                                    Var1, "_shared_genes_with_velocity.h5ad"),
                sce = paste0(topdir, "/output/sce/sce_", Var1, ".rds"))
stopifnot(all(file.exists(dfiles$veloemb)))
stopifnot(all(file.exists(dfiles$veloall)))
stopifnot(all(file.exists(dfiles$veloshared)))
stopifnot(all(file.exists(dfiles$sce)))

## Calculate the length of the velocity projection for each cell
## The idea here is that if this projection is long, one is relatively 
## sure of which direction a cell is heading. 
res <- do.call(dplyr::bind_rows, lapply(seq_len(nrow(dfiles)), function(i) {
  read.csv(dfiles$veloemb[i], header = TRUE, as.is = TRUE) %>%
    dplyr::mutate(velocity_length = sqrt(X0 ^ 2 + X1 ^ 2),
                  method = dfiles$Var1[i], 
                  dimred = dfiles$Var2[i]) %>%
    dplyr::left_join(methods_short, by = "method") %>%
    dplyr::select(method_short, mtype, rtype, dimred, velocity_length)
}))

pdf(gsub("\\.rds", "_velocity_length.pdf", outrds), width = 9, height = 9)
ggplot(res, aes(x = dimred, y = velocity_length, fill = mtype)) + 
  geom_boxplot() + 
  facet_wrap(~ method_short) + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = "", y = "Embedded velocity length",
       title = "Length of embedded velocity vector")
dev.off()

## Calculate the average dot product between the velocity projection of 
## a cell and those of its nearest neighbors. 
## The idea is that if the low-dimensional representation is good, the 
## dynamics should be visible in there (and the dot product should be high)
## Then do the same for the high-dimensional velocities, checking whether
## points mapping close together in the low-dim space also have similar velocities
## Distances are always the same for a given dimension reduction, so pre-calculate these
sce <- readRDS(dfiles$sce[1])
tmp <- unique(dfiles$Var2)
names(tmp) <- tmp
distsdr <- lapply(tmp, function(dr) {
  as.matrix(dist(reducedDim(sce, dr)))
})
resdot <- do.call(dplyr::bind_rows, lapply(seq_len(nrow(dfiles)), function(i) {
  print(dfiles[i, ])
  
  ## distances in low-dimensional representation
  dists <- distsdr[[dfiles$Var2[i]]]
  
  set.seed(1234)
  s <- sample(seq_len(ncol(dists)), ncol(dists))
  
  ## low-dimensional embedding of velocities
  a <- read.csv(dfiles$veloemb[i], header = TRUE, as.is = TRUE) %>%
    tibble::column_to_rownames("index") %>%
    as.matrix()
  ## dot product of velocity projections
  dots <- a %*% t(a)
  stopifnot(all(rownames(dots) == rownames(dists)))
  stopifnot(all(colnames(dots) == colnames(dists)))
  d1 <- data.frame(cell = rownames(dists),
                   method = dfiles$Var1[i],
                   dimred = dfiles$Var2[i],
                   dtype = "Embedding",
                   nn = nn,
                   ave_dot = rowSums(dots * t(apply(dists, 1, function(x) {
                     seq_len(ncol(dists)) %in% order(x)[2:(nn + 1)]
                   }))),
                   ave_dot_shuffled = rowSums(dots * t(apply(dists, 1, function(x) {
                     seq_len(ncol(dists)) %in% order(x)[s][2:(nn + 1)]
                   }))),
                   stringsAsFactors = FALSE) %>%
    dplyr::left_join(methods_short, by = "method") %>%
    dplyr::select(cell, method_short, mtype, rtype, dimred, dtype, ave_dot, ave_dot_shuffled)
  
  # ## high-dimensional velocities, individual genes
  # a <- Seurat::ReadH5AD(dfiles$veloall[i])
  # a <- as.matrix(GetAssayData(GetAssay(a, "velocity")))
  # a <- a[rowSums(is.na(a)) == 0, ]
  # ## scale, to get more manageable values
  # a <- a/sqrt(nrow(a))
  # ## dot product of velocities
  # dots <- t(a) %*% a
  # stopifnot(all(rownames(dots) == rownames(dists)))
  # stopifnot(all(colnames(dots) == colnames(dists)))
  # d2 <- data.frame(cell = rownames(dists),
  #                  method = dfiles$Var1[i],
  #                  dimred = dfiles$Var2[i],
  #                  dtype = "Velocities, all genes",
  #                  nn = nn,
  #                  ave_dot = rowSums(dots * t(apply(dists, 1, function(x) {
  #                    seq_len(ncol(dists)) %in% order(x)[2:(nn + 1)]
  #                  }))),
  #                  ave_dot_shuffled = rowSums(dots * t(apply(dists, 1, function(x) {
  #                    seq_len(ncol(dists)) %in% order(x)[s][2:(nn + 1)]
  #                  }))),
  #                  stringsAsFactors = FALSE) %>%
  #   dplyr::left_join(methods_short, by = "method") %>%
  #   dplyr::select(cell, method_short, mtype, rtype, dimred, dtype, ave_dot, ave_dot_shuffled)
  
  # ## high-dimensional velocities, shared genes
  # a <- Seurat::ReadH5AD(dfiles$veloshared[i])
  # a <- as.matrix(GetAssayData(GetAssay(a, "velocity")))
  # a <- a[rowSums(is.na(a)) == 0, ]
  # ## scale, to get more manageable values
  # a <- a/sqrt(nrow(a))
  # ## dot product of velocities
  # dots <- t(a) %*% a
  # stopifnot(all(rownames(dots) == rownames(dists)))
  # stopifnot(all(colnames(dots) == colnames(dists)))
  # d3 <- data.frame(cell = rownames(dists),
  #                  method = dfiles$Var1[i],
  #                  dimred = dfiles$Var2[i],
  #                  dtype = "Velocities, shared genes",
  #                  nn = nn,
  #                  ave_dot = rowSums(dots * t(apply(dists, 1, function(x) {
  #                    seq_len(ncol(dists)) %in% order(x)[2:(nn + 1)]
  #                  }))),
  #                  ave_dot_shuffled = rowSums(dots * t(apply(dists, 1, function(x) {
  #                    seq_len(ncol(dists)) %in% order(x)[s][2:(nn + 1)]
  #                  }))),
  #                  stringsAsFactors = FALSE) %>%
  #   dplyr::left_join(methods_short, by = "method") %>%
  #   dplyr::select(cell, method_short, mtype, rtype, dimred, dtype, ave_dot, ave_dot_shuffled)
  
  # dplyr::bind_rows(d1, d2, d3)
  d1
}))

pdf(gsub("\\.rds", "_ave_dot.pdf", outrds), width = 9, height = 9)
ggplot(resdot %>% tidyr::gather(key = "drtype", value = "ave_dot", 
                                ave_dot, ave_dot_shuffled) %>%
         dplyr::filter(drtype == "ave_dot") %>%
         tidyr::unite("dimred", dimred, drtype, sep = "_") %>%
         dplyr::mutate(dimred = gsub("_ave_dot", "", dimred)), 
       aes(x = dimred, y = ave_dot, fill = mtype)) + 
  geom_boxplot() + 
  facet_wrap(~ method_short) + 
  # facet_wrap(method_short ~ dtype, scales = "free_y", ncol = 3) + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = "", y = "Average dot product",
       title = "Average dot product with velocities of neighboring cells")

dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
