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

print(topdir)
print(plothelperscript)
print(methods)
print(outrds)

## Number of neighbors
nn <- 10

source(plothelperscript)

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods
methods_short <- shorten_methods(methods)

dimreds <- c("PCA", "TSNE", "UMAP")
inputs <- c("alevin_spliced", "starsolo", "starsolo_unspliced",
            "starsolo_summed", "starsolo_concatenated")

dfiles <- expand.grid(dimreds, inputs) %>%
  tidyr::unite(col = "dfile", Var1, Var2, sep = "_") %>% 
  dplyr::pull(dfile)
dfiles <- expand.grid(methods, dfiles) %>%
  dplyr::mutate(file = paste0(topdir, "/plots/velocity/anndata_", Var1, 
                              "/anndata_", Var1, "_velocity_", Var2, ".csv"),
                sce = paste0(topdir, "/output/sce/sce_", Var1, ".rds"))
stopifnot(all(file.exists(dfiles$file)))
stopifnot(all(file.exists(dfiles$sce)))

## Calculate the length of the velocity projection for each cell
## The idea here is that if this projection is long, one is relatively 
## sure of which direction a cell is heading. 
res <- do.call(dplyr::bind_rows, lapply(seq_len(nrow(dfiles)), function(i) {
  read.csv(dfiles$file[i], header = TRUE, as.is = TRUE) %>%
    dplyr::mutate(velocity_length = sqrt(X0 ^ 2 + X1 ^ 2),
                  method = dfiles$Var1[i], 
                  dimred = dfiles$Var2[i]) %>%
    dplyr::left_join(methods_short, by = "method") %>%
    dplyr::select(method_short, mtype, rtype, dimred, velocity_length)
}))

pdf(gsub("rds", "pdf", outrds), width = 10, height = 10)

ggplot(res, aes(x = dimred, y = velocity_length, fill = mtype)) + 
  geom_boxplot() + 
  facet_wrap(~ method_short) + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = "", y = "Projected velocity length",
       title = "Length of projected velocity vector")

## Calculate the average dot product between the velocity projection of 
## a cell and those of its nearest neighbors. 
## The idea is that if the low-dimensional representation is good, the 
## dynamics should be visible in there (and the dot product should be high)
resdot <- do.call(dplyr::bind_rows, lapply(seq_len(nrow(dfiles)), function(i) {
  print(dfiles[i, ])
  a <- read.csv(dfiles$file[i], header = TRUE, as.is = TRUE) %>%
    tibble::column_to_rownames("index") %>%
    as.matrix()
  sce <- readRDS(dfiles$sce[i])
  ## dot product of velocity projections
  dots <- a %*% t(a)
  ## distances in low-dimensional representation
  dists <- as.matrix(dist(reducedDim(sce, dfiles$Var2[i])))
  stopifnot(all(rownames(dots) == rownames(dists)))
  stopifnot(all(colnames(dots) == colnames(dists)))
  data.frame(cell = rownames(dists),
             method = dfiles$Var1[i],
             dimred = dfiles$Var2[i],
             nn = nn,
             ave_dot = rowSums(dots * t(apply(dists, 1, function(x) {
               seq_len(ncol(dists)) %in% order(x)[2:(nn + 1)]
             }))),
             stringsAsFactors = FALSE) %>%
    dplyr::left_join(methods_short, by = "method") %>%
    dplyr::select(method_short, mtype, rtype, dimred, ave_dot)
}))

ggplot(resdot, aes(x = dimred, y = ave_dot, fill = mtype)) + 
  geom_boxplot() + 
  facet_wrap(~ method_short) + 
  scale_fill_manual(values = base_method_colors, name = "") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = "", y = "Average dot product",
       title = "Average dot product with projected velocities of neighboring cells")

dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
