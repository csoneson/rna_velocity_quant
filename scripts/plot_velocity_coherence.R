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
  library(Seurat)
})
source(plothelperscript)

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

## Plot total length of velocity vectors for different methods
## Check number of NAs for different methods
## Potentially reduce to genes with non-NA value in all methods?

print(plothelperscript)
print(topdir)
print(dataset)
print(methods)
print(genetxt)
print(outrds)

methods_short <- shorten_methods(methods)

genes <- read.delim(genetxt, header = FALSE, as.is = TRUE)[, 1]

## Read counts, subset to genes selected in all scVelo runs
sces <- lapply(methods, function(nm) {
  readRDS(file.path(topdir, paste0("output/sce/sce_", nm, ".rds")))
})
sces <- lapply(sces, function(w) w[match(genes, rownames(w)), ])

## Read velocities from scVelo run on shared genes
seurats <- lapply(methods, function(nm) {
  w <- ReadH5AD(file.path(topdir, paste0("output/anndata_with_velocity/anndata_", 
                                         nm, "_shared_genes_with_velocity.h5ad")))
  stopifnot(all(rownames(w) == genes))
  w
})
velocities <- lapply(seurats, function(w) as.matrix(GetAssayData(GetAssay(w, "velocity"))))
cells <- colnames(velocities[[1]])
genes <- rownames(velocities[[1]])
stopifnot(all(sapply(velocities, function(w) all(colnames(w) == cells))))
stopifnot(all(sapply(velocities, function(w) all(rownames(w) == genes))))

## Check number of genes with NA velocities for each method
na_vel <- do.call(cbind, lapply(velocities, function(w) is.na(rowSums(w)))) + 0
pdf(gsub("\\.rds$", "_genes_with_na_velocity.pdf", outrds), height = 5, width = 10)
UpSetR::upset(data.frame(na_vel), nsets = ncol(na_vel), nintersects = 1000, 
              order.by = "freq", decreasing = TRUE)
dev.off()

## For each method, check whether genes with NA velocities have higher/lower counts
na_vel_vs_count <- do.call(dplyr::bind_rows, lapply(methods, function(m) {
  dplyr::inner_join(data.frame(gene = rownames(velocities[[m]]),
                               na_vel = is.na(rowSums(velocities[[m]])),
                               stringsAsFactors = FALSE),
                    data.frame(gene = rownames(sces[[m]]),
                               spliced_count = rowSums(assay(sces[[m]], "spliced")),
                               unspliced_count = rowSums(assay(sces[[m]], "unspliced")),
                               total_count = rowSums(assay(sces[[m]], "spliced")) + 
                                 rowSums(assay(sces[[m]], "unspliced")),
                               spliced_fraczero = rowMeans(assay(sces[[m]], "spliced") == 0),
                               unspliced_fraczero = rowMeans(assay(sces[[m]], "unspliced") == 0),
                               stringsAsFactors = FALSE),
                    by = "gene"
  ) %>% dplyr::mutate(method = m)
})) %>% dplyr::left_join(methods_short, by = "method")

pdf(gsub("\\.rds$", "_na_vel_vs_count.pdf", outrds), width = 12, height = 8)
cowplot::plot_grid(
  ggplot(na_vel_vs_count, aes(x = na_vel, y = spliced_count)) + 
    geom_boxplot(aes(fill = mtype), alpha = 0.5, outlier.size = -1) + 
    geom_point(alpha = 0.5, position = position_jitter(width = 0.15)) + 
    facet_wrap(~ method, ncol = 3) + 
    theme_bw() + scale_y_log10() + 
    scale_fill_manual(values = base_method_colors) + 
    theme(legend.position = "none") + 
    labs(x = "Velocity is N/A", y = "UMI count across cells", title = "Spliced count"),
  ggplot(na_vel_vs_count, aes(x = na_vel, y = unspliced_count)) + 
    geom_boxplot(aes(fill = mtype), alpha = 0.5, outlier.size = -1) + 
    geom_point(alpha = 0.5, position = position_jitter(width = 0.15)) + 
    facet_wrap(~ method, ncol = 3) + 
    theme_bw() + scale_y_log10() + 
    scale_fill_manual(values = base_method_colors) + 
    theme(legend.position = "none") + 
    labs(x = "Velocity is N/A", y = "UMI count across cells", title = "Unspliced count"),
  nrow = 1
)
dev.off()

pdf(gsub("\\.rds$", "_na_vel_vs_fraczero.pdf", outrds), width = 12, height = 8)
cowplot::plot_grid(
  ggplot(na_vel_vs_count, aes(x = na_vel, y = spliced_fraczero)) + 
    geom_boxplot(aes(fill = mtype), alpha = 0.5, outlier.size = -1) + 
    geom_point(alpha = 0.5, position = position_jitter(width = 0.15)) + 
    facet_wrap(~ method, ncol = 3) + 
    theme_bw() + 
    scale_fill_manual(values = base_method_colors) + 
    theme(legend.position = "none") + 
    labs(x = "Velocity is N/A", y = "UMI count across cells", title = "Spliced count"),
  ggplot(na_vel_vs_count, aes(x = na_vel, y = unspliced_fraczero)) + 
    geom_boxplot(aes(fill = mtype), alpha = 0.5, outlier.size = -1) + 
    geom_point(alpha = 0.5, position = position_jitter(width = 0.15)) + 
    facet_wrap(~ method, ncol = 3) + 
    theme_bw() + 
    scale_fill_manual(values = base_method_colors) + 
    theme(legend.position = "none") + 
    labs(x = "Velocity is N/A", y = "UMI count across cells", title = "Unspliced count"),
  nrow = 1
)
dev.off()

## Find genes with largest variability of velocities
vs <- Reduce('+', velocities)
vs <- apply(vs, 1, function(w) sqrt(sum(w ^ 2)))
sn <- lapply(velocities, function(w) apply(w, 1, function(w) sqrt(sum(w ^ 2))))
sn <- Reduce('+', sn)
concordance <- vs/sn



lapply(velocities, function(w) {
  table(rowSums(is.na(w)))
})
valid_in_all <- Reduce(
  intersect, lapply(velocities, function(w) rownames(w)[rowSums(is.na(w)) == 0])
)
velocities_sub <- lapply(velocities, function(w) w[match(valid_in_all, rownames(w)), ])

velocities_total_length <- data.frame(sapply(velocities, function(w) {
  sqrt(sum(w ^ 2, na.rm = TRUE))
})) %>%
  tibble::rownames_to_column("method") %>%
  setNames(c("method", "total_length"))
velocities_sub_total_length <- data.frame(sapply(velocities_sub, function(w) {
  sqrt(sum(w ^ 2, na.rm = TRUE))
})) %>%
  tibble::rownames_to_column("method") %>%
  setNames(c("method", "total_length_allvalid"))

## Total count vs total contribution to velocity
tcvs <- lapply(names(sces), function(i) {
  data.frame(method = i,
             gene = rownames(sces[[i]]), 
             total_count = rowSums(assay(sces[[i]], "spliced")) + 
               rowSums(assay(sces[[i]], "unspliced")),
             stringsAsFactors = FALSE) %>%
    dplyr::left_join(
      data.frame(gene = rownames(velocities[[i]]),
                 total_velocity = sqrt(rowSums(velocities[[i]] ^ 2)),
                 stringsAsFactors = FALSE),
      by = "gene"
    )
})
cowplot::plot_grid(plotlist = lapply(tcvs, function(df) {
  m <- unique(df$method)
  ggplot(df, aes(x = total_count, y = total_velocity)) + 
    geom_point(alpha = 0.5) + 
    scale_x_log10() + scale_y_log10() + 
    labs(title = m)
}))

## Normalize velocity vector for each cell 
velocities_norm <- lapply(velocities, function(w) {
  scale(w, center = FALSE, scale = TRUE)
})
velocities_sub_norm <- lapply(velocities_sub, function(w) {
  scale(w, center = FALSE, scale = TRUE)
})


## Mean velocity across all methods
velocities_means <- Reduce(`+`, velocities)/length(velocities)
velocities_norm_means <- Reduce(`+`, velocities_norm)/length(velocities_norm)
velocities_sub_means <- Reduce(`+`, velocities_sub)/length(velocities_sub)
velocities_sub_norm_means <- Reduce(`+`, velocities_sub_norm)/length(velocities_sub_norm)

## Length of mean velocity vector for each cell
velocities_means_length <- sqrt(colSums(velocities_means ^ 2, na.rm = TRUE))
velocities_norm_means_length <- sqrt(colSums(velocities_norm_means ^ 2, na.rm = TRUE))
velocities_sub_means_length <- sqrt(colSums(velocities_sub_means ^ 2, na.rm = TRUE))
velocities_sub_norm_means_length <- sqrt(colSums(velocities_sub_norm_means ^ 2, na.rm = TRUE))

## Mean length of velocity vector for each cell across methods
velocities_lengths_mean <- 
  Reduce(`+`, lapply(velocities, function(w) sqrt(colSums(w ^ 2, na.rm = TRUE))))/
  length(velocities)
velocities_norm_lengths_mean <- 
  Reduce(`+`, lapply(velocities_norm, function(w) sqrt(colSums(w ^ 2, na.rm = TRUE))))/
  length(velocities_norm)
velocities_sub_lengths_mean <- 
  Reduce(`+`, lapply(velocities_sub, function(w) sqrt(colSums(w ^ 2, na.rm = TRUE))))/
  length(velocities_sub)
velocities_sub_norm_lengths_mean <- 
  Reduce(`+`, lapply(velocities_sub_norm, function(w) sqrt(colSums(w ^ 2, na.rm = TRUE))))/
  length(velocities_sub_norm)

stopifnot(all(names(velocities_lengths_mean) == names(velocities_means_length)))
stopifnot(all(names(velocities_norm_lengths_mean) == names(velocities_norm_means_length)))
stopifnot(all(names(velocities_lengths_mean) == colnames(sces[[1]])))
stopifnot(all(names(velocities_sub_lengths_mean) == names(velocities_sub_means_length)))
stopifnot(all(names(velocities_sub_norm_lengths_mean) == names(velocities_sub_norm_means_length)))
stopifnot(all(names(velocities_sub_lengths_mean) == colnames(sces[[1]])))

sce <- sces[[1]]
sce$velocities_lengths_mean <- velocities_lengths_mean
sce$velocities_norm_lengths_mean <- velocities_norm_lengths_mean
sce$velocities_means_length <- velocities_means_length
sce$velocities_norm_means_length <- velocities_norm_means_length
sce$velocities_sub_lengths_mean <- velocities_sub_lengths_mean
sce$velocities_sub_norm_lengths_mean <- velocities_sub_norm_lengths_mean
sce$velocities_sub_means_length <- velocities_sub_means_length
sce$velocities_sub_norm_means_length <- velocities_sub_norm_means_length

ggplot(as.data.frame(colData(sce)), 
       aes(x = velocities_lengths_mean, y = velocities_means_length)) + 
  geom_point(alpha = 0.5) + 
  theme_bw()
ggplot(as.data.frame(colData(sce)), 
       aes(x = velocities_norm_lengths_mean, y = velocities_norm_means_length)) + 
  geom_point(alpha = 0.5) + 
  theme_bw()
ggplot(as.data.frame(colData(sce)), 
       aes(x = velocities_sub_lengths_mean, y = velocities_sub_means_length)) + 
  geom_point(alpha = 0.5) + 
  theme_bw()
ggplot(as.data.frame(colData(sce)), 
       aes(x = velocities_sub_norm_lengths_mean, y = velocities_sub_norm_means_length)) + 
  geom_point(alpha = 0.5) + 
  theme_bw()

plotdf <- data.frame(reducedDim(sce, "UMAP_alevin_spliced"),
                     colData(sce))
pdf(gsub("\\.rds$", "_umap_concordance.pdf", outrds))
ggplot(plotdf, aes(x = X1, y = X2, color = velocities_norm_means_length)) + 
  geom_point(alpha = 0.5) + 
  labs(x = "UMAP 1", y = "UMAP 2", title = dataset) + 
  scale_color_gradientn(colors = viridis::viridis(21), na.value = "grey50",
                        name = "Length of\nmean\nnormalized\nvelocity\nvector") + 
  theme_bw()
dev.off()

pdf(gsub("\\.rds$", "_umap_concordance_sub.pdf", outrds))
ggplot(plotdf, aes(x = X1, y = X2, color = velocities_sub_norm_means_length)) + 
  geom_point(alpha = 0.5) + 
  labs(x = "UMAP 1", y = "UMAP 2", title = dataset) + 
  scale_color_gradientn(colors = viridis::viridis(21), na.value = "grey50",
                        name = "Length of\nmean\nnormalized\nvelocity\nvector") + 
  theme_bw()
dev.off()


