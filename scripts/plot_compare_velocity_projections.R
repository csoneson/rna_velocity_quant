args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
})

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(topdir)
print(plothelperscript)
print(dataset)
print(methods)
print(outrds)

source(plothelperscript)
methods_short <- shorten_methods(methods)

## Read one data set (only to get the UMAP_alevin_spliced representation)
sce <- readRDS(file.path(topdir, paste0("output/sce/sce_", methods[1], ".rds")))

## ------------------------------------------------------------------------- ##
## Compare velocity projections with "all" genes to those with only shared 
## genes, for each method
## ------------------------------------------------------------------------- ##
plotlist <- lapply(methods, function(m) {
  sm <- methods_short$method_short[match(m, methods_short$method)]
  v_all <- read.delim(
    file.path(topdir, paste0("plots/velocity/anndata_", m, "/anndata_", 
                             m, "_velocity_UMAP_alevin_spliced_gentrome.csv")), 
    header = TRUE, as.is = TRUE, sep = ",")
  v_shared <- read.delim(
    file.path(topdir, paste0("plots/velocity/anndata_", m, "_shared_genes/anndata_", 
                             m, "_shared_genes_velocity_UMAP_alevin_spliced_gentrome.csv")), 
    header = TRUE, as.is = TRUE, sep = ",")
  
  tmp <- v_all %>% dplyr::rename(X0_all = X0, X1_all = X1) %>%
    dplyr::inner_join(v_shared %>% dplyr::rename(X0_shared = X0, X1_shared = X1),
                      by = "index") %>%
    dplyr::mutate(X0sum = X0_all + X0_shared,
                  X1sum = X1_all + X1_shared) %>%
    dplyr::mutate(lengthOfSum = sqrt(X0sum ^ 2 + X1sum ^ 2),
                  sumOfLength = sqrt(X0_all ^ 2 + X1_all ^ 2) + sqrt(X0_shared ^ 2 + X1_shared ^ 2)) %>%
    dplyr::mutate(los_over_sol = lengthOfSum/sumOfLength)
  
  stopifnot(all(tmp$index == colnames(sce)))
  df <- data.frame(reducedDim(sce, "UMAP_alevin_spliced_gentrome"), los_over_sol = tmp$los_over_sol,
                   stringsAsFactors = FALSE)
  ggplot(df, aes(x = X1, y = X2, color = los_over_sol)) + 
    geom_point(size = 0.8, alpha = 0.5) +
    scale_color_gradientn(colors = viridis::viridis(21), na.value = "grey50", limits = c(0, 1),
                          name = "Concordance between velocity projections") + 
    labs(x = "UMAP_alevin_spliced 1", y = "UMAP_alevin_spliced 2",
         title = sm) + 
    theme_void()
})

png(gsub("\\.rds", "_all_vs_shared_genes_bymethod.png", outrds), width = 8, height = 12,
    unit = "in", res = 200)
print(cowplot::plot_grid(
  cowplot::plot_grid(plotlist = lapply(plotlist, function(w) w + theme(legend.position = "none")), ncol = 3),
  cowplot::get_legend(plotlist[[1]] + theme(legend.position = "bottom")),
  ncol = 1, rel_heights = c(1, 0.05)
))
dev.off()

## ------------------------------------------------------------------------- ##
## Compare velocity projections across methods, using "all"/shared genes
## ------------------------------------------------------------------------- ##
## Help function
summarize_res <- function(res) {
  ## For each cell, calculate the sum of the individual (method-wise) 
  ## velocity vector lengths
  res1 <- res %>% dplyr::group_by(index, method) %>%
    dplyr::summarize(speed = sqrt(sum(value ^ 2))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(index) %>%
    dplyr::summarize(sum_speed = sum(speed))
  
  ## For each cell, calculate the sum of the velocity vectors across methods, 
  ## and then the length of that sum
  res2 <- res %>% dplyr::group_by(index, coord) %>%
    dplyr::summarize(value = sum(value)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(index) %>%
    dplyr::summarize(speed_sum = sqrt(sum(value ^ 2)))
  
  ## Combine
  speeds <- res1 %>% dplyr::inner_join(res2, by = "index") %>%
    dplyr::mutate(speed_ratio = speed_sum/sum_speed)
  speeds
}

## ------------------------------------------------------------------------- ##
## "All" (individually selected) genes
## ------------------------------------------------------------------------- ##
res <- do.call(dplyr::bind_rows, lapply(methods, function(m) {
  read.delim(
    file.path(topdir, paste0("plots/velocity/anndata_", m, "/anndata_", 
                             m, "_velocity_UMAP_alevin_spliced_gentrome.csv")), 
    header = TRUE, as.is = TRUE, sep = ",") %>%
    tidyr::gather(key = "coord", value = "value", X0, X1) %>%
    dplyr::mutate(method = m)
}))

speeds <- summarize_res(res)
stopifnot(length(intersect(speeds$index, colnames(sce))) == ncol(sce))
speeds <- speeds[match(colnames(sce), speeds$index), ]
stopifnot(all(speeds$index == colnames(sce)))
  
df <- data.frame(reducedDim(sce, "UMAP_alevin_spliced_gentrome"), speed_ratio = speeds$speed_ratio,
                 stringsAsFactors = FALSE)
g1 <- ggplot(df, aes(x = X1, y = X2, color = speed_ratio)) + 
  geom_point(alpha = 0.5) +
  scale_color_gradientn(colors = viridis::viridis(21), na.value = "grey50", limits = c(0, 1),
                        name = "Concordance between\nvelocity projections") + 
  labs(x = "UMAP_alevin_spliced 1", y = "UMAP_alevin_spliced 2",
       title = paste0(dataset, ", individually selected genes")) + 
  theme_bw() + theme(legend.position = "bottom")

## ------------------------------------------------------------------------- ##
## Shared genes
## ------------------------------------------------------------------------- ##
res <- do.call(dplyr::bind_rows, lapply(methods, function(m) {
  read.delim(
    file.path(topdir, paste0("plots/velocity/anndata_", m, "_shared_genes/anndata_", 
                             m, "_shared_genes_velocity_UMAP_alevin_spliced_gentrome.csv")), 
    header = TRUE, as.is = TRUE, sep = ",") %>%
    tidyr::gather(key = "coord", value = "value", X0, X1) %>%
    dplyr::mutate(method = m)
}))

speeds <- summarize_res(res)
stopifnot(length(intersect(speeds$index, colnames(sce))) == ncol(sce))
speeds <- speeds[match(colnames(sce), speeds$index), ]
stopifnot(all(speeds$index == colnames(sce)))

df <- data.frame(reducedDim(sce, "UMAP_alevin_spliced_gentrome"), speed_ratio = speeds$speed_ratio,
                 stringsAsFactors = FALSE)
g2 <- ggplot(df, aes(x = X1, y = X2, color = speed_ratio)) + 
  geom_point(alpha = 0.5) +
  scale_color_gradientn(colors = viridis::viridis(21), na.value = "grey50", limits = c(0, 1),
                        name = "Concordance between\nvelocity projections") + 
  labs(x = "UMAP_alevin_spliced 1", y = "UMAP_alevin_spliced 2",
       title = paste0(dataset, ", shared genes")) + 
  theme_bw() + theme(legend.position = "bottom")


png(gsub("\\.rds", "_across_methods.png", outrds), width = 4.5, height = 8,
    unit = "in", res = 200)
print(cowplot::plot_grid(
  cowplot::plot_grid(
    g1 + theme(legend.position = "none"), 
    g2 + theme(legend.position = "none"),
    ncol = 1, rel_heights = c(1, 1)),
  cowplot::get_legend(g1), ncol = 1, rel_heights = c(1, 0.2)
  ))
dev.off()


saveRDS(NULL, file = outrds)
date()
sessionInfo()
