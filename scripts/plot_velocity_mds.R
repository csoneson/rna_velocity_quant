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
  library(ggrepel)
})
source(plothelperscript)

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(plothelperscript)
print(topdir)
print(dataset)
print(methods)
print(genetxt)
print(outrds)

methods_short <- shorten_methods(methods)

genes <- read.delim(genetxt, header = FALSE, as.is = TRUE)[, 1]

## Read velocities from scVelo run on shared genes
seurats <- lapply(methods, function(nm) {
  w <- ReadH5AD(file.path(topdir, paste0("output/anndata_with_velocity/anndata_", 
                                         nm, "_shared_genes_with_velocity.h5ad")))
  stopifnot(all(rownames(w) == genes))
  w
})
velocities_shared <- lapply(seurats, function(w) as.matrix(GetAssayData(GetAssay(w, "velocity"))))
keep_genes <- Reduce(intersect, lapply(velocities_shared, function(w) rownames(w)[!is.na(rowSums(w))]))
velocities_shared <- lapply(velocities_shared, function(w) w[keep_genes, ])
cells <- colnames(velocities_shared[[1]])
genes <- rownames(velocities_shared[[1]])
stopifnot(all(sapply(velocities_shared, function(w) all(colnames(w) == cells))))
stopifnot(all(sapply(velocities_shared, function(w) all(rownames(w) == genes))))

## Scale velocities for each gene
velocities_shared_scaled <- lapply(velocities_shared, function(w) {
  t(scale(t(w), center = FALSE, scale = TRUE))
})

## Read velocity projections from scVelo run on all genes
projs_all <- lapply(methods, function(nm) {
  w <- read.csv(file.path(topdir, paste0("plots/velocity/anndata_", nm, "/anndata_", 
                                         nm, "_velocity_UMAP_alevin_spliced.csv")), 
                header = TRUE, as.is = TRUE) %>%
    tibble::column_to_rownames("index")
})
cells <- rownames(projs_all[[1]])
stopifnot(all(sapply(projs_all, function(w) all(rownames(w) == cells))))

## Read velocity projections from scVelo run on shared genes
projs_shared <- lapply(methods, function(nm) {
  w <- read.csv(file.path(topdir, paste0("plots/velocity/anndata_", nm, "_shared_genes/anndata_", 
                                         nm, "_shared_genes_velocity_UMAP_alevin_spliced.csv")), 
                header = TRUE, as.is = TRUE) %>%
    tibble::column_to_rownames("index")
})
cells <- rownames(projs_shared[[1]])
stopifnot(all(sapply(projs_shared, function(w) all(rownames(w) == cells))))

## Calculate distances between methods
d_vel_shared <- do.call(dplyr::bind_rows, lapply(names(velocities_shared), function(s1) {
  do.call(dplyr::bind_rows, lapply(names(velocities_shared), function(s2) {
    stopifnot(all(rownames(velocities_shared[[s1]]) == rownames(velocities_shared[[s2]])))
    stopifnot(all(colnames(velocities_shared[[s1]]) == colnames(velocities_shared[[s2]])))
    data.frame(m1 = s1, m2 = s2, 
               dist = sqrt(mean((velocities_shared[[s1]] - 
                                   velocities_shared[[s2]])^2)),
               stringsAsFactors = FALSE)
  }))
}))

d_projs_all <- do.call(dplyr::bind_rows, lapply(names(projs_all), function(s1) {
  do.call(dplyr::bind_rows, lapply(names(projs_all), function(s2) {
    stopifnot(all(rownames(projs_all[[s1]]) == rownames(projs_all[[s2]])))
    stopifnot(all(colnames(projs_all[[s1]]) == colnames(projs_all[[s2]])))
    data.frame(m1 = s1, m2 = s2, 
               dist = sqrt(mean((as.matrix(projs_all[[s1]]) - 
                                   as.matrix(projs_all[[s2]]))^2)),
               stringsAsFactors = FALSE)
  }))
}))

d_projs_shared <- do.call(dplyr::bind_rows, lapply(names(projs_shared), function(s1) {
  do.call(dplyr::bind_rows, lapply(names(projs_shared), function(s2) {
    stopifnot(all(rownames(projs_shared[[s1]]) == rownames(projs_shared[[s2]])))
    stopifnot(all(colnames(projs_shared[[s1]]) == colnames(projs_shared[[s2]])))
    data.frame(m1 = s1, m2 = s2, 
               dist = sqrt(mean((as.matrix(projs_shared[[s1]]) - 
                                   as.matrix(projs_shared[[s2]]))^2)),
               stringsAsFactors = FALSE)
  }))
}))

d_vel_shared <- d_vel_shared %>% dplyr::select(m1, m2, dist) %>% 
  tidyr::spread(key = m2, value = dist) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("m1") %>%
  as.matrix()
stopifnot(all(rownames(d_vel_shared) == colnames(d_vel_shared)))

d_projs_all <- d_projs_all %>% dplyr::select(m1, m2, dist) %>% 
  tidyr::spread(key = m2, value = dist) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("m1") %>%
  as.matrix()
stopifnot(all(rownames(d_projs_all) == colnames(d_projs_all)))

d_projs_shared <- d_projs_shared %>% dplyr::select(m1, m2, dist) %>% 
  tidyr::spread(key = m2, value = dist) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("m1") %>%
  as.matrix()
stopifnot(all(rownames(d_projs_shared) == colnames(d_projs_shared)))

cmd_vel_shared <- as.data.frame(cmdscale(d = d_vel_shared, k = 2)) %>%
  tibble::rownames_to_column("method") %>%
  dplyr::left_join(methods_short, by = "method")
cmd_projs_all <- as.data.frame(cmdscale(d = d_projs_all, k = 2)) %>%
  tibble::rownames_to_column("method") %>%
  dplyr::left_join(methods_short, by = "method")
cmd_projs_shared <- as.data.frame(cmdscale(d = d_projs_shared, k = 2)) %>%
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

pdf(gsub("\\.rds$", "_mds.pdf", outrds), width = 8, height = 8)
cowplot::plot_grid(
  cowplot::plot_grid(
    ggplot(cmd_projs_all) + 
      plset + ggtitle("MDS, projections, all genes") + 
      theme(legend.position = "none"),
    ggplot(cmd_projs_shared) + 
      plset + ggtitle("MDS, projections, shared genes") + 
      theme(legend.position = "none"),
    nrow = 1, labels = c("A", "B"), rel_widths = c(1, 1)
  ),
  cowplot::plot_grid(
    ggplot(cmd_vel_shared) + 
      plset + ggtitle("MDS, velocities, shared genes") + 
      theme(legend.position = "none"),
    cowplot::get_legend(ggplot(cmd_projs_all) + 
                          plset + theme(legend.position = "right")),
    nrow = 1, labels = c("C"), rel_widths = c(1, 1)
  ),
  ncol = 1, labels = "", rel_heights = c(1, 1)
)
dev.off()

## Plot distribution of velocity projection lengths between methods
tmp1 <- dplyr::bind_rows(
  do.call(dplyr::bind_rows, lapply(names(projs_all), function(nm) {
    data.frame(m1 = nm, m2 = nm,
               dtype = "All genes", 
               vtype = "Length of\nprojections",
               value = sqrt(rowSums(projs_all[[nm]] ^ 2)),
               stringsAsFactors = FALSE)
  })),
  do.call(dplyr::bind_rows, lapply(names(projs_shared), function(nm) {
    data.frame(m1 = nm, m2 = nm,
               dtype = "Shared genes", 
               vtype = "Length of\nprojections",
               value = sqrt(rowSums(projs_shared[[nm]] ^ 2)),
               stringsAsFactors = FALSE)
  }))
) %>% 
  dplyr::mutate(m1 = factor(methods_short$method_short[match(m1, methods_short$method)]),
                m2 = factor(methods_short$method_short[match(m2, methods_short$method)]))
  
## Plot distribution of distances
tmp2 <- dplyr::bind_rows(
  do.call(dplyr::bind_rows, lapply(names(projs_all), function(m1) {
    do.call(dplyr::bind_rows, lapply(names(projs_all), function(m2) {
      data.frame(m1 = m1, m2 = m2,
                 dtype = "All genes", 
                 vtype = "Distance\nbetween\nprojections",
                 value = sqrt(rowSums((projs_all[[m1]] - projs_all[[m2]]) ^ 2)),
                 stringsAsFactors = FALSE)
      
    }))
  })),
  do.call(dplyr::bind_rows, lapply(names(projs_shared), function(m1) {
    do.call(dplyr::bind_rows, lapply(names(projs_shared), function(m2) {
      data.frame(m1 = m1, m2 = m2,
                 dtype = "Shared genes", 
                 vtype = "Distance\nbetween\nprojections",
                 value = sqrt(rowSums((projs_shared[[m1]] - projs_shared[[m2]]) ^ 2)),
                 stringsAsFactors = FALSE)
      
    }))
  })),
) %>% dplyr::filter(m1 != m2) %>%
  dplyr::mutate(m1 = factor(methods_short$method_short[match(m1, methods_short$method)]),
                m2 = factor(methods_short$method_short[match(m2, methods_short$method)])) %>%
  dplyr::filter(as.numeric(m1) > as.numeric(m2))

## Plot distribution of distances in high-dim space
tmp3 <- do.call(dplyr::bind_rows, lapply(names(velocities_shared), function(m1) {
  do.call(dplyr::bind_rows, lapply(names(velocities_shared), function(m2) {
    data.frame(m1 = m1, m2 = m2,
               dtype = "Shared genes", 
               vtype = "Distance\nbetween\nvelocities/10",
               value = sqrt(colSums((velocities_shared[[m1]] - velocities_shared[[m2]]) ^ 2))/10,
               stringsAsFactors = FALSE)
    
  }))
})) %>% dplyr::filter(m1 != m2) %>%
  dplyr::mutate(m1 = factor(methods_short$method_short[match(m1, methods_short$method)]),
                m2 = factor(methods_short$method_short[match(m2, methods_short$method)])) %>%
  dplyr::filter(as.numeric(m1) < as.numeric(m2))

png(gsub("\\.rds", "_velocity_dists.png", outrds), height = 13, width = 13,
    unit = "in", res = 400)
tmp <- dplyr::bind_rows(tmp1, tmp2, tmp3)
ggplot(tmp, 
       aes(x = dtype, y = value, fill = vtype)) + 
  geom_boxplot(alpha = 0.5) + facet_grid(m1 ~ m2) + 
  theme_bw() + 
  labs(x = "", y = "Euclidean distance between/length of velocity vectors", title = dataset) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 5),
        legend.key.size = unit(3, 'lines')) + 
  scale_fill_manual(values = c("red", "blue", "green"), name = "")
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
