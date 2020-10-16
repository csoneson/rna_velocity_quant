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
dataset <- gsub("_", " ", dataset)

print(plothelperscript)
print(velosuffix)
print(topdir)
print(dataset)
print(methods)
print(outrds)

methods_short <- shorten_methods(methods)

## Read velocities from scVelo
seurats <- lapply(methods, function(nm) {
  w <- ReadH5AD(file.path(topdir, paste0(
    "output/anndata_with_velocity", velosuffix, "/anndata_", 
    nm, "_with_velocity.h5ad")))
  w
})
velocities_shared <- lapply(seurats, function(w) as.matrix(GetAssayData(GetAssay(w, "velocity"))))
keep_genes <- Reduce(intersect, lapply(velocities_shared, function(w) rownames(w)[!is.na(rowSums(w))]))
velocities_shared <- lapply(velocities_shared, function(w) w[keep_genes, ])
cells <- colnames(velocities_shared[[1]])
genes <- rownames(velocities_shared[[1]])
stopifnot(all(sapply(velocities_shared, function(w) all(colnames(w) == cells))))
stopifnot(all(sapply(velocities_shared, function(w) all(rownames(w) == genes))))

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

## Reshape to distance matrix
d_vel_shared <- d_vel_shared %>% dplyr::select(m1, m2, dist) %>% 
  tidyr::spread(key = m2, value = dist) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("m1") %>%
  as.matrix()
stopifnot(all(rownames(d_vel_shared) == colnames(d_vel_shared)))

## Apply classical MDS
cmd_vel_shared <- as.data.frame(cmdscale(d = d_vel_shared, k = 2)) %>%
  tibble::rownames_to_column("method") %>%
  dplyr::left_join(methods_short, by = "method")

## Define plot settings
plset <- list(
  aes(x = V1, y = V2, shape = rtype, label = method_short),
  geom_point(aes(color = mtype), size = 7, alpha = 0.7),
  geom_label_repel(size = 3),
  theme_minimal(),
  scale_color_manual(values = base_method_colors, name = ""),
  scale_shape_discrete(name = ""),
  labs(x = "MDS1", y = "MDS2")
)

pdf(gsub("\\.rds$", "_mds.pdf", outrds), width = 7, height = 0.85 * 7)
print(ggplot(cmd_vel_shared) + 
        plset + ggtitle(paste0(dataset, ", MDS, velocities, shared genes")) + 
        theme(legend.position = "bottom") + 
        guides(color = guide_legend(nrow = 2, byrow = TRUE)))
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
