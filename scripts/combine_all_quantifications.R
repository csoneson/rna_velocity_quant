args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods
print(scedir)
print(methods)
print(outfile)

sces <- lapply(methods, function(m) {
    readRDS(file.path(scedir, paste0("sce_", m, ".rds")))
})

## Check that rownames and column names are identical, as well as the dataset
rn <- rownames(sces[[1]])
cn <- colnames(sces[[1]])
ds <- metadata(sces[[1]])$dataset
stopifnot(all(vapply(sces, function(s) {
    all(rownames(s) == rn) && all(colnames(s) == cn) && 
        metadata(s)$dataset == ds
}, FALSE)))

## colData columns to retain
cdcols <- setdiff(grep("cluster_", colnames(colData(sces[[1]])), value = TRUE, invert = TRUE), "cluster")
print(cdcols)

## reduced dimension representation to retain
reddims <- c("PCA_alevin_spliced_gentrome", "TSNE_alevin_spliced_gentrome",
             "UMAP_alevin_spliced_gentrome")

## Build the object
sce <- SingleCellExperiment(
    assays = list(),
    colData = colData(sces[[1]])[, cdcols],
    rowData = rowData(sces[[1]]),
    metadata = list(dataset = ds),
    reducedDims = reducedDims(sces[[1]])[reddims]
)
for (m in methods) {
    if (all(c("spliced", "unspliced") %in% assayNames(sces[[m]]))) {
        assays(sce)[[paste0(m, "_spliced")]] <- assays(sces[[m]])[["spliced"]]
        assays(sce)[[paste0(m, "_unspliced")]] <- assays(sces[[m]])[["unspliced"]]
    } else if ("counts" %in% assayNames(sces[[m]])) {
        assays(sce)[[paste0(m, "_counts")]] <- assays(sces[[m]])[["counts"]]
    } else {
        print(paste0("No valid assays for ", m))
    }
}
print(assayNames(sce))
print(sce)

saveRDS(sce, file = outfile)

date()
sessionInfo()
