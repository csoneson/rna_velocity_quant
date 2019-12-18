args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(scater)
  library(readr)
})

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(topdir)
print(methods)
print(outtxt)

geneinfo <- lapply(methods, function(nm) {
  readr::read_csv(paste0(topdir, "/plots/velocity/anndata_", nm, "/anndata_", nm, 
                         "_gene_info.csv")) %>%
    dplyr::mutate(method = nm)
})

## Extract all genes that are selected with all methods, and have valid fits with all methods
allgenes <- as.character(Reduce(union, lapply(geneinfo, function(w) w$index[!is.na(w$fit_alpha)])))
selgenes <- do.call(cbind, lapply(geneinfo, function(w) as.numeric(allgenes %in% w$index[!is.na(w$fit_alpha)])))
rownames(selgenes) <- allgenes
selgenes <- data.frame(selgenes)

n_methods <- ncol(selgenes)
genes_in_all <- selgenes %>% 
  tibble::rownames_to_column("gene") %>%
  tidyr::gather(key = "method", value = "selected", -gene) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(n = sum(selected)) %>% 
  dplyr::filter(n == n_methods) %>%
  dplyr::pull(gene)

write.table(genes_in_all, file = outtxt, row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = "\t")

date()
sessionInfo()
