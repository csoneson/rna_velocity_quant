args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

## Various diagnostics

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(SingleCellExperiment)
})

print(solodir)
print(outrds)

## ------------------------------------------------------------------------- ##
## Number and fraction of genes for which the STAR/GeneFull count 
## is lower than the STAR/Gene count
## Consider total count across all cells
## ------------------------------------------------------------------------- ##
## Gene body
fgenes <- readr::read_delim(file.path(solodir, "GeneFull/raw/features.tsv"), col_names = FALSE, delim = "\t")
fcells <- readr::read_delim(file.path(solodir, "GeneFull/raw/barcodes.tsv"), col_names = FALSE, delim = "\t")
fn_genes <- nrow(fgenes)
fn_cells <- nrow(fcells)
## Read matrix. col 1 = gene, col 2 = cell, col 3 = spliced, col 4 = unspliced, col 5 = ambiguous
fmatr <- readr::read_delim(file.path(solodir, "GeneFull/raw/matrix.mtx"), skip = 3, delim = " ", col_names = FALSE)
fcounts <- Matrix::sparseMatrix(
  i = fmatr$X1,
  j = fmatr$X2,
  x = fmatr$X3,
  dims = c(fn_genes, fn_cells)
)
rownames(fcounts) <- fgenes$X2
colnames(fcounts) <- fcells$X1

## Exons
sgenes <- readr::read_delim(file.path(solodir, "Gene/raw/features.tsv"), col_names = FALSE, delim = "\t")
scells <- readr::read_delim(file.path(solodir, "Gene/raw/barcodes.tsv"), col_names = FALSE, delim = "\t")
sn_genes <- nrow(sgenes)
sn_cells <- nrow(scells)
## Read matrix. col 1 = gene, col 2 = cell, col 3 = spliced, col 4 = unspliced, col 5 = ambiguous
smatr <- readr::read_delim(file.path(solodir, "Gene/raw/matrix.mtx"), skip = 3, delim = " ", col_names = FALSE)
scounts <- Matrix::sparseMatrix(
  i = smatr$X1,
  j = smatr$X2,
  x = smatr$X3,
  dims = c(sn_genes, sn_cells)
)
rownames(scounts) <- sgenes$X2
colnames(scounts) <- scells$X1

stopifnot(all(rownames(scounts) == rownames(fcounts)))
stopifnot(all(colnames(scounts) == colnames(fcounts)))

## Subset to only cells and genes used for the analysis
sce <- readRDS(file.path("output/sce/sce_starsolo_subtr.rds"))

scounts <- scounts[rownames(sce), colnames(sce)]
fcounts <- fcounts[rownames(sce), colnames(sce)]

ftot <- rowSums(fcounts)
stot <- rowSums(scounts)

print(sum(ftot < stot))
print(mean(ftot < stot))

stopifnot(all(names(ftot) == names(stot)))
print(data.frame(gene = names(ftot),
                 ftot = ftot, 
                 stot = stot,
                 stringsAsFactors = FALSE) %>%
        dplyr::mutate(stot_minus_ftot = stot - ftot) %>%
        dplyr::arrange(desc(stot_minus_ftot)) %>%
        head(10)) 

saveRDS(NULL, file = outrds)
date()
sessionInfo()
