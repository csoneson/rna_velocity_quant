args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

## Use SingleR to assign cell types to the neuron data set
suppressPackageStartupMessages({
  library(SingleR)
  library(SingleCellExperiment)
  library(dplyr)
  library(scater)
})

print(topdir)
print(tx2gene)
print(helperscript)
print(outcsv)

source(helperscript)
tx2gene <- readRDS(tx2gene)

sce <- read_alevin_spliced(
  alevindir = file.path(topdir, "quants/alevin_spliced/alevin"),
  sampleid = "neuron", tx2gene = tx2gene
)
sce <- scater::logNormCounts(sce)

ref <- SingleR::MouseRNAseqData()
pred <- SingleR::SingleR(test = sce, ref = ref, labels = ref$label.fine, prune = FALSE)

df <- data.frame(index = rownames(pred), 
                 clusters = pred$labels,
                 clusters_coarse = pred$labels,
                 stringsAsFactors = FALSE)
write.table(df, file = outcsv, row.names = FALSE, 
            col.names = TRUE, sep = ",", quote = FALSE)

date()
sessionInfo()
