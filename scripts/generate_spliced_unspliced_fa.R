args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(Biostrings)
  library(rtracklayer)
  library(dplyr)
  library(GenomicFeatures)
  library(BSgenome)
})

print(scriptdir)
print(gtf)
print(genome)
print(tx2gene)
print(outfa)
print(outt2g)

source(file.path(scriptdir, "extractTxSeqs.R"))

genome <- Biostrings::readDNAStringSet(genome)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)

## Extract transcript and pre-mRNA sequences
tx <- extractTxSeqs(gtf = gtf, genome = genome, type = "spliced")
premrna <- extractTxSeqs(gtf = gtf, genome = genome, type = "unspliced")
names(premrna) <- paste0(names(premrna), "_unspliced")

## Generate tx2gene
t2g <- readRDS(tx2gene)
t2gpre <- t2g %>% dplyr::mutate(transcript_id = paste0(transcript_id, "_unspliced"),
                                gene_id = paste0(gene_id, "_unspliced"),
                                gene_name = paste0(gene_name, "_unspliced"))
t2g <- rbind(t2g, t2gpre)

## Combine
mrna <- c(tx, premrna)

## Save
write.table(t2g %>% dplyr::select(transcript_id, gene_id), 
            file = outt2g, row.names = FALSE, col.names = FALSE, 
            sep = "\t", quote = FALSE)
Biostrings::writeXStringSet(mrna, file = outfa)

date()
sessionInfo()

