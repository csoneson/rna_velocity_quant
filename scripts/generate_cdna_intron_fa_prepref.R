args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(Biostrings)
  library(rtracklayer)
  library(dplyr)
  library(GenomicFeatures)
  library(BiocGenerics)
  library(BSgenome)
  library(GenomicRanges)
})

print(scriptdir)
print(gtf)
print(genome)
print(isoform_action)
print(flanklength)
print(outdir)

source(file.path(scriptdir, "extractIntronSeqs.R"))
source(file.path(scriptdir, "extractTxSeqs.R"))

## Extract intronic sequences flanked by L-1 bases 
## of exonic sequences where L is the biological read length
genome <- Biostrings::readDNAStringSet(genome)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
gtfdf <- as.data.frame(rtracklayer::import(gtf))

## Extract transcript and intron sequences
tx <- extractTxSeqs(gtf = gtf, genome = genome, type = "spliced")
intr <- extractIntronSeqs(gtf = gtf, genome = genome, type = isoform_action, 
                          flanklength = flanklength,
                          joinOverlappingIntrons = FALSE)

## Generate transcript/intron-to-gene mapping
t2gtx <- gtfdf %>% dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, gene_id) %>%
  dplyr::distinct()
if (isoform_action == "collapse") {
  ## Intron names already contain gene name
  t2gin <- data.frame(intr = names(intr),
                      gene = gsub("\\-I[0-9]*$", "", names(intr)),
                      stringsAsFactors = FALSE)
} else if (isoform_action == "separate") {
  ## Intron names contain transcript name
  t2gin <- data.frame(intr = names(intr),
                      transcript_id = gsub("\\-I[0-9]*$", "", names(intr)),
                      stringsAsFactors = FALSE) %>%
    dplyr::left_join(t2gtx, by = "transcript_id") %>%
    dplyr::select(intr, gene_id)
} else {
  stop("Unknown isoform_action")
}
colnames(t2gin) <- colnames(t2gtx)
t2g <- rbind(t2gtx, t2gin)

Biostrings::writeXStringSet(c(tx, intr), file.path(outdir, "cDNA_introns.fa.gz"), 
                            compress = TRUE)
write.table(names(tx), file = file.path(outdir, "cDNA_tx_to_capture.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(names(intr), file = file.path(outdir, "introns_tx_to_capture.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(t2g, file = file.path(outdir, "tr2g.tsv"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

date()
sessionInfo()
