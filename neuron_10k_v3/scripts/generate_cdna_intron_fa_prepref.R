args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(Biostrings)
  library(devtools)
  library(rtracklayer)
  library(dplyr)
})

devtools::load_all("/tungstenfs/groups/gbioinfo/sonechar/Projects/alevin_velocity/Rlib/prepref")

print(gtf)
print(genome)
print(isoform_action)
print(outdir)

## Extract intronic sequences flanked by L-1 bases 
## of exonic sequences where L is the biological read length of 
## the single cell technology of interest
genome <- Biostrings::readDNAStringSet(genome)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
gtfdf <- as.data.frame(rtracklayer::import(gtf))

tx <- extractTxSeqs(gtf = gtf, genome = genome, type = "spliced")
intr <- extractIntronSeqs(gtf = gtf, genome = genome, type = isoform_action, 
                          flanklength = 90)

t2gtx <- gtfdf %>% dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, gene_id)
if (isoform_action == "collapse") {
  t2gin <- data.frame(intr = names(intr),
                      gene = gsub("\\-I[0-9]*$", "", names(intr)),
                      stringsAsFactors = FALSE)
} else if (isoform_action == "separate") {
  t2gin <- data.frame(intr = names(intr),
                      transcript_id = gsub("\\-I[0-9]*$", "", names(intr)),
                      stringsAsFactors = FALSE) %>%
    dplyr::left_join(t2gtx, by = "transcript_id") %>%
    dplyr::select(intr, gene_id)
} else {
  stop("Don't know what to do with this isoform_action")
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
