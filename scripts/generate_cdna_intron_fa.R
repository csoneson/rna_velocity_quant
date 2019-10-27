args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(Biostrings)
  library(BUSpaRse)
})

print(gtf)
print(genome)
print(isoform_action)
print(outdir)

## This function extracts intronic sequences flanked by L-1 bases 
## of exonic sequences where L is the biological read length of 
## the single cell technology of interest
genome <- Biostrings::readDNAStringSet(genome)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
BUSpaRse::get_velocity_files(
  X = gtf, L = 91, Genome = genome, Transcriptome = NULL, 
  out_path = outdir, isoform_action = isoform_action, 
  exon_option = "full", compress_fa = TRUE, 
  transcript_id = "transcript_id", gene_id = "gene_id", 
  transcript_version = "transcript_version", 
  gene_version = "gene_version", version_sep = "."
)

date()
sessionInfo()
