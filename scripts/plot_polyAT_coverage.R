args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(BiocGenerics)
  library(Biostrings)
  library(IRanges)
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
  library(rtracklayer)
})

print(genome)
print(gtf)
print(bwpos)
print(bwneg)
print(outrds)

## Read genome sequence
dss <- Biostrings::readDNAStringSet(genome)
names(dss) <- sapply(strsplit(names(dss), " "), .subset, 1)

## Find polyA/T stretches, at least 15nt, max 1 mismatch per 15 nt
polya <- IRanges::reduce(
  as(Biostrings::vmatchPattern("AAAAAAAAAAAAAAA", dss, max.mismatch = 1), "GRanges")
)
polyt <- IRanges::reduce(
  as(Biostrings::vmatchPattern("TTTTTTTTTTTTTTT", dss, max.mismatch = 1), "GRanges")
)

## Get unambiguously intronic regions and split by gene strand
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf")
grl <- GenomicFeatures::exonsBy(txdb, by = "gene")
grl <- GenomicRanges::reduce(grl)
gri <- BiocGenerics::setdiff(range(grl), grl)
introns <- BiocGenerics::unlist(gri)
exons <- BiocGenerics::unlist(grl)
introns <- BiocGenerics::setdiff(introns, exons) ## remove exons from other genes
introns_fwd <- subset(introns, strand == "+")
introns_rev <- subset(introns, strand == "-")

## Keep only polyA/T stretches fully in introns
## Note that the "polyA" is always extracted from the forward strand of the genome
## So a polyT here in a gene on the negative strand represents a polyA in the actual transcript
polyai_fwd <- IRanges::subsetByOverlaps(x = polya, ranges = introns_fwd, type = "within")
polyai_rev <- IRanges::subsetByOverlaps(x = polya, ranges = introns_rev, type = "within")
polyti_fwd <- IRanges::subsetByOverlaps(x = polyt, ranges = introns_fwd, type = "within")
polyti_rev <- IRanges::subsetByOverlaps(x = polyt, ranges = introns_rev, type = "within")

## polyA, genes on forward strand, reads on forward strand
## expand each polyA range to 600bp upstream and 600bp downstream
polyai_gfwd_rfwd <- colMeans(as.matrix(rtracklayer::import(
  bwpos, selection = rtracklayer::BigWigSelection(
    ranges = GenomicRanges::resize(
      GenomicRanges::resize(polyai_fwd, width = 601, fix = "end"),
      width = 1201, fix = "start"
    )
  ), as = "NumericList"
)))

## polyT, genes on forward strand, reads on reverse strand
## expand each polyT range to 600bp upstream and 600bp downstream
polyti_gfwd_rrev <- colMeans(as.matrix(rtracklayer::import(
  bwneg, selection = rtracklayer::BigWigSelection(
    ranges = GenomicRanges::resize(
      GenomicRanges::resize(polyti_fwd, width = 601, fix = "start"),
      width = 1201, fix = "end"
    )
  ), as = "NumericList"
)))

## polyT, genes on reverse strand, reads on reverse strand
## expand each polyT range to 600bp upstream and 600bp downstream
polyti_grev_rrev <- colMeans(as.matrix(rtracklayer::import(
  bwneg, selection = rtracklayer::BigWigSelection(
    ranges = GenomicRanges::resize(
      GenomicRanges::resize(polyti_rev, width = 601, fix = "start"),
      width = 1201, fix = "end"
    )
  ), as = "NumericList"
)))

## polyA, genes on reverse strand, reads on forward strand
## expand each polyA range to 600bp upstream and 600bp downstream
polyai_grev_rfwd <- colMeans(as.matrix(rtracklayer::import(
  bwpos, selection = rtracklayer::BigWigSelection(
    ranges = GenomicRanges::resize(
      GenomicRanges::resize(polyai_rev, width = 601, fix = "end"),
      width = 1201, fix = "start"
    )
  ), as = "NumericList"
)))



df <- dplyr::bind_rows(
  data.frame(gtype = "polyA, + genes", 
             position = (-600):600,
             coverage = polyai_gfwd_rfwd,
             ctype = "concordant",
             stringsAsFactors = FALSE),
  data.frame(gtype = "polyT, + genes", 
             position = (-600):600,
             coverage = polyti_gfwd_rrev,
             ctype = "discordant",
             stringsAsFactors = FALSE),
  data.frame(gtype = "polyT, - genes", 
             position = (-600):600,
             coverage = polyti_grev_rrev,
             ctype = "concordant",
             stringsAsFactors = FALSE),
  data.frame(gtype = "polyA, - genes", 
             position = (-600):600,
             coverage = polyai_grev_rfwd,
             ctype = "discordant",
             stringsAsFactors = FALSE)
)

pdf(gsub("rds", "pdf", outrds), height = 4, width = 6)
ggplot(df, aes(x = position, y = coverage, color = gtype)) + 
  geom_line(aes(size = ctype)) + 
  scale_size_manual(values = c(discordant = 1.5, concordant = 3),
                    name = "") + 
  scale_color_manual(values = c(`polyA, + genes` = "forestgreen",
                                `polyA, - genes` = "olivedrab3",
                                `polyT, - genes` = "darkred",
                                `polyT, + genes` = "orange"),
                     name = "") + 
  theme_bw() + 
  labs(x = "Position relative to end of polyA/T stretch",
       y = "Average coverage")
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()

