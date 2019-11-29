args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(tidyr)
  library(GenomicFeatures)
  library(ggplot2)
  library(cowplot)
})

print(txomefa)
print(gtffile)
print(refdir)
print(flanklength)
print(outrds)

## ----------------------------------------------------------------------------
## Read reference files
## ----------------------------------------------------------------------------
## Transcriptome
txome <- Biostrings::readDNAStringSet(txomefa)
names(txome) <- sapply(strsplit(names(txome), "\\|"), .subset, 1)

## gtf file
gtf <- rtracklayer::import(gtffile)

## get the strand and the number of exons for each transcript
gtftx <- as.data.frame(gtf) %>% 
  dplyr::filter(type == "exon") %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::mutate(strand = as.character(strand)) %>% 
  dplyr::summarize(nbr_exons = length(type), 
                   strand = unique(strand))

## get the number of transcripts for each gene
gtfg <- as.data.frame(gtf) %>%
  dplyr::filter(type == "transcript") %>% 
  dplyr::group_by(gene_id) %>%
  dplyr::summarize(nbr_tx = length(type))

## transcript-to-gene mapping for genes with a single isoform
gtfg1 <- as.data.frame(subset(gtf, gene_id %in% gtfg$gene_id[gtfg$nbr_tx == 1] & type == "transcript")) %>%
  dplyr::select(gene_id, transcript_id)

## construct txdb and group exons and introns by transcript/gene
txdb <- GenomicFeatures::makeTxDbFromGFF(gtffile)
ebt <- GenomicFeatures::exonsBy(txdb, "tx", use.names = TRUE)
ebg <- GenomicFeatures::exonsBy(txdb, "gene")
ibt <- GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)

## ----------------------------------------------------------------------------
## Read extracted transcript and intron sequences
## ----------------------------------------------------------------------------
busparse_separate <- Biostrings::readDNAStringSet(
  file.path(refdir, "busparse_isoseparate/cDNA_introns.fa.gz")
)
names(busparse_separate) <- gsub("\\.-I", "-I", gsub("\\.$", "", names(busparse_separate)))
prepref_separate <- Biostrings::readDNAStringSet(
  file.path(refdir, "prepref_isoseparate/cDNA_introns.fa.gz")
)
busparse_collapse <- Biostrings::readDNAStringSet(
  file.path(refdir, "busparse_isocollapse/cDNA_introns.fa.gz")
)
names(busparse_collapse) <- gsub("\\.-I", "-I", gsub("\\.$", "", names(busparse_collapse)))
prepref_collapse <- Biostrings::readDNAStringSet(
  file.path(refdir, "prepref_isocollapse/cDNA_introns.fa.gz")
)

## ----------------------------------------------------------------------------
## Compare extracted transcripts to annotated ones
## ----------------------------------------------------------------------------
gtftx <- gtftx[match(names(txome), gtftx$transcript_id), ]
busparse_separate_txome <- busparse_separate[match(names(txome), names(busparse_separate))]
prepref_separate_txome <- prepref_separate[match(names(txome), names(prepref_separate))]
busparse_collapse_txome <- busparse_collapse[match(names(txome), names(busparse_collapse))]
prepref_collapse_txome <- prepref_collapse[match(names(txome), names(prepref_collapse))]

busparse_separate_tb <- 
  table(correct_seq = factor(busparse_separate_txome == txome, levels = c(FALSE, TRUE)), 
        single_exon = gtftx$nbr_exons == 1, strand = gtftx$strand)
prepref_separate_tb <- 
  table(correct_seq = factor(prepref_separate_txome == txome, levels = c(FALSE, TRUE)),
        single_exon = gtftx$nbr_exons == 1, strand = gtftx$strand)
busparse_collapse_tb <- 
  table(correct_seq = factor(busparse_collapse_txome == txome, levels = c(FALSE, TRUE)), 
        single_exon = gtftx$nbr_exons == 1, strand = gtftx$strand)
prepref_collapse_tb <- 
  table(correct_seq = factor(prepref_collapse_txome == txome, levels = c(FALSE, TRUE)),
        single_exon = gtftx$nbr_exons == 1, strand = gtftx$strand)

dftx <- dplyr::bind_rows(
  as.data.frame(busparse_separate_tb) %>% dplyr::mutate(type = "separate", method = "BUSpaRse"),
  as.data.frame(prepref_separate_tb) %>% dplyr::mutate(type = "separate", method = "prepref"),
  as.data.frame(busparse_collapse_tb) %>% dplyr::mutate(type = "collapse", method = "BUSpaRse"),
  as.data.frame(prepref_collapse_tb) %>% dplyr::mutate(type = "collapse", method = "prepref")
)

g1 <- ggplot(
  dftx %>% dplyr::mutate(single_exon = ifelse(single_exon == TRUE, "Single exon", "Multiple exons")) %>%
    dplyr::mutate(strand = ifelse(strand == "+", "Positive strand", "Negative strand")) %>%
    dplyr::mutate(correct_seq = ifelse(correct_seq == TRUE, "Correct", "Incorrect")), 
  aes(x = method, y = Freq, fill = correct_seq)) + geom_bar(stat = "identity") + 
  facet_grid(type ~ paste(single_exon, strand, sep = ", ")) + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  labs(y = "Number of transcripts", title = "Transcripts") + 
  scale_fill_manual(values = c(Incorrect = "pink", Correct = "lightblue"), name = "Inferred sequence")

## ----------------------------------------------------------------------------
## Compare lengths of extracted introns to annotated ones 
## (all transcripts for separate, genes with single transcripts for collapse)
## ----------------------------------------------------------------------------
## Separate
## ----------------------------------------------------------------------------
intrs_separate <- grep("-I", names(busparse_separate), value = TRUE)
stopifnot(all(intrs_separate %in% names(prepref_separate)))
busparse_separate_intrs <- busparse_separate[match(intrs_separate, names(busparse_separate))]
prepref_separate_intrs <- prepref_separate[match(intrs_separate, names(prepref_separate))]

## Example of problematic transcript
ebt$ENSMUST00000195885.1
busparse_separate_intrs[grep("ENSMUST00000195885.1", names(busparse_separate_intrs))]
prepref_separate_intrs[grep("ENSMUST00000195885.1", names(prepref_separate_intrs))]

## Get introns per transcript and extract their lengths
ann_intr <- sapply(ibt, function(w) paste(sort(width(w) + 2 * flanklength), collapse = ","))
busparse_separate_widths <- sapply(
  split(width(busparse_separate_intrs), 
        f = sapply(strsplit(names(busparse_separate_intrs), "-"), .subset, 1)), 
  function(w) paste(sort(w), collapse = ","))
prepref_separate_widths <- sapply(
  split(width(prepref_separate_intrs), 
        f = sapply(strsplit(names(prepref_separate_intrs), "-"), .subset, 1)), 
  function(w) paste(sort(w), collapse = ","))
busparse_separate_widths <- busparse_separate_widths[match(names(ann_intr), names(busparse_separate_widths))]
prepref_separate_widths <- prepref_separate_widths[match(names(ann_intr), names(prepref_separate_widths))]
gtfintr <- gtftx[match(names(ann_intr), gtftx$transcript_id), ]

busparse_separate_tb <- 
  table(correct_length = factor(busparse_separate_widths == ann_intr, levels = c(FALSE, TRUE)), 
        single_exon = gtfintr$nbr_exons == 1, strand = gtfintr$strand)
prepref_separate_tb <- 
  table(correct_length = factor(prepref_separate_widths == ann_intr, levels = c(FALSE, TRUE)),
        single_exon = gtfintr$nbr_exons == 1, strand = gtfintr$strand)

## Collapse
## ----------------------------------------------------------------------------
intrs_collapse <- grep("-I", names(busparse_collapse), value = TRUE)
stopifnot(all(intrs_collapse %in% names(prepref_collapse)))
busparse_collapse_intrs <- busparse_collapse[match(intrs_collapse, names(busparse_collapse))]
prepref_collapse_intrs <- prepref_collapse[match(intrs_collapse, names(prepref_collapse))]

## Example of problematic gene
ebg$ENSMUSG00000103757.1
busparse_collapse_intrs[grep("ENSMUSG00000103757.1", names(busparse_collapse_intrs))]
prepref_collapse_intrs[grep("ENSMUSG00000103757.1", names(prepref_collapse_intrs))]

## Get introns per gene and extract their lengths
busparse_collapse_widths <- sapply(
  split(width(busparse_collapse_intrs), 
        f = sapply(strsplit(names(busparse_collapse_intrs), "-"), .subset, 1)), 
  function(w) paste(sort(w), collapse = ","))
prepref_collapse_widths <- sapply(
  split(width(prepref_collapse_intrs), 
        f = sapply(strsplit(names(prepref_collapse_intrs), "-"), .subset, 1)), 
  function(w) paste(sort(w), collapse = ","))

## Transcripts from genes with a single transcript
ints <- intersect(gtfg1$transcript_id, names(ann_intr))
ann_intr1 <- ann_intr[match(ints, names(ann_intr))]
single_tx_genes <- gtfg1$gene_id[match(names(ann_intr1), gtfg1$transcript_id)]
busparse_collapse_widths <- busparse_collapse_widths[match(single_tx_genes, names(busparse_collapse_widths))]
busparse_collapse_widths[is.na(busparse_collapse_widths)] <- ""
prepref_collapse_widths <- prepref_collapse_widths[match(single_tx_genes, names(prepref_collapse_widths))]
prepref_collapse_widths[is.na(prepref_collapse_widths)] <- ""

gtfintr <- gtftx[match(names(ann_intr1), gtftx$transcript_id), ]

busparse_collapse_tb <- 
  table(correct_length = factor(busparse_collapse_widths == ann_intr1, levels = c(FALSE, TRUE)), 
        single_exon = gtfintr$nbr_exons == 1, strand = gtfintr$strand)
prepref_collapse_tb <- 
  table(correct_length = factor(prepref_collapse_widths == ann_intr1, levels = c(FALSE, TRUE)),
        single_exon = gtfintr$nbr_exons == 1, strand = gtfintr$strand)


dfintr <- dplyr::bind_rows(
  as.data.frame(busparse_separate_tb) %>% dplyr::mutate(type = "separate", method = "BUSpaRse"),
  as.data.frame(prepref_separate_tb) %>% dplyr::mutate(type = "separate", method = "prepref"),
  as.data.frame(busparse_collapse_tb) %>% dplyr::mutate(type = "collapse", method = "BUSpaRse"),
  as.data.frame(prepref_collapse_tb) %>% dplyr::mutate(type = "collapse", method = "prepref")
)

g2 <- ggplot(dfintr %>% dplyr::filter(type == "separate") %>% 
         dplyr::mutate(single_exon = ifelse(single_exon == TRUE, "Single exon", "Multiple exons")) %>%
         dplyr::mutate(strand = ifelse(strand == "+", "Positive strand", "Negative strand")) %>%
         dplyr::mutate(correct_length = ifelse(correct_length == TRUE, "Correct", "Incorrect")) %>%
         dplyr::filter(single_exon == "Multiple exons"), 
       aes(x = method, y = Freq, fill = correct_length)) + geom_bar(stat = "identity") + 
  facet_grid(type ~ paste(single_exon, strand, sep = ", ")) + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  labs(y = "Number of transcripts", title = "Introns") + 
  scale_fill_manual(values = c(Incorrect = "pink", Correct = "lightblue"), name = "Inferred intron lengths")
g3 <- ggplot(dfintr %>% dplyr::filter(type == "collapse") %>% 
               dplyr::mutate(single_exon = ifelse(single_exon == TRUE, "Single exon", "Multiple exons")) %>%
               dplyr::mutate(strand = ifelse(strand == "+", "Positive strand", "Negative strand")) %>%
               dplyr::mutate(correct_length = ifelse(correct_length == TRUE, "Correct", "Incorrect")) %>%
               dplyr::filter(single_exon == "Multiple exons"), 
             aes(x = method, y = Freq, fill = correct_length)) + geom_bar(stat = "identity") + 
  facet_grid(type ~ paste(single_exon, strand, sep = ", ")) + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  labs(y = "Number of genes", title = "Introns") + 
  scale_fill_manual(values = c(Incorrect = "pink", Correct = "lightblue"), name = "Inferred intron lengths")

pdf(gsub("rds$", "pdf", outrds), width = 11, height = 9)
cowplot::plot_grid(g1, cowplot::plot_grid(g2, g3, nrow = 1, labels = c("B", "C")), ncol = 1, labels = c("A", ""), rel_heights = c(0.8, 0.55))
dev.off()

saveRDS(list(dftx = dftx, dfintr = dfintr), file = outrds)

date()
sessionInfo()
