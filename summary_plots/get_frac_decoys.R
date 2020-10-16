args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggbeeswarm)
  library(tibble)
})

topdir <- ".."
datasets <- c("spermatogenesis_mouse", "dentate_gyrus_mouse", 
              "pancreas_mouse", "pfc_mouse", "oldbrain_mouse")

dirs <- do.call(c, lapply(datasets, function(ds) {
  list.files(file.path(topdir, ds, "quants"), pattern = "gentrome", 
             full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
}))

dirs <- dirs[grep("asdecoy_gentrome", dirs, invert = TRUE)]
files <- paste0(dirs, "/logs/salmon_quant.log")
frac_decoys <- sapply(files, function(f) {
  x <- readLines(f)
  ndecoys <- as.numeric(gsub(
    ",", "", strsplit(
      grep("Number of fragments discarded because they are best-mapped to decoys : ", 
           x, value = TRUE), 
      "decoys : ")[[1]][2]
  ))
  neqreads <- as.numeric(gsub(
    ",", "", strsplit(strsplit(grep("total reads in the equivalence classes", 
                                    x, value = TRUE), " Counted ")[[1]][2], 
                      " total reads in the equivalence classes")[[1]][1]
    ))
  frac_discarded <- ndecoys/(ndecoys + neqreads)
  frac_discarded
})
print(summary(frac_decoys))

df <- as.data.frame(frac_decoys) %>% 
  tibble::rownames_to_column("run")

pdf("fraction_decoys.pdf", height = 8, width = 5)
ggplot(df %>% dplyr::mutate(method = sapply(run, function(w) {
  a <- strsplit(w, "/")[[1]]; a[grep("alevin", a)]
})),
aes(x = method, y = frac_decoys)) + ggbeeswarm::geom_beeswarm(size = 2) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(title = "Fraction of the reads discarded\nfor mapping better to decoys", 
       x = "", y = "")
dev.off()

write.table(df,
            file = "fraction_decoys.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

date()
sessionInfo()
