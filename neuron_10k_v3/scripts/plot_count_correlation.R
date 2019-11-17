args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(tximeta)
  library(dplyr)
  library(ggplot2)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(scater)
})

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(topdir)
print(methods)
print(genetxt)
print(outrds)

sces <- lapply(methods, function(nm) {
  readRDS(file.path(topdir, paste0("output/sce/sce_", nm, ".rds")))
})
for (nm in names(sces)) {
  metadata(sces[[nm]]) <- list(dataset = nm)
}

genes <- read.delim(genetxt, header = FALSE, as.is = TRUE)[, 1]
sces <- lapply(sces, function(w) w[match(genes, rownames(w)), ])

corrs <- do.call(
  dplyr::bind_rows, 
  lapply(seq_len(length(methods) - 1), function(j) {
    jj <- methods[j]
    do.call(
      dplyr::bind_rows,  
      lapply((j + 1):(length(methods)), function(k) {
        kk <- methods[k]
        message(jj, " - ", kk)
        a <- scale(t(as.matrix(assay(sces[[jj]], "spliced"))), center = TRUE, scale = TRUE)
        b <- scale(t(as.matrix(assay(sces[[kk]], "spliced"))), center = TRUE, scale = TRUE)
        corrs_spliced <- colSums(a * b)/(nrow(a) - 1)
        a <- scale(t(as.matrix(assay(sces[[jj]], "unspliced"))), center = TRUE, scale = TRUE)
        b <- scale(t(as.matrix(assay(sces[[kk]], "unspliced"))), center = TRUE, scale = TRUE)
        corrs_unspliced <- colSums(a * b)/(nrow(a) - 1)
        dplyr::bind_rows(
          data.frame(method1 = jj,
                     method2 = kk, 
                     gene = names(corrs_spliced), 
                     ctype = "spliced",
                     corrs = corrs_spliced, 
                     stringsAsFactors = FALSE),
          data.frame(method1 = jj,
                     method2 = kk,
                     gene = names(corrs_unspliced), 
                     ctype = "unspliced",
                     corrs = corrs_unspliced, 
                     stringsAsFactors = FALSE)
        )
      }))
  }))

corrs$method1 <- factor(corrs$method1, levels = methods[methods %in% corrs$method1])
corrs$method2 <- factor(corrs$method2, levels = methods[methods %in% corrs$method2])

pdf(gsub("rds$", "pdf", outrds), width = 30, height = 30)

ggplot(corrs, aes(x = ctype, y = corrs, fill = ctype)) + 
  geom_boxplot() + 
  facet_grid(method1 ~ method2) + 
  theme_bw() + 
  theme(strip.text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        title = element_text(size = 20)) + 
  labs(title = "Correlation, spliced and unspliced counts",
       subtitle = "Over genes selected by all methods")

ggplot(corrs, aes(x = paste(method1, method2, sep = "-"), y = corrs, fill = ctype)) + 
  geom_boxplot() + 
  theme_bw() + 
  labs(title = "Correlation, spliced and unspliced counts",
       subtitle = "Over genes selected by all methods", 
       x = "", y = "Correlation") + 
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5))

dev.off()

saveRDS(NULL, file = outrds)

date()
sessionInfo()
