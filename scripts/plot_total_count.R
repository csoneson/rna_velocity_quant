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
  library(cowplot)
  library(pheatmap)
  library(ggrepel)
})
source(plothelperscript)

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(plothelperscript)
print(topdir)
print(tx2gene)
print(methods)
print(outrds)

sces <- lapply(methods, function(nm) {
  readRDS(file.path(topdir, paste0("output/sce/sce_", nm, ".rds")))
})

tx2gene <- readRDS(tx2gene)

methods_short <- shorten_methods(methods)

## ------------------------------------------------------------------------- ##
## Total number of assigned reads
## ------------------------------------------------------------------------- ##
sumdf_bygene <- do.call(dplyr::bind_rows, lapply(sces, function(w) {
  data.frame(gene = rownames(w),
             method = metadata(w)$method, 
             spliced = rowSums(assay(w, "spliced")),
             unspliced = rowSums(assay(w, "unspliced")),
             total = rowSums(assay(w, "unspliced")) + rowSums(assay(w, "spliced")),
             frac_unspliced = rowSums(assay(w, "unspliced"))/(rowSums(assay(w, "unspliced")) + rowSums(assay(w, "spliced"))),
             stringsAsFactors = FALSE
  )
}))

## Add info about uniqueness
uniq <- dplyr::bind_rows(
  read.delim(
    file.path(topdir, "reference/prepref_isoseparate_uniqueness.txt"),
    header = TRUE, as.is = TRUE
  ) %>% 
    dplyr::mutate(ctype = c("exonic", "intronic")[grepl("I\\.", gene) + 1]) %>%
    dplyr::mutate(gene = gsub("I\\.", "", gene)) %>%
    dplyr::mutate(frac_unique = unique/total) %>%
    dplyr::select(gene, ctype, frac_unique) %>%
    dplyr::mutate(atype = "separate"),
  read.delim(
    file.path(topdir, "reference/prepref_isocollapse_uniqueness.txt"),
    header = TRUE, as.is = TRUE
  ) %>% 
    dplyr::mutate(ctype = c("exonic", "intronic")[grepl("I\\.", gene) + 1]) %>%
    dplyr::mutate(gene = gsub("I\\.", "", gene)) %>%
    dplyr::mutate(frac_unique = unique/total) %>%
    dplyr::select(gene, ctype, frac_unique) %>%
    dplyr::mutate(atype = "collapse"),
  read.delim(
    file.path(topdir, "reference/prepref_isoseparate_uniqueness_overall.txt"),
    header = TRUE, as.is = TRUE
  ) %>% 
    dplyr::mutate(ctype = "overall") %>%
    dplyr::mutate(gene = gsub("I\\.", "", gene)) %>%
    dplyr::mutate(frac_unique = unique/total) %>%
    dplyr::select(gene, ctype, frac_unique) %>%
    dplyr::mutate(atype = "separate"),
  read.delim(
    file.path(topdir, "reference/prepref_isocollapse_uniqueness_overall.txt"),
    header = TRUE, as.is = TRUE
  ) %>% 
    dplyr::mutate(ctype = "overall") %>%
    dplyr::mutate(gene = gsub("I\\.", "", gene)) %>%
    dplyr::mutate(frac_unique = unique/total) %>%
    dplyr::select(gene, ctype, frac_unique) %>%
    dplyr::mutate(atype = "collapse")
) %>%
  dplyr::mutate(frac_unique_bin = Hmisc::cut2(frac_unique, 
                                              cuts = c(0, 0.001, 0.5, 0.999, 1))) %>% 
  dplyr::mutate(gene = tx2gene$gene_name[match(gene, tx2gene$gene_id)]) %>%
  dplyr::filter(gene %in% sumdf_bygene$gene) %>%
  dplyr::group_by(ctype, atype, frac_unique_bin) %>%
  dplyr::mutate(nbr_genes = length(gene)) %>% 
  dplyr::ungroup() %>%
  tidyr::unite("frac_unique_bin", frac_unique_bin, nbr_genes, sep = ", n = ")

## ------------------------------------------------------------------------- ##
## Plot uniqueness
## ------------------------------------------------------------------------- ##
pdf(gsub("\\.rds$", "_uniqueness.pdf", outrds), width = 8, height = 5)
ggplot(uniq %>% dplyr::select(-frac_unique_bin) %>% dplyr::group_by(atype) %>% 
         tidyr::spread(key = "ctype", value = "frac_unique"),
       aes(x = exonic, y = intronic)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point(alpha = 0.25) + 
  facet_wrap(~ atype) + 
  theme_bw() + 
  labs(x = "Fraction unique k-mers, exonic features",
       y = "Fraction unique k-mers, intronic features",
       title = "Fraction unique k-mers")

ggplot(uniq %>% dplyr::select(-frac_unique_bin) %>% dplyr::group_by(ctype) %>% 
         tidyr::spread(key = "atype", value = "frac_unique"),
       aes(x = separate, y = collapse)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point(alpha = 0.25) + 
  facet_wrap(~ ctype) + 
  theme_bw() + 
  labs(x = "Fraction unique k-mers, separate",
       y = "Fraction unique k-mers, collapse",
       title = "Fraction unique k-mers")

ggplot(uniq %>% dplyr::select(-frac_unique_bin), aes(x = frac_unique)) + 
  geom_histogram(bins = 100) + 
  facet_grid(atype ~ ctype) + 
  theme_bw() + 
  labs(x = "Fraction unique k-mers",
       y = "count",
       title = "Fraction unique k-mers")
dev.off()

## ------------------------------------------------------------------------- ##
## Plot counts
## ------------------------------------------------------------------------- ##
pdf(gsub("\\.rds$", "_total_count.pdf", outrds), width = 7, height = 12)
for (ct in c("exonic", "intronic", "overall")) {
  for (at in c("collapse", "separate")) {
    print(ggplot(dplyr::bind_rows(
      sumdf_bygene %>% 
        dplyr::left_join(uniq %>% dplyr::filter(ctype == ct & 
                                                  atype == at)) %>% 
        dplyr::mutate(frac_unique_bin = as.character(frac_unique_bin)) %>% 
        dplyr::group_by(method, frac_unique_bin) %>% 
        dplyr::summarize(spliced = sum(spliced), 
                         unspliced = sum(unspliced), 
                         total = sum(total)) %>% 
        tidyr::gather(key = "ctype", value = "count", spliced, unspliced, total) %>%
        dplyr::mutate(ctype = factor(ctype, 
                                     levels = c("total", "spliced", "unspliced"))) %>%
        dplyr::left_join(methods_short, by = "method"), 
      sumdf_bygene %>% group_by(method) %>% 
        dplyr::summarize(spliced = sum(spliced), 
                         unspliced = sum(unspliced), 
                         total = sum(total)) %>% 
        tidyr::gather(key = "ctype", value = "count", spliced, unspliced, total) %>%
        dplyr::mutate(frac_unique_bin = "overall") %>%
        dplyr::mutate(ctype = factor(ctype, 
                                     levels = c("total", "spliced", "unspliced"))) %>%
        dplyr::left_join(methods_short, by = "method")
    ) %>%
      dplyr::mutate(frac_unique_bin = relevel(factor(frac_unique_bin), ref = "overall")),
    aes(x = method_short, y = count, fill = mtype)) + 
      geom_bar(stat = "identity") + 
      facet_grid(frac_unique_bin ~ ctype, scale = "free_y") + 
      theme_bw() + 
      scale_fill_manual(values = base_method_colors, name = "") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "none") + 
      labs(x = "",
           y = "Total UMI count across genes in bin",
           title = paste0("Total count, stratified by ", ct, " uniqueness (", at, ")")))
  }
}
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
