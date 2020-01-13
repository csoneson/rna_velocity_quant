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
  library(GGally)
})
source(plothelperscript)

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(plothelperscript)
print(topdir)
print(refdir)  ## directory where uniqueness files are
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
uniq <- merge_uniq(refdir = refdir, tx2gene = tx2gene,
                   keepgenes = sumdf_bygene$gene)

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

## ------------------------------------------------------------------------- ##
## Plot scatter plot of total counts
## ------------------------------------------------------------------------- ##
lowerfun <- function(data, mapping) {
  ggplot(data = data, mapping = mapping) +
    geom_point(size = 0.3, alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "blue")
}  

for (i in c("spliced", "unspliced")) {
  png(gsub("\\.rds$", paste0("_total_count_scatter_", i, ".png"), outrds), 
      width = 11, height = 11, unit = "in", res = 400)
  print(GGally::ggpairs(
    sumdf_bygene %>% dplyr::left_join(methods_short, by = "method") %>%
      dplyr::select(c("gene", "method_short", i)) %>%
      tidyr::spread(key = "method_short", value = i) %>% 
      tibble::column_to_rownames("gene") %>%
      dplyr::mutate_all(.funs = log1p),
    lower = list(continuous = wrap(lowerfun)),
    title = i, xlab = "log(total gene count + 1)", 
    ylab = "log(total gene count + 1)") + 
      theme(strip.text = element_text(size = 5))
  )
  dev.off()
}

## Check number of genes where the count is different between methods
res <- do.call(dplyr::bind_rows, lapply(c("spliced", "unspliced"), function(tp) {
  tmp <- sumdf_bygene %>% dplyr::left_join(methods_short, by = "method") %>%
    dplyr::select(c("gene", "method_short", tp)) %>%
    tidyr::spread(key = "method_short", value = tp) %>% 
    tibble::column_to_rownames("gene")
  do.call(dplyr::bind_rows, lapply(seq_len(ncol(tmp)), function(i) {
    do.call(dplyr::bind_rows, lapply(seq_len(i - 1), function(j) {
      data.frame(
        m1 = ifelse(tp == "spliced", colnames(tmp)[i], colnames(tmp[j])),
        m2 = ifelse(tp == "spliced", colnames(tmp)[j], colnames(tmp[i])),
        ctype = tp,
        tot_genes = nrow(tmp),
        nbr_equal = sum(tmp[, i] == tmp[, j]),
        nbr_equal_not0 = sum(tmp[, i] == tmp[, j] & tmp[, i] != 0),
        nbr_notequal = sum(tmp[, i] != tmp[, j]),
        nbr_morethan5pdiff = sum(abs(tmp[, i] - tmp[, j]) > 
                                   0.05 * ((tmp[, i] + tmp[, j])/2)),
        nbr_morethan10pdiff = sum(abs(tmp[, i] - tmp[, j]) > 
                                    0.1 * ((tmp[, i] + tmp[, j])/2)),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::mutate(
          frac_equal = nbr_equal/tot_genes,
          frac_equal_not0 = nbr_equal_not0/tot_genes,
          frac_notequal = nbr_notequal/tot_genes,
          frac_morethan5pdiff = nbr_morethan5pdiff/tot_genes,
          frac_morethan10pdiff = nbr_morethan10pdiff/tot_genes
        )
    }))
  }))
}))

pdf(gsub("\\.rds$", "_frac_diff_btw_methods.pdf", outrds), width = 13, height = 13)
print(res %>% dplyr::select(c("m1", "m2", "ctype", contains("frac"))) %>%
        tidyr::gather(key = "ftype", value = "fraction", -m1, -m2, -ctype) %>%
        dplyr::mutate(ftype = replace(ftype, ftype == "frac_equal", "Equal"),
                      ftype = replace(ftype, ftype == "frac_equal_not0", "Equal (not 0)"),
                      ftype = replace(ftype, ftype == "frac_notequal", "Difference > 0"),
                      ftype = replace(ftype, ftype == "frac_morethan5pdiff", "Difference > 5%"),
                      ftype = replace(ftype, ftype == "frac_morethan10pdiff", "Difference > 10%")) %>%
        dplyr::mutate(ftype = factor(ftype, levels = c("Equal", "Equal (not 0)",
                                                       "Difference > 0",
                                                       "Difference > 5%",
                                                       "Difference > 10%"))) %>%
        ggplot(aes(x = ftype, y = fraction)) + 
        geom_bar(stat = "identity", alpha = 0.7, aes(fill = ctype)) + 
        facet_grid(m1 ~ m2) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              strip.text = element_text(size = 5)) + 
        scale_fill_manual(values = c(spliced = "red", unspliced = "blue"),
                          name = "") + 
        labs(title = "Differences in total gene count across cells",
             x = "", y = "Fraction of genes")
)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
