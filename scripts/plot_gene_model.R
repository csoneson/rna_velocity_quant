args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(swissknife)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(magick)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
})

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(topdir)
print(plothelperscript)
print(gtf)
print(methods)
print(bigwigfile)
print(samplename)
print(showgene)
print(outpdf)

source(plothelperscript)

gr <- prepareGTF(gtf)
bwf <- bigwigfile
names(bwf) <- samplename

sces <- lapply(methods, function(nm) {
  readRDS(file.path(topdir, paste0("output/sce/sce_", nm, ".rds")))
})

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

methods_short <- shorten_methods(methods)

plot_gene_model <- function(gr, bwf, showgene, sumdfg, methodsdf) {
  
  rn <- round(1e7 * runif(1))
  tmpdir <- tempdir()
  grDevices::png(paste0(tmpdir, "/gviz", rn, ".png"), width = 10.5,
                 height = 5.25, unit = "in", res = 400)
  swissknife::plotGeneRegion(granges = gr, bigwigFiles = bwf, 
                             showgene = showgene, geneTrackTitle = showgene)
  grDevices::dev.off()
  
  p2 <- ggplot(sumdfg %>% dplyr::filter(gene == showgene) %>%
                 dplyr::select(gene, method, spliced, unspliced) %>%
                 tidyr::gather(key = "ctype", value = "count", spliced, unspliced) %>%
                 dplyr::left_join(methods_short, by = "method"),
               aes(x = method_short, y = count, fill = mtype)) + 
    geom_bar(stat = "identity") + facet_wrap(~ ctype, ncol = 1) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none") + 
    scale_fill_manual(values = base_method_colors) + 
    labs(x = "", title = "Total UMI count", y = "")
  
  plt <- cowplot::plot_grid(
    cowplot::ggdraw() +
      cowplot::draw_image(paste0(tmpdir, "/gviz", rn, ".png")),
    p2,
    nrow = 1, rel_widths = c(1, 0.3), labels = ""
  )
  unlink(paste0(tmpdir, "/gviz", rn, ".png"))
  plt
}

pdf(outpdf, width = 10, height = 4.5)
print(plot_gene_model(gr = gr, bwf = bwf, showgene = showgene,
                      sumdfg = sumdf_bygene, methodsdf = methods_short))
dev.off()

date()
sessionInfo()
