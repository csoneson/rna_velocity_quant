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
print(refdir)  ## directory where uniqueness files are
print(plothelperscript)
print(gtf)
print(tx2gene)
print(methods)
print(bigwigfileplus)
print(bigwigfileminus)
print(samplename)
print(showgene)
print(outpdf)

source(plothelperscript)

gr <- prepareGTF(gtf)
bwf <- c(bigwigfileplus, bigwigfileminus)
names(bwf) <- c("pos", "neg")
bwc <- c("pos", "neg")
names(bwc) <- c("pos", "neg")
tx2gene <- readRDS(tx2gene)

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

uniq <- merge_uniq(refdir = refdir, tx2gene = tx2gene,
                   keepgenes = sumdf_bygene$gene)

methods_short <- shorten_methods(methods)

plot_gene_model <- function(gr, bwf, bwc, showgene, sumdfg, methodsdf, uniq) {
  
  rn <- round(1e7 * runif(1))
  tmpdir <- tempdir()
  grDevices::png(paste0(tmpdir, "/gviz", rn, ".png"), width = 10.5,
                 height = 5.25, unit = "in", res = 400)
  swissknife::plotGeneRegion(granges = gr, bigwigFiles = bwf, 
                             bigwigCond = bwc, 
                             showgene = showgene, geneTrackTitle = showgene,
                             colorByStrand = TRUE, 
                             condColors = c(pos = "blue", neg = "red"),
                             scaleDataTracks = TRUE)
  grDevices::dev.off()
  
  plotdfg <- sumdfg %>% dplyr::filter(gene == showgene) %>%
    dplyr::select(gene, method, spliced, unspliced, frac_unspliced) %>%
    tidyr::gather(key = "ctype", value = "count", spliced, unspliced, frac_unspliced) %>%
    dplyr::left_join(methodsdf, by = "method") %>%
    dplyr::mutate(ctype = replace(ctype, ctype == "frac_unspliced", "frac unspl")) %>%
    dplyr::mutate(ctype = factor(ctype, levels = c("spliced", "unspliced", "frac unspl")))
  dummy <- data.frame(method_short = plotdfg$method_short[1], 
                      ctype = c("spliced", "unspliced", "frac unspl"),
                      mtype = plotdfg$mtype[1],
                      count = c(rep(max(plotdfg$count), 2), 1),
                      stringsAsFactors = FALSE) %>%
    dplyr::mutate(ctype = factor(ctype, levels = c("spliced", "unspliced", "frac unspl")))
  
  p2 <- ggplot(plotdfg,
               aes(x = method_short, y = count, fill = mtype)) + 
    geom_blank(data = dummy) + 
    geom_bar(stat = "identity") + 
    facet_wrap(~ ctype, nrow = 1, scales = "free_x") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none") + 
    scale_fill_manual(values = base_method_colors) + 
    labs(x = "", title = "Total UMI count/fraction unspliced count", y = "") + 
    scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
    coord_flip()
  
  ## Change relative widths of facets
  ## See https://stackoverflow.com/questions/52341385/how-to-automatically-adjust-the-width-of-each-facet-for-facet-wrap
  gp <- ggplotGrob(p2)
  # get gtable columns corresponding to the facets
  facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]
  gp$widths[facet.columns] <- unit(c(1, 1, 0.5), "null")
  p2 <- gp

  p3 <- ggplot(uniq %>% dplyr::filter(gene == showgene) %>%
                 tidyr::unite(col = "categ", ctype, atype, sep = "."),
               aes(x = categ, y = frac_unique)) + 
    geom_bar(stat = "identity") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none") + 
    labs(x = "", title = "Uniqueness", y = "") + 
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0, 0.05, 0)) + 
    coord_flip()
  
  plt <- cowplot::plot_grid(
    cowplot::ggdraw() +
      cowplot::draw_image(paste0(tmpdir, "/gviz", rn, ".png")),
    cowplot::plot_grid(p2, p3, nrow = 1, rel_widths = c(1, 0.45),
                       align = "h", axis = "b"),
    ncol = 1, rel_heights = c(1, 0.7), labels = ""
  )
  unlink(paste0(tmpdir, "/gviz", rn, ".png"))
  plt
}

pdf(outpdf, width = 8.5, height = 6)
print(plot_gene_model(gr = gr, bwf = bwf, bwc = bwc, showgene = showgene,
                      sumdfg = sumdf_bygene, methodsdf = methods_short,
                      uniq = uniq))
dev.off()

date()
sessionInfo()
