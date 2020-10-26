shorten_methods <- function(methods) {
  data.frame(method = methods, 
             stringsAsFactors = FALSE) %>%
    dplyr::mutate(method_short = 
                    gsub("ude", "", 
                         gsub("kallisto_bustools", "kallisto|bus", 
                              gsub("collapse", "coll", 
                                   gsub("separate", "sep",
                                        gsub("iso", "", 
                                             gsub("prepref_", "",
                                                  gsub("_cdna_introns", "", method)))))))) %>%
    dplyr::mutate(method_short = replace(method_short, method_short == "starsolo_subtr", 
                                         "starsolo_diff")) %>%
    dplyr::mutate(method_short = replace(method_short, method_short == "kb_python_lamanno", 
                                         "kb_python")) %>%
    dplyr::mutate(method_short = gsub("gentrome", "gtr", method_short)) %>%
    dplyr::mutate(
      mtype = stringr::str_extract(
        method_short, "alevin|kallisto\\|bus|starsolo|velocyto|dropest"
      ),
      rtype = stringr::str_extract(
        method, "separate|collapse"
      )
    ) %>%
    dplyr::mutate(mtype = replace(mtype, method_short == "kb_python", "kallisto|bus")) %>%
    dplyr::mutate(rtype = replace(rtype, method_short == "kb_python", "separate")) %>% 
    dplyr::mutate(rtype = replace(rtype, is.na(rtype), "N/A")) %>%
    dplyr::mutate(rtype = factor(rtype, levels = c("collapse", "separate", "N/A")))
}

base_method_colors <- c(alevin = "#999999", `kallisto|bus` = "#009E73",
                        starsolo = "#0072B2", velocyto = "#CC79A7", dropest = "#CFAE2B")

merge_uniq <- function(refdir, tx2gene, keepgenes) {
  uniq <- dplyr::bind_rows(
    read.delim(
      file.path(refdir, "prepref_isoseparate_uniqueness.txt"),
      header = TRUE, as.is = TRUE
    ) %>% 
      dplyr::mutate(ctype = c("exonic", "intronic")[grepl("I\\.", gene) + 1]) %>%
      dplyr::mutate(gene = gsub("I\\.", "", gene)) %>%
      dplyr::mutate(frac_unique = unique/total) %>%
      dplyr::select(gene, ctype, frac_unique) %>%
      dplyr::mutate(atype = "separate"),
    read.delim(
      file.path(refdir, "prepref_isocollapse_uniqueness.txt"),
      header = TRUE, as.is = TRUE
    ) %>% 
      dplyr::mutate(ctype = c("exonic", "intronic")[grepl("I\\.", gene) + 1]) %>%
      dplyr::mutate(gene = gsub("I\\.", "", gene)) %>%
      dplyr::mutate(frac_unique = unique/total) %>%
      dplyr::select(gene, ctype, frac_unique) %>%
      dplyr::mutate(atype = "collapse"),
    read.delim(
      file.path(refdir, "prepref_isoseparate_uniqueness_overall.txt"),
      header = TRUE, as.is = TRUE
    ) %>% 
      dplyr::mutate(ctype = "overall") %>%
      dplyr::mutate(gene = gsub("I\\.", "", gene)) %>%
      dplyr::mutate(frac_unique = unique/total) %>%
      dplyr::select(gene, ctype, frac_unique) %>%
      dplyr::mutate(atype = "separate"),
    read.delim(
      file.path(refdir, "prepref_isocollapse_uniqueness_overall.txt"),
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
    dplyr::filter(gene %in% keepgenes) %>%
    dplyr::group_by(ctype, atype, frac_unique_bin) %>%
    dplyr::mutate(nbr_genes = length(gene)) %>% 
    dplyr::ungroup() %>%
    tidyr::unite("frac_unique_bin", frac_unique_bin, nbr_genes, sep = ", n = ")
  
  uniq
}

cluster_levels <- list(
  Dentate_gyrus = c("nIPC", "Neuroblast", "Mossy", "Cck-Tox", 
                    "Granule immature", "Granule mature", "Microglia", 
                    "Endothelial", "Radial Glia-like", "Astrocytes",
                    "OPC", "OL", "GABA", "Cajal Retzius"),
  Pancreas = c("Ductal", "Ngn3 low EP", "Ngn3 high EP", 
               "Pre-endocrine", "Epsilon", "Delta", "Alpha", "Beta"),
  PFC = c("Astro", "Endo", "Excitatory", "Inhibitory", "Microglia",
          "NF Oligo", "Oligo", "OPC"),
  Spermatogenesis = c("A3-A4-In-B Differentiating spermatogonia", 
                      "Leptotene/Zygotene spermatocytes", 
                      "Pachytene spermatocytes", "DIplotene/Secondary spermatocytes",
                      "Early Round spermatids", "Mid Round spermatids",
                      "Late Round spermatids"),
  Neuron = NULL
)
