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
    dplyr::mutate(
      mtype = stringr::str_extract(
        method_short, "alevin|kallisto\\|bus|starsolo|velocyto"
      ),
      rtype = stringr::str_extract(
        method, "separate|collapse"
      )
    ) %>%
    dplyr::mutate(rtype = replace(rtype, is.na(rtype), "N/A")) %>%
    dplyr::mutate(rtype = factor(rtype, levels = c("collapse", "separate", "N/A")))
}

base_method_colors <- c(alevin = "#999999", `kallisto|bus` = "#009E73",
                        starsolo = "#0072B2", velocyto = "#CC79A7")
