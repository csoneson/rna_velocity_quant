args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
source(plothelperscript)

methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods
dataset <- gsub("_", " ", dataset)

print(plothelperscript)
print(velosuffix)
print(topdir)
print(dataset)
print(methods)
print(outrds)

methods_short <- shorten_methods(methods) %>%
  dplyr::filter(method %in% methods) %>%
  dplyr::arrange(as.factor(method_short))
methods <- methods_short$method
names(methods) <- methods

## ----------------------------------------------------------------------------
## Read cell info for each shuffled/perturbed instance
shufres_cell <- do.call(dplyr::bind_rows, lapply(methods, function(nm) {
  idxs <- readr::read_csv(file.path(
    topdir, paste0("plots/velocity_shuffle", velosuffix, 
                   "/anndata_", nm, "/velocity_shuffle.out")
  ), col_names = FALSE, col_types = "d")
  idxp <- readr::read_csv(file.path(
    topdir, paste0("plots/velocity_perturb", velosuffix, 
                   "/anndata_", nm, "/velocity_perturb.out")
  ), col_names = FALSE, col_types = "d")
  
  if (length(idxs) > 0) {
    idxs <- idxs[[1]]
    shuf <- lapply(idxs, function(i) {
      readr::read_csv(file.path(
        topdir, paste0("plots/velocity_shuffle", velosuffix, "/anndata_", nm, "/", 
                       "anndata_", nm, "_cell_info_", i, ".csv"))) %>%
        dplyr::select(index, clusters, max_cosine_corr, velocity_length, 
                      velocity_confidence, contains("latent_time")) %>%
        dplyr::mutate(iter = paste0("S", i), type = "shuffle", 
                      dataset = dataset)
    })
  } else {
    shuf <- NULL
  }
  if (length(idxp) > 0) {
    idxp <- idxp[[1]]
    pert <- lapply(idxp, function(i) {
      readr::read_csv(file.path(
        topdir, paste0("plots/velocity_perturb", velosuffix, "/anndata_", nm, "/", 
                       "anndata_", nm, "_cell_info_", i, ".csv"))) %>%
        dplyr::select(index, clusters, max_cosine_corr, velocity_length, 
                      velocity_confidence, contains("latent_time")) %>%
        dplyr::mutate(iter = paste0("P", i), type = "perturb",
                      dataset = dataset)
    })
  } else {
    pert <- NULL
  }
  do.call(dplyr::bind_rows, c(shuf, pert)) %>%
    dplyr::mutate(method = methods_short$method_short[match(nm, methods_short$method)])
}))

## ----------------------------------------------------------------------------
## Read gene info for each shuffled/perturbed instance
shufres_gene <- do.call(dplyr::bind_rows, lapply(methods, function(nm) {
  idxs <- readr::read_csv(file.path(
    topdir, paste0("plots/velocity_shuffle", velosuffix, 
                   "/anndata_", nm, "/velocity_shuffle.out")
  ), col_names = FALSE, col_types = "d")
  idxp <- readr::read_csv(file.path(
    topdir, paste0("plots/velocity_perturb", velosuffix, 
                   "/anndata_", nm, "/velocity_perturb.out")
  ), col_names = FALSE, col_types = "d")
  
  if (length(idxs) > 0) {
    idxs <- idxs[[1]]
    shuf <- lapply(idxs, function(i) {
      readr::read_csv(file.path(
        topdir, paste0("plots/velocity_shuffle", velosuffix, "/anndata_", nm, "/", 
                       "anndata_", nm, "_gene_info_", i, ".csv"))) %>%
        dplyr::select(index, fit_alpha, fit_beta, fit_gamma, fit_likelihood,
                      velocity_genes) %>%
        dplyr::mutate(iter = paste0("S", i), type = "shuffle",
                      dataset = dataset)
    })
  } else {
    shuf <- NULL
  }
  if (length(idxp) > 0) {
    idxp <- idxp[[1]]
    pert <- lapply(idxp, function(i) {
      readr::read_csv(file.path(
        topdir, paste0("plots/velocity_perturb", velosuffix, "/anndata_", nm, "/", 
                       "anndata_", nm, "_gene_info_", i, ".csv"))) %>%
        dplyr::select(index, fit_alpha, fit_beta, fit_gamma, fit_likelihood,
                      velocity_genes) %>%
        dplyr::mutate(iter = paste0("P", i), type = "perturb",
                      dataset = dataset)
    })
  } else {
    pert <- NULL
  }
  do.call(dplyr::bind_rows, c(shuf, pert)) %>%
    dplyr::mutate(method = methods_short$method_short[match(nm, methods_short$method)])
}))

## ----------------------------------------------------------------------------
## Read projected velocity information for each shuffled/perturbed instance
shufres_velo <- do.call(dplyr::bind_rows, lapply(methods, function(nm) {
  idxs <- readr::read_csv(file.path(
    topdir, paste0("plots/velocity_shuffle", velosuffix, 
                   "/anndata_", nm, "/velocity_shuffle.out")
  ), col_names = FALSE, col_types = "d")
  idxp <- readr::read_csv(file.path(
    topdir, paste0("plots/velocity_perturb", velosuffix, 
                   "/anndata_", nm, "/velocity_perturb.out")
  ), col_names = FALSE, col_types = "d")
  
  if (length(idxs) > 0) {
    idxs <- idxs[[1]]
    shuf <- lapply(idxs, function(i) {
      readr::read_csv(file.path(
        topdir, paste0("plots/velocity_shuffle", velosuffix, "/anndata_", nm, "/", 
                       "anndata_", nm, 
                       "_velocity_UMAP_alevin_spliced_gentrome_", i, ".csv"))) %>%
        setNames(c("index", "xvel", "yvel")) %>%
        dplyr::mutate(iter = paste0("S", i), type = "shuffle",
                      dataset = dataset)
    })
  } else {
    shuf <- NULL
  }
  if (length(idxp) > 0) {
    idxp <- idxp[[1]]
    pert <- lapply(idxp, function(i) {
      readr::read_csv(file.path(
        topdir, paste0("plots/velocity_perturb", velosuffix, "/anndata_", nm, "/", 
                       "anndata_", nm, 
                       "_velocity_UMAP_alevin_spliced_gentrome_", i, ".csv"))) %>%
        setNames(c("index", "xvel", "yvel")) %>%
        dplyr::mutate(iter = paste0("P", i), type = "perturb",
                      dataset = dataset)
    })
  } else {
    pert <- NULL
  }
  do.call(dplyr::bind_rows, c(shuf, pert)) %>%
    dplyr::mutate(method = methods_short$method_short[match(nm, methods_short$method)])
})) %>%
  dplyr::mutate(ds_m_tp = paste0(dataset, "__", method, "__", type))

shufres_velo_list <- split(shufres_velo, f = shufres_velo$ds_m_tp)
velores <- do.call(dplyr::bind_rows, lapply(shufres_velo_list, function(ds) {
  d <- unique(ds$dataset)
  m <- unique(ds$method)
  tp <- unique(ds$type)
  do.call(dplyr::bind_rows, lapply(unique(ds$iter), function(s1) {
    do.call(dplyr::bind_rows, lapply(unique(ds$iter), function(s2) {
      if (s1 < s2) {
        d1 <- ds %>% 
          dplyr::filter(dataset == d & 
                          method == m & 
                          type == tp & 
                          iter == s1) %>%
          dplyr::select(index, xvel, yvel)
        d2 <- ds %>% 
          dplyr::filter(dataset == d & 
                          method == m & 
                          type == tp & 
                          iter == s2) %>%
          dplyr::select(index, xvel, yvel)
        dplyr::left_join(d1, d2, by = "index") %>%
          dplyr::mutate(costheta = (xvel.x * xvel.y + 
                                      yvel.x * yvel.y) /
                          (sqrt(xvel.x^2 + yvel.x^2) * 
                             sqrt(xvel.y^2 + yvel.y^2))) %>%
          dplyr::mutate(dataset = d, method = m, type = tp, 
                        iter1 = s1, iter2 = s2)
      } else {
        NULL
      }
    }))
  }))
}))

pdf(gsub("\\.rds", ".pdf", outrds), width = 8, height = 8)

## ----------------------------------------------------------------------------
## Plot number of valid instances
ggplot(shufres_cell %>% 
         dplyr::group_by(method, type) %>%
         dplyr::summarize(n_iters = length(unique(iter))),
       aes(x = method, y = n_iters, fill = type)) + 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = "", y = "Number of valid iterations")

## ----------------------------------------------------------------------------
## Plot number of velocity genes
ggplot(shufres_gene %>% dplyr::filter(velocity_genes) %>%
         dplyr::group_by(dataset, iter, type, method) %>%
         dplyr::tally(),
       aes(x = method, y = n, fill = type)) + 
  geom_boxplot(position = position_dodge(preserve = "single")) + 
  scale_y_log10() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = "", y = "Number of velocity genes")

## ----------------------------------------------------------------------------
## Plot distributions of cell-wise characteristics
for (char in c("max_cosine_corr", "velocity_length", "velocity_confidence")) {
  g <- ggplot(shufres_cell,
              aes(x = method, y = !!rlang::sym(char), fill = type)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    labs(x = "", y = char)
  if (char == "velocity_length") {
    g <- g + scale_y_sqrt()
  }
  print(g)
}

## ----------------------------------------------------------------------------
## Plot distributions of correlations of cell-wise characteristics among iterations
shufres_cell_long <- shufres_cell %>% 
  tidyr::gather(key = "measure", value = "value", 
                max_cosine_corr, velocity_length, 
                velocity_confidence, latent_time) %>%
  dplyr::mutate(method_measure = paste0(method, "__", measure, "__", type),
                iter_measure = paste0(iter, "__", measure, "__", type))
corrs_within_method_cell <- do.call(
  dplyr::bind_rows, 
  lapply(split(shufres_cell_long, shufres_cell_long$method_measure), function(sr) {
    cor(sr %>% dplyr::select(index, value, iter) %>%
          tidyr::spread(key = iter, value = value) %>%
          dplyr::select(-index), use = "pairwise.complete.obs") %>%
      as.data.frame() %>% 
      tibble::rownames_to_column("iter1") %>%
      tidyr::gather(key = "iter2", value = "correlation", -iter1) %>%
      dplyr::filter(iter1 < iter2) %>%
      dplyr::mutate(method = unique(sr$method), measure = unique(sr$measure), 
                    type = unique(sr$type), dataset = unique(sr$dataset))
  })
)
for (char in c("max_cosine_corr", "velocity_length",
               "velocity_confidence", "latent_time")) {
  print(ggplot(corrs_within_method_cell %>% dplyr::filter(measure == char),
               aes(x = method, y = correlation, fill = type)) + 
          geom_boxplot(position = position_dodge(preserve = "single")) + 
          theme_bw() + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
          labs(x = "", y = "Correlation among iterations", title = char))
}

## ----------------------------------------------------------------------------
## Plot distributions of correlations of gene-wise characteristics among iterations
shufres_gene_long <- shufres_gene %>% 
  tidyr::gather(key = "measure", value = "value", 
                fit_alpha, fit_beta, 
                fit_gamma, fit_likelihood) %>%
  dplyr::mutate(method_measure = paste0(method, "__", measure, "__", type),
                iter_measure = paste0(iter, "__", measure, "__", type))
corrs_within_method_gene <- do.call(
  dplyr::bind_rows, 
  lapply(split(shufres_gene_long, shufres_gene_long$method_measure), function(sr) {
    cor(sr %>% dplyr::select(index, value, iter) %>%
          tidyr::spread(key = iter, value = value) %>%
          dplyr::select(-index), use = "pairwise.complete.obs") %>%
      as.data.frame() %>% 
      tibble::rownames_to_column("iter1") %>%
      tidyr::gather(key = "iter2", value = "correlation", -iter1) %>%
      dplyr::filter(iter1 < iter2) %>%
      dplyr::mutate(method = unique(sr$method), measure = unique(sr$measure), 
                    type = unique(sr$type), dataset = unique(sr$dataset))
  })
)
## Not clear if it makes sense to look at correlations in the shuffled data 
## since the number of retained genes is so low
for (char in c("fit_alpha", "fit_beta", "fit_gamma", "fit_likelihood")) {
  print(ggplot(corrs_within_method_gene %>% dplyr::filter(type == "perturb") %>% 
                 dplyr::filter(measure == char),
               aes(x = method, y = correlation, fill = type)) + 
          geom_boxplot(position = position_dodge(preserve = "single")) + 
          theme_bw() + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
          labs(x = "", y = "Correlation among iterations", title = char))
}

## ----------------------------------------------------------------------------
## Plot distribution of cos(theta), where theta is the angle between the 
## projected velocity vector for a cell in two different instances
ggplot(velores,
       aes(x = method, y = costheta, fill = type)) + 
  geom_violin(scale = "width") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = "", y = "cos(theta)")


dev.off()
saveRDS(list(shufres_cell = shufres_cell, shufres_gene = shufres_gene, 
             corrs_within_method_cell = corrs_within_method_cell, 
             corrs_within_method_gene = corrs_within_method_gene,
             velores = velores), file = outrds)

date()
sessionInfo()
