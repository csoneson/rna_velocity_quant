args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(BiocParallel)
  library(BiocSingular)
})

print(topdir)
print(dataset)
print(helperscript)
print(tx2gene)
print(cellfile)
print(samplename)
print(outrds)

source(helperscript)

tx2gene <- readRDS(tx2gene)
if (cellfile != "") {
  cells <- read.csv(cellfile)
} else {
  cells <- NULL
}

sces <- list()

## ========================================================================= ##
## Read quantifications and create SingleCellExperiment objects
## ========================================================================= ##
## CellRanger + velocyto
sces$velocyto <- read_velocyto(
  loomfile = file.path(topdir, paste0("quants/cellranger/", samplename, "/velocyto/", samplename, ".loom")), 
  sampleid = samplename
)

## STARsolo
sces$starsolo <- read_starsolo(
  solodir = file.path(topdir, paste0("quants/starsolo/Solo.out/Velocyto/raw")),
  sampleid = samplename
)

sces$starsolo_subtr <- read_starsolo_subtract(
  solodir = file.path(topdir, paste0("quants/starsolo/Solo.out")),
  sampleid = samplename
)

## cDNA/introns separately (with decoys)
for (m in c("prepref")) {
  for (v in c("separate", "collapse")) {
    sces[[paste0("alevin_", m, "_iso", v, "_cdna_introns_decoy")]] <- 
      read_alevin_with_decoys(
        spliceddir = file.path(topdir, paste0("quants/alevin_", m, "_iso", v, "_cdna_intronsasdecoy/alevin")),
        unspliceddir = file.path(topdir, paste0("quants/alevin_", m, "_iso", v, "_introns_cdnaasdecoy/alevin")),
        sampleid = samplename, tx2gene = tx2gene
      )
    
    sces[[paste0("alevin_", m, "_iso", v, "_cdna_introns_decoy_gentrome")]] <- 
      read_alevin_with_decoys(
        spliceddir = file.path(topdir, paste0("quants/alevin_", m, "_iso", v, "_cdna_intronsasdecoy_gentrome/alevin")),
        unspliceddir = file.path(topdir, paste0("quants/alevin_", m, "_iso", v, "_introns_cdnaasdecoy_gentrome/alevin")),
        sampleid = samplename, tx2gene = tx2gene
      )
  }
}

## cDNA/introns quantified jointly
for (m in c("prepref")) {
  for (v in c("separate", "collapse")) {
    sces[[paste0("alevin_", m, "_iso", v, "_cdna_introns")]] <- 
      read_alevin_cdna_introns(
        alevindir = file.path(topdir, paste0("quants/alevin_", m, "_iso", v, "_cdna_introns/alevin")),
        sampleid = samplename, tx2gene = tx2gene)
    
    sces[[paste0("alevin_", m, "_iso", v, "_cdna_introns_gentrome")]] <- 
      read_alevin_cdna_introns(
        alevindir = file.path(topdir, paste0("quants/alevin_", m, "_iso", v, "_cdna_introns_gentrome/alevin")),
        sampleid = samplename, tx2gene = tx2gene)
    
    if (v == "separate") {
      sces[[paste0("alevin_", m, "_iso", v, "_cdna_introns_gentrome_unstranded")]] <- 
        read_alevin_cdna_introns(
          alevindir = file.path(topdir, paste0("quants/alevin_", m, "_iso", v, "_cdna_introns_gentrome_unstranded/alevin")),
          sampleid = samplename, tx2gene = tx2gene)
      
      for (u in c("_flankL20", "_flankL40")) {
        sces[[paste0("alevin_", m, "_iso", v, u, "_cdna_introns_gentrome")]] <- 
          read_alevin_cdna_introns(
            alevindir = file.path(topdir, paste0("quants/alevin_", m, "_iso", v, u, "_cdna_introns_gentrome/alevin")),
            sampleid = samplename, tx2gene = tx2gene)
      }
    }
  }
}

## alevin spliced/unspliced
sces$alevin_spliced_unspliced <- 
  read_alevin_spliced_unspliced(
    alevindir = file.path(topdir, "quants/alevin_spliced_unspliced/alevin"),
    sampleid = samplename, tx2gene = tx2gene
  )
sces$alevin_spliced_unspliced_gentrome <- 
  read_alevin_spliced_unspliced(
    alevindir = file.path(topdir, "quants/alevin_spliced_unspliced_gentrome/alevin"),
    sampleid = samplename, tx2gene = tx2gene
  )

## alevin spliced
sces$alevin_spliced <- 
  read_alevin_spliced(
    alevindir = file.path(topdir, "quants/alevin_spliced/alevin"),
    sampleid = samplename, tx2gene = tx2gene
  )
sces$alevin_spliced_gentrome <- 
  read_alevin_spliced(
    alevindir = file.path(topdir, "quants/alevin_spliced_gentrome/alevin"),
    sampleid = samplename, tx2gene = tx2gene
  )

## kallisto/bustools
for (m in c("prepref")) {
  for (v in c("separate", "collapse")) {
    for (k in c("exclude", "include")) {
      sces[[paste0("kallisto_bustools_", m, "_iso", v, "_", k)]] <- 
        read_kallisto_bustools(
          kallistodir = file.path(topdir, paste0("quants/kallisto_bustools_", 
                                                 m, "_iso", v, "_cdna_introns")),
          splicedname = paste0("spliced.", k),
          unsplicedname = paste0("unspliced.", k)
        )
    }
  }
}

## kb-python
sces[["kb_python_lamanno"]] <- 
  read_kallisto_bustools(
    kallistodir = file.path(topdir, "quants/kb_python_lamanno/counts_unfiltered"),
    splicedname = "spliced",
    unsplicedname = "unspliced"
  )

## ========================================================================= ##
## subset to shared cells/genes
## ========================================================================= ##
shared_cells <- as.character(Reduce(intersect, lapply(sces, colnames)))
if (!is.null(cells)) {
  shared_cells <- intersect(shared_cells, as.character(cells$index))
}
shared_genes <- as.character(Reduce(intersect, lapply(sces, rownames)))

sces <- lapply(sces, function(w) w[shared_genes, shared_cells])

if (!is.null(cells)) {
  cells <- cells[match(shared_cells, cells$index), ] %>%
    dplyr::rename(cell_index = index)
  cells$cell_index <- as.character(cells$cell_index)
  cells$clusters <- as.character(cells$clusters)
  cells$clusters_coarse <- as.character(cells$clusters_coarse)
  
  sces <- lapply(sces, function(w) {
    colData(w) <- cbind(colData(w), DataFrame(cells))
    w
  })
}

## ========================================================================= ##
## Add reduced dimension representation + clusters
## ========================================================================= ##
do_dimred <- function(sce) {
  sce <- scater::logNormCounts(sce)
  set.seed(1)
  sce <- scater::runPCA(sce, exprs_values = "logcounts", 
                        ncomponents = 30,
                        BSPARAM = BiocSingular::IrlbaParam())
  sce <- scater::runTSNE(sce, dimred = "PCA", 
                         ncomponents = 2)
  sce <- scater::runUMAP(sce, dimred = "PCA", 
                         ncomponents = 2)
  
  snn.gr <- scran::buildSNNGraph(sce, use.dimred = "PCA")
  clusters <- igraph::cluster_walktrap(snn.gr)
  sce$cluster <- factor(clusters$membership)
  
  sce
}

## Calculate logcounts and reduced dimensions
for (nm in names(sces)) {
  message(nm)
  sces[[nm]] <- do_dimred(sces[[nm]])
}

## Add common representations to all data sets
## From alevin_spliced_gentrome
sces <- lapply(sces, function(w) {
  reducedDim(w, "PCA_alevin_spliced_gentrome") <- reducedDim(sces[["alevin_spliced_gentrome"]], "PCA")
  reducedDim(w, "TSNE_alevin_spliced_gentrome") <- reducedDim(sces[["alevin_spliced_gentrome"]], "TSNE")
  reducedDim(w, "UMAP_alevin_spliced_gentrome") <- reducedDim(sces[["alevin_spliced_gentrome"]], "UMAP")
  w$cluster_alevin_spliced_gentrome <- sces[["alevin_spliced_gentrome"]]$cluster
  w
})

## From starsolo
sces <- lapply(sces, function(w) {
  reducedDim(w, "PCA_starsolo") <- reducedDim(sces[["starsolo"]], "PCA")
  reducedDim(w, "TSNE_starsolo") <- reducedDim(sces[["starsolo"]], "TSNE")
  reducedDim(w, "UMAP_starsolo") <- reducedDim(sces[["starsolo"]], "UMAP")
  w$cluster_starsolo <- sces[["starsolo"]]$cluster
  w
})

for (m in c("starsolo")) {
  message(m)
  ## Concatenated spliced and unspliced
  tmp <- SingleCellExperiment(
    assays = list(counts = rbind(assay(sces[[m]], "spliced"), 
                                 assay(sces[[m]], "unspliced"))))
  tmp <- do_dimred(tmp)
  sces <- lapply(sces, function(w) {
    reducedDim(w, paste0("PCA_", m, "_concatenated")) <- reducedDim(tmp, "PCA")
    reducedDim(w, paste0("TSNE_", m, "_concatenated")) <- reducedDim(tmp, "TSNE")
    reducedDim(w, paste0("UMAP_", m, "_concatenated")) <- reducedDim(tmp, "UMAP")
    colData(w)[[paste0("cluster_", m, "_concatenated")]] <- tmp$cluster
    w
  })
  
  ## Summed spliced and unspliced
  tmp <- SingleCellExperiment(
    assays = list(counts = assay(sces[[m]], "spliced") +  
                    assay(sces[[m]], "unspliced")))
  tmp <- do_dimred(tmp)
  sces <- lapply(sces, function(w) {
    reducedDim(w, paste0("PCA_", m, "_summed")) <- reducedDim(tmp, "PCA")
    reducedDim(w, paste0("TSNE_", m, "_summed")) <- reducedDim(tmp, "TSNE")
    reducedDim(w, paste0("UMAP_", m, "_summed")) <- reducedDim(tmp, "UMAP")
    colData(w)[[paste0("cluster_", m, "_summed")]] <- tmp$cluster
    w
  })

  ## Only unspliced
  tmp <- SingleCellExperiment(
    assays = list(counts = assay(sces[[m]], "unspliced")))
  tmp <- do_dimred(tmp)
  sces <- lapply(sces, function(w) {
    reducedDim(w, paste0("PCA_", m, "_unspliced")) <- reducedDim(tmp, "PCA")
    reducedDim(w, paste0("TSNE_", m, "_unspliced")) <- reducedDim(tmp, "TSNE")
    reducedDim(w, paste0("UMAP_", m, "_unspliced")) <- reducedDim(tmp, "UMAP")
    colData(w)[[paste0("cluster_", m, "_unspliced")]] <- tmp$cluster
    w
  })
  
}

## ========================================================================= ##
## Save
## ========================================================================= ##
for (nm in names(sces)) {
  metadata(sces[[nm]]) <- list(method = nm, dataset = dataset)
  saveRDS(sces[[nm]], file = file.path(dirname(outrds), paste0("sce_", nm, ".rds")))
}
write.table(paste0(shared_cells, "-1"), 
            file = file.path(dirname(outrds), 
                             paste0("retained_cell_barcodes.csv")),
            row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
saveRDS(NULL, file = outrds)

date()
sessionInfo()
