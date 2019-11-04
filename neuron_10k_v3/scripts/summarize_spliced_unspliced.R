args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

source("scripts/sce_helpers.R")

print(topdir)
print(tx2gene)
print(outrds)

tx2gene <- readRDS(tx2gene)

sces <- list()

## ========================================================================= ##
## Read quantifications and create SingleCellExperiment objects
## ========================================================================= ##
## CellRanger + velocyto
sces$velocyto <- read_velocyto(loomfile = file.path(topdir, "quants/cellranger/neuron_10k_v3/velocyto/neuron_10k_v3.loom"), 
                               sampleid = "neuron_10k_v3")


## cDNA/introns separately (with decoys)
for (m in c("busparse", "prepref")) {
  for (v in c("separate", "collapse")) {
    sces[[paste0("alevin_", m, "_iso", v, "_cdna_introns_decoy")]] <- 
      read_alevin_with_decoys(
        spliceddir = file.path(topdir, paste0("quants/alevin", m, "_iso", v, "_cdna_intronsasdecoy/alevin")),
        unspliceddir = file.path(topdir, paste0("quants/alevin", m, "_iso", v, "_introns_cdnaasdecoy/alevin")),
        sampleid = "neuron_10k_v3", tx2gene = tx2gene
      )
  }
}

## cDNA/introns quantified jointly
for (m in c("busparse", "prepref")) {
  for (v in c("separate", "collapse")) {
    sces[[paste0("alevin_", m, "_iso", v, "_cdna_introns_decoy")]] <- 
      read_alevin_cdna_introns(
        alevindir = file.path(topdir, paste0("quants/alevin_", m, "_iso", v, "_cdna_introns/alevin")),
        sampleid = "neuron_10k_v3", tx2gene = tx2gene)
  }
}

sces$alevin_spliced_unspliced <- 
  read_alevin_spliced_unspliced(alevindir = file.path(topdir, "quants/alevin_spliced_unspliced/alevin"),
                                sampleid = "neuron_10k_v3", tx2gene = tx2gene)

sces$alevin_spliced <- 
  read_alevin_spliced(alevindir = file.path(topdir, "quants/alevin_spliced/alevin"),
                      sampleid = "neuron_10k_v3", tx2gene = tx2gene)


for (m in c("busparse", "prepref")) {
  for (v in c("separate", "collapse")) {
    for (k in c("exclude", "include")) {
      sces[[paste0("kallisto_bustools_", m, "_iso", v, "_", k)]] <- 
        read_kallisto_bustools(kallistodir = file.path(topdir, paste0("quants/kallisto_bustools_", 
                                                                      m, "_iso", v, "_cdna_introns")),
                               splicedname = paste0("spliced.", k),
                               unsplicedname = paste0("unspliced.", k))
    }
  }
}

## ========================================================================= ##
## subset to shared cells/genes
## ========================================================================= ##
shared_cells <- Reduce(intersect, lapply(sces, colnames))
shared_genes <- Reduce(intersect, lapply(sces, rownames))

sces <- lapply(sces, function(w) w[shared_genes, shared_cells])

## ========================================================================= ##
## Save
## ========================================================================= ##
for (nm in names(sces)) {
  saveRDS(sces[[nm]], file = file.path(dirname(outrds), paste0("sce_", nm, ".rds")))
}
saveRDS(NULL, file = outrds)

date()
sessionInfo()
