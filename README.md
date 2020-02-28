# Preprocessing choices affect RNA velocity results for droplet-based scRNA-seq data

This repository contains the scripts used to generate the results in 

* C Soneson, A Srivastava, R Patro, MB Stadler: Preprocessing choices affect RNA velocity results for droplet-based scRNA-seq data. bioRxiv (2020).

## Repository structure

In this manuscript, four real scRNA-seq data sets generated by 10x Genomics technology were analyzed:

* `Dentate gyrus` ([Hochgerner _et al_. 2018](https://www.ncbi.nlm.nih.gov/pubmed/29335606)) - GEO accession [GSE95315](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95315)
* `Pancreas` ([Bastidas-Ponce _et al_. 2019](https://www.ncbi.nlm.nih.gov/pubmed/31160421)) - GEO accession [GSE132188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132188)
* `PFC` ([Bhattacherjee _et al_. 2019](https://www.ncbi.nlm.nih.gov/pubmed/31519873)) - GEO accession [GSE124952](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124952) (sample accession [GSM3559979](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3559979))
* `Spermatogenesis` ([Hermann _et al_. 2018](https://www.ncbi.nlm.nih.gov/pubmed/30404016)) - GEO accession [GSE109033](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109033) (sample accession [GSM2928341](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2928341))

For each of these data sets, the corresponding subfolder of this repository contains the Makefile that was used to run the analysis. The R and Python scripts called in these Makefiles are provided in the `scripts/` folder. 

## Data preprocessing

The downloaded data were preprocessed as follows before being used in this evaluation:

### Dentate gyrus

* FASTQ files for individual cells from the p12 and p35 time points were downloaded (accession date November 23, 2019) and merged using the [GSE95315_download.sh](dentate_gyrus_mouse/data_preprocessing/GSE95315_download.sh) script.
* The sample annotation file ([`cells_kept_in_scvelo_exampledata.csv`](dentate_gyrus_mouse/data_preprocessing/cells_kept_in_scvelo_exampledata.csv)) was obtained via the [scVelo](https://scvelo.readthedocs.io/) Python package (v0.1.24):

```
import scvelo as scv
adata = scv.datasets.dentategyrus()
adata.obs.to_csv("cells_kept_in_scvelo_exampledata.csv")
```

### Pancreas

* FASTQ files for the E15.5 sample were downloaded from [ENA](https://www.ebi.ac.uk/ena/data/view/SRR9201794) (accession date November 24, 2019)
* The FASTQ files were subsequently renamed, replacing `_X` with `S1_L001_RX_001` (for `X` = 1,2)
* The sample annotation file ([`cells_kept_in_scvel_exampledata.csv`](pancreas_mouse/data_preprocessing/cells_kept_in_scvelo_exampledata.csv)) was obtained via the [scVelo](https://scvelo.readthedocs.io/) Python package (v0.1.24):

```
import scvelo as scv
aadata = scv.datasets.pancreatic_endocrinogenesis()
aadata.obs.to_csv("cells_kept_in_scvelo_exampledata.csv")
```

### PFC

* The BAM file with the reads was downloaded from [SRA](https://sra-pub-src-2.s3.amazonaws.com/SRR8433692/PFC_Sample2.bam.1)
* FASTQ files were extracted from the BAM file using [bamtofastq](https://github.com/10XGenomics/bamtofastq) v1.1.2:

```
bamtofastq_1.1.2/bamtofastq --reads-per-fastq=500000000 PFC_Sample2.bam FASTQtmp

## Concatenate FASTQ files from all flowcells and lanes
mkdir -p FASTQ
cat FASTQtmp/PFC_sample2_MissingLibrary_1_HJW35BCXY/bamtofastq_S1_L001_I1_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HJW35BCXY/bamtofastq_S1_L002_I1_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HMNNFBCXY/bamtofastq_S1_L001_I1_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HMNNFBCXY/bamtofastq_S1_L002_I1_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HMNTLBCXY/bamtofastq_S1_L001_I1_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HMNTLBCXY/bamtofastq_S1_L002_I1_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HNGGMBCXY/bamtofastq_S1_L001_I1_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HNGGMBCXY/bamtofastq_S1_L002_I1_001.fastq.gz > \
FASTQ/PFCsample2_S1_L001_I1_001.fastq.gz

cat FASTQtmp/PFC_sample2_MissingLibrary_1_HJW35BCXY/bamtofastq_S1_L001_R1_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HJW35BCXY/bamtofastq_S1_L002_R1_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HMNNFBCXY/bamtofastq_S1_L001_R1_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HMNNFBCXY/bamtofastq_S1_L002_R1_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HMNTLBCXY/bamtofastq_S1_L001_R1_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HMNTLBCXY/bamtofastq_S1_L002_R1_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HNGGMBCXY/bamtofastq_S1_L001_R1_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HNGGMBCXY/bamtofastq_S1_L002_R1_001.fastq.gz > \
FASTQ/PFCsample2_S1_L001_R1_001.fastq.gz

cat FASTQtmp/PFC_sample2_MissingLibrary_1_HJW35BCXY/bamtofastq_S1_L001_R2_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HJW35BCXY/bamtofastq_S1_L002_R2_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HMNNFBCXY/bamtofastq_S1_L001_R2_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HMNNFBCXY/bamtofastq_S1_L002_R2_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HMNTLBCXY/bamtofastq_S1_L001_R2_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HMNTLBCXY/bamtofastq_S1_L002_R2_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HNGGMBCXY/bamtofastq_S1_L001_R2_001.fastq.gz \
FASTQtmp/PFC_sample2_MissingLibrary_1_HNGGMBCXY/bamtofastq_S1_L002_R2_001.fastq.gz > \
FASTQ/PFCsample2_S1_L001_R2_001.fastq.gz
```

* The sample annotation file (`GSE124952_meta_data.csv`) was downloaded from the [GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124952). It was processed to retain only cells from PFCsample2 using the following R code, generating the final [`cells_in_pfcsample2.csv`](pfc_mouse/data_preprocessing/cells_in_pfcsample2.csv) file:

```
suppressPackageStartupMessages({
  library(dplyr)
})

csv <- read.delim("GSE124952_meta_data.csv", header = TRUE, as.is = TRUE, sep = ",")

## Only Sample2
csv <- csv %>% dplyr::filter(Sample == "PFCSample2") %>%
  dplyr::mutate(Sample = gsub("Sample", "sample", Sample)) %>%
  dplyr::mutate(index = gsub("PFCSample2_", "", X)) %>%
  dplyr::rename(clusters = CellType,
                clusters_coarse = L2_clusters) %>%
  dplyr::select(index, clusters_coarse, clusters, everything())

write.table(csv, file = "cells_in_pfcsample2.csv",
            row.names = FALSE, col.names = TRUE,
            sep = ",", quote = FALSE)
```

### Spermatogenesis

* The BAM file from the AdultMouse_Rep3 sample was downloaded (accession date January 23, 2020) from [SRA](https://sra-pub-src-1.s3.amazonaws.com/SRR6459157/AdultMouse_Rep3_possorted_genome_bam.bam.1).
* FASTQ files were extracted from the BAM file using [bamtofastq](https://github.com/10XGenomics/bamtofastq) v1.1.2:

```
bamtofastq_1.1.2/bamtofastq --reads-per-fastq=500000000 AdultMouse_Rep3_possorted_genome_bam.bam FASTQtmp
mkdir -p FASTQ
mv FASTQtmp/Ad-Ms-Total-Sorted_20k_count_MissingLibrary_1_HK2GNBBXX/bamtofastq_S1_L006_I1_001.fastq.gz FASTQ/AdultMouseRep3_S1_L001_I1_001.fastq.gz
mv FASTQtmp/Ad-Ms-Total-Sorted_20k_count_MissingLibrary_1_HK2GNBBXX/bamtofastq_S1_L006_R1_001.fastq.gz FASTQ/AdultMouseRep3_S1_L001_R1_001.fastq.gz
mv FASTQtmp/Ad-Ms-Total-Sorted_20k_count_MissingLibrary_1_HK2GNBBXX/bamtofastq_S1_L006_R2_001.fastq.gz FASTQ/AdultMouseRep3_S1_L001_R2_001.fastq.gz
```

* The sample annotation file ([`spermatogenesis_loupe_celltypes_repl3.csv`](spermatogenesis_mouse/data_preprocessing/spermatogenesis_loupe_celltypes_repl3.csv)) was obtained by downloading the loupe file `Mouse Unselected Spermatogenic cells.cloupe` from [here](https://data.mendeley.com/datasets/kxd5f8vpt4/1#file-fe79c10b-c42e-472e-9c7e-9a9873d9b3d8), loading it into the 10x Genomics Loupe browser and exporting the "Cell Type" labels. The resulting file was subset to only include replicate 3 and cell types with at least 5 cells, using the following R code:

```
suppressPackageStartupMessages({
  library(dplyr)
})

csv <- read.delim("Spermatogenesis_Loupe_CellTypes.csv", 
                  header = TRUE, as.is = TRUE, sep = ",")

## Only Replicate 3
csv <- csv %>% dplyr::filter(grepl("-3", Barcode)) %>%
  dplyr::mutate(index = gsub("-3", "", Barcode)) %>%
  dplyr::mutate(clusters = Cell.Types,
                clusters_coarse = Cell.Types) %>%
  dplyr::select(index, clusters_coarse, clusters)

## Remove rare cell types (<5 cells)
tbl <- table(csv$clusters)
kp <- names(tbl)[tbl >= 5]
csv <- csv %>% dplyr::filter(clusters %in% kp)

write.table(csv, file = "spermatogenesis_loupe_celltypes_repl3.csv",
            row.names = FALSE, col.names = TRUE,
            sep = ",", quote = FALSE)
```