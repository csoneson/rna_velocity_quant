SHELL := /bin/bash

comma := ,
empty :=
space := $(empty) $(empty)

genome := /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/GRCm38.primary_assembly.genome.fa
gtf := /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/gencode.vM21.annotation.gtf
ref10x := /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/ref10x_spliced/cellranger3_0_2_GRCm38.primary_assembly_gencodeM21_spliced
tx2genefull := /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/tx2gene/gencode.vM21.annotation.tx2gene_full.rds
tx2gene := /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/tx2gene/gencode.vM21.annotation.tx2gene.txt
txome := /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/gencode.vM21.transcripts.fa
chrlengths := /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/refstar_spliced/star2_7_0c_GRCm38.primary_assembly_gencodeM21_spliced_sjdb100/chrNameLength.txt
starindex := /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/refstar_spliced/star2_7_3a_GRCm38.primary_assembly_gencodeM21_spliced_sjdb150

bustools := /tungstenfs/groups/gbioinfo/sonechar/Software/bustools_0.39.3/bustools
subsetbam := /tungstenfs/groups/gbioinfo/sonechar/Software/subset-bam-1.0-x86_64-linux/subset-bam

datadir := /tungstenfs/groups/gbioinfo/sonechar/Data/GSE124952_pfc_mouse/GSM3559979/FASTQ
topdir := /tungstenfs/groups/gbioinfo/sonechar/Projects/alevin_velocity/pfc_mouse
## Use reference files built with the same flanklength
refdir := /tungstenfs/groups/gbioinfo/sonechar/Projects/alevin_velocity/dentate_gyrus_mouse/reference
scriptdir := /tungstenfs/groups/gbioinfo/sonechar/Projects/alevin_velocity/scripts

flanklength := 97

methods := velocyto $(foreach k,exclude include,$(foreach m,prepref,$(foreach v,collapse separate,kallisto_bustools_$(m)_iso$(v)_$(k)))) \
$(foreach m,prepref,$(foreach v,collapse separate,alevin_$(m)_iso$(v)_cdna_introns)) \
$(foreach m,prepref,$(foreach v,collapse separate,alevin_$(m)_iso$(v)_cdna_introns_decoy)) \
alevin_spliced_unspliced starsolo starsolo_subtr
methodsc := $(subst $(space),$(comma),$(methods))

all: quantific veloc plot

tmp: $(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_introns/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_intronsasdecoy/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_introns_cdnaasdecoy/lib_format_counts.json)) \
quants/alevin_spliced_unspliced/lib_format_counts.json \
quants/alevin_spliced/lib_format_counts.json

quantific: $(foreach m,prepref,$(foreach v,collapse separate,quants/kallisto_bustools_$(m)_iso$(v)_cdna_introns/run_info.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/kallisto_bustools_$(m)_iso$(v)_cdna_introns/output.correct.sort.bus)) \
$(foreach s,spliced unspliced,$(foreach d,exclude include,$(foreach m,prepref,$(foreach v,collapse separate,quants/kallisto_bustools_$(m)_iso$(v)_cdna_introns/$(s).$(d).bus)))) \
$(foreach s,spliced unspliced,$(foreach d,exclude include,$(foreach m,prepref,$(foreach v,collapse separate,quants/kallisto_bustools_$(m)_iso$(v)_cdna_introns/$(s).$(d).mtx)))) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_introns/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_intronsasdecoy/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_introns_cdnaasdecoy/lib_format_counts.json)) \
quants/alevin_spliced_unspliced/lib_format_counts.json \
quants/alevin_spliced/lib_format_counts.json \
quants/cellranger/PFCsample2/outs/web_summary.html \
quants/cellranger/PFCsample2/velocyto/PFCsample2.loom \
quants/cellranger/PFCsample2/outs/possorted_genome_bam_retcells.bam \
quants/cellranger/PFCsample2/outs/possorted_genome_bam_retcells.bw \
quants/cellranger/PFCsample2/outs/possorted_genome_bam_retcells_plusstrand.bw \
quants/cellranger/PFCsample2/outs/possorted_genome_bam_retcells_minusstrand.bw \
quants/starsolo/Solo.out/Gene/Summary.csv \
output/sce/sces_spliced_unspliced.rds \
$(foreach x,$(methods) alevin_spliced,output/sce/sce_$(x).rds) \
output/sce/retained_cell_barcodes.csv

veloc: $(foreach x,$(methods),output/anndata/anndata_$(x).h5ad) \
$(foreach x,$(methods),output/anndata_with_velocity/anndata_$(x)_with_velocity.h5ad) \
output/genesets/shared_genes.txt \
$(foreach x,$(methods),output/anndata_with_velocity/anndata_$(x)_shared_genes_with_velocity.h5ad)

plot: plots/velocity_confidence/pfc_plot_velocity_confidence.rds \
plots/count_correlation/pfc_plot_count_correlation.rds \
plots/total_count/pfc_plot_total_count.rds \
plots/spliced_unspliced_mds/pfc_plot_spliced_unspliced_mds.rds \
plots/frac_unspliced/pfc_plot_frac_unspliced.rds \
plots/merge_umaps/pfc_plot_merge_umaps.rds \
plots/compare_velocity_projections/pfc_plot_compare_velocity_projections.rds \
plots/shared_genes_scvelo/pfc_plot_shared_genes_scvelo.rds
#$(foreach g,$(genestoplot),plots/genes/pfc_$(g).pdf) \

## -----------------------------------------------------------------------------------------------------
## Standard CellRanger + velocyto workflow
## -----------------------------------------------------------------------------------------------------
## Run CellRanger
quants/cellranger/PFCsample2/outs/web_summary.html: $(ref10x)/reference.json \
$(datadir)/PFCsample2_S1_L001_R1_001.fastq.gz $(datadir)/PFCsample2_S1_L001_R2_001.fastq.gz
	mkdir -p quants/cellranger
	cd quants/cellranger && module load CellRanger/3.0.2 && \
	cellranger count --id=PFCsample2 --fastqs=$(datadir) \
	--sample=PFCsample2 --nosecondary --transcriptome=$(ref10x) \
	--localcores=24 --localmem=64 --lanes=1

## Run velocyto
quants/cellranger/PFCsample2/velocyto/PFCsample2.loom: \
quants/cellranger/PFCsample2/outs/web_summary.html $(gtf)
	module purge && module load Anaconda3/5.3.0 && \
	source activate velocyto_0.17 && \
	velocyto run10x --samtools-threads 16 quants/cellranger/PFCsample2 $(gtf) && \
	source deactivate

## Subset bam file to only reads from retained cell barcodes
quants/cellranger/PFCsample2/outs/possorted_genome_bam_retcells.bam: \
quants/cellranger/PFCsample2/outs/web_summary.html \
output/sce/retained_cell_barcodes.csv
	$(subsetbam) --bam quants/cellranger/PFCsample2/outs/possorted_genome_bam.bam \
	--cell-barcodes output/sce/retained_cell_barcodes.csv \
	--out-bam $@ --cores 24 --log-level info --bam-tag CB

## Generate bigwig file
quants/cellranger/PFCsample2/outs/possorted_genome_bam_retcells.bw: \
quants/cellranger/PFCsample2/outs/possorted_genome_bam_retcells.bam $(chrlengths)
	module load BEDTools/2.27.1-foss-2018b && \
	bedtools genomecov -split -ibam $< -bg | \
	LC_COLLATE=C sort -k1,1 -k2,2n > PFCsample2.tmp.bedGraph
	module load Kent_tools/20190212-linux.x86_64 && \
	bedGraphToBigWig PFCsample2.tmp.bedGraph $(chrlengths) $@
	rm -f PFCsample2.tmp.bedGraph

## Separate strands
quants/cellranger/PFCsample2/outs/possorted_genome_bam_retcells_plusstrand.bw: \
quants/cellranger/PFCsample2/outs/possorted_genome_bam_retcells.bam $(chrlengths)
	module load BEDTools/2.27.1-foss-2018b && \
	bedtools genomecov -split -ibam $< -bg -strand '+' | \
	LC_COLLATE=C sort -k1,1 -k2,2n > PFCsample2.pos.tmp.bedGraph
	module load Kent_tools/20190212-linux.x86_64 && \
	bedGraphToBigWig PFCsample2.pos.tmp.bedGraph $(chrlengths) $@
	rm -f PFCsample2.pos.tmp.bedGraph

quants/cellranger/PFCsample2/outs/possorted_genome_bam_retcells_minusstrand.bw: \
quants/cellranger/PFCsample2/outs/possorted_genome_bam_retcells.bam $(chrlengths)
	module load BEDTools/2.27.1-foss-2018b && \
	bedtools genomecov -split -ibam $< -bg -strand '-' | \
	LC_COLLATE=C sort -k1,1 -k2,2n > PFCsample2.neg.tmp.bedGraph
	module load Kent_tools/20190212-linux.x86_64 && \
	bedGraphToBigWig PFCsample2.neg.tmp.bedGraph $(chrlengths) $@
	rm -f PFCsample2.neg.tmp.bedGraph

## -----------------------------------------------------------------------------------------------------
## STARsolo
## -----------------------------------------------------------------------------------------------------
quants/starsolo/Solo.out/Gene/Summary.csv: \
$(datadir)/PFCsample2_S1_L001_R1_001.fastq.gz $(datadir)/PFCsample2_S1_L001_R2_001.fastq.gz \
reference/cellranger2.2.0-737K-august-2016.txt
	mkdir -p quants/starsolo
	module load STAR/2.7.3a-foss-2019b && \
	STAR --genomeDir $(starindex) --readFilesIn $(datadir)/PFCsample2_S1_L001_R2_001.fastq.gz $(datadir)/PFCsample2_S1_L001_R1_001.fastq.gz \
	--soloFeatures Gene GeneFull Velocyto --soloType CB_UMI_Simple --soloCBwhitelist reference/cellranger2.2.0-737K-august-2016.txt \
	--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --soloStrand Forward \
	--outFileNamePrefix quants/starsolo/ --runThreadN 48 --readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI CR UR CB UB GX GN

## -----------------------------------------------------------------------------------------------------
## kallisto/bustools
## -----------------------------------------------------------------------------------------------------
## Run kallisto bus
define kallistobusrule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/run_info.json: \
$(refdir)/$(2)_iso$(1)_cdna_introns.kidx \
$(datadir)/PFCsample2_S1_L001_R1_001.fastq.gz $(datadir)/PFCsample2_S1_L001_R2_001.fastq.gz
	mkdir -p $$(@D)
	module load kallisto/0.46.0-foss-2019b && \
	kallisto bus -i $$< -o $$(@D) -x 10xv2 -t 24 \
	$(datadir)/PFCsample2_S1_L001_R1_001.fastq.gz $(datadir)/PFCsample2_S1_L001_R2_001.fastq.gz
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call kallistobusrule,$(v),$(m)))))

## Correct and sort output
define kallistosortrule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/output.correct.sort.bus: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/run_info.json \
reference/cellranger2.2.0-737K-august-2016.txt
	$(bustools) correct -w reference/cellranger2.2.0-737K-august-2016.txt -p $$(<D)/output.bus | \
	$(bustools) sort -o $$@ -t 24 -
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call kallistosortrule,$(v),$(m)))))

define kallistosplicedexcluderule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/spliced.exclude.bus: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/output.correct.sort.bus
	$(bustools) capture -s -x -o $$@ \
	-c $(refdir)/$(2)_iso$(1)/introns_tx_to_capture.txt \
	-e $$(<D)/matrix.ec -t $$(<D)/transcripts.txt $$<
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call kallistosplicedexcluderule,$(v),$(m)))))

define kallistounsplicedexcluderule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/unspliced.exclude.bus: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/output.correct.sort.bus
	$(bustools) capture -s -x -o $$@ \
	-c $(refdir)/$(2)_iso$(1)/cDNA_tx_to_capture.txt \
	-e $$(<D)/matrix.ec -t $$(<D)/transcripts.txt $$<
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call kallistounsplicedexcluderule,$(v),$(m)))))

define kallistosplicedincluderule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/spliced.include.bus: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/output.correct.sort.bus
	$(bustools) capture -s -o $$@ \
	-c $(refdir)/$(2)_iso$(1)/cDNA_tx_to_capture.txt \
	-e $$(<D)/matrix.ec -t $$(<D)/transcripts.txt $$<
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call kallistosplicedincluderule,$(v),$(m)))))

define kallistounsplicedincluderule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/unspliced.include.bus: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/output.correct.sort.bus
	$(bustools) capture -s -o $$@ \
	-c $(refdir)/$(2)_iso$(1)/introns_tx_to_capture.txt \
	-e $$(<D)/matrix.ec -t $$(<D)/transcripts.txt $$<
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call kallistounsplicedincluderule,$(v),$(m)))))

## Generate spliced and unspliced count matrices
define kallistocountrule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/$(4).$(3).mtx: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/$(4).$(3).bus
	$(bustools) count -o $$(<D)/$(4).$(3) -g $(refdir)/$(2)_iso$(1)/tr2g.tsv \
	-e $$(<D)/matrix.ec -t $$(<D)/transcripts.txt --genecounts $$<
endef
$(foreach s,spliced unspliced,$(foreach d,exclude include,$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call kallistocountrule,$(v),$(m),$(d),$(s)))))))

## -----------------------------------------------------------------------------------------------------
## Salmon/alevin
## -----------------------------------------------------------------------------------------------------
## cDNA + introns
define alevincdnaintronrule
quants/alevin_$(2)_iso$(1)_cdna_introns/lib_format_counts.json: \
$(refdir)/$(2)_iso$(1)_cdna_introns.sidx/seq.bin \
$(refdir)/$(2)_iso$(1)/tr2g_modified.tsv \
$(datadir)/PFCsample2_S1_L001_R1_001.fastq.gz $(datadir)/PFCsample2_S1_L001_R2_001.fastq.gz
	module load Salmon/1.0.0 && \
	salmon alevin -l ISR -i $$(<D) \
	-1 $(datadir)/PFCsample2_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/PFCsample2_S1_L001_R2_001.fastq.gz \
	-o $$(@D) -p 24 --tgMap $$(word 2,$$^) \
	--chromium --dumpFeatures --dumpBfh --expectCells 1300
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call alevincdnaintronrule,$(v),$(m)))))

define alevincdnaintrondecoyrule
quants/alevin_$(2)_iso$(1)_cdna_intronsasdecoy/lib_format_counts.json: \
$(refdir)/$(2)_iso$(1)_cdna_intronsasdecoy.sidx/seq.bin \
$(refdir)/$(2)_iso$(1)/tr2g_modified.tsv \
$(datadir)/PFCsample2_S1_L001_R1_001.fastq.gz $(datadir)/PFCsample2_S1_L001_R2_001.fastq.gz
	module load Salmon/1.0.0 && \
	salmon alevin -l ISR -i $$(<D) \
	-1 $(datadir)/PFCsample2_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/PFCsample2_S1_L001_R2_001.fastq.gz \
	-o $$(@D) -p 24 --tgMap $$(word 2,$$^) \
	--chromium --dumpFeatures --dumpBfh --expectCells 1300
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call alevincdnaintrondecoyrule,$(v),$(m)))))

define alevinintroncdnadecoyrule
quants/alevin_$(2)_iso$(1)_introns_cdnaasdecoy/lib_format_counts.json: \
$(refdir)/$(2)_iso$(1)_introns_cdnaasdecoy.sidx/seq.bin \
$(refdir)/$(2)_iso$(1)/tr2g_modified.tsv \
$(datadir)/PFCsample2_S1_L001_R1_001.fastq.gz $(datadir)/PFCsample2_S1_L001_R2_001.fastq.gz
	module load Salmon/1.0.0 && \
	salmon alevin -l ISR -i $$(<D) \
	-1 $(datadir)/PFCsample2_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/PFCsample2_S1_L001_R2_001.fastq.gz \
	-o $$(@D) -p 24 --tgMap $$(word 2,$$^) \
	--chromium --dumpFeatures --dumpBfh --expectCells 1300
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call alevinintroncdnadecoyrule,$(v),$(m)))))

## spliced + unspliced transcripts
quants/alevin_spliced_unspliced/lib_format_counts.json: \
$(refdir)/spliced_unspliced.sidx/seq.bin \
$(refdir)/gencode.vM21.spliced.unspliced.tx2gene.txt \
$(datadir)/PFCsample2_S1_L001_R1_001.fastq.gz $(datadir)/PFCsample2_S1_L001_R2_001.fastq.gz
	module load Salmon/1.0.0 && \
	salmon alevin -l ISR -i $(<D) \
	-1 $(datadir)/PFCsample2_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/PFCsample2_S1_L001_R2_001.fastq.gz \
	-o $(@D) -p 24 --tgMap $(word 2,$^) \
	--chromium --dumpFeatures --dumpBfh --expectCells 1300

## spliced transcripts only (regular annotation)
quants/alevin_spliced/lib_format_counts.json: \
$(refdir)/spliced.sidx/seq.bin \
$(tx2gene) \
$(datadir)/PFCsample2_S1_L001_R1_001.fastq.gz $(datadir)/PFCsample2_S1_L001_R2_001.fastq.gz
	module load Salmon/1.0.0 && \
	salmon alevin -l ISR -i $(<D) \
	-1 $(datadir)/PFCsample2_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/PFCsample2_S1_L001_R2_001.fastq.gz \
	-o $(@D) -p 16 --tgMap $(word 2,$^) \
	--chromium --dumpFeatures --dumpBfh --expectCells 1300

## -----------------------------------------------------------------------------------------------------
## Analysis
## -----------------------------------------------------------------------------------------------------
## Create SCEs with shared genes/cells
output/sce/sces_spliced_unspliced.rds: $(scriptdir)/sce_helpers.R $(scriptdir)/summarize_spliced_unspliced.R \
$(foreach s,spliced unspliced,$(foreach d,exclude include,$(foreach m,prepref,$(foreach v,collapse separate,quants/kallisto_bustools_$(m)_iso$(v)_cdna_introns/$(s).$(d).mtx)))) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_introns/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_intronsasdecoy/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_introns_cdnaasdecoy/lib_format_counts.json)) \
quants/alevin_spliced_unspliced/lib_format_counts.json \
quants/alevin_spliced/lib_format_counts.json \
quants/cellranger/PFCsample2/velocyto/PFCsample2.loom \
quants/starsolo/Solo.out/Gene/Summary.csv
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' dataset='PFC' helperscript='$(scriptdir)/sce_helpers.R' tx2gene='$(tx2genefull)' cellfile='/work/gbioinfo/sonechar/Data/GSE124952_pfc_mouse/GSM3559979/cells_in_pfcsample2.csv' samplename='PFCsample2' outrds='$@'" $(scriptdir)/summarize_spliced_unspliced.R Rout/summarize_spliced_unspliced.Rout

define scerule
output/sce/sce_$(1).rds: output/sce/sces_spliced_unspliced.rds
	touch $$@
endef
$(foreach x,$(methods) alevin_spliced,$(eval $(call scerule,$(x))))

output/sce/retained_cell_barcodes.csv: output/sce/sces_spliced_unspliced.rds
	touch $@

plots/total_count/pfc_plot_total_count.rds: output/sce/sces_spliced_unspliced.rds \
$(scriptdir)/plot_total_count.R $(tx2genefull) $(scriptdir)/plot_helpers.R \
$(foreach v,separate collapse,$(refdir)/prepref_iso$(v)_uniqueness.txt) \
$(foreach v,separate collapse,$(refdir)/prepref_iso$(v)_uniqueness_overall.txt)
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' refdir='$(refdir)' plothelperscript='$(scriptdir)/plot_helpers.R' tx2gene='$(tx2genefull)' methods='$(methodsc)' outrds='$@'" $(scriptdir)/plot_total_count.R Rout/plot_total_count.Rout

plots/spliced_unspliced_mds/pfc_plot_spliced_unspliced_mds.rds: output/sce/sces_spliced_unspliced.rds \
$(scriptdir)/plot_spliced_unspliced_mds.R $(scriptdir)/plot_helpers.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' plothelperscript='$(scriptdir)/plot_helpers.R' methods='$(methodsc)' outrds='$@'" $(scriptdir)/plot_spliced_unspliced_mds.R Rout/plot_spliced_unspliced_mds.Rout

plots/frac_unspliced/pfc_plot_frac_unspliced.rds: output/sce/sces_spliced_unspliced.rds \
$(scriptdir)/plot_frac_unspliced.R $(scriptdir)/plot_helpers.R $(tx2genefull) \
$(foreach v,separate collapse,$(refdir)/prepref_iso$(v)_uniqueness.txt) \
$(foreach v,separate collapse,$(refdir)/prepref_iso$(v)_uniqueness_overall.txt)
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' refdir='$(refdir)' plothelperscript='$(scriptdir)/plot_helpers.R' tx2gene='$(tx2genefull)' methods='$(methodsc)' outrds='$@'" $(scriptdir)/plot_frac_unspliced.R Rout/plot_frac_unspliced.Rout

plots/velocity_confidence/pfc_plot_velocity_confidence.rds: $(foreach x,$(methods),output/anndata_with_velocity/anndata_$(x)_with_velocity.h5ad) \
$(scriptdir)/plot_velocity_confidence.R $(scriptdir)/plot_helpers.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' plothelperscript='$(scriptdir)/plot_helpers.R' dataset='PFC' methods='$(methodsc)' outrds='$@'" $(scriptdir)/plot_velocity_confidence.R Rout/plot_velocity_confidence.Rout

plots/count_correlation/pfc_plot_count_correlation.rds: output/genesets/shared_genes.txt output/sce/sces_spliced_unspliced.rds \
$(scriptdir)/plot_count_correlation.R 
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' methods='$(methodsc)' genetxt='$<' outrds='$@'" $(scriptdir)/plot_count_correlation.R Rout/plot_count_correlation.Rout

plots/merge_umaps/pfc_plot_merge_umaps.rds: $(foreach m,$(methods),output/anndata_with_velocity/anndata_$(m)_with_velocity.h5ad) \
$(foreach m,$(methods),output/anndata_with_velocity/anndata_$(m)_shared_genes_with_velocity.h5ad) $(scriptdir)/plot_merge_umaps.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load ImageMagick/7.0.9-8-GCCcore-8.3.0 && \
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' methods='$(methodsc)' outrds='$@'" $(scriptdir)/plot_merge_umaps.R Rout/plot_merge_umaps.Rout

plots/compare_velocity_projections/pfc_plot_compare_velocity_projections.rds: \
$(foreach m,$(methods),output/anndata_with_velocity/anndata_$(m)_with_velocity.h5ad) \
$(foreach m,$(methods),output/anndata_with_velocity/anndata_$(m)_shared_genes_with_velocity.h5ad) \
$(scriptdir)/plot_compare_velocity_projections.R $(scriptdir)/plot_helpers.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' methods='$(methodsc)' plothelperscript='$(scriptdir)/plot_helpers.R' outrds='$@'" $(scriptdir)/plot_compare_velocity_projections.R Rout/plot_compare_velocity_projections.Rout

plots/shared_genes_scvelo/pfc_plot_shared_genes_scvelo.rds: \
$(foreach m,$(methods),output/anndata_with_velocity/anndata_$(m)_with_velocity.h5ad) \
$(scriptdir)/plot_shared_genes_scvelo.R $(scriptdir)/plot_helpers.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' dataset='PFC' methods='$(methodsc)' plothelperscript='$(scriptdir)/plot_helpers.R' outrds='$@'" $(scriptdir)/plot_shared_genes_scvelo.R Rout/plot_shared_genes_scvelo.Rout

define plotgenerule
plots/genes/pfc_$(1).pdf: quants/cellranger/PFCsample2/outs/possorted_genome_bam_retcells_plusstrand.bw \
quants/cellranger/PFCsample2/outs/possorted_genome_bam_retcells_minusstrand.bw \
$(gtf) $(scriptdir)/plot_helpers.R $(scriptdir)/plot_gene_model.R $(tx2genefull) \
$(foreach v,separate collapse,$(refdir)/prepref_iso$(v)_uniqueness.txt) \
$(foreach v,separate collapse,$(refdir)/prepref_iso$(v)_uniqueness_overall.txt)
	mkdir -p $$(@D)
	mkdir -p Rout
	module load ImageMagick/7.0.9-8-GCCcore-8.3.0 && \
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' refdir='$(refdir)' plothelperscript='$(scriptdir)/plot_helpers.R' tx2gene='$(tx2genefull)' methods='$(methodsc)' gtf='$(gtf)' bigwigfileplus='$$(word 1,$$^)' bigwigfileminus='$$(word 2,$$^)' samplename='PFC' showgene='$(1)' outpdf='$$@'" $(scriptdir)/plot_gene_model.R Rout/plot_gene_model_$(1).Rout
endef
$(foreach g,$(genestoplot),$(eval $(call plotgenerule,$(g))))

## Convert SCE to anndata
define anndatarule
output/anndata/anndata_$(1).h5ad: output/sce/sce_$(1).rds \
$(scriptdir)/convert_sce_to_anndata.py
	mkdir -p $$(@D)
	mkdir -p pylog
	module purge && module load Anaconda3/5.3.0 && \
	set +u && source activate scvelo_0.1.24 && set -u && \
	python $(scriptdir)/convert_sce_to_anndata.py $$< $$@ >pylog/$(1).anndata.log 2>&1 && \
	set +u && source deactivate && set -u && module purge
endef
$(foreach x,$(methods),$(eval $(call anndatarule,$(x))))

## Calculate velocity with scVelo
define scvelorule
output/anndata_with_velocity/anndata_$(1)_with_velocity.h5ad: output/anndata/anndata_$(1).h5ad \
$(scriptdir)/run_scVelo_pfc.py
	mkdir -p plots/velocity/anndata_$(1)
	mkdir -p $$(@D)
	mkdir -p pylog
	module purge && module load Anaconda3/5.3.0 && \
	set +u && source activate scvelo_0.1.24 && set -u && \
	python $(scriptdir)/run_scVelo_pfc.py $$< plots/velocity $$@ "None" >pylog/$(1).scvelo.log 2>&1 && \
	set +u && source deactivate && set -u && module purge
endef
$(foreach x,$(methods),$(eval $(call scvelorule,$(x))))

## Make text file with genes retained by all methods
output/genesets/shared_genes.txt: $(scriptdir)/make_gene_subsets.R \
$(foreach x,$(methods),output/anndata_with_velocity/anndata_$(x)_with_velocity.h5ad)
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' methods='$(methodsc)' outtxt='$@'" $(scriptdir)/make_gene_subsets.R Rout/make_gene_subsets.Rout

## Velocity, based on genes selected by all methods
define scveloruleshared
output/anndata_with_velocity/anndata_$(1)_shared_genes$(2)_with_velocity.h5ad: output/anndata/anndata_$(1).h5ad \
$(scriptdir)/run_scVelo_pfc.py output/genesets/shared_genes$(2).txt
	mkdir -p plots/velocity/anndata_$(1)_shared_genes$(2)
	mkdir -p $$(@D)
	mkdir -p pylog
	module purge && module load Anaconda3/5.3.0 && \
	set +u && source activate scvelo_0.1.24 && set -u && \
	python $(scriptdir)/run_scVelo_pfc.py $$< plots/velocity $$@ output/genesets/shared_genes$(2).txt >pylog/$(1).scvelo_shared_genes$(2).log 2>&1 && \
	set +u && source deactivate && set -u && module purge
endef
$(foreach x,$(methods),$(eval $(call scveloruleshared,$(x),)))
