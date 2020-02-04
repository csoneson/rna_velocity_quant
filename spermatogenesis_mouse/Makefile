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
chrnames := /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/refstar_spliced/star2_7_0c_GRCm38.primary_assembly_gencodeM21_spliced_sjdb100/chrName.txt
starindex := /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/refstar_spliced/star2_7_3a_GRCm38.primary_assembly_gencodeM21_spliced_sjdb150

bustools := /tungstenfs/groups/gbioinfo/sonechar/Software/bustools_0.39.3/bustools
subsetbam := /tungstenfs/groups/gbioinfo/sonechar/Software/subset-bam-1.0-x86_64-linux/subset-bam

datadir := /tungstenfs/groups/gbioinfo/sonechar/Data/GSE109033_sperm_mouse/GSM2928341/FASTQ
topdir := /tungstenfs/groups/gbioinfo/sonechar/Projects/alevin_velocity/spermatogenesis_mouse
scriptdir := /tungstenfs/groups/gbioinfo/sonechar/Projects/alevin_velocity/scripts
## file with assigned cell types
cellfile := /tungstenfs/groups/gbioinfo/sonechar/Data/GSE109033_sperm_mouse/GSM2928341/spermatogenesis_loupe_celltypes_repl3.csv
samplename := AdultMouseRep3
## prefix for plot output files. Should also be the suffix of the scVelo.py file
outprefix := spermatogenesis
## dataset label for plots
datalabel := Spermatogenesis

flanklength := 99

ncores := 24

genestoplot := Id4 Gfra1 Zbtb16 Stra8 Magea4 Rhox13 Sycp3 Dmc1 Piwil1 Pgk2 Acr Gapdhs Prm1 \
Gm498 Spaca9 1700008N11Rik Gm11978 Gm5893

methods := velocyto $(foreach k,exclude include,$(foreach m,prepref,$(foreach v,collapse separate,kallisto_bustools_$(m)_iso$(v)_$(k)))) \
$(foreach m,prepref,$(foreach v,collapse separate,alevin_$(m)_iso$(v)_cdna_introns)) \
$(foreach m,prepref,$(foreach v,collapse separate,alevin_$(m)_iso$(v)_cdna_introns_decoy)) \
alevin_spliced_unspliced starsolo starsolo_subtr \
$(foreach m,prepref,$(foreach v,collapse separate,alevin_$(m)_iso$(v)_cdna_introns_gentrome)) \
$(foreach m,prepref,$(foreach v,collapse separate,alevin_$(m)_iso$(v)_cdna_introns_decoy_gentrome)) \
alevin_spliced_unspliced_gentrome kb_python_lamanno
methodsc := $(subst $(space),$(comma),$(methods))

methodssel := velocyto $(foreach k,exclude include,$(foreach m,prepref,$(foreach v,collapse separate,kallisto_bustools_$(m)_iso$(v)_$(k)))) \
starsolo starsolo_subtr \
$(foreach m,prepref,$(foreach v,collapse separate,alevin_$(m)_iso$(v)_cdna_introns_gentrome)) \
$(foreach m,prepref,$(foreach v,collapse separate,alevin_$(m)_iso$(v)_cdna_introns_decoy_gentrome)) \
alevin_spliced_unspliced_gentrome kb_python_lamanno
methodsselc := $(subst $(space),$(comma),$(methodssel))

all: ref quantific veloc plot

ref: $(foreach m,prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)/cDNA_introns.fa.gz)) \
$(foreach m,prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)/introns_cDNA.fa.gz)) \
$(foreach m,prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)/tr2g_modified.tsv)) \
$(foreach m,prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)_cdna_introns.kidx)) \
reference/kb_python_lamanno/kb_python_lamanno.kbidx \
$(foreach m,prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)_cdna_introns.sidx/seq.bin)) \
$(foreach m,prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)_cdna_introns_gentrome.sidx/seq.bin)) \
$(foreach m,prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)_cdna_intronsasdecoy.sidx/seq.bin)) \
$(foreach m,prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)_cdna_intronsasdecoy_gentrome.sidx/seq.bin)) \
$(foreach m,prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)_introns_cdnaasdecoy.sidx/seq.bin)) \
$(foreach m,prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)_introns_cdnaasdecoy_gentrome.sidx/seq.bin)) \
reference/spliced_unspliced.sidx/seq.bin \
reference/spliced_unspliced_gentrome.sidx/seq.bin \
reference/spliced.sidx/seq.bin \
reference/spliced_gentrome.sidx/seq.bin \
$(foreach m,prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)_uniqueness.txt)) \
$(foreach m,prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)_uniqueness_overall.txt))

quantific: $(foreach m,prepref,$(foreach v,collapse separate,quants/kallisto_bustools_$(m)_iso$(v)_cdna_introns/run_info.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/kallisto_bustools_$(m)_iso$(v)_cdna_introns/output.correct.sort.bus)) \
$(foreach s,spliced unspliced,$(foreach d,exclude include,$(foreach m,prepref,$(foreach v,collapse separate,quants/kallisto_bustools_$(m)_iso$(v)_cdna_introns/$(s).$(d).bus)))) \
$(foreach s,spliced unspliced,$(foreach d,exclude include,$(foreach m,prepref,$(foreach v,collapse separate,quants/kallisto_bustools_$(m)_iso$(v)_cdna_introns/$(s).$(d).mtx)))) \
quants/kb_python_lamanno/inspect.unspliced.json \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_introns/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_introns_gentrome/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_intronsasdecoy/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_intronsasdecoy_gentrome/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_introns_cdnaasdecoy/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_introns_cdnaasdecoy_gentrome/lib_format_counts.json)) \
quants/alevin_spliced_unspliced/lib_format_counts.json \
quants/alevin_spliced_unspliced_gentrome/lib_format_counts.json \
quants/alevin_spliced/lib_format_counts.json \
quants/alevin_spliced_gentrome/lib_format_counts.json \
quants/cellranger/$(samplename)/outs/web_summary.html \
quants/cellranger/$(samplename)/velocyto/$(samplename).loom \
quants/cellranger/$(samplename)/outs/possorted_genome_bam_retcells.bam \
quants/cellranger/$(samplename)/outs/possorted_genome_bam_retcells.bw \
quants/cellranger/$(samplename)/outs/possorted_genome_bam_retcells_plusstrand.bw \
quants/cellranger/$(samplename)/outs/possorted_genome_bam_retcells_minusstrand.bw \
quants/starsolo/Solo.out/Gene/Summary.csv \
output/sce/sces_spliced_unspliced.rds \
$(foreach x,$(methods) alevin_spliced alevin_spliced_gentrome,output/sce/sce_$(x).rds) \
output/sce/retained_cell_barcodes.csv

veloc: $(foreach x,$(methods),output/anndata/anndata_$(x).h5ad) \
$(foreach x,$(methods),output/anndata_with_velocity/anndata_$(x)_with_velocity.h5ad) \
output/genesets/shared_genes.txt \
$(foreach x,$(methodssel),output/anndata_with_velocity/anndata_$(x)_shared_genes_with_velocity.h5ad)

plot: plots/velocity_confidence/$(outprefix)_plot_velocity_confidence.rds \
plots/count_correlation/$(outprefix)_plot_count_correlation.rds \
plots/total_count/$(outprefix)_plot_total_count.rds \
plots/spliced_unspliced_mds/$(outprefix)_plot_spliced_unspliced_mds.rds \
plots/frac_unspliced/$(outprefix)_plot_frac_unspliced.rds \
$(foreach g,$(genestoplot),plots/genes/$(outprefix)_$(g).pdf) \
plots/merge_umaps/$(outprefix)_plot_merge_umaps.rds \
plots/compare_velocity_projections/$(outprefix)_plot_compare_velocity_projections.rds \
plots/shared_genes_scvelo/$(outprefix)_plot_shared_genes_scvelo.rds \
plots/velocity_mds/$(outprefix)_plot_velocity_mds.rds \
$(foreach g,$(genestoplot),plots/fracunspliced_velocity_singlegene/$(outprefix)_$(g)_velocity.png) \
plots/velocity_confidence_embedding/$(outprefix)_plot_velocity_confidence_embedding.rds
#plots/compare_dimred_velocities/$(outprefix)_plot_compare_dimred_velocities.rds \

## -----------------------------------------------------------------------------------------------------
## Standard CellRanger + velocyto workflow
## -----------------------------------------------------------------------------------------------------
## Run CellRanger
quants/cellranger/$(samplename)/outs/web_summary.html: $(ref10x)/reference.json \
$(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz
	mkdir -p quants/cellranger
	cd quants/cellranger && module load CellRanger/3.0.2 && \
	cellranger count --id=$(samplename) --fastqs=$(datadir) \
	--sample=$(samplename) --nosecondary --transcriptome=$(ref10x) \
	--localcores=$(ncores) --localmem=64 --lanes=1

## Run velocyto
quants/cellranger/$(samplename)/velocyto/$(samplename).loom: \
quants/cellranger/$(samplename)/outs/web_summary.html $(gtf)
	module purge && module load Anaconda3/5.3.0 && \
	source activate velocyto_0.17 && \
	velocyto run10x --samtools-threads $(ncores) quants/cellranger/$(samplename) $(gtf) && \
	source deactivate

## Subset bam file to only reads from retained cell barcodes
quants/cellranger/$(samplename)/outs/possorted_genome_bam_retcells.bam: \
quants/cellranger/$(samplename)/outs/web_summary.html \
output/sce/retained_cell_barcodes.csv
	$(subsetbam) --bam quants/cellranger/$(samplename)/outs/possorted_genome_bam.bam \
	--cell-barcodes output/sce/retained_cell_barcodes.csv \
	--out-bam $@ --cores $(ncores) --log-level info --bam-tag CB

## Generate bigwig file
quants/cellranger/$(samplename)/outs/possorted_genome_bam_retcells.bw: \
quants/cellranger/$(samplename)/outs/possorted_genome_bam_retcells.bam $(chrlengths)
	module load BEDTools/2.27.1-foss-2018b && \
	bedtools genomecov -split -ibam $< -bg | \
	LC_COLLATE=C sort -k1,1 -k2,2n > $(samplename).tmp.bedGraph
	module load Kent_tools/20190212-linux.x86_64 && \
	bedGraphToBigWig $(samplename).tmp.bedGraph $(chrlengths) $@
	rm -f $(samplename).tmp.bedGraph

## Separate strands
quants/cellranger/$(samplename)/outs/possorted_genome_bam_retcells_plusstrand.bw: \
quants/cellranger/$(samplename)/outs/possorted_genome_bam_retcells.bam $(chrlengths)
	module load BEDTools/2.27.1-foss-2018b && \
	bedtools genomecov -split -ibam $< -bg -strand '+' | \
	LC_COLLATE=C sort -k1,1 -k2,2n > $(samplename).pos.tmp.bedGraph
	module load Kent_tools/20190212-linux.x86_64 && \
	bedGraphToBigWig $(samplename).pos.tmp.bedGraph $(chrlengths) $@
	rm -f $(samplename).pos.tmp.bedGraph

quants/cellranger/$(samplename)/outs/possorted_genome_bam_retcells_minusstrand.bw: \
quants/cellranger/$(samplename)/outs/possorted_genome_bam_retcells.bam $(chrlengths)
	module load BEDTools/2.27.1-foss-2018b && \
	bedtools genomecov -split -ibam $< -bg -strand '-' | \
	LC_COLLATE=C sort -k1,1 -k2,2n > $(samplename).neg.tmp.bedGraph
	module load Kent_tools/20190212-linux.x86_64 && \
	bedGraphToBigWig $(samplename).neg.tmp.bedGraph $(chrlengths) $@
	rm -f $(samplename).neg.tmp.bedGraph

## -----------------------------------------------------------------------------------------------------
## STARsolo
## -----------------------------------------------------------------------------------------------------
quants/starsolo/Solo.out/Gene/Summary.csv: \
$(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz \
reference/cellranger2.2.0-737K-august-2016.txt
	mkdir -p quants/starsolo
	module load STAR/2.7.3a-foss-2019b && \
	STAR --genomeDir $(starindex) --readFilesIn $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz $(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz \
	--soloFeatures Gene GeneFull Velocyto --soloType CB_UMI_Simple --soloCBwhitelist reference/cellranger2.2.0-737K-august-2016.txt \
	--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --soloStrand Forward \
	--outFileNamePrefix quants/starsolo/ --runThreadN $(ncores) --readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI CR UR CB UB GX GN \
	--limitBAMsortRAM 58140431943


## -----------------------------------------------------------------------------------------------------
## Generate fasta file with transcripts + introns (with some flanking sequence)
## -----------------------------------------------------------------------------------------------------
## Extracts intronic sequences flanked by L-1 bases of exonic sequences where L is the biological 
## read length (for this data set, 91)
define farefrule
reference/$(2)_iso$(1)/cDNA_introns.fa.gz: $(genome) $(gtf) \
$(scriptdir)/generate_cdna_intron_fa_$(2).R $(3) $(4)
	mkdir -p Rout
	mkdir -p $$(@D)
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args scriptdir='$(scriptdir)' gtf='$(gtf)' genome='$(genome)' isoform_action='$(1)' flanklength=$(flanklength) outdir='$$(@D)'" $(scriptdir)/generate_cdna_intron_fa_$(2).R Rout/generate_cdna_intron_fa_$(2)_$(1).Rout
endef
$(foreach v,collapse separate,$(eval $(call farefrule,$(v),prepref,$(scriptdir)/extractIntronSeqs.R,$(scriptdir)/extractTxSeqs.R)))

## Swap order of cDNAs and introns, to allow quantifying introns with cDNAs as decoy sequences
define swaprule
reference/$(2)_iso$(1)/introns_cDNA.fa.gz: \
reference/$(2)_iso$(1)/cDNA_introns.fa.gz \
$(scriptdir)/swap_order_fasta_records.R
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args infa='$$<' outfa='$$@'" $(scriptdir)/swap_order_fasta_records.R Rout/swap_order_fasta_records_$(2)_$(1).Rout
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call swaprule,$(v),$(m)))))

## Add suffix to gene ID for introns in transcript-to-gene mapping (to avoid both 
## spliced and unspliced reads being assigned to the same gene ID and thus be indistinguishable)
define suffixrule
reference/$(2)_iso$(1)/tr2g_modified.tsv: \
reference/$(2)_iso$(1)/cDNA_introns.fa.gz \
$(scriptdir)/identify_unspliced_genes_in_tx2gene.R
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args intsv='$$(<D)/tr2g.tsv' outtsv='$$@'" $(scriptdir)/identify_unspliced_genes_in_tx2gene.R Rout/identify_unspliced_genes_in_tx2gene_$(2)_$(1).Rout
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call suffixrule,$(v),$(m)))))

## Generate fasta file with spliced + unspliced transcripts
reference/gencode.vM21.transcripts.spliced.and.unspliced.fa: \
$(gtf) $(genome) $(tx2genefull) $(scriptdir)/generate_spliced_unspliced_fa.R
	mkdir -p Rout
	mkdir -p reference
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args scriptdir='$(scriptdir)' gtf='$(gtf)' genome='$(genome)' tx2gene='$(tx2genefull)' outfa='$@' outt2g='$(@D)/gencode.vM21.spliced.unspliced.tx2gene.txt'" $(scriptdir)/generate_spliced_unspliced_fa.R Rout/generate_spliced_unspliced_fa.Rout

reference/gencode.vM21.spliced.unspliced.tx2gene.txt: \
reference/gencode.vM21.transcripts.spliced.and.unspliced.fa
	touch $@

## Get uniqueness for each (spliced and unspliced) gene
define uniquenessrule
reference/$(2)_iso$(1)_uniqueness.txt: reference/$(2)_iso$(1)/cDNA_introns.fa.gz reference/$(2)_iso$(1)/tr2g_modified.tsv \
$(scriptdir)/frac_unique.py
	module purge && module load Anaconda3/5.3.0 && \
	source activate velocyto_0.17 && \
	python $(scriptdir)/frac_unique.py <(zcat reference/$(2)_iso$(1)/cDNA_introns.fa.gz) reference/$(2)_iso$(1)/tr2g_modified.tsv $(flanklength) > $$@ && \
	source deactivate && module purge
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call uniquenessrule,$(v),$(m)))))

define uniquenessruleoverall
reference/$(2)_iso$(1)_uniqueness_overall.txt: reference/$(2)_iso$(1)/cDNA_introns.fa.gz reference/$(2)_iso$(1)/cDNA_introns.fa.gz \
$(scriptdir)/frac_unique.py
	module purge && module load Anaconda3/5.3.0 && \
	source activate velocyto_0.17 && \
	python $(scriptdir)/frac_unique.py <(zcat reference/$(2)_iso$(1)/cDNA_introns.fa.gz) reference/$(2)_iso$(1)/tr2g.tsv $(flanklength) > $$@ && \
	source deactivate && module purge
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call uniquenessruleoverall,$(v),$(m)))))

## -----------------------------------------------------------------------------------------------------
## Build indexes
## -----------------------------------------------------------------------------------------------------
## Build kallisto index
define kidxrule
reference/$(2)_iso$(1)_cdna_introns.kidx: \
reference/$(2)_iso$(1)/cDNA_introns.fa.gz
	module load kallisto/0.46.0-foss-2019b && \
	kallisto index -k 31 -i $$@ $$<
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call kidxrule,$(v),$(m)))))

## Build kallisto index with kb-python
reference/kb_python_lamanno/kb_python_lamanno.kbidx: \
$(genome) $(gtf)
	mkdir -p $(@D)
	module purge && module load Anaconda3/5.3.0 && \
	source activate kb-python-0.24.4 && \
	kb ref --workflow lamanno -i $@ -g $(@D)/t2g.tsv \
	-f1 $(@D)/kb_python_lamanno_cdna.fa -f2 $(@D)/kb_python_lamanno_intron.fa \
	-c1 $(@D)/cdna_tx_to_capture.txt -c2 $(@D)/intron_tx_to_capture.txt $(genome) $(gtf) && \
	source deactivate && module purge

## Build Salmon index
define sidxrule1
reference/$(2)_iso$(1)_cdna_introns.sidx/seq.bin: \
reference/$(2)_iso$(1)/cDNA_introns.fa.gz
	module load Salmon/1.0.0 && \
	salmon index -k 31 -i $$(@D) -t $$< --gencode -p $(ncores) --type puff 
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call sidxrule1,$(v),$(m)))))

## Build Salmon index with genome decoy
define sidxrule1genome
reference/$(2)_iso$(1)_cdna_introns_gentrome.sidx/seq.bin: \
reference/$(2)_iso$(1)/cDNA_introns.fa.gz $(genome) $(chrnames)
	module load Salmon/1.0.0 && \
	salmon index -k 31 -i $$(@D) -t <(cat <(zcat $$<) $(genome)) \
	--gencode -p $(ncores) --type puff -d $(chrnames) 
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call sidxrule1genome,$(v),$(m)))))

## cDNAs, introns as decoys
define sidxrule2
reference/$(2)_iso$(1)_cdna_intronsasdecoy.sidx/seq.bin: \
reference/$(2)_iso$(1)/cDNA_introns.fa.gz
	module load Salmon/1.0.0 && \
	salmon index -k 31 -i $$(@D) -t $$< --gencode -p $(ncores) --type puff \
	-d $$(<D)/introns_tx_to_capture.txt
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call sidxrule2,$(v),$(m)))))

## cDNAs, introns + genome as decoys
define sidxrule2genome
reference/$(2)_iso$(1)_cdna_intronsasdecoy_gentrome.sidx/seq.bin: \
reference/$(2)_iso$(1)/cDNA_introns.fa.gz $(genome) $(chrnames)
	module load Salmon/1.0.0 && \
	salmon index -k 31 -i $$(@D) -t <(cat <(zcat $$<) $(genome)) \
	--gencode -p $(ncores) --type puff \
	-d <(cat $$(<D)/introns_tx_to_capture.txt $(chrnames))
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call sidxrule2genome,$(v),$(m)))))

## introns, cDNAs as decoys
define sidxrule3
reference/$(2)_iso$(1)_introns_cdnaasdecoy.sidx/seq.bin: \
reference/$(2)_iso$(1)/introns_cDNA.fa.gz
	module load Salmon/1.0.0 && \
	salmon index -k 31 -i $$(@D) -t $$< --gencode -p $(ncores) --type puff \
	-d $$(<D)/cDNA_tx_to_capture.txt
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call sidxrule3,$(v),$(m)))))

## introns, cDNAs + genome as decoys
define sidxrule3genome
reference/$(2)_iso$(1)_introns_cdnaasdecoy_gentrome.sidx/seq.bin: \
reference/$(2)_iso$(1)/introns_cDNA.fa.gz $(genome) $(chrnames)
	module load Salmon/1.0.0 && \
	salmon index -k 31 -i $$(@D) -t <(cat <(zcat $$<) $(genome)) \
	--gencode -p $(ncores) --type puff \
	-d <(cat $$(<D)/cDNA_tx_to_capture.txt $(chrnames))
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call sidxrule3genome,$(v),$(m)))))

## From spliced + unspliced annotated transcripts
reference/spliced_unspliced.sidx/seq.bin: \
reference/gencode.vM21.transcripts.spliced.and.unspliced.fa
	module load Salmon/1.0.0 && \
	salmon index -k 31 -i $(@D) -t $< --gencode -p $(ncores) --type puff 

## From spliced + unspliced annotated transcripts, genome as decoys
reference/spliced_unspliced_gentrome.sidx/seq.bin: \
reference/gencode.vM21.transcripts.spliced.and.unspliced.fa $(genome) $(chrnames)
	module load Salmon/1.0.0 && \
	salmon index -k 31 -i $(@D) -t <(cat $< $(genome)) \
	--gencode -p $(ncores) --type puff -d $(chrnames) 

## From only the regular (spliced) transcripts, for comparison
reference/spliced.sidx/seq.bin: $(txome)
	module load Salmon/1.0.0 && \
	salmon index -k 31 -i $(@D) -t $< --gencode -p $(ncores) --type puff 

## From only the regular (spliced) transcripts, for comparison, genome as decoys
reference/spliced_gentrome.sidx/seq.bin: $(txome)
	module load Salmon/1.0.0 && \
	salmon index -k 31 -i $(@D) -t <(cat $< $(genome)) \
	--gencode -p $(ncores) --type puff  -d $(chrnames)

## -----------------------------------------------------------------------------------------------------
## kallisto/bustools
## -----------------------------------------------------------------------------------------------------
## Run kallisto bus
define kallistobusrule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/run_info.json: \
reference/$(2)_iso$(1)_cdna_introns.kidx \
$(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz
	mkdir -p $$(@D)
	module load kallisto/0.46.0-foss-2019b && \
	kallisto bus -i $$< -o $$(@D) -x 10xv2 -t $(ncores) \
	$(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call kallistobusrule,$(v),$(m)))))

## Correct and sort output
define kallistosortrule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/output.correct.sort.bus: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/run_info.json \
reference/cellranger2.2.0-737K-august-2016.txt
	$(bustools) correct -w reference/cellranger2.2.0-737K-august-2016.txt -p $$(<D)/output.bus | \
	$(bustools) sort -o $$@ -t $(ncores) -
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call kallistosortrule,$(v),$(m)))))

define kallistosplicedexcluderule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/spliced.exclude.bus: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/output.correct.sort.bus
	$(bustools) capture -s -x -o $$@ \
	-c reference/$(2)_iso$(1)/introns_tx_to_capture.txt \
	-e $$(<D)/matrix.ec -t $$(<D)/transcripts.txt $$<
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call kallistosplicedexcluderule,$(v),$(m)))))

define kallistounsplicedexcluderule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/unspliced.exclude.bus: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/output.correct.sort.bus
	$(bustools) capture -s -x -o $$@ \
	-c reference/$(2)_iso$(1)/cDNA_tx_to_capture.txt \
	-e $$(<D)/matrix.ec -t $$(<D)/transcripts.txt $$<
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call kallistounsplicedexcluderule,$(v),$(m)))))

define kallistosplicedincluderule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/spliced.include.bus: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/output.correct.sort.bus
	$(bustools) capture -s -o $$@ \
	-c reference/$(2)_iso$(1)/cDNA_tx_to_capture.txt \
	-e $$(<D)/matrix.ec -t $$(<D)/transcripts.txt $$<
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call kallistosplicedincluderule,$(v),$(m)))))

define kallistounsplicedincluderule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/unspliced.include.bus: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/output.correct.sort.bus
	$(bustools) capture -s -o $$@ \
	-c reference/$(2)_iso$(1)/introns_tx_to_capture.txt \
	-e $$(<D)/matrix.ec -t $$(<D)/transcripts.txt $$<
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call kallistounsplicedincluderule,$(v),$(m)))))

## Generate spliced and unspliced count matrices
define kallistocountrule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/$(4).$(3).mtx: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/$(4).$(3).bus
	$(bustools) count -o $$(<D)/$(4).$(3) -g reference/$(2)_iso$(1)/tr2g.tsv \
	-e $$(<D)/matrix.ec -t $$(<D)/transcripts.txt --genecounts $$<
endef
$(foreach s,spliced unspliced,$(foreach d,exclude include,$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call kallistocountrule,$(v),$(m),$(d),$(s)))))))

## kb-python
quants/kb_python_lamanno/inspect.unspliced.json: \
reference/kb_python_lamanno/kb_python_lamanno.kbidx \
$(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz
	mkdir -p $(@D)
	module purge && module load Anaconda3/5.3.0 && \
	source activate kb-python-0.24.4 && \
	kb count -o $(@D) -w reference/cellranger2.2.0-737K-august-2016.txt -t $(ncores) --workflow lamanno \
	-i $< -g $(<D)/t2g.tsv -x 10XV2 -c1 $(<D)/cdna_tx_to_capture.txt -c2 $(<D)/intron_tx_to_capture.txt \
	$(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz && \
	source deactivate && module purge

## -----------------------------------------------------------------------------------------------------
## Salmon/alevin
## -----------------------------------------------------------------------------------------------------
## cDNA + introns
define alevincdnaintronrule
quants/alevin_$(2)_iso$(1)_cdna_introns$(3)/lib_format_counts.json: \
reference/$(2)_iso$(1)_cdna_introns$(3).sidx/seq.bin \
reference/$(2)_iso$(1)/tr2g_modified.tsv \
$(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz
	module load Salmon/1.0.0 && \
	salmon alevin -l ISR -i $$(<D) \
	-1 $(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz \
	-o $$(@D) -p $(ncores) --tgMap $$(word 2,$$^) \
	--chromium --dumpFeatures --dumpBfh --expectCells 1850
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call alevincdnaintronrule,$(v),$(m),))))
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call alevincdnaintronrule,$(v),$(m),_gentrome))))

define alevincdnaintrondecoyrule
quants/alevin_$(2)_iso$(1)_cdna_intronsasdecoy$(3)/lib_format_counts.json: \
reference/$(2)_iso$(1)_cdna_intronsasdecoy$(3).sidx/seq.bin \
reference/$(2)_iso$(1)/tr2g_modified.tsv \
$(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz
	module load Salmon/1.0.0 && \
	salmon alevin -l ISR -i $$(<D) \
	-1 $(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz \
	-o $$(@D) -p $(ncores) --tgMap $$(word 2,$$^) \
	--chromium --dumpFeatures --dumpBfh --expectCells 1850
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call alevincdnaintrondecoyrule,$(v),$(m),))))
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call alevincdnaintrondecoyrule,$(v),$(m),_gentrome))))

define alevinintroncdnadecoyrule
quants/alevin_$(2)_iso$(1)_introns_cdnaasdecoy$(3)/lib_format_counts.json: \
reference/$(2)_iso$(1)_introns_cdnaasdecoy$(3).sidx/seq.bin \
reference/$(2)_iso$(1)/tr2g_modified.tsv \
$(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz
	module load Salmon/1.0.0 && \
	salmon alevin -l ISR -i $$(<D) \
	-1 $(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz \
	-o $$(@D) -p $(ncores) --tgMap $$(word 2,$$^) \
	--chromium --dumpFeatures --dumpBfh --expectCells 1850
endef
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call alevinintroncdnadecoyrule,$(v),$(m),))))
$(foreach m,prepref,$(foreach v,collapse separate,$(eval $(call alevinintroncdnadecoyrule,$(v),$(m),_gentrome))))

## spliced + unspliced transcripts
define alevinsplicedunsplicedrule
quants/alevin_spliced_unspliced$(1)/lib_format_counts.json: \
reference/spliced_unspliced$(1).sidx/seq.bin \
reference/gencode.vM21.spliced.unspliced.tx2gene.txt \
$(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz
	module load Salmon/1.0.0 && \
	salmon alevin -l ISR -i $$(<D) \
	-1 $(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz \
	-o $$(@D) -p $(ncores) --tgMap $$(word 2,$$^) \
	--chromium --dumpFeatures --dumpBfh --expectCells 1850
endef
$(eval $(call alevinsplicedunsplicedrule,))
$(eval $(call alevinsplicedunsplicedrule,_gentrome))

## spliced transcripts only (regular annotation)
define alevinsplicedrule
quants/alevin_spliced$(1)/lib_format_counts.json: \
reference/spliced$(1).sidx/seq.bin \
$(tx2gene) \
$(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz
	module load Salmon/1.0.0 && \
	salmon alevin -l ISR -i $$(<D) \
	-1 $(datadir)/$(samplename)_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/$(samplename)_S1_L001_R2_001.fastq.gz \
	-o $$(@D) -p $(ncores) --tgMap $$(word 2,$$^) \
	--chromium --dumpFeatures --dumpBfh --expectCells 1850
endef
$(eval $(call alevinsplicedrule,))
$(eval $(call alevinsplicedrule,_gentrome))

## -----------------------------------------------------------------------------------------------------
## Analysis
## -----------------------------------------------------------------------------------------------------
## Create SCEs with shared genes/cells
output/sce/sces_spliced_unspliced.rds: $(scriptdir)/sce_helpers.R $(scriptdir)/summarize_spliced_unspliced.R \
$(foreach s,spliced unspliced,$(foreach d,exclude include,$(foreach m,prepref,$(foreach v,collapse separate,quants/kallisto_bustools_$(m)_iso$(v)_cdna_introns/$(s).$(d).mtx)))) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_introns/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_introns_gentrome/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_intronsasdecoy/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_intronsasdecoy_gentrome/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_introns_cdnaasdecoy/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_introns_cdnaasdecoy_gentrome/lib_format_counts.json)) \
quants/alevin_spliced_unspliced/lib_format_counts.json \
quants/alevin_spliced_unspliced_gentrome/lib_format_counts.json \
quants/alevin_spliced/lib_format_counts.json \
quants/alevin_spliced_gentrome/lib_format_counts.json \
quants/cellranger/$(samplename)/velocyto/$(samplename).loom \
quants/starsolo/Solo.out/Gene/Summary.csv \
quants/kb_python_lamanno/inspect.unspliced.json \
$(cellfile)
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' dataset='$(datalabel)' helperscript='$(scriptdir)/sce_helpers.R' tx2gene='$(tx2genefull)' cellfile='$(cellfile)' samplename='$(samplename)' outrds='$@'" $(scriptdir)/summarize_spliced_unspliced.R Rout/summarize_spliced_unspliced.Rout

define scerule
output/sce/sce_$(1).rds: output/sce/sces_spliced_unspliced.rds
	touch $$@
endef
$(foreach x,$(methods) alevin_spliced alevin_spliced_gentrome,$(eval $(call scerule,$(x))))

output/sce/retained_cell_barcodes.csv: output/sce/sces_spliced_unspliced.rds
	touch $@

plots/total_count/$(outprefix)_plot_total_count.rds: output/sce/sces_spliced_unspliced.rds \
$(scriptdir)/plot_total_count.R $(tx2genefull) $(scriptdir)/plot_helpers.R \
$(foreach v,separate collapse,reference/prepref_iso$(v)_uniqueness.txt) \
$(foreach v,separate collapse,reference/prepref_iso$(v)_uniqueness_overall.txt)
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' refdir='reference' plothelperscript='$(scriptdir)/plot_helpers.R' tx2gene='$(tx2genefull)' methods='$(methodsselc)' outrds='$@'" $(scriptdir)/plot_total_count.R Rout/plot_total_count.Rout

plots/spliced_unspliced_mds/$(outprefix)_plot_spliced_unspliced_mds.rds: output/sce/sces_spliced_unspliced.rds \
$(scriptdir)/plot_spliced_unspliced_mds.R $(scriptdir)/plot_helpers.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' plothelperscript='$(scriptdir)/plot_helpers.R' methods='$(methodsc)' outrds='$@'" $(scriptdir)/plot_spliced_unspliced_mds.R Rout/plot_spliced_unspliced_mds.Rout

plots/frac_unspliced/$(outprefix)_plot_frac_unspliced.rds: output/sce/sces_spliced_unspliced.rds \
$(scriptdir)/plot_frac_unspliced.R $(scriptdir)/plot_helpers.R $(tx2genefull) \
$(foreach v,separate collapse,reference/prepref_iso$(v)_uniqueness.txt) \
$(foreach v,separate collapse,reference/prepref_iso$(v)_uniqueness_overall.txt)
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' refdir='reference' plothelperscript='$(scriptdir)/plot_helpers.R' tx2gene='$(tx2genefull)' methods='$(methodsselc)' outrds='$@'" $(scriptdir)/plot_frac_unspliced.R Rout/plot_frac_unspliced.Rout

plots/velocity_confidence/$(outprefix)_plot_velocity_confidence.rds: $(foreach x,$(methodssel),output/anndata_with_velocity/anndata_$(x)_with_velocity.h5ad) \
$(scriptdir)/plot_velocity_confidence.R $(scriptdir)/plot_helpers.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' plothelperscript='$(scriptdir)/plot_helpers.R' dataset='$(datalabel)' methods='$(methodsselc)' outrds='$@'" $(scriptdir)/plot_velocity_confidence.R Rout/plot_velocity_confidence.Rout

plots/velocity_confidence_embedding/$(outprefix)_plot_velocity_confidence_embedding.rds: output/sce/sces_spliced_unspliced.rds \
$(foreach m,$(methodssel),output/anndata_with_velocity/anndata_$(m)_with_velocity.h5ad) \
$(foreach m,$(methodssel),output/anndata_with_velocity/anndata_$(m)_shared_genes_with_velocity.h5ad) \
$(scriptdir)/plot_velocity_confidence_embedding.R $(scriptdir)/plot_helpers.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' plothelperscript='$(scriptdir)/plot_helpers.R' methods='$(methodsselc)' outrds='$@'" $(scriptdir)/plot_velocity_confidence_embedding.R Rout/plot_velocity_confidence_embedding.Rout

plots/count_correlation/$(outprefix)_plot_count_correlation.rds: output/genesets/shared_genes.txt output/sce/sces_spliced_unspliced.rds \
$(foreach m,$(methodssel),output/anndata_with_velocity/anndata_$(m)_shared_genes_with_velocity.h5ad) \
$(scriptdir)/plot_count_correlation.R $(scriptdir)/plot_helpers.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' plothelperscript='$(scriptdir)/plot_helpers.R' methods='$(methodsselc)' genetxt='$<' outrds='$@'" $(scriptdir)/plot_count_correlation.R Rout/plot_count_correlation.Rout

plots/merge_umaps/$(outprefix)_plot_merge_umaps.rds: $(foreach m,$(methodssel),output/anndata_with_velocity/anndata_$(m)_with_velocity.h5ad) \
$(foreach m,$(methodssel),output/anndata_with_velocity/anndata_$(m)_shared_genes_with_velocity.h5ad) $(scriptdir)/plot_merge_umaps.R \
$(scriptdir)/plot_helpers.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load ImageMagick/7.0.9-8-GCCcore-8.3.0 && \
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' plothelperscript='$(scriptdir)/plot_helpers.R' dimred='UMAP' methods='$(methodsselc)' outrds='$@'" $(scriptdir)/plot_merge_umaps.R Rout/plot_merge_umaps.Rout

plots/compare_velocity_projections/$(outprefix)_plot_compare_velocity_projections.rds: \
$(foreach m,$(methodssel),output/anndata_with_velocity/anndata_$(m)_with_velocity.h5ad) \
$(foreach m,$(methodssel),output/anndata_with_velocity/anndata_$(m)_shared_genes_with_velocity.h5ad) \
$(scriptdir)/plot_compare_velocity_projections.R $(scriptdir)/plot_helpers.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' methods='$(methodsselc)' plothelperscript='$(scriptdir)/plot_helpers.R' outrds='$@'" $(scriptdir)/plot_compare_velocity_projections.R Rout/plot_compare_velocity_projections.Rout

plots/compare_dimred_velocities/$(outprefix)_plot_compare_dimred_velocities.rds: \
$(foreach m,$(methodssel),output/anndata_with_velocity/anndata_$(m)_with_velocity.h5ad) \
$(scriptdir)/plot_compare_dimred_velocities.R $(scriptdir)/plot_helpers.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' methods='$(methodsselc)' plothelperscript='$(scriptdir)/plot_helpers.R' outrds='$@'" $(scriptdir)/plot_compare_dimred_velocities.R Rout/plot_compare_dimred_velocities.Rout

plots/velocity_mds/$(outprefix)_plot_velocity_mds.rds: output/genesets/shared_genes.txt \
$(foreach m,$(methodssel),output/anndata_with_velocity/anndata_$(m)_shared_genes_with_velocity.h5ad) \
$(scriptdir)/plot_velocity_mds.R $(scriptdir)/plot_helpers.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' dataset='$(datalabel)' methods='$(methodsselc)' genetxt='$<' plothelperscript='$(scriptdir)/plot_helpers.R' outrds='$@'" $(scriptdir)/plot_velocity_mds.R Rout/plot_velocity_mds.Rout

plots/shared_genes_scvelo/$(outprefix)_plot_shared_genes_scvelo.rds: \
$(foreach m,$(methodssel),output/anndata_with_velocity/anndata_$(m)_with_velocity.h5ad) \
$(scriptdir)/plot_shared_genes_scvelo.R $(scriptdir)/plot_helpers.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' dataset='$(datalabel)' methods='$(methodsselc)' plothelperscript='$(scriptdir)/plot_helpers.R' outrds='$@'" $(scriptdir)/plot_shared_genes_scvelo.R Rout/plot_shared_genes_scvelo.Rout

define plotgenerule
plots/genes/$(outprefix)_$(1).pdf: quants/cellranger/$(samplename)/outs/possorted_genome_bam_retcells_plusstrand.bw \
quants/cellranger/$(samplename)/outs/possorted_genome_bam_retcells_minusstrand.bw \
$(gtf) $(scriptdir)/plot_helpers.R $(scriptdir)/plot_gene_model.R $(tx2genefull) \
$(foreach v,separate collapse,reference/prepref_iso$(v)_uniqueness.txt) \
$(foreach v,separate collapse,reference/prepref_iso$(v)_uniqueness_overall.txt)
	mkdir -p $$(@D)
	mkdir -p Rout
	module load ImageMagick/7.0.9-8-GCCcore-8.3.0 && \
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' refdir='reference' plothelperscript='$(scriptdir)/plot_helpers.R' tx2gene='$(tx2genefull)' methods='$(methodsselc)' gtf='$(gtf)' bigwigfileplus='$$(word 1,$$^)' bigwigfileminus='$$(word 2,$$^)' samplename='$(datalabel)' showgene='$(1)' outpdf='$$@'" $(scriptdir)/plot_gene_model.R Rout/plot_gene_model_$(1).Rout
endef
$(foreach g,$(genestoplot),$(eval $(call plotgenerule,$(g))))

define plotfuvgenerule
plots/fracunspliced_velocity_singlegene/$(outprefix)_$(1)_velocity.png: \
output/sce/sces_spliced_unspliced.rds \
$(foreach m,$(methodssel),output/anndata_with_velocity/anndata_$(m)_with_velocity.h5ad) \
$(scriptdir)/plot_fracunspliced_velocity_singlegene.R $(scriptdir)/plot_helpers.R
	mkdir -p $$(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' plothelperscript='$(scriptdir)/plot_helpers.R' methods='$(methodsselc)' showgene='$(1)' outvelopng='$$@'" $(scriptdir)/plot_fracunspliced_velocity_singlegene.R Rout/plot_fracunspliced_velocity_singlegene_$(1).Rout
endef
$(foreach g,$(genestoplot),$(eval $(call plotfuvgenerule,$(g))))

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
$(scriptdir)/run_scVelo_$(outprefix).py
	mkdir -p plots/velocity/anndata_$(1)
	mkdir -p $$(@D)
	mkdir -p pylog
	module purge && module load Anaconda3/5.3.0 && \
	set +u && source activate scvelo_0.1.24 && set -u && \
	python $(scriptdir)/run_scVelo_$(outprefix).py $$< plots/velocity $$@ "None" >pylog/$(1).scvelo.log 2>&1 && \
	set +u && source deactivate && set -u && module purge
endef
$(foreach x,$(methods),$(eval $(call scvelorule,$(x))))

## Make text file with genes retained by all methods
output/genesets/shared_genes.txt: $(scriptdir)/make_gene_subsets.R \
$(foreach x,$(methodssel),output/anndata_with_velocity/anndata_$(x)_with_velocity.h5ad)
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' methods='$(methodsselc)' outtxt='$@'" $(scriptdir)/make_gene_subsets.R Rout/make_gene_subsets.Rout

## Velocity, based on genes selected by all methods
define scveloruleshared
output/anndata_with_velocity/anndata_$(1)_shared_genes$(2)_with_velocity.h5ad: output/anndata/anndata_$(1).h5ad \
$(scriptdir)/run_scVelo_$(outprefix).py output/genesets/shared_genes$(2).txt
	mkdir -p plots/velocity/anndata_$(1)_shared_genes$(2)
	mkdir -p $$(@D)
	mkdir -p pylog
	module purge && module load Anaconda3/5.3.0 && \
	set +u && source activate scvelo_0.1.24 && set -u && \
	python $(scriptdir)/run_scVelo_$(outprefix).py $$< plots/velocity $$@ output/genesets/shared_genes$(2).txt >pylog/$(1).scvelo_shared_genes$(2).log 2>&1 && \
	set +u && source deactivate && set -u && module purge
endef
$(foreach x,$(methodssel),$(eval $(call scveloruleshared,$(x),)))
