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

datadir := /tungstenfs/groups/gbioinfo/sonechar/Data/10xGenomics_neuron_10k_v3/neuron_10k_v3_fastqs
topdir := /tungstenfs/groups/gbioinfo/sonechar/Projects/alevin_velocity/neuron_10k_v3_mouse
scriptdir := /tungstenfs/groups/gbioinfo/sonechar/Projects/alevin_velocity/scripts

flanklength := 90

methodsall := velocyto $(foreach k,exclude include,$(foreach m,busparse prepref,$(foreach v,collapse separate,kallisto_bustools_$(m)_iso$(v)_$(k)))) \
$(foreach m,busparse prepref,$(foreach v,collapse separate,alevin_$(m)_iso$(v)_cdna_introns)) \
$(foreach m,busparse prepref,$(foreach v,collapse separate,alevin_$(m)_iso$(v)_cdna_introns_decoy)) \
alevin_spliced_unspliced starsolo starsolo_subtr
methodsallc := $(subst $(space),$(comma),$(methodsall))

methods := velocyto $(foreach k,exclude include,$(foreach m,prepref,$(foreach v,collapse separate,kallisto_bustools_$(m)_iso$(v)_$(k)))) \
$(foreach m,prepref,$(foreach v,collapse separate,alevin_$(m)_iso$(v)_cdna_introns)) \
$(foreach m,prepref,$(foreach v,collapse separate,alevin_$(m)_iso$(v)_cdna_introns_decoy)) \
alevin_spliced_unspliced starsolo starsolo_subtr
methodsc := $(subst $(space),$(comma),$(methods))

all: ref quantific veloc plot

ref: $(foreach m,busparse prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)/cDNA_introns.fa.gz)) \
$(foreach m,busparse prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)/introns_cDNA.fa.gz)) \
$(foreach m,busparse prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)/tr2g_modified.tsv)) \
$(foreach m,busparse prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)_cdna_introns.kidx)) \
$(foreach m,busparse prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)_cdna_introns.sidx/seq.bin)) \
$(foreach m,busparse prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)_cdna_intronsasdecoy.sidx/seq.bin)) \
$(foreach m,busparse prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)_introns_cdnaasdecoy.sidx/seq.bin)) \
reference/spliced_unspliced.sidx/seq.bin \
reference/spliced.sidx/seq.bin \
$(foreach m,prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)_uniqueness.txt)) \
$(foreach m,prepref,$(foreach v,collapse separate,reference/$(m)_iso$(v)_uniqueness_overall.txt))

quantific: $(foreach m,prepref,$(foreach v,collapse separate,quants/kallisto_bustools_$(m)_iso$(v)_cdna_introns/run_info.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/kallisto_bustools_$(m)_iso$(v)_cdna_introns/output.correct.sort.bus)) \
$(foreach s,spliced unspliced,$(foreach d,exclude include,$(foreach m,prepref,$(foreach v,collapse separate,quants/kallisto_bustools_$(m)_iso$(v)_cdna_introns/$(s).$(d).bus)))) \
$(foreach s,spliced unspliced,$(foreach d,exclude include,$(foreach m,prepref,$(foreach v,collapse separate,quants/kallisto_bustools_$(m)_iso$(v)_cdna_introns/$(s).$(d).mtx)))) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_introns/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_cdna_intronsasdecoy/lib_format_counts.json)) \
$(foreach m,prepref,$(foreach v,collapse separate,quants/alevin_$(m)_iso$(v)_introns_cdnaasdecoy/lib_format_counts.json)) \
quants/alevin_spliced_unspliced/lib_format_counts.json \
quants/alevin_spliced/lib_format_counts.json \
quants/cellranger/neuron_10k_v3/outs/web_summary.html \
quants/cellranger/neuron_10k_v3/velocyto/neuron_10k_v3.loom \
quants/cellranger/neuron_10k_v3/outs/possorted_genome_bam_retcells.bam \
quants/cellranger/neuron_10k_v3/outs/possorted_genome_bam_retcells.bw \
quants/starsolo/Solo.out/Gene/Summary.csv \
output/sce/sces_spliced_unspliced.rds \
$(foreach x,$(methods) alevin_spliced,output/sce/sce_$(x).rds) \
output/sce/retained_cell_barcodes.csv

veloc: $(foreach x,$(methods),output/anndata/anndata_$(x).h5ad) \
$(foreach x,$(methods),output/anndata_with_velocity/anndata_$(x)_with_velocity.h5ad) \
output/genesets/shared_genes.txt \
$(foreach x,$(methods),output/anndata_with_velocity/anndata_$(x)_shared_genes_with_velocity.h5ad)

plot: plots/misc/plot_velocity_confidence.rds \
plots/misc/plot_count_correlation.rds \
plots/misc/plot_busparse_prepref.rds \
plots/misc/plot_total_count.rds \
plots/misc/plot_spliced_unspliced_mds.rds \
plots/misc/plot_frac_unspliced.rds

## -----------------------------------------------------------------------------------------------------
## Standard CellRanger + velocyto workflow
## -----------------------------------------------------------------------------------------------------
## Run CellRanger
quants/cellranger/neuron_10k_v3/outs/web_summary.html: $(ref10x)/reference.json \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
	mkdir -p quants/cellranger
	cd quants/cellranger && module load CellRanger/3.0.2 && \
	cellranger count --id=neuron_10k_v3 --fastqs=$(datadir) \
	--sample=neuron_10k_v3 --nosecondary --transcriptome=$(ref10x) \
	--localcores=24 --localmem=64 --lanes=1,2

## Run velocyto
quants/cellranger/neuron_10k_v3/velocyto/neuron_10k_v3.loom: \
quants/cellranger/neuron_10k_v3/outs/web_summary.html $(gtf)
	module purge && module load Anaconda3/5.3.0 && \
	source activate velocyto_0.17 && \
	velocyto run10x --samtools-threads 16 quants/cellranger/neuron_10k_v3 $(gtf) && \
	source deactivate

## Subset bam file to only reads from retained cell barcodes
quants/cellranger/neuron_10k_v3/outs/possorted_genome_bam_retcells.bam: \
quants/cellranger/neuron_10k_v3/outs/web_summary.html \
output/sce/retained_cell_barcodes.csv
	$(subsetbam) --bam quants/cellranger/neuron_10k_v3/outs/possorted_genome_bam.bam \
	--cell-barcodes output/sce/retained_cell_barcodes.csv \
	--out-bam $@ --cores 24 --log-level info --bam-tag CB

## Generate bigwig file
quants/cellranger/neuron_10k_v3/outs/possorted_genome_bam_retcells.bw: \
quants/cellranger/neuron_10k_v3/outs/possorted_genome_bam_retcells.bam $(chrlengths)
	module load BEDTools/2.27.1-foss-2018b && \
	bedtools genomecov -split -ibam $< -bg | \
	LC_COLLATE=C sort -k1,1 -k2,2n > neuron_10k_v3.tmp.bedGraph
	module load Kent_tools/20190212-linux.x86_64 && \
	bedGraphToBigWig neuron_10k_v3.tmp.bedGraph $(chrlengths) $@
	rm -f neuron_10k_v3.tmp.bedGraph

## -----------------------------------------------------------------------------------------------------
## STARsolo
## -----------------------------------------------------------------------------------------------------
quants/starsolo/Solo.out/Gene/Summary.csv: \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz \
reference/cellranger3.1.0-3M-february-2018.txt
	mkdir -p quants/starsolo
	module load STAR/2.7.3a-foss-2019b && \
	STAR --genomeDir $(starindex) --readFilesIn $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz,$(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
	$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz,$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz \
	--soloFeatures Gene GeneFull Velocyto --soloType CB_UMI_Simple --soloCBwhitelist reference/cellranger3.1.0-3M-february-2018.txt \
	--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloStrand Forward \
	--outFileNamePrefix quants/starsolo/ --runThreadN 48 --readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI CR UR CB UB GX GN

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
$(foreach v,collapse separate,$(eval $(call farefrule,$(v),busparse,,)))
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
$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call swaprule,$(v),$(m)))))

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
$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call suffixrule,$(v),$(m)))))

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

plots/misc/plot_busparse_prepref.rds: $(txome) $(gtf) $(scriptdir)/compare_busparse_prepref.R \
$(foreach v,collapse separate,$(foreach m,busparse prepref,reference/$(m)_iso$(v)/cDNA_introns.fa.gz)) 
	mkdir -p Rout
	mkdir -p $(@D)
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args txomefa='$(txome)' gtffile='$(gtf)' refdir='reference' flanklength=$(flanklength) outrds='$(@)'" $(scriptdir)/compare_busparse_prepref.R Rout/compare_busparse_prepref.Rout

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
$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call kidxrule,$(v),$(m)))))

## Build Salmon index
define sidxrule1
reference/$(2)_iso$(1)_cdna_introns.sidx/seq.bin: \
reference/$(2)_iso$(1)/cDNA_introns.fa.gz
	module load Salmon/1.0.0 && \
	salmon index -k 31 -i $$(@D) -t $$< --gencode -p 16 --type puff 
endef
$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call sidxrule1,$(v),$(m)))))

## cDNAs, introns as decoys
define sidxrule2
reference/$(2)_iso$(1)_cdna_intronsasdecoy.sidx/seq.bin: \
reference/$(2)_iso$(1)/cDNA_introns.fa.gz
	module load Salmon/1.0.0 && \
	salmon index -k 31 -i $$(@D) -t $$< --gencode -p 16 --type puff \
	-d $$(<D)/introns_tx_to_capture.txt
endef
$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call sidxrule2,$(v),$(m)))))

## introns, cDNAs as decoys
define sidxrule3
reference/$(2)_iso$(1)_introns_cdnaasdecoy.sidx/seq.bin: \
reference/$(2)_iso$(1)/introns_cDNA.fa.gz
	module load Salmon/1.0.0 && \
	salmon index -k 31 -i $$(@D) -t $$< --gencode -p 16 --type puff \
	-d $$(<D)/cDNA_tx_to_capture.txt
endef
$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call sidxrule3,$(v),$(m)))))

## From spliced + unspliced annotated transcripts
reference/spliced_unspliced.sidx/seq.bin: \
reference/gencode.vM21.transcripts.spliced.and.unspliced.fa
	module load Salmon/1.0.0 && \
	salmon index -k 31 -i $(@D) -t $< --gencode -p 16 --type puff 

## From only the regular (spliced) transcripts, for comparison
reference/spliced.sidx/seq.bin: $(txome)
	module load Salmon/1.0.0 && \
	salmon index -k 31 -i $(@D) -t $< --gencode -p 16 --type puff 

## -----------------------------------------------------------------------------------------------------
## kallisto/bustools
## -----------------------------------------------------------------------------------------------------
## Run kallisto bus
define kallistobusrule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/run_info.json: \
reference/$(2)_iso$(1)_cdna_introns.kidx \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
	mkdir -p $$(@D)
	module load kallisto/0.46.0-foss-2019b && \
	kallisto bus -i $$< -o $$(@D) -x 10xv3 -t 24 \
	$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
	$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
endef
$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call kallistobusrule,$(v),$(m)))))

## Correct and sort output
define kallistosortrule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/output.correct.sort.bus: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/run_info.json \
reference/cellranger3.1.0-3M-february-2018.txt
	$(bustools) correct -w reference/cellranger3.1.0-3M-february-2018.txt -p $$(<D)/output.bus | \
	$(bustools) sort -o $$@ -t 24 -
endef
$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call kallistosortrule,$(v),$(m)))))

define kallistosplicedexcluderule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/spliced.exclude.bus: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/output.correct.sort.bus
	$(bustools) capture -s -x -o $$@ \
	-c reference/$(2)_iso$(1)/introns_tx_to_capture.txt \
	-e $$(<D)/matrix.ec -t $$(<D)/transcripts.txt $$<
endef
$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call kallistosplicedexcluderule,$(v),$(m)))))

define kallistounsplicedexcluderule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/unspliced.exclude.bus: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/output.correct.sort.bus
	$(bustools) capture -s -x -o $$@ \
	-c reference/$(2)_iso$(1)/cDNA_tx_to_capture.txt \
	-e $$(<D)/matrix.ec -t $$(<D)/transcripts.txt $$<
endef
$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call kallistounsplicedexcluderule,$(v),$(m)))))

define kallistosplicedincluderule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/spliced.include.bus: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/output.correct.sort.bus
	$(bustools) capture -s -o $$@ \
	-c reference/$(2)_iso$(1)/cDNA_tx_to_capture.txt \
	-e $$(<D)/matrix.ec -t $$(<D)/transcripts.txt $$<
endef
$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call kallistosplicedincluderule,$(v),$(m)))))

define kallistounsplicedincluderule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/unspliced.include.bus: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/output.correct.sort.bus
	$(bustools) capture -s -o $$@ \
	-c reference/$(2)_iso$(1)/introns_tx_to_capture.txt \
	-e $$(<D)/matrix.ec -t $$(<D)/transcripts.txt $$<
endef
$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call kallistounsplicedincluderule,$(v),$(m)))))

## Generate spliced and unspliced count matrices
define kallistocountrule
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/$(4).$(3).mtx: \
quants/kallisto_bustools_$(2)_iso$(1)_cdna_introns/$(4).$(3).bus
	$(bustools) count -o $$(<D)/$(4).$(3) -g reference/$(2)_iso$(1)/tr2g.tsv \
	-e $$(<D)/matrix.ec -t $$(<D)/transcripts.txt --genecounts $$<
endef
$(foreach s,spliced unspliced,$(foreach d,exclude include,$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call kallistocountrule,$(v),$(m),$(d),$(s)))))))

## -----------------------------------------------------------------------------------------------------
## Salmon/alevin
## -----------------------------------------------------------------------------------------------------
## cDNA + introns
define alevincdnaintronrule
quants/alevin_$(2)_iso$(1)_cdna_introns/lib_format_counts.json: \
reference/$(2)_iso$(1)_cdna_introns.sidx/seq.bin \
reference/$(2)_iso$(1)/tr2g_modified.tsv \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
	module load Salmon/1.0.0 && \
	salmon alevin -l ISR -i $$(<D) \
	-1 $(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz \
	-o $$(@D) -p 24 --tgMap $$(word 2,$$^) \
	--chromiumV3 --dumpFeatures --dumpBfh
endef
$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call alevincdnaintronrule,$(v),$(m)))))

define alevincdnaintrondecoyrule
quants/alevin_$(2)_iso$(1)_cdna_intronsasdecoy/lib_format_counts.json: \
reference/$(2)_iso$(1)_cdna_intronsasdecoy.sidx/seq.bin \
reference/$(2)_iso$(1)/tr2g_modified.tsv \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
	module load Salmon/1.0.0 && \
	salmon alevin -l ISR -i $$(<D) \
	-1 $(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz \
	-o $$(@D) -p 24 --tgMap $$(word 2,$$^) \
	--chromiumV3 --dumpFeatures --dumpBfh
endef
$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call alevincdnaintrondecoyrule,$(v),$(m)))))

define alevinintroncdnadecoyrule
quants/alevin_$(2)_iso$(1)_introns_cdnaasdecoy/lib_format_counts.json: \
reference/$(2)_iso$(1)_introns_cdnaasdecoy.sidx/seq.bin \
reference/$(2)_iso$(1)/tr2g_modified.tsv \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
	module load Salmon/1.0.0 && \
	salmon alevin -l ISR -i $$(<D) \
	-1 $(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz \
	-o $$(@D) -p 24 --tgMap $$(word 2,$$^) \
	--chromiumV3 --dumpFeatures --dumpBfh
endef
$(foreach m,busparse prepref,$(foreach v,collapse separate,$(eval $(call alevinintroncdnadecoyrule,$(v),$(m)))))

## spliced + unspliced transcripts
quants/alevin_spliced_unspliced/lib_format_counts.json: \
reference/spliced_unspliced.sidx/seq.bin \
reference/gencode.vM21.spliced.unspliced.tx2gene.txt \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
	module load Salmon/1.0.0 && \
	salmon alevin -l ISR -i $(<D) \
	-1 $(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz \
	-o $(@D) -p 24 --tgMap $(word 2,$^) \
	--chromiumV3 --dumpFeatures --dumpBfh

## spliced transcripts only (regular annotation)
quants/alevin_spliced/lib_format_counts.json: \
reference/spliced.sidx/seq.bin \
$(tx2gene) \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
	module load Salmon/1.0.0 && \
	salmon alevin -l ISR -i $(<D) \
	-1 $(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz \
	-o $(@D) -p 16 --tgMap $(word 2,$^) \
	--chromiumV3 --dumpFeatures --dumpBfh

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
quants/cellranger/neuron_10k_v3/velocyto/neuron_10k_v3.loom \
quants/starsolo/Solo.out/Gene/Summary.csv
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' dataset='Neuron' helperscript='$(scriptdir)/sce_helpers.R' tx2gene='$(tx2genefull)' cellfile='' samplename='neuron_10k_v3' outrds='$@'" $(scriptdir)/summarize_spliced_unspliced.R Rout/summarize_spliced_unspliced.Rout

define scerule
output/sce/sce_$(1).rds: output/sce/sces_spliced_unspliced.rds
	touch $$@
endef
$(foreach x,$(methods) alevin_spliced,$(eval $(call scerule,$(x))))

output/sce/retained_cell_barcodes.csv: output/sce/sces_spliced_unspliced.rds
	touch $@

plots/misc/plot_total_count.rds: output/sce/sces_spliced_unspliced.rds \
$(scriptdir)/plot_total_count.R $(tx2genefull) $(scriptdir)/plot_helpers.R \
$(foreach v,separate collapse,reference/prepref_iso$(v)_uniqueness.txt) \
$(foreach v,separate collapse,reference/prepref_iso$(v)_uniqueness_overall.txt)
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' plothelperscript='$(scriptdir)/plot_helpers.R' tx2gene='$(tx2genefull)' methods='$(methodsc)' outrds='$@'" $(scriptdir)/plot_total_count.R Rout/plot_total_count.Rout

## Plot summaries of spliced/unspliced counts
plots/misc/plot_spliced_unspliced_mds.rds: output/sce/sces_spliced_unspliced.rds \
$(scriptdir)/plot_spliced_unspliced_mds.R $(scriptdir)/plot_helpers.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' plothelperscript='$(scriptdir)/plot_helpers.R' methods='$(methodsc)' outrds='$@'" $(scriptdir)/plot_spliced_unspliced_mds.R Rout/plot_spliced_unspliced_mds.Rout

plots/misc/plot_frac_unspliced.rds: output/sce/sces_spliced_unspliced.rds \
$(scriptdir)/plot_frac_unspliced.R $(scriptdir)/plot_helpers.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' plothelperscript='$(scriptdir)/plot_helpers.R' methods='$(methodsc)' outrds='$@'" $(scriptdir)/plot_frac_unspliced.R Rout/plot_frac_unspliced.Rout

plots/misc/plot_velocity_confidence.rds: $(foreach x,$(methods),output/anndata_with_velocity/anndata_$(x)_with_velocity.h5ad) \
$(scriptdir)/plot_velocity_confidence_neuron.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' methods='$(methodsc)' outrds='$@'" $(scriptdir)/plot_velocity_confidence_neuron.R Rout/plot_velocity_confidence.Rout

plots/misc/plot_count_correlation.rds: output/genesets/shared_genes.txt output/sce/sces_spliced_unspliced.rds \
$(scriptdir)/plot_count_correlation_neuron.R 
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' methods='$(methodsc)' genetxt='$<' outrds='$@'" $(scriptdir)/plot_count_correlation_neuron.R Rout/plot_count_correlation.Rout

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
$(scriptdir)/run_scVelo_neuron.py
	mkdir -p plots/velocity/anndata_$(1)
	mkdir -p $$(@D)
	mkdir -p pylog
	module purge && module load Anaconda3/5.3.0 && \
	set +u && source activate scvelo_0.1.24 && set -u && \
	python $(scriptdir)/run_scVelo_neuron.py $$< plots/velocity $$@ "None" >pylog/$(1).scvelo.log 2>&1 && \
	set +u && source deactivate && set -u && module purge
endef
$(foreach x,$(methods),$(eval $(call scvelorule,$(x))))

## Make text file with genes retained by all methods (except busparse-based)
output/genesets/shared_genes.txt: $(scriptdir)/make_gene_subsets.R \
$(foreach x,$(methods),output/anndata_with_velocity/anndata_$(x)_with_velocity.h5ad)
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' methods='$(methodsc)' outtxt='$@'" $(scriptdir)/make_gene_subsets.R Rout/make_gene_subsets.Rout

## Velocity, based on genes selected by all methods
define scveloruleshared
output/anndata_with_velocity/anndata_$(1)_shared_genes_with_velocity.h5ad: output/anndata/anndata_$(1).h5ad \
$(scriptdir)/run_scVelo_neuron.py output/genesets/shared_genes.txt
	mkdir -p plots/velocity/anndata_$(1)_shared_genes
	mkdir -p $$(@D)
	mkdir -p pylog
	module purge && module load Anaconda3/5.3.0 && \
	set +u && source activate scvelo_0.1.24 && set -u && \
	python $(scriptdir)/run_scVelo_neuron.py $$< plots/velocity $$@ output/genesets/shared_genes.txt >pylog/$(1).scvelo_shared_genes.log 2>&1 && \
	set +u && source deactivate && set -u && module purge
endef
$(foreach x,$(methods),$(eval $(call scveloruleshared,$(x))))
