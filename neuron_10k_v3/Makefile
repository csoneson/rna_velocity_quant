genome := /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/GRCm38.primary_assembly.genome.fa
gtf := /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/gencode.vM21.annotation.gtf
ref10x := /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/ref10x_spliced/cellranger3_0_2_GRCm38.primary_assembly_gencodeM21_spliced
tx2genefull := /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/tx2gene/gencode.vM21.annotation.tx2gene_full.rds

bustools := /tungstenfs/groups/gbioinfo/sonechar/Software/bustools_0.39.3/bustools
salmon_new := /tungstenfs/groups/gbioinfo/sonechar/Software/salmon-latest_linux_x86_64_20191027/bin/salmon

datadir := /tungstenfs/groups/gbioinfo/sonechar/Data/10xGenomics_neuron_10k_v3/neuron_10k_v3_fastqs
topdir := /tungstenfs/groups/gbioinfo/sonechar/Projects/alevin_velocity/neuron_10k_v3

methods := velocyto cdna_introns_decoy cdna_introns cdna_intronscollapsed spliced_unspliced kallistobus_separate \
kallistobus_collapse kallistobus_sep_permissive kallistobus_coll_permissive

all: quantific veloc

quantific: reference/busparse_gencodevM21_isoseparate_exonfull/cDNA_introns.fa.gz \
reference/busparse_gencodevM21_isoseparate_exonfull_cDNA_introns.kidx \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/run_info.json \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/output.correct.sort.bus \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/spliced.bus \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/spliced.txt \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/unspliced.bus \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/unspliced.txt \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/unspliced.mtx \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/spliced.mtx \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/spliced.permissive.bus \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/unspliced.permissive.bus \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/unspliced.permissive.mtx \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/spliced.permissive.mtx \
reference/salmon_gencodevM21_cDNA_introns.sidx/seq.bin \
reference/salmon_gencodevM21_spliced_unspliced.sidx/seq.bin \
reference/salmon_gencodevM21_spliced.sidx/seq.bin \
reference/busparse_gencodevM21_isoseparate_exonfull/tr2g_modified.tsv \
reference/salmon_gencodevM21_cDNA_intronsasdecoy.sidx/seq.bin \
reference/salmon_gencodevM21_introns_cDNAasdecoy.sidx/seq.bin \
quants/salmon_gencodevM21_cDNA_introns/lib_format_counts.json \
quants/salmon_gencodevM21_spliced_unspliced/lib_format_counts.json \
quants/salmon_gencodevM21_spliced/lib_format_counts.json \
quants/salmon_gencodevM21_cDNA_intronsasdecoy/lib_format_counts.json \
quants/salmon_gencodevM21_introns_cDNAasdecoy/lib_format_counts.json \
quants/cellranger/neuron_10k_v3/outs/web_summary.html \
quants/cellranger/neuron_10k_v3/velocyto/neuron_10k_v3.loom \
output/summarize_spliced_unspliced.rds \
reference/busparse_gencodevM21_isocollapse_exonfull/cDNA_introns.fa.gz \
reference/busparse_gencodevM21_isocollapse_exonfull/introns_cDNA.fa.gz \
reference/busparse_gencodevM21_isocollapse_exonfull/tr2g_modified.tsv \
reference/busparse_gencodevM21_isocollapse_exonfull_cDNA_introns.kidx \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/run_info.json \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/output.correct.sort.bus \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/spliced.bus \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/spliced.txt \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/unspliced.bus \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/unspliced.txt \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/unspliced.mtx \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/spliced.mtx \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/spliced.permissive.bus \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/unspliced.permissive.bus \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/unspliced.permissive.mtx \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/spliced.permissive.mtx \
reference/salmon_gencodevM21_cDNA_intronscollapsed.sidx/seq.bin \
quants/salmon_gencodevM21_cDNA_intronscollapsed/lib_format_counts.json \
$(foreach x,$(methods) spliced,output/sce_$(x).rds)

veloc: $(foreach x,$(methods),velocity/anndata_$(x).h5ad)

plot: output/plot_spliced_unspliced.rds

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
	module load Anaconda3/5.3.0 && \
	source activate velocyto_0.17 && \
	velocyto run10x --samtools-threads 16 quants/cellranger/neuron_10k_v3 $(gtf) && \
	source deactivate

## -----------------------------------------------------------------------------------------------------
## Generate fasta file with transcripts + introns (with some flanking sequence)
## -----------------------------------------------------------------------------------------------------
## Extracts intronic sequences flanked by L-1 bases of exonic sequences where L is the biological 
## read length (for this data set, 91)
## Isoforms are considered separately when introns are extracted
reference/busparse_gencodevM21_isoseparate_exonfull/cDNA_introns.fa.gz: $(genome) $(gtf) \
scripts/generate_cdna_intron_fa.R
	mkdir -p Rout
	mkdir -p $(@D)
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args gtf='$(gtf)' genome='$(genome)' isoform_action='separate' outdir='$(@D)'" scripts/generate_cdna_intron_fa.R Rout/generate_cdna_intron_fa_GRCm38_gencodevM21.Rout

## Swap order of cDNAs and introns, to allow quantifying introns with cDNAs as decoy sequences
reference/busparse_gencodevM21_isoseparate_exonfull/introns_cDNA.fa.gz: \
reference/busparse_gencodevM21_isoseparate_exonfull/cDNA_introns.fa.gz \
scripts/swap_order_fasta_records.R
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args infa='$<' outfa='$@'" scripts/swap_order_fasta_records.R Rout/swap_order_fasta_records.Rout

## Add suffix to gene ID for introns in transcript-to-gene mapping (to avoid both 
## spliced and unspliced reads being assigned to the same gene ID and thus be indistinguishable)
reference/busparse_gencodevM21_isoseparate_exonfull/tr2g_modified.tsv: \
reference/busparse_gencodevM21_isoseparate_exonfull/cDNA_introns.fa.gz \
scripts/identify_unspliced_genes_in_tx2gene.R
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args intsv='$(<D)/tr2g.tsv' outtsv='$@'" scripts/identify_unspliced_genes_in_tx2gene.R Rout/identify_unspliced_genes_in_tx2gene.Rout


## As above, but isoforms are considered separately when introns are extracted
reference/busparse_gencodevM21_isocollapse_exonfull/cDNA_introns.fa.gz: $(genome) $(gtf) \
scripts/generate_cdna_intron_fa.R
	mkdir -p Rout
	mkdir -p $(@D)
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args gtf='$(gtf)' genome='$(genome)' isoform_action='collapse' outdir='$(@D)'" scripts/generate_cdna_intron_fa.R Rout/generate_cdna_intron_fa_GRCm38_gencodevM21_collapse.Rout

## Swap order of cDNAs and introns, to allow quantifying introns with cDNAs as decoy sequences
reference/busparse_gencodevM21_isocollapse_exonfull/introns_cDNA.fa.gz: \
reference/busparse_gencodevM21_isocollapse_exonfull/cDNA_introns.fa.gz \
scripts/swap_order_fasta_records.R
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args infa='$<' outfa='$@'" scripts/swap_order_fasta_records.R Rout/swap_order_fasta_records_isocollapse.Rout

## Add suffix to gene ID for introns in transcript-to-gene mapping (to avoid both 
## spliced and unspliced reads being assigned to the same gene ID and thus be indistinguishable)
reference/busparse_gencodevM21_isocollapse_exonfull/tr2g_modified.tsv: \
reference/busparse_gencodevM21_isocollapse_exonfull/cDNA_introns.fa.gz \
scripts/identify_unspliced_genes_in_tx2gene.R
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args intsv='$(<D)/tr2g.tsv' outtsv='$@'" scripts/identify_unspliced_genes_in_tx2gene.R Rout/identify_unspliced_genes_in_tx2gene_isocollapse.Rout

## -----------------------------------------------------------------------------------------------------
## kallisto/bustools
## -----------------------------------------------------------------------------------------------------
## Build kallisto index
reference/busparse_gencodevM21_isoseparate_exonfull_cDNA_introns.kidx: \
reference/busparse_gencodevM21_isoseparate_exonfull/cDNA_introns.fa.gz
	module load kallisto/0.46.0-foss-2019b && \
	kallisto index -k 31 -i $@ $<

## Run kallisto bus
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/run_info.json: \
reference/busparse_gencodevM21_isoseparate_exonfull_cDNA_introns.kidx \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
	mkdir -p $(@D)
	module load kallisto/0.46.0-foss-2019b && \
	kallisto bus -i $< -o $(@D) -x 10xv3 -t 16 \
	$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
	$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz

## Correct and sort output
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/output.correct.sort.bus: \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/run_info.json \
reference/cellranger3.1.0-3M-february-2018.txt
	$(bustools) correct -w reference/cellranger3.1.0-3M-february-2018.txt -p $(<D)/output.bus | \
	$(bustools) sort -o $@ -t 16 -

## Generate spliced and unspliced bus files
## Here we're taking the complement of the introns as a capture list. 
## In other words, exclude any equivalence class that contains any of the introns
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/spliced.bus: \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/output.correct.sort.bus
	$(bustools) capture -s -x -o $@ \
	-c reference/busparse_gencodevM21_isoseparate_exonfull/introns_tx_to_capture.txt \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt $<

quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/unspliced.bus: \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/output.correct.sort.bus
	$(bustools) capture -s -x -o $@ \
	-c reference/busparse_gencodevM21_isoseparate_exonfull/cDNA_tx_to_capture.txt \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt $<

## Also try a permissive approach, where we keep all equivalence classes that contain any of the provided features
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/spliced.permissive.bus: \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/output.correct.sort.bus
	$(bustools) capture -s -o $@ \
	-c reference/busparse_gencodevM21_isoseparate_exonfull/cDNA_tx_to_capture.txt \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt $<

quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/unspliced.permissive.bus: \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/output.correct.sort.bus
	$(bustools) capture -s -o $@ \
	-c reference/busparse_gencodevM21_isoseparate_exonfull/introns_tx_to_capture.txt \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt $<

## Generate spliced and unspliced count matrices
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/unspliced.mtx: \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/unspliced.bus
	$(bustools) count -o $(<D)/unspliced -g reference/busparse_gencodevM21_isoseparate_exonfull/tr2g.tsv \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt --genecounts $<

quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/spliced.mtx: \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/spliced.bus
	$(bustools) count -o $(<D)/spliced -g reference/busparse_gencodevM21_isoseparate_exonfull/tr2g.tsv \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt --genecounts $<

## Permissive approach
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/unspliced.permissive.mtx: \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/unspliced.permissive.bus
	$(bustools) count -o $(<D)/unspliced.permissive -g reference/busparse_gencodevM21_isoseparate_exonfull/tr2g.tsv \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt --genecounts $<

quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/spliced.permissive.mtx: \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/spliced.permissive.bus
	$(bustools) count -o $(<D)/spliced.permissive -g reference/busparse_gencodevM21_isoseparate_exonfull/tr2g.tsv \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt --genecounts $<

## convert bus files to text files
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/unspliced.txt: \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/unspliced.bus
	$(bustools) text -o $@ $<

quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/spliced.txt: \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/spliced.bus
	$(bustools) text -o $@ $<

## -----------------------------------------------------------------------------------------------------
## kallisto/bustools, collapsed introns
## -----------------------------------------------------------------------------------------------------
## Build kallisto index
reference/busparse_gencodevM21_isocollapse_exonfull_cDNA_introns.kidx: \
reference/busparse_gencodevM21_isocollapse_exonfull/cDNA_introns.fa.gz
	module load kallisto/0.46.0-foss-2019b && \
	kallisto index -k 31 -i $@ $<

## Run kallisto bus
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/run_info.json: \
reference/busparse_gencodevM21_isocollapse_exonfull_cDNA_introns.kidx \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
	mkdir -p $(@D)
	module load kallisto/0.46.0-foss-2019b && \
	kallisto bus -i $< -o $(@D) -x 10xv3 -t 16 \
	$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
	$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz

## Correct and sort output
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/output.correct.sort.bus: \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/run_info.json \
reference/cellranger3.1.0-3M-february-2018.txt
	$(bustools) correct -w reference/cellranger3.1.0-3M-february-2018.txt -p $(<D)/output.bus | \
	$(bustools) sort -o $@ -t 16 -

## Generate spliced and unspliced bus files
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/spliced.bus: \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/output.correct.sort.bus
	$(bustools) capture -s -x -o $@ \
	-c reference/busparse_gencodevM21_isocollapse_exonfull/introns_tx_to_capture.txt \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt $<

quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/unspliced.bus: \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/output.correct.sort.bus
	$(bustools) capture -s -x -o $@ \
	-c reference/busparse_gencodevM21_isocollapse_exonfull/cDNA_tx_to_capture.txt \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt $<

## Also try a permissive approach, where we keep all equivalence classes that contain any of the provided features
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/spliced.permissive.bus: \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/output.correct.sort.bus
	$(bustools) capture -s -o $@ \
	-c reference/busparse_gencodevM21_isocollapse_exonfull/cDNA_tx_to_capture.txt \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt $<

quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/unspliced.permissive.bus: \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/output.correct.sort.bus
	$(bustools) capture -s -o $@ \
	-c reference/busparse_gencodevM21_isocollapse_exonfull/introns_tx_to_capture.txt \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt $<

## Generate spliced and unspliced count matrices
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/unspliced.mtx: \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/unspliced.bus
	$(bustools) count -o $(<D)/unspliced -g reference/busparse_gencodevM21_isocollapse_exonfull/tr2g.tsv \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt --genecounts $<

quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/spliced.mtx: \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/spliced.bus
	$(bustools) count -o $(<D)/spliced -g reference/busparse_gencodevM21_isocollapse_exonfull/tr2g.tsv \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt --genecounts $<

## Permissive approach
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/unspliced.permissive.mtx: \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/unspliced.permissive.bus
	$(bustools) count -o $(<D)/unspliced.permissive -g reference/busparse_gencodevM21_isocollapse_exonfull/tr2g.tsv \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt --genecounts $<

quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/spliced.permissive.mtx: \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/spliced.permissive.bus
	$(bustools) count -o $(<D)/spliced.permissive -g reference/busparse_gencodevM21_isocollapse_exonfull/tr2g.tsv \
	-e $(<D)/matrix.ec -t $(<D)/transcripts.txt --genecounts $<

## convert bus files to text files
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/unspliced.txt: \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/unspliced.bus
	$(bustools) text -o $@ $<

quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/spliced.txt: \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/spliced.bus
	$(bustools) text -o $@ $<


## -----------------------------------------------------------------------------------------------------
## Salmon/alevin
## -----------------------------------------------------------------------------------------------------
## Build Salmon indexes
## From cDNAs + introns as generated by BUSpaRse
reference/salmon_gencodevM21_cDNA_introns.sidx/seq.bin: \
reference/busparse_gencodevM21_isoseparate_exonfull/cDNA_introns.fa.gz
	module load Salmon/0.99.0-beta2 && \
	salmon index -k 31 -i $(@D) -t $< --gencode -p 16 --type puff 

## From cDNAs + intronscollapsed as generated by BUSpaRse
reference/salmon_gencodevM21_cDNA_intronscollapsed.sidx/seq.bin: \
reference/busparse_gencodevM21_isocollapse_exonfull/cDNA_introns.fa.gz
	module load Salmon/0.99.0-beta2 && \
	salmon index -k 31 -i $(@D) -t $< --gencode -p 16 --type puff 

## From cDNAs only, using introns as decoy sequences
reference/salmon_gencodevM21_cDNA_intronsasdecoy.sidx/seq.bin: \
reference/busparse_gencodevM21_isoseparate_exonfull/cDNA_introns.fa.gz
	module load Salmon/0.99.0-beta2 && \
	salmon index -k 31 -i $(@D) -t $< --gencode -p 16 --type puff -d $(<D)/introns_tx_to_capture.txt

## From introns only, using cDNAs as decoy sequences
reference/salmon_gencodevM21_introns_cDNAasdecoy.sidx/seq.bin: \
reference/busparse_gencodevM21_isoseparate_exonfull/introns_cDNA.fa.gz
	module load Salmon/0.99.0-beta2 && \
	salmon index -k 31 -i $(@D) -t $< --gencode -p 16 --type puff -d $(<D)/cDNA_tx_to_capture.txt

## From spliced + unspliced annotated transcripts
reference/salmon_gencodevM21_spliced_unspliced.sidx/seq.bin: \
/tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/gencode.vM21.transcripts.spliced.and.unspliced.min23.fa
	module load Salmon/0.99.0-beta2 && \
	salmon index -k 31 -i $(@D) -t $< --gencode -p 16 --type puff 

## From only the regular (spliced) transcripts, for comparison
reference/salmon_gencodevM21_spliced.sidx/seq.bin: \
/tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/gencode.vM21.transcripts.fa
	module load Salmon/0.99.0-beta2 && \
	salmon index -k 31 -i $(@D) -t $< --gencode -p 16 --type puff 

## Quantification
## cDNA + introns
quants/salmon_gencodevM21_cDNA_introns/lib_format_counts.json: \
reference/salmon_gencodevM21_cDNA_introns.sidx/seq.bin \
reference/busparse_gencodevM21_isoseparate_exonfull/tr2g_modified.tsv \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
	$(salmon_new) alevin -l ISR -i $(<D) \
	-1 $(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz \
	-o $(@D) -p 16 --tgMap $(word 2,$^) \
	--chromiumV3 --dumpFeatures --dumpBfh

## cDNA + intronscollapsed
quants/salmon_gencodevM21_cDNA_intronscollapsed/lib_format_counts.json: \
reference/salmon_gencodevM21_cDNA_intronscollapsed.sidx/seq.bin \
reference/busparse_gencodevM21_isocollapse_exonfull/tr2g_modified.tsv \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
	$(salmon_new) alevin -l ISR -i $(<D) \
	-1 $(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz \
	-o $(@D) -p 16 --tgMap $(word 2,$^) \
	--chromiumV3 --dumpFeatures --dumpBfh

## cDNA, introns as decoy sequences
quants/salmon_gencodevM21_cDNA_intronsasdecoy/lib_format_counts.json: \
reference/salmon_gencodevM21_cDNA_intronsasdecoy.sidx/seq.bin \
reference/busparse_gencodevM21_isoseparate_exonfull/tr2g_modified.tsv \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
	$(salmon_new) alevin -l ISR -i $(<D) \
	-1 $(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz \
	-o $(@D) -p 16 --tgMap $(word 2,$^) \
	--chromiumV3 --dumpFeatures --dumpBfh

## introns, cDNA as decoy sequences
quants/salmon_gencodevM21_introns_cDNAasdecoy/lib_format_counts.json: \
reference/salmon_gencodevM21_introns_cDNAasdecoy.sidx/seq.bin \
reference/busparse_gencodevM21_isoseparate_exonfull/tr2g_modified.tsv \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
	$(salmon_new) alevin -l ISR -i $(<D) \
	-1 $(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz \
	-o $(@D) -p 16 --tgMap $(word 2,$^) \
	--chromiumV3 --dumpFeatures --dumpBfh

## spliced + unspliced transcripts
quants/salmon_gencodevM21_spliced_unspliced/lib_format_counts.json: \
reference/salmon_gencodevM21_spliced_unspliced.sidx/seq.bin \
/tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/tx2gene/gencode.vM21.annotation.spliced.and.unspliced.tx2gene.txt \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
	$(salmon_new) alevin -l ISR -i $(<D) \
	-1 $(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz \
	-o $(@D) -p 16 --tgMap $(word 2,$^) \
	--chromiumV3 --dumpFeatures --dumpBfh

## spliced transcripts only (regular annotation)
quants/salmon_gencodevM21_spliced/lib_format_counts.json: \
reference/salmon_gencodevM21_spliced.sidx/seq.bin \
/tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M21/tx2gene/gencode.vM21.annotation.tx2gene.txt \
$(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
$(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz
	$(salmon_new) alevin -l ISR -i $(<D) \
	-1 $(datadir)/neuron_10k_v3_S1_L002_R1_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R1_001.fastq.gz \
	-2 $(datadir)/neuron_10k_v3_S1_L002_R2_001.fastq.gz $(datadir)/neuron_10k_v3_S1_L001_R2_001.fastq.gz \
	-o $(@D) -p 16 --tgMap $(word 2,$^) \
	--chromiumV3 --dumpFeatures --dumpBfh

## -----------------------------------------------------------------------------------------------------
## Analysis
## -----------------------------------------------------------------------------------------------------
output/summarize_spliced_unspliced.rds: quants/cellranger/neuron_10k_v3/velocyto/neuron_10k_v3.loom \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/unspliced.mtx \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/spliced.mtx \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/unspliced.mtx \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/spliced.mtx \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/unspliced.permissive.mtx \
quants/kallisto_bustools_gencodevM21_isoseparate_exonfull_cDNA_introns/spliced.permissive.mtx \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/unspliced.permissive.mtx \
quants/kallisto_bustools_gencodevM21_isocollapse_exonfull_cDNA_introns/spliced.permissive.mtx \
quants/salmon_gencodevM21_cDNA_introns/lib_format_counts.json \
quants/salmon_gencodevM21_cDNA_intronscollapsed/lib_format_counts.json \
quants/salmon_gencodevM21_spliced_unspliced/lib_format_counts.json \
quants/salmon_gencodevM21_spliced/lib_format_counts.json \
quants/salmon_gencodevM21_cDNA_intronsasdecoy/lib_format_counts.json \
quants/salmon_gencodevM21_introns_cDNAasdecoy/lib_format_counts.json \
scripts/summarize_spliced_unspliced.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' tx2gene='$(tx2genefull)' outrds='$@'" scripts/summarize_spliced_unspliced.R Rout/summarize_spliced_unspliced.Rout

output/plot_spliced_unspliced.rds: output/summarize_spliced_unspliced.rds \
scripts/plot_spliced_unspliced.R
	mkdir -p $(@D)
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore "--args topdir='$(topdir)' outrds='$@'" scripts/plot_spliced_unspliced.R Rout/plot_spliced_unspliced.Rout

define scerule
output/sce_$(1).rds: output/summarize_spliced_unspliced.rds
	touch $$@
endef
$(foreach x,$(methods) spliced,$(eval $(call scerule,$(x))))

## Convert SCE to anndata
define anndatarule
velocity/anndata_$(1).h5ad: output/sce_$(1).rds
	mkdir -p velocity
	mkdir -p pylog
	module purge && module load Anaconda3/5.3.0 && \
	set +u && source activate anndata2ri_1.0 && set -u && \
	python scripts/convert_sce_to_anndata.py $$< $$@ >pylog/$(1).anndata.log 2>&1 && \
	set +u && source deactivate && set -u && module purge
endef
$(foreach x,$(methods),$(eval $(call anndatarule,$(x))))





