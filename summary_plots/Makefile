SHELL := /bin/bash

all: summary_velocity.pdf fraction_decoys.txt

summary_velocity.pdf: ../spermatogenesis_mouse/plots/velocity_confidence/spermatogenesis_plot_velocity_confidence.rds \
../pancreas_mouse/plots/velocity_confidence/pancreas_plot_velocity_confidence.rds \
../dentate_gyrus_mouse/plots/velocity_confidence/dentate_gyrus_plot_velocity_confidence.rds \
../pfc_mouse/plots/velocity_confidence/pfc_plot_velocity_confidence.rds \
summary_velocity.R
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore summary_velocity.R Rout/summary_velocity.Rout

fraction_decoys.txt: ../spermatogenesis_mouse/output/sce/sces_spliced_unspliced.rds \
../pancreas_mouse/output/sce/sces_spliced_unspliced.rds \
../dentate_gyrus_mouse/output/sce/sces_spliced_unspliced.rds \
../pfc_mouse/output/sce/sces_spliced_unspliced.rds \
get_frac_decoys.R
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore get_frac_decoys.R Rout/get_frac_decoys.Rout