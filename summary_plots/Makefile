SHELL := /bin/bash

all: summary_velocity.pdf

summary_velocity.pdf: ../neuron_10k_v3_mouse/plots/velocity_confidence/neuron_plot_velocity_confidence.rds \
../pancreas_mouse/plots/velocity_confidence/pancreas_plot_velocity_confidence.rds \
../dentate_gyrus_mouse/plots/velocity_confidence/dentate_gyrus_plot_velocity_confidence.rds \
summary_velocity.R
	mkdir -p Rout
	module load R-BioC/3.6-3.10-foss-2019b && \
	R CMD BATCH --no-save --no-restore summary_velocity.R Rout/summary_velocity.Rout
