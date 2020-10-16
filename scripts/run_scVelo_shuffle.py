import sys
import anndata
import scvelo as scv
import os
import matplotlib
import pandas as pd
import numpy as np
import random

matplotlib.use('AGG')

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.set_figure_params('scvelo')  # for beautified visualization

# get input file from command line
n = len(sys.argv)
if len(sys.argv) >= 4:
	adfile = sys.argv[1]
	plotdir = sys.argv[2]
	genesetfile = sys.argv[3]
else:
	raise ValueError("must have >=three arguments (adfile, plotdir, genesetfile)")

print("input adfile: ", adfile)
print("output directory: ", plotdir)
print("gene set file: ", genesetfile)

base = os.path.basename(adfile)
base = os.path.splitext(base)[0]
mname = base.replace("anndata_", "")

if genesetfile != "None":
	basegs = os.path.basename(genesetfile)
	basegs = "_" + os.path.splitext(basegs)[0]
else:
	basegs = ""

# read input
adata = anndata.read(adfile)
scv.utils.show_proportions(adata)

# normalize the data
if genesetfile == "None":
	scv.pp.filter_genes(adata, min_shared_counts = 20)
	scv.pp.normalize_per_cell(adata, enforce = True)
	scv.pp.filter_genes_dispersion(adata, n_top_genes = 2000)
	scv.pp.log1p(adata)
else:
	geneset = [f.rstrip('\n') for f in open(genesetfile, 'r').readlines()]
	scv.pp.filter_genes(adata, min_shared_counts = 20)
	scv.pp.normalize_per_cell(adata, enforce = True)
	adata = adata[:, geneset]

## Convert unspliced matrix to dense format. 
adata.layers['unspliced'] = adata.layers['unspliced'].todense()

# Repeatedly shuffle the unspliced counts for each gene and estimate velocities
ix = np.arange(adata.layers['unspliced'].shape[0])
random.seed(30)
outfile = open(plotdir + "/" + base + basegs + "/velocity_shuffle.out", "w")
for i in range(15):
	try:
		adata2 = adata.copy()
		for j in range(adata2.shape[1]):
			adata2.layers['unspliced'][:, j] = adata2.layers['unspliced'][np.random.permutation(ix), j]
		scv.pp.moments(adata2, n_pcs = 30, n_neighbors = 30)
	
		# compute velocity and velocity graph
		scv.tl.recover_dynamics(adata2)
		scv.tl.velocity(adata2, mode = 'dynamical')
		scv.tl.velocity_graph(adata2)
		## Maximal cosine correlation for each cell (used by scVelo to calculate self-transition probabilities)
		adata2.obs['max_cosine_corr'] = adata2.uns['velocity_graph'].max(1).A.flatten()
	
		scv.settings.figdir = plotdir + '/' + base + basegs + '/'
		scv.settings.plot_prefix = base + basegs + '_scvelo_'
		scv.settings.set_figure_params(dpi_save = 300, vector_friendly = True)
	
		scv.pl.velocity_embedding_stream(adata2, basis='PCA_alevin_spliced_gentrome', save="PCA_alevin_spliced_gentrome_stream_" + str(i) + ".png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='clusters', title='')
		scv.pl.velocity_embedding_stream(adata2, basis='TSNE_alevin_spliced_gentrome', save="TSNE_alevin_spliced_gentrome_stream_" + str(i) + ".png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='clusters', title='')
		scv.pl.velocity_embedding_stream(adata2, basis='UMAP_alevin_spliced_gentrome', save="UMAP_alevin_spliced_gentrome_stream_" + str(i) + ".png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='clusters', title='')
	
		scv.tl.velocity_embedding(adata2, basis = 'UMAP_alevin_spliced_gentrome', all_comps = False, autoscale = False)
		pd.DataFrame(adata2.obsm['velocity_UMAP_alevin_spliced_gentrome'], index = adata2.obs.index).to_csv(plotdir + "/" + base + basegs + "/" + base + basegs + "_velocity_UMAP_alevin_spliced_gentrome_" + str(i) + ".csv")

		try:
			scv.tl.recover_latent_time(adata2)
		except:
			print('Latent time/top genes could not be extracted')
	
		scv.tl.velocity_confidence(adata2)
		adata2.obs.to_csv(plotdir + "/" + base + basegs + "/" + base + basegs + "_cell_info_" + str(i) + ".csv")
		adata2.var.to_csv(plotdir + "/" + base + basegs + "/" + base + basegs + "_gene_info_" + str(i) + ".csv")
		outfile.write(str(i) + "\n") 
	except:
		print('Velocity calculation failed')

# session info
scv.logging.print_version()

outfile.close()
