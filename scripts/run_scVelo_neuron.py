import sys
import anndata
import scvelo as scv
import os
import matplotlib
import pandas as pd

matplotlib.use('AGG')

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.set_figure_params('scvelo')  # for beautified visualization

# get input file from command line
n = len(sys.argv)
if len(sys.argv) >= 5:
	adfile = sys.argv[1]
	plotdir = sys.argv[2]
	adfileout = sys.argv[3]
	genesetfile = sys.argv[4]
else:
	raise ValueError("must have >=four arguments (adfile, plotdir, adfileout, genesetfile)")

print("input adfile: ", adfile)
print("output directory: ", plotdir)
print("output adfile: ", adfileout)
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

# preprocess the data
# set log = False if data is already log-transformed
# runs the following steps:
# scv.pp.filter_genes(adata)
# scv.pp.normalize_per_cell(adata)
# scv.pp.filter_genes_dispersion(adata)
# scv.pp.log1p(adata)
# it will determine whether to log-transform X:
# log_advised = np.allclose(adata.X[:10].sum(), adata.layers['spliced'][:10].sum())
# if log and log_advised: log1p(adata)
# elif log and not log_advised: logg.warn('Did not modify X as it looks preprocessed already.')
# in other words, only if X is close to spliced will X be log-transformed (spliced and unspliced will not be transformed anyway)
# defaults: min_counts=None, min_counts_u=None, min_cells=None, min_cells_u=None,
#           min_shared_cells=None, flavor='seurat'
if genesetfile == "None":
	adata2 = scv.pp.filter_and_normalize(adata, min_shared_counts = 20, n_top_genes = 2000, copy = True, log = True)
else:
	geneset = [f.rstrip('\n') for f in open(genesetfile, 'r').readlines()]
	scv.pp.normalize_per_cell(adata)
	adata2 = adata[:, geneset]

scv.pp.moments(adata2, n_pcs = 30, n_neighbors = 30)

# compute velocity and velocity graph
#scv.tl.velocity(adata2, mode = 'steady_state')
#scv.tl.velocity(adata2, mode = 'stochastic')
scv.tl.recover_dynamics(adata2)
scv.tl.velocity(adata2, mode = 'dynamical')
scv.tl.velocity_graph(adata2)

scv.settings.figdir = plotdir + '/' + base + basegs + '/'
scv.settings.plot_prefix = base + basegs + '_scvelo_'
scv.settings.set_figure_params(dpi_save = 300, vector_friendly = True)

scv.pl.velocity_embedding_stream(adata2, basis='X_pca', save="PCA_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')
scv.pl.velocity_embedding_stream(adata2, basis='X_tsne', save="TSNE_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')
scv.pl.velocity_embedding_stream(adata2, basis='X_umap', save="UMAP_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')

scv.pl.velocity_embedding_stream(adata2, basis='PCA_starsolo_concatenated', save="PCA_starsolo_concatenated_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')
scv.pl.velocity_embedding_stream(adata2, basis='TSNE_starsolo_concatenated', save="TSNE_starsolo_concatenated_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')
scv.pl.velocity_embedding_stream(adata2, basis='UMAP_starsolo_concatenated', save="UMAP_starsolo_concatenated_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')

scv.pl.velocity_embedding_stream(adata2, basis='PCA_starsolo_summed', save="PCA_starsolo_summed_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')
scv.pl.velocity_embedding_stream(adata2, basis='TSNE_starsolo_summed', save="TSNE_starsolo_summed_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')
scv.pl.velocity_embedding_stream(adata2, basis='UMAP_starsolo_summed', save="UMAP_starsolo_summed_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')

scv.pl.velocity_embedding_stream(adata2, basis='PCA_starsolo_unspliced', save="PCA_starsolo_unspliced_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')
scv.pl.velocity_embedding_stream(adata2, basis='TSNE_starsolo_unspliced', save="TSNE_starsolo_unspliced_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')
scv.pl.velocity_embedding_stream(adata2, basis='UMAP_starsolo_unspliced', save="UMAP_starsolo_unspliced_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')

scv.pl.velocity_embedding_stream(adata2, basis='PCA_starsolo', save="PCA_starsolo_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')
scv.pl.velocity_embedding_stream(adata2, basis='TSNE_starsolo', save="TSNE_starsolo_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')
scv.pl.velocity_embedding_stream(adata2, basis='UMAP_starsolo', save="UMAP_starsolo_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')

scv.pl.velocity_embedding_stream(adata2, basis='PCA_alevin_spliced', save="PCA_alevin_spliced_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')
scv.pl.velocity_embedding_stream(adata2, basis='TSNE_alevin_spliced', save="TSNE_alevin_spliced_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')
scv.pl.velocity_embedding_stream(adata2, basis='UMAP_alevin_spliced', save="UMAP_alevin_spliced_stream.png", figsize=(7,5), size = 50, legend_fontsize = 12, show=False, color='cluster', title='')

scv.pl.velocity_graph(adata2, basis='UMAP_alevin_spliced', save='UMAP_alevin_spliced_velocitygraph.png', figsize=(12,9), show=False, color='cluster_alevin_spliced', title=mname)

scv.tl.velocity_embedding(adata2, basis = 'UMAP_alevin_spliced', all_comps = False, autoscale = False)
pd.DataFrame(adata2.obsm['velocity_UMAP_alevin_spliced'], index = adata2.obs.index).to_csv(plotdir + "/" + base + basegs + "/" + base + basegs + "_velocity_UMAP_alevin_spliced.csv")

try:
	scv.tl.recover_latent_time(adata2)
	top_genes = adata2.var_names[adata2.var.fit_likelihood.argsort()[::-1]][:300]
	scv.pl.heatmap(adata2, var_names=top_genes, tkey='latent_time', n_convolve=100, save="top_genes_heatmap.png", show=False)
	scv.pl.scatter(adata2, basis=top_genes[:10], legend_loc='none', size=80, frameon=False, ncols=5, fontsize=20, save="top_genes_scatter.png", show=False)
	scv.pl.velocity_embedding_stream(adata2, basis='X_umap', save="UMAP_stream_latent_time.png", figsize=(12,9), show=False, color='latent_time', title=mname)
	scv.pl.velocity_embedding_stream(adata2, basis='UMAP_alevin_spliced', save="UMAP_alevin_spliced_stream_latent_time.png", figsize=(12,9), show=False, color='latent_time', title=mname)
	scv.pl.scatter(adata2, basis='UMAP_alevin_spliced', save="UMAP_alevin_spliced_latent_time.png", figsize=(12,9), size=100, show=False, color='latent_time', color_map='gnuplot', perc=[2,98], rescale_color=[0,1], title=mname)
except:
	print('Latent time/top genes could not be extracted')

try:
	scv.tl.rank_velocity_genes(adata2, match_with = "cluster_alevin_spliced", n_genes = 25)
	generank = pd.DataFrame(adata2.uns['rank_velocity_genes']['names']).head(25)
	genescore = pd.DataFrame(adata2.uns['rank_velocity_genes']['scores']).head(25)
	adata2.var.to_csv(plotdir + "/" + base + basegs + "/" + base + basegs + "_gene_info.csv")
	generank.to_csv(plotdir + "/" + base + basegs + "/" + base + basegs + "_gene_rank_velocity.csv")
	genescore.to_csv(plotdir + "/" + base + basegs + "/" + base + basegs + "_gene_rank_velocity_score.csv")
	
	genes = generank.iloc[0]
	scv.pl.velocity(adata2, var_names = genes, save="velocity_genes.png", color='cluster_alevin_spliced', basis='UMAP_alevin_spliced', show=False)
except:
	print('Velocity genes could not be extracted')

scv.tl.velocity_confidence(adata2)
adata2.obs.to_csv(plotdir + "/" + base + basegs + "/" + base + basegs + "_cell_info.csv")

# session info
scv.logging.print_version()

# save 
adata2.write(filename = adfileout)