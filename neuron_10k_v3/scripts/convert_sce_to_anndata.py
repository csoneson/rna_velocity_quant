import sys
import scanpy as sc
import numpy as np
import anndata2ri

sc.settings.verbosity = 3

from rpy2.robjects import r

anndata2ri.activate()

# get input output files from command line
scefiles = None
adfiles = None

n = len(sys.argv)
if n > 2 and (n % 2) == 1:
    n1 = int((n - 1) / 2) + 1
    scefiles = sys.argv[1:n1]
    adfiles = sys.argv[n1:n]
else:
    raise ValueError("must have an even number of arguments")

print("input sce: ", scefiles, flush = True)
print("output h5ad: ", adfiles, flush = True)


# loop over inputs/output
for i in range(n1 - 1):
    scefile = scefiles[i]
    adfile = adfiles[i]

    print(f'  converting {scefile} to {adfile}...', flush = True)

    # read sce and convert it to annData object
    adata = r(f'sce <- readRDS("{scefile}")')

    # remark: assay(sce, 'counts') goes into adata.X
    #         reducedDim(sce, 'PCA') goes into adata.obsm['X_pca']
    #         reducedDim(sce, 'TSNE') goes into adata.obsm['X_tsne']

    # save annData object
    adata.write(f'{adfile}')

    print('  done\n\n', flush = True)

print('all finished', flush = True)
sc.logging.print_versions()
