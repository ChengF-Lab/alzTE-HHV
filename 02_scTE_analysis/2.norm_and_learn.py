import logging, matplotlib, os, sys
import anndata
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=300)

adata = sc.read('../01_pack/raw_data.h5ad')

mito_genes=adata.var_names.str.startswith('MT-')
adata.obs['percent_mito']=np.sum(adata[:,mito_genes].X,axis=1).A1/np.sum(adata.X,axis=1).A1
adata.obs['n_counts']=adata.X.sum(axis=1).A1
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True)
plt.savefig("QC_violin.pdf")
sc.pl.scatter(adata, x='n_counts', y='percent_mito')
plt.savefig("QC_dot.pdf")
sc.pl.scatter(adata,x='n_counts',y='n_genes')
plt.savefig("QC_dot1.pdf")
adata = adata[adata.obs.percent_mito < 0.05, :]


sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
print(adata)

print('Number of cells: {:d}'.format(adata.n_obs))

sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=2000)
sc.pl.highly_variable_genes(adata, show=False, save='highly_variable.pdf')

# Calculate the visualizations
sc.pp.pca(adata, n_comps=30, use_highly_variable=True, svd_solver='arpack') # PC=20 from Nature paper
sc.pp.neighbors(adata)
sc.tl.tsne(adata, n_jobs=3)
sc.tl.umap(adata, min_dist=0.6)
sc.tl.diffmap(adata)

sc.pl.pca_variance_ratio(adata, log=True, show=False, save='pca_variance.pdf')

# Perform clustering - using highly variable genes
sc.tl.leiden(adata, resolution=1.0, key_added='leiden_r1')
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_r0.5')
sc.tl.leiden(adata, resolution=0.4, key_added='leiden_r0.4')
sc.tl.leiden(adata, resolution=0.35, key_added='leiden_r0.35')
sc.tl.leiden(adata, resolution=0.3, key_added='leiden_r0.3')
sc.tl.leiden(adata, resolution=0.25, key_added='leiden_r0.25')
sc.tl.leiden(adata, resolution=0.2, key_added='leiden_r0.2')
sc.tl.leiden(adata, resolution=0.1, key_added='leiden_r0.1')

adata.write('./learned.h5ad')

todraw = ['leiden_r1', 'leiden_r0.5', 'leiden_r0.4', 'leiden_r0.35', 'leiden_r0.3', 'leiden_r0.25', 'leiden_r0.2', 'leiden_r0.1', 'replicate']

#Visualize the clustering and how this is reflected by different technical covariates
sc.pl.tsne(adata, color=todraw, size=20, show=False, save='tsne.pdf')
sc.pl.umap(adata, color=todraw, size=20, show=False, save='umap.pdf')

finaldraw = ['leiden_r0.25', 'replicate']
sc.pl.tsne(adata, color=finaldraw, size=20, show=False, save='final_tsne.pdf')
sc.pl.umap(adata, color=finaldraw, size=20, show=False, save='final_umap.pdf')


