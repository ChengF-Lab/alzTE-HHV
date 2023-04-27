import logging, matplotlib, os, sys
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import pandas as pd
from glbase3 import genelist
plt.rcParams['figure.figsize']=(8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=300, dpi_save=300)
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['font.size']=10

sc.settings.figdir = 'diffexp'

adata = sc.read('../02_norm/learned.h5ad')

adata.obs['new_clusters'] = (adata.obs['leiden_r0.25'].map(lambda x: {"0":"Microglia","1": "Astrocyte","2":"Astrocyte","3":"OLs","4":"Mural_cells","5":"Endothelial","6":"Microglia","7":"Microglia"}.get(x, x)).astype("category"))

marker_genes = ["GFAP","SLC1A2","AQP4","ATP1B2","CSF1R","C3","CIITA","CX3CR1","MAP1B","RBFOX3","PLP1","MOBP","OPALIN","CLDN5"]


######### diff genes and TEs for new_clusters
adata.uns['log1p']["base"] = None
sc.tl.rank_genes_groups(adata, 'new_clusters', method='wilcoxon')

sc.pl.tsne(adata, color=["ZNF146"],color_map="YlGnBu",save="_OTC_ZNF146_exp_tsne.pdf",size = 2)
sc.pl.violin(adata,["ZNF146"], groupby='new_clusters',save="_OTC_ZNF146_exp_violin.pdf")


