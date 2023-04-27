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

sc.pl.umap(adata, color='new_clusters', legend_loc='on data', title='', frameon=True,save="_celltype_cluster_label.pdf")

sc.pl.umap(adata, color='new_clusters', legend_loc='right margin', title='', frameon=True,save="_celltype_cluster_nolabel.pdf")

sc.pl.tsne(adata, color='new_clusters', legend_loc='on data',size=2,show=False,save="_celltype_cluster_label.pdf")

sc.pl.tsne(adata, color='new_clusters', legend_loc='right margin',size=2, show=False,save="_celltype_cluster_nolabel.pdf")


######### diff genes and TEs for new_clusters
adata.uns['log1p']["base"] = None
sc.tl.rank_genes_groups(adata, 'new_clusters', method='wilcoxon', n_genes=2000)
adata.write('./de.h5ad')

adata = sc.read('./de.h5ad')

sc.pl.rank_genes_groups(adata, n_genes=50, sharey=True, show=False, save='genes-top25.pdf')
sc.pl.rank_genes_groups(adata, key='rank_genes_groups', show=False, save='genes.pdf')
sc.pl.rank_genes_groups_dotplot(adata, key='rank_genes_groups', show=False, save='genes-top25.pdf')

#print(pd.DataFrame(adata.uns['rank_genes_groups']))

print(pd.DataFrame(adata.uns['rank_genes_groups']['names']))

print()

pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
adata.uns['rank_genes_groups'].keys()
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'pvals']}).head(5)
res = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'pvals','logfoldchanges','pvals_adj','scores']})
res.to_csv("dif_celltype.csv")


topall = pd.DataFrame(adata.uns['rank_genes_groups']['names']) # get all;
fcs = pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges'])
padj = pd.DataFrame(adata.uns['rank_genes_groups']['pvals_adj'])
topall.to_csv('top100_celltype.csv')


sc.pl.violin(adata, ["GFAP","SLC1A2","LTR16","L1P4a","CLDN5","LTR7B","L1MEa","CSF1R","C3","CIITA","L1MCa","L1M6B","PDGFRB","MAP1B","RBFOX3","L1MCa","L1PA12","L1M3d","PLP1","MOBP","L1MA6","L1M3","LTR16A","LTR8"], groupby='new_clusters',save="_OC_marker_exp_Vio.pdf")


sc.pl.tsne(adata, color=["GFAP","SLC1A2","LTR16","L1P4a","CLDN5","LTR7B","L1MEa","CSF1R","C3","CIITA","L1MCa","L1M6B","PDGFRB","MAP1B","RBFOX3","L1MCa","L1PA12","L1M3d","PLP1","MOBP","L1MA6","L1M3","LTR16A","LTR8"],save="_OC_marker_exp_tsne.pdf")


sc.pl.dotplot(adata, ["GFAP","SLC1A2","LTR16","L1P4a","CLDN5","LTR7B","L1MEa","CSF1R","C3","CIITA","L1MCa","L1M6B","PDGFRB","MAP1B","RBFOX3","L1MCa","L1PA12","L1M3d","PLP1","MOBP","L1MA6","L1M3","LTR16A","LTR8"], groupby='new_clusters',save="_OC_marker_exp_dotplot.pdf")


sc.pl.stacked_violin(adata, ["GFAP","SLC1A2","LTR16","L1P4a","CLDN5","LTR7B","L1MEa","CSF1R","C3","CIITA","L1MCa","L1M6B","PDGFRB","MAP1B","RBFOX3","L1MCa","L1PA12","L1M3d","PLP1","MOBP","L1MA6","L1M3","LTR16A","LTR8"], groupby='new_clusters',save="_OC_marker_exp_stacked_violin.pdf")



#adata.obs['new_clusters'] = (adata.obs['leiden_r0.4'].map(lambda x: {"0": "Astrocyte", "1": "Astrocyte","2": "Astrocyte","3": "Microglia","4": "Neuron","5": "Microglia","6": "OLs","7": "Astrocyte","8": "Endothelial","9": "OLs","10": "Neuron"}.get(x, x)).astype("category"))


# Go through and trim the TEs:
#TEs = set(genelist(filename='../../TE_genes_id.hg38.txt', format={'name': 0, 'force_tsv': True})['name'])

#newcols = {}


