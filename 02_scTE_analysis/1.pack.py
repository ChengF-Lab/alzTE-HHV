"""

Pack the scRNA-seq data using scanpy, prep for scran normalisation

"""

import logging, matplotlib, os, sys
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from anndata import AnnData
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
from rpy2.robjects.packages import importr
plt.rcParams['figure.figsize'] = (8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 10
sc.settings.autoshow = False

def sparsify(filename):
    data = pd.read_csv(filename, index_col=0, header=0)
    genes = data.columns
    cells = data.index
    data = sp.sparse.csr_matrix(data.to_numpy())
    data.astype('float32')

    '''
    oh = open('gene_names.{0}.tsv'.format(os.path.split(filename)[1]), 'w')
    for g in genes:
        oh.write('%s\n' % g)
    oh.close()
    '''

    print('Loaded {0}'.format(filename))
    ad = AnnData(data, obs={'obs_names': cells}, var={'var_names': genes})
    del data
    return ad

######### Change the input files according your analysis #########

sam1 = sparsify("../../00_scTE/Donor10OTC.csv.gz")    ; sam1.obs['stage'] = "virus.OTC.AD"   ; sam1.obs['replicate'] = "virus.OTC.AD-1"
sam2 = sparsify("../../00_scTE/Donor11OTC.csv.gz")    ; sam2.obs['stage'] = "OTC.AD"   ; sam2.obs['replicate'] = "OTC.AD-1"
sam3 = sparsify("../../00_scTE/Donor12OTC.csv.gz")   ; sam3.obs['stage'] = "OTC.AD"   ; sam3.obs['replicate'] = "OTC.AD-2"
sam4 = sparsify("../../00_scTE/Donor13OTC.csv.gz")   ; sam4.obs['stage'] = "OTC.AD"  ; sam4.obs['replicate'] = "OTC.AD-3"
sam5 = sparsify("../../00_scTE/Donor14OTC.csv.gz")   ; sam5.obs['stage'] = "OTC.AD"   ; sam5.obs['replicate'] = "OTC.AD-4"
sam6 = sparsify("../../00_scTE/Donor15OTC.csv.gz")   ; sam6.obs['stage'] = "OTC.AD"   ; sam6.obs['replicate'] = "OTC.AD-5"
sam7 = sparsify("../../00_scTE/Donor16OTC.csv.gz")   ; sam7.obs['stage'] = "OTC.AD"   ; sam7.obs['replicate'] = "OTC.AD-6"
sam8 = sparsify("../../00_scTE/Donor17OTC.csv.gz")   ; sam8.obs['stage'] = "OTC.AD"   ; sam8.obs['replicate'] = "OTC.AD-7"
sam9 = sparsify("../../00_scTE/Donor18OTC.csv.gz")   ; sam9.obs['stage'] = "OTC.AD"   ; sam9.obs['replicate'] = "OTC.AD-8"
sam10 = sparsify("../../00_scTE/Donor1OTC.csv.gz") ; sam10.obs['stage'] = "OTC.CTR" ; sam10.obs['replicate'] = "OTC.CTR-1"
sam11 = sparsify("../../00_scTE/Donor2OTC.csv.gz") ; sam11.obs['stage'] = "OTC.CTR" ; sam11.obs['replicate'] = "OTC.CTR-2"
sam12 = sparsify("../../00_scTE/Donor3OTC.csv.gz") ; sam12.obs['stage'] = "OTC.CTR" ; sam12.obs['replicate'] = "OTC.CTR-3"
sam13 = sparsify("../../00_scTE/Donor5OTC.csv.gz")   ; sam13.obs['stage'] = "OTC.CTR"  ; sam13.obs['replicate'] = "OTC.CTR-4"
sam14 = sparsify("../../00_scTE/Donor6OTC.csv.gz")   ; sam14.obs['stage'] = "OTC.CTR"  ; sam14.obs['replicate'] = "OTC.CTR-5"
sam15 = sparsify("../../00_scTE/Donor7OTC.csv.gz")   ; sam15.obs['stage'] = "OTC.CTR"  ; sam15.obs['replicate'] = "OTC.CTR-6"
sam16 = sparsify("../../00_scTE/Donor8OTC.csv.gz")   ; sam16.obs['stage'] = "OTC.CTR"  ; sam16.obs['replicate'] = "OTC.CTR-7"
sam17 = sparsify("../../00_scTE/Donor9OTC.csv.gz")  ; sam17.obs['stage'] = "OTC.AD"  ; sam17.obs['replicate'] = "OTC.AD-9"

print('Loaded Samples...')

# Do very simple prefiltering:
samples = [sam1, sam2, sam3, sam4,sam5, sam6,
           sam7, sam8, sam9, sam10,sam11,sam12, 
		   sam13, sam14, sam15,sam16, sam17]

# Quick pre-filtering, these should be low, otherwise it can mess up downstream analysis, but also can get rid of trivial uninteresting things
[sc.pp.filter_cells(sam, min_genes=200) for sam in samples]
[sc.pp.filter_cells(sam, max_counts=100000) for sam in samples]
[sc.pp.filter_cells(sam, min_counts=1000) for sam in samples]
# Do not filter gene here; concatenate joins on the union, so if a gene fails in a single sample, it will also be deleted from all other samples;

print('Concatenating')
adata = sam1.concatenate(samples[1:])

del samples

adata.X = adata.X.astype('float32')

print(adata)

sc.pl.violin(adata, ['n_genes', 'n_counts'], groupby='replicate', size=0, log=False, cut=0, show=False, save='qc1-pre-norm-replicates.pdf')

# Base filtering for trivial QC failures:
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_cells(adata, min_counts=1500)
sc.pp.filter_cells(adata, max_counts=100000)
sc.pp.filter_genes(adata, min_cells=50) # Only filter genes here;

print('Number of cells after gene filter: {:d}'.format(adata.n_obs))

#sc.pl.violin(adata, ['n_genes','n_counts'], groupby='stage', size=0, log=False, cut=0, show=False, save='qc1.pdf')
sc.pl.violin(adata, ['n_genes','n_counts'], groupby='replicate', size=0, log=False, cut=0, show=False, save='qc1-replicates.pdf')

p = sb.distplot(adata.obs['n_counts'], kde=False)
p.get_figure().savefig('figures/distplot_ncounts1.pdf')
p = sb.distplot(adata.obs['n_counts'][adata.obs['n_counts']<4000], kde=False, bins=60)
p.get_figure().savefig('figures/distplot_ncounts2.pdf')
p = sb.distplot(adata.obs['n_counts'][adata.obs['n_counts']>10000], kde=False, bins=60)
p.get_figure().savefig('figures/distplot_ncounts3.pdf')
#Thresholding decision: genes
p = sb.distplot(adata.obs['n_genes'], kde=False, bins=60)
p.get_figure().savefig('figures/distplot_ngenes1.pdf')
p = sb.distplot(adata.obs['n_genes'][adata.obs['n_genes']<2000], kde=False, bins=60)
p.get_figure().savefig('figures/distplot_ngenes2.pdf')

print('Total number of cells: {:d}'.format(adata.n_obs))
print('Total number of genes: {:d}'.format(adata.n_vars))

adata.write('./raw_data.h5ad')
