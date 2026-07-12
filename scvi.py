import os
import tempfile

import scanpy as sc
import scvi
import seaborn as sns
import torch
from rich import print
from scib_metrics.benchmark import Benchmarker
import numpy as np

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()

import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import pandas as pd

adata = sc.read('MyoFtest.h5ad')
print(adata)
print(adata.obs.head())
print(adata.var.head())

print(adata.raw.X.shape)
print(adata.X.shape)

print(np.max(adata.X))
print(np.min(adata.X))
print(np.min(adata.raw.X))
print(np.max(adata.raw.X))

adata.layers["counts"] = adata.raw.X
adata.raw = adata  # keep full dimension safe
print(adata)
print(adata.obs.head())
print(adata.var.head())
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2000,
    layer="counts",
    batch_key=None,
    subset=False,
)
print(adata)
print(adata.obs.head())
print(adata.var.head())

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="orig.ident")
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train()

print(adata)
print(adata.obs.head())
print(adata.var.head())


SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

print(adata)

sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
print(adata)

sc.tl.umap(adata)
print(adata)
sc.pl.umap(
    adata,
    color=["percent.mt", "nFeature_RNA"],
    wspace=0.5,
    ncols=2,
    legend_loc="on data",
)
plt.savefig("MyoFtest.umap.MTperc_nFeatureRNA.pdf")

sc.tl.embedding_density(adata, groupby='Group2', basis='umap') #computing density on 'umap'
print(adata)
with rc_context({'figure.figsize': (20, 15)}):
    sc.pl.embedding_density(adata, key='umap'+'_density_'+'Group2', basis='umap')
plt.savefig("MyoFtest.scVI.DensityPlot."+'Group2'+".sc."+'X_umap'+".pdf")


for res in [0.1, 0.2, 0.3, 0.5]:
    print(res)
    resolution_key = f"leiden_res_{res:4.2f}"
    sc.tl.leiden(adata, key_added=resolution_key, resolution=res)#Please install the leiden algorithm: `conda install -c conda-forge leidenalg` or `pip3 install leidenalg`.
    # Obtain cluster-specific differentially expressed genes
    sc.tl.rank_genes_groups(adata, groupby=resolution_key, method="wilcoxon")
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
    plt.savefig("pl_top20markers.leiden_res_"+str(res)+".pdf")
    print(adata)
    print(adata.obs.head())
    print(adata.var.head())
    # save DEG file
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    df = pd.DataFrame(
        {group + '_' + key: result[key][group]
         for group in groups for key in ['names', 'pvals', 'pvals_adj', 'logfoldchanges']}
    )
    df.to_csv(f"rank_genes_groups_{resolution_key}.csv",sep='\t')

print(adata)
print(adata.obs.head())
print(adata.var.head())

sc.pl.umap(
    adata,
    color=["leiden_res_0.10", "leiden_res_0.20", "leiden_res_0.30", "leiden_res_0.50"],
    wspace=0.5,
    ncols=2,
    legend_loc="on data",
)
plt.savefig("MyoFtest.umap.cluster_res.pdf")


metadatanew2=adata.obs
metadatanew2.to_csv("MyoFtest.scVI.h5ad"+'.metadata.xls', index=True, sep='\t')
adata.write_h5ad("MyoFtest.scVI.h5ad")
