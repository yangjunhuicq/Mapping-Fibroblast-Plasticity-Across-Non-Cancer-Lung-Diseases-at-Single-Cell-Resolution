# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import anndata as ad
from matplotlib.pyplot import rc_context
import os
import argparse
import getpass
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description="using scanpy to show density plot",formatter_class=RawTextHelpFormatter)
parser.add_argument('--infile',help="h5ad file",required=True)
parser.add_argument('--outdir',help="output path(default: ./)",default=os.path.abspath('./'))
parser.add_argument('--outpre',help="the prefix for output files",default='out.')
argv = vars(parser.parse_args())


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

os.chdir(argv['outdir'])

adata_combined = sc.read_h5ad(argv['infile'])
print(adata_combined)
print(adata_combined.obs.head())
print(adata_combined.var.head())
print(adata_combined.var_names)
adata_combined.var['features'] = adata_combined.var_names
print(adata_combined.var.head())

metadata = pd.read_csv('clinicalinfo.all.xls',delimiter='\t')
metadata.set_index('orig.ident', inplace=True)
print(metadata.head())

adata_combined.obs = adata_combined.obs.join(metadata, on='orig.ident')
metadatanew=adata_combined.obs
metadatanew.to_csv(argv['infile']+'.metadata.new.csv', index=True, sep='\t')
print(adata_combined.obs.head())

print(adata_combined.var.head())

obs_df = adata_combined.obs
X = adata_combined.X

results = {}
for gse_alias, group in obs_df.groupby('gse_alias'):
    indices = group.index
    mask = adata_combined.obs_names.isin(indices)
    X_subset = X[mask, :]
    max_value = X_subset.max()
    min_value = X_subset.min()
    results[gse_alias] = {'max_value': max_value, 'min_value': min_value}

for gse_alias, values in results.items():
    print(f"{gse_alias}: max_value = {values['max_value']}, min_value = {values['min_value']}")
X = False
metadataori = False
metadatanew = False

print(np.max(adata_combined.X))
print(np.max(adata_combined.raw.X))
print(adata_combined.X.shape)
print(adata_combined.raw.X.shape)
adata_combined.X = adata_combined.raw.X
print(np.max(adata_combined.X))
print(np.max(adata_combined.raw.X))
print(adata_combined.obs.head())
print(adata_combined.var.head())

adata_combined.var['mt'] = adata_combined.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'. Add metadata for each feature, i.e. set feature startswith MT- True.
print(adata_combined.var.loc[adata_combined.var['mt'] == True, 'features'])
print(adata_combined.var.head())

sc.pp.calculate_qc_metrics(adata_combined, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)
#scanpy.pp.calculate_qc_metrics(adata, *, expr_type='counts', var_type='genes', qc_vars=(), percent_top=(50, 100, 200, 500), layer=None, use_raw=False, inplace=False, log1p=True, parallel=None). log1p=True时，会对obs中的'n_genes_by_counts', 'total_counts', 'total_counts_mt'进行Loge计算，对var中的'mean_counts'和'total_counts'进行Loge计算。
print(adata_combined)
print(adata_combined.obs.head())
print(adata_combined.var.head())

sc.pl.violin(adata_combined, keys = ["n_genes_by_counts", "total_counts", "pct_counts_mt"], groupby= "gse_alias", jitter=0.4, multi_panel=True, rotation=90)
plt.savefig(argv['outpre']+".QC.by.gse_alias.png")

plt.subplots(2, 2)
sc.pl.scatter(adata_combined, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata_combined, x='total_counts', y='n_genes_by_counts')
plt.savefig(argv['outpre']+".QC.scatter.pdf")
# 计算 Pearson 相关系数
x = adata_combined.obs['total_counts']
y = adata_combined.obs['pct_counts_mt']
from scipy.stats import pearsonr
corr, _ = pearsonr(x, y)
print(corr)

# Saving count data
adata_combined.layers["counts"] = adata_combined.X.copy()
print(np.max(adata_combined.X))

#same to Seurat NormalizeData()
sc.pp.normalize_total(adata_combined, target_sum=1e4) #Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells.
print(np.max(adata_combined.X))
# Logarithmize the data
sc.pp.log1p(adata_combined)
print(np.max(adata_combined.X))

sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000, batch_key=None)
sc.pl.highly_variable_genes(adata_combined)
plt.savefig(argv['outpre']+".highly_variable_genes.pdf")

print(adata_combined)
print(adata_combined.obs.head())
print(adata_combined.var.head())





#Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. This simply freezes the state of the AnnData object.
#You can get back an AnnData of the object in .raw by calling .raw.to_adata().
adata_combined.raw = adata_combined #此步操作后，adata_combined.raw包含的是所有features的normalized data；因为下一步要对adata_combined筛选仅高可变基因。后面sc.pl.umap()函数中，当 use_raw=False 时，Scanpy 将从 .X 中获取数据，当 use_raw=True 时，是从 .raw.X 中获取。

#Actually do the filtering
adata_combined = adata_combined[:, adata_combined.var.highly_variable]
print(adata_combined)

#same to Seurat ScaleData()
#Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.
#sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
#Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata_combined, max_value=10) #此时adata_combined.X是scaled.data，adata_combined.raw是normalized data
print(adata_combined)
print(np.max(adata_combined.X))
print(np.min(adata_combined.X))
print(np.max(adata_combined.raw.X))
print(np.min(adata_combined.raw.X))


#RunPCA
sc.tl.pca(adata_combined, svd_solver='arpack', use_highly_variable=True)

sc.pl.pca(adata_combined, color=['gse_alias', 'gse_alias', "n_genes_by_counts", "n_genes_by_counts"], dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)], ncols=2, size=2)
plt.savefig(argv['outpre']+".pca.pdf")
sc.pl.pca_variance_ratio(adata_combined, log=True)
plt.savefig(argv['outpre']+".pca_variance_ratio.pdf")
print(adata_combined)

sc.external.pp.harmony_integrate(adata_combined, key='orig.ident')
#adata.write(results_file)
print(adata_combined)
print(adata_combined.obs.head())
print(adata_combined.var.head())

sc.pp.neighbors(adata_combined, use_rep='X_pca_harmony')
print(adata_combined)

sc.tl.umap(adata_combined)
print(adata_combined)
#sc.pl.umap(pbmc, color=['CD79A', 'MS4A1', 'IGJ', 'CD3D', 'FCER1A', 'FCGR3A', 'n_counts', 'bulk_labels'], s=50, frameon=False, ncols=4, vmax='p99')
sc.pl.umap(adata_combined, color=['EPCAM', 'COL1A1', 'PECAM1','PTPRC']) #same to sc.pl.umap(adata_combined, color=['EPCAM', 'COL1A1', 'PECAM1','PTPRC'], use_raw=True). Use the .raw.X for coloring with gene expression. For the scatter plots, the value to plot is given as the color argument. This can be any gene or any column in .obs, where .obs is a DataFrame containing the annotations per observation/cell, see AnnData for more information.
plt.savefig(argv['outpre']+".umap.EPCAM_COL1A1_PECAM_PTPRC.pdf")
#sc.pl.umap(adata_combined, color=['EPCAM', 'COL1A1', 'PECAM1','PTPRC'], use_raw=False) #当 use_raw=False 时，Scanpy 将从 .X 中获取数据，而不是从 .raw.X 中获取。确保你感兴趣的基因存在于当前使用的数据表示中。如果前面筛选了highly_variable基因，那么当指定的基因不在highly_variable基因中时，就会报错。
#plt.savefig("umap.EPCAM_COL1A1_PECAM_PTPRC.not_raw.pdf")



for res in [0.1, 0.2, 0.3, 0.5, 1.0, 1.5, 2.0]:
    print(res)
    resolution_key = f"leiden_res_{res:4.2f}"
    sc.tl.leiden(adata_combined, key_added=resolution_key, resolution=res)#Please install the leiden algorithm: `conda install -c conda-forge leidenalg` or `pip3 install leidenalg`.
    # Obtain cluster-specific differentially expressed genes
    sc.tl.rank_genes_groups(adata_combined, groupby=resolution_key, method="wilcoxon")
    sc.pl.rank_genes_groups(adata_combined, n_genes=20, sharey=False)
    plt.savefig("pl_top20markers.leiden_res_"+str(res)+".pdf")
    print(adata_combined)
    print(adata_combined.obs.head())
    print(adata_combined.var.head())
    # 保存差异基因到文件
    result = adata_combined.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    df = pd.DataFrame(
        {group + '_' + key: result[key][group]
         for group in groups for key in ['names', 'pvals', 'pvals_adj', 'logfoldchanges']}
    )
    df.to_csv(f"rank_genes_groups_{resolution_key}.csv",sep='\t')
    #
    # 获取top 20 markers用于violin图
    #top20_genes = {group: result['names'][group][:20] for group in groups}
    #top20_genes_flat = [gene for sublist in top20_genes.values() for gene in sublist]
    #
    # 绘制violin图
    #sc.pl.violin(adata_combined, top20_genes_flat, groupby=resolution_key)
    #plt.savefig("top20markers_violin.leiden_res_"+str(res)+".pdf")
    #sc.pl.rank_genes_groups_dotplot(adata_combined, groupby=resolution_key, standard_scale="var", n_genes=5)
    #plt.savefig("top5markers_rankdotplot.leiden_res_"+str(res)+".pdf")
    #sc.pl.rank_genes_groups_violin(adata_combined, n_genes=20, jitter=False)
    #plt.savefig("top20markers_rankviolin.leiden_res_"+str(res)+".pdf")
print(adata_combined)
print(adata_combined.obs.head())
print(adata_combined.var.head())


sc.pl.umap(
    adata_combined,
    color=["leiden_res_0.10", "pct_counts_mt", "n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
    legend_loc="on data",
)
plt.savefig(argv['outpre']+".umap.cluster_MTperc_nFeatureRNA.pdf")


sc.pl.umap(
    adata_combined,
    color=["leiden_res_0.10", "leiden_res_0.20", "leiden_res_0.30", "leiden_res_0.50", "leiden_res_1.00", "leiden_res_1.50", "leiden_res_2.00"],
    legend_loc="on data",
)
plt.savefig(argv['outpre']+".umap.clusters_res.pdf")

#marker_genes = {
#    "CD14+ Mono": ["FCN1", "CD14"],
#    "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
#}
marker_genes = {}
with open('markers3.use.add', 'r') as f:
    for line in f:
        if line.strip():  # 确保不是空行
            gene, cell_type = line.strip().split('\t')
            gene = gene.upper()
            if gene in adata_combined.raw.var_names:
                if cell_type not in marker_genes:
                    marker_genes[cell_type] = []
                marker_genes[cell_type].append(gene)
sc.pl.dotplot(adata_combined, marker_genes, groupby="leiden_res_0.30", standard_scale="var")
plt.savefig(argv['outpre']+".celltype_markers.leiden_res_0.30.pdf")

metadatanew2=adata_combined.obs
metadatanew2.to_csv(argv['outpre']+".ana.h5ad"+'.metadata.xls', index=True, sep='\t')
adata_combined.write_h5ad(argv['outpre']+".ana.h5ad")

