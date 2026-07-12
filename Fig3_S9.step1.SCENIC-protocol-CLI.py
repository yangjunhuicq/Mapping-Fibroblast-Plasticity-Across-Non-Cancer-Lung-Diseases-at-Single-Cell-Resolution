# import dependencies
#coding=gb18030
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import sys

# set variables for file paths to read from and write to:

# set a working directory
wdir = sys.argv[1]
os.chdir( wdir )

# path to unfiltered loom file (this will be created in the optional steps below)
f_loom_path_unfilt = sys.argv[4]+"_unfiltered.loom" # test dataset, n=500 cells

# # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
f_loom_path_scenic = sys.argv[4]+"_filtered_scenic.loom"

# path to anndata object, which will be updated to store Scanpy results as they are generated below
f_anndata_path = sys.argv[4]+".anndata.h5ad"

# path to pyscenic output
f_pyscenic_output = sys.argv[4]+".pyscenic_output.loom"

# loom output, generated from a combination of Scanpy and pySCENIC results:
f_final_loom = sys.argv[4]+'_scenic_integrated-output.loom'

sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.set_figure_params(dpi=150, fontsize=10, dpi_save=600)

# Set maximum number of jobs for Scanpy.
sc.settings.njobs = 20

#Expression data import
adata = sc.read_10x_mtx(path=sys.argv[2],prefix=sys.argv[3])
#adata2=sc.read_text('testcountMat.tsv', delimiter='\t', first_column_names=True)
#f_mtx_dir = '/home/yang_junhui/test_softwares/pySCENIC/filtered_feature_bc_matrix/'
#adata = sc.read_10x_mtx(f_mtx_dir ,                 # the directory with the `.mtx` file
#    var_names='gene_symbols',   # use gene symbols for the variable names (variables-axis index)
#    cache=True) 

#write to an unfiltered loom file
#Here, we use the loompy functions directly
row_attrs = { "Gene": np.array(adata.var.index) ,}
col_attrs = { 
    "CellID":  np.array(adata.obs.index) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}

lp.create( f_loom_path_unfilt, adata.X.transpose(), row_attrs, col_attrs )

#Initial/basic filtering
import seaborn as sns
import matplotlib.pyplot as plt
# read unfiltered data from a loom file
adata = sc.read_loom( f_loom_path_unfilt )
nCountsPerGene = np.sum(adata.X, axis=0)
nCellsPerGene = np.sum(adata.X>0, axis=0)
# Show info
print("Number of counts (in the dataset units) per gene:", nCountsPerGene.min(), " - " ,nCountsPerGene.max())
print("Number of cells in which each gene is detected:", nCellsPerGene.min(), " - " ,nCellsPerGene.max())

nCells=adata.X.shape[0]

# pySCENIC thresholds
minCountsPerGene=3*.01*nCells # 3 counts in 1% of cells
print("minCountsPerGene: ", minCountsPerGene)
minSamples=.01*nCells # 1% of cells
print("minSamples: ", minSamples)

# simply compute the number of genes per cell (computers 'n_genes' column)
sc.pp.filter_cells(adata, min_genes=0)
# mito and genes/counts cuts
mito_genes = adata.var_names.str.startswith('MT-')
# for each cell compute fraction of counts in mito genes vs. all genes
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

#Diagnostic plots, pre-filtering
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=True)

x = adata.obs['n_genes']
x_lowerbound = 1500
x_upperbound = 2000
nbins=100

sns.distplot(x, ax=ax1, norm_hist=True, bins=nbins)
sns.distplot(x, ax=ax2, norm_hist=True, bins=nbins)
sns.distplot(x, ax=ax3, norm_hist=True, bins=nbins)

ax2.set_xlim(0,x_lowerbound)
ax3.set_xlim(x_upperbound, adata.obs['n_genes'].max() )

for ax in (ax1,ax2,ax3): 
  ax.set_xlabel('')

ax1.title.set_text('n_genes')
ax2.title.set_text('n_genes, lower bound')
ax3.title.set_text('n_genes, upper bound')

fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='x-large')
fig.text(0.5, 0.0, 'Genes expressed per cell', ha='center', va='center', size='x-large')

fig.tight_layout()
fig.savefig(sys.argv[4]+'.nFeature.pdf', dpi=600, bbox_inches='tight')

#Percentage of mitochondrial reads per cell
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=True)

x = adata.obs['percent_mito']
x_lowerbound = [0.0, 0.07 ]
x_upperbound = [ 0.10, 0.3 ]
nbins=100

sns.distplot(x, ax=ax1, norm_hist=True, bins=nbins)
sns.distplot(x, ax=ax2, norm_hist=True, bins=int(nbins/(x_lowerbound[1]-x_lowerbound[0])) )
sns.distplot(x, ax=ax3, norm_hist=True, bins=int(nbins/(x_upperbound[1]-x_upperbound[0])) )

ax2.set_xlim(x_lowerbound[0], x_lowerbound[1])
ax3.set_xlim(x_upperbound[0], x_upperbound[1] )
for ax in (ax1,ax2,ax3): 
  ax.set_xlabel('')

ax1.title.set_text('percent_mito')
ax2.title.set_text('percent_mito, lower bound')
ax3.title.set_text('percent_mito, upper bound')

fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='x-large')
fig.text(0.5, 0.0, 'Mitochondrial read fraction per cell', ha='center', va='center', size='x-large')

fig.tight_layout()
fig.savefig(sys.argv[4]+'.percent_mito.pdf', dpi=600, bbox_inches='tight')

#Three-panel summary plots
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=False)

sns.distplot( adata.obs['n_genes'], ax=ax1, norm_hist=True, bins=100)
sns.distplot( adata.obs['n_counts'], ax=ax2, norm_hist=True, bins=100)
sns.distplot( adata.obs['percent_mito'], ax=ax3, norm_hist=True, bins=100)

ax1.title.set_text('Number of genes expressed per cell')
ax2.title.set_text('Counts per cell')
ax3.title.set_text('Mitochondrial read fraction per cell')

fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='x-large')

fig.tight_layout()

fig.savefig(sys.argv[4]+'.filtering_panel_prefilter.pdf', dpi=600, bbox_inches='tight')

sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
    jitter=0.4, multi_panel=True )
plt.savefig(sys.argv[4]+".violin.raw.pdf")
#Scatter plot, n_genes by n_counts
sc.pl.scatter(adata, x='n_counts', y='n_genes', color='percent_mito')
plt.savefig(sys.argv[4]+".scatter.nFeature_nGenes.raw.pdf")


#Carry out the filtering steps:
# initial cuts
#sc.pp.filter_cells(adata, min_genes=200 )
#sc.pp.filter_genes(adata, min_cells=3 )
#adata = adata[adata.obs['n_genes'] < 4000, :]
#adata = adata[adata.obs['percent_mito'] < 0.15, :]

#Diagnostic plots, post filtering
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=False)
adata.obs['n_genes']
sns.distplot( adata.obs['n_genes'], ax=ax1, norm_hist=True, bins=100)
sns.distplot( adata.obs['n_counts'], ax=ax2, norm_hist=True, bins=100)
sns.distplot( adata.obs['percent_mito'], ax=ax3, norm_hist=True, bins=100)
ax1.title.set_text('Number of genes expressed per cell')
ax2.title.set_text('Counts per cell')
ax3.title.set_text('Mitochondrial read fraction per cell')
fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='x-large')
fig.tight_layout()
fig.savefig(sys.argv[4]+'.filtering_panel_postfilter.pdf', dpi=600, bbox_inches='tight')

sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
    jitter=0.4, multi_panel=True )
plt.savefig(sys.argv[4]+".violin.after_QC.pdf")
sc.pl.scatter(adata, x='n_counts', y='n_genes', color='percent_mito')
plt.savefig(sys.argv[4]+".scatter.nFeature_nGenes.after_QC.pdf")

#Update the anndata file, to be used for further processing, clustering, visualization, etc..
adata.write( f_anndata_path )
#Output the basic filtered expression matrix to a loom file.
# create basic row and column attributes for the loom file:
row_attrs = {
    "Gene": np.array(adata.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.obs_names) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}
lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs)


#SCENIC steps
#STEP 1: Gene regulatory network inference, and generation of co-expression modules
#Phase Ia: GRN inference using the GRNBoost2 algorithm  [Co-expression modules between transcription factors and candidate target genes are inferred using GENIE3 or GRNBoost.]
#For this step the CLI version of SCENIC is used. This step can be deployed on an High Performance Computing system. We use the counts matrix (without log transformation or further processing) from the loom file we wrote earlier. Output: List of adjacencies between a TF and its targets stored in ADJACENCIES_FNAME.
# transcription factors list
f_tfs = "/home/yang_junhui/test_softwares/pySCENIC/database/allTFs_hg38.txt"
#f_tfs = "/ddn1/vol1/staging/leuven/stg_00002/lcb/cflerin/resources/allTFs_hg38.txt" # human
# f_tfs = "/ddn1/vol1/staging/leuven/stg_00002/lcb/cflerin/resources/allTFs_dmel.txt" # drosophila
# f_tfs = "/ddn1/vol1/staging/leuven/stg_00002/lcb/cflerin/resources/allTFs_mm.txt"   # mouse
# tf_names = load_tf_names( f_tfs )

#ipython中的用法：   !pyscenic grn {f_loom_path_scenic} {f_tfs} -o adj.csv --num_workers 20
os.system('/home/yang_junhui/softwares/anaconda/Anaconda3-5.3.1/envs/scenic_protocol/bin/pyscenic grn %s %s  -o %s --num_workers 20 >pyscenic_grn.log.o 2>pyscenic_grn.log.e' % (f_loom_path_scenic, f_tfs, sys.argv[4]+'.adj.csv'))
#os.system('/PUBLIC/software/public/System/Perl-5.18.2/bin/perl %s -pwd %s -list %s -date %s -single %s -rawopts "%s" 2>/dev/null' % (generate_qc, analydir, infile, times, qc_opts, argv['rawopts']))
#read in the adjacencies matrix:
adjacencies = pd.read_csv(sys.argv[4]+".adj.csv", index_col=False, sep=',')
print(adjacencies.head())

#STEP 2-3: Regulon prediction aka cisTarget from CLI [RcisTarget identifies those modules where the regulator’s binding motif is significantly enriched across the target genes; and creates regulons with only direct targets. 对第一步找到的共表达module进一步分析，识别出那些“靶基因上显著富集转录因子结合motif”的TF-target module，只对存在直接作用关系的module构建regulon。]
#For this step the CLI version of SCENIC is used. This step can be deployed on an High Performance Computing system.
#Output: List of adjacencies between a TF and its targets stored in MOTIFS_FNAME.
#locations for ranking databases, and motif annotations:
import glob
# ranking databases
f_db_glob = "/home/yang_junhui/test_softwares/pySCENIC/database/*feather"
#f_db_glob = "/ddn1/vol1/staging/leuven/res_00001/databases/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/*feather"
f_db_names = ' '.join( glob.glob(f_db_glob) )

# motif databases
f_motif_path = "/home/yang_junhui/test_softwares/pySCENIC/database/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
#f_motif_path = "/ddn1/vol1/staging/leuven/res_00001/databases/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
#Here, we use the --mask_dropouts option, which affects how the correlation between TF and target genes is calculated during module creation. It is important to note that prior to pySCENIC v0.9.18, the default behavior was to mask dropouts, while in v0.9.18 and later, the correlation is performed using the entire set of cells (including those with zero expression). When using the modules_from_adjacencies function directly in python instead of via the command line, the rho_mask_dropouts option can be used to control this.
# !pyscenic ctx adj.tsv \
#    {f_db_names} \
#    --annotations_fname {f_motif_path} \
#    --expression_mtx_fname {f_loom_path_scenic} \
#    --output reg.csv \
#    --mask_dropouts \
#    --num_workers 20
os.system('/home/yang_junhui/softwares/anaconda/Anaconda3-5.3.1/envs/scenic_protocol/bin/pyscenic ctx %s %s --annotations_fname %s --expression_mtx_fname %s --output %s --mask_dropouts --num_workers 20 >pyscenic_ctx.log.o 2>pyscenic_ctx.log.e' % (sys.argv[4]+'.adj.csv', f_db_names, f_motif_path, f_loom_path_scenic, sys.argv[4]+'.reg.csv'))

#STEP 4: Cellular enrichment (aka AUCell) from CLI  [AUCell scores the activity of each regulon in each cell, yielding a binarized activity matrix. Cell states are based on the shared activity of regulatory subnetworks. 在每个细胞中，对每个regulon打分，就可以得到一个regulon activity matrix；基于这个matrix就可以识别出share regulatory subnetworks的cell states。]
#It is important to check that most cells have a substantial fraction of expressed/detected genes in the calculation of the AUC. The following histogram gives an idea of the distribution and allows selection of an appropriate threshold. In this plot, a few thresholds are highlighted, with the number of genes selected shown in red text and the corresponding percentile in parentheses). See the relevant section in the R tutorial for more information.
#By using the default setting for --auc_threshold of 0.05, we see that 1192 genes are selected for the rankings based on the plot below.
'''
nGenesDetectedPerCellbefore = np.sum(adata.X>0, axis=1)
nGenesDetectedPerCell = pd.Series(nGenesDetectedPerCellbefore)
percentiles = nGenesDetectedPerCell.quantile([.01, .05, .10, .50, 1])
print(percentiles)
sns.distplot(nGenesDetectedPerCell, norm_hist=False, kde=False, bins='fd')
for i,x in enumerate(percentiles):
    fig.gca().axvline(x=x, ymin=0,ymax=1, color='red')
    ax.text(x=x, y=ax.get_ylim()[1], s=f'{int(x)} ({percentiles.index.values[i]*100}%)', color='red', rotation=30, size='x-small',rotation_mode='anchor' )
ax.set_xlabel('# of genes')
ax.set_ylabel('# of cells')
fig.tight_layout()
fig.savefig(sys.argv[4]+'.nGenesDetectedPerCell.pdf', dpi=600, bbox_inches='tight')
'''

#!pyscenic aucell \
#    {f_loom_path_scenic} \
#    reg.csv \
#    --output {f_pyscenic_output} \
#    --num_workers 20
os.system('/home/yang_junhui/softwares/anaconda/Anaconda3-5.3.1/envs/scenic_protocol/bin/pyscenic aucell %s %s --output %s --num_workers 20 >pyscenic_aucell.log.o 2>pyscenic_aucell.log.e' % (f_loom_path_scenic, sys.argv[4]+'.reg.csv', f_pyscenic_output))

#Visualization of SCENIC's AUC matrix
#First, load the relevant data from the loom we just created
import json
import zlib
import base64

# collect SCENIC AUCell output
lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()
auc_mtx.to_csv( sys.argv[4]+".scenic_aucell_mtx.xls", sep='\t')

import umap

# UMAP
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap( auc_mtx )
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv( sys.argv[4]+".scenic_umap.txt", sep='\t')
# tSNE
tsne = TSNE( n_jobs=20 )
dr_tsne = tsne.fit_transform( auc_mtx )
pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv( sys.argv[4]+".scenic_tsne.txt", sep='\t')

#Integrate the output
#Here, we combine the results from SCENIC and the Scanpy analysis into a SCope-compatible loom file
