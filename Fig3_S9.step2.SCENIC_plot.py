import os, glob, re, pickle
from functools import partial
from collections import OrderedDict
import operator as op
from cytoolz import compose

import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss

from IPython.display import HTML, display
import sys
import math

# Set maximum number of jobs for Scanpy.
sc.settings.njobs = 32


AUCELL_MTX_FNAME = sys.argv[2]
BIN_MTX_FNAME = sys.argv[1]+'/'+sys.argv[3]+'.bin.csv'
THR_FNAME = sys.argv[1]+'/'+sys.argv[3]+'.thresholds.csv'
CELL_ANNOTATIONS_FNAME = sys.argv[4]
groupfac = sys.argv[5]

df_annotations = pd.read_csv(CELL_ANNOTATIONS_FNAME, sep='\t', index_col=0)
print(df_annotations.head())

auc_mtx = pd.read_csv(AUCELL_MTX_FNAME, index_col=0, sep='\t')
bin_mtx, thresholds = binarize(auc_mtx, num_workers=20)
bin_mtx.to_csv(BIN_MTX_FNAME, sep='\t')
thresholds.to_frame().rename(columns={0:'threshold'}).to_csv(THR_FNAME, sep='\t')

bin_mtx = pd.read_csv(BIN_MTX_FNAME, index_col=0, sep='\t')
thresholds = pd.read_csv(THR_FNAME, index_col=0, sep='\t').threshold


#Create heatmap with binarized regulon activity.

df_annotations_use=df_annotations.loc[auc_mtx.index.values, : ]


#Save clustered binarized heatmap to Excel for further inspection.
bin_mtx_clustered = bin_mtx.T.copy()
bin_mtx_clustered.rename(columns=df_annotations_use.loc[:,groupfac].to_dict(), inplace=True)
bin_mtx_clustered.to_csv( sys.argv[1]+'/'+sys.argv[3]+'.bin.clustered.csv', sep='\t')

rss = regulon_specificity_scores( auc_mtx, df_annotations_use.loc[:,groupfac] )
print(rss.head())
rss.to_csv( sys.argv[1]+'/'+sys.argv[3]+'.rss.csv', sep='\t')
#Select the top 5 regulons from each cell type
cats = sorted(list(set(df_annotations_use.loc[:,groupfac])))
topreg = []
for i,c in enumerate(cats):
    topreg.extend(
        list(rss.T[c].sort_values(ascending=False)[:5].index)
    )
print(cats)
print(topreg)
print(len(topreg))
topreg = list(set(topreg))
print(len(topreg))

nr=math.ceil(len(cats)/5)
fig = plt.figure(figsize=(15, 8))
for c,num in zip(cats, range(1,len(cats)+1)):
    x=rss.T[c]
    ax = fig.add_subplot(nr,5,num)
    plot_rss(rss, c, top_n=5, max_n=None, ax=ax)
    ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    #adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )
 
fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
        'figure.titlesize': 'large' ,
        'axes.labelsize': 'medium',
        'axes.titlesize':'large',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'
        })
plt.savefig(sys.argv[1]+'/'+sys.argv[3]+".plots_rss_top5.pdf", dpi=600, bbox_inches = "tight")

#Generate a Z-score for each regulon to enable comparison between regulons
auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
for col in list(auc_mtx.columns):
    auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)
#auc_mtx_Z.sort_index(inplace=True)
auc_mtx_Z.to_csv(sys.argv[1]+'/'+sys.argv[3]+".auc_mtx_Zcore.csv", sep='\t')


nr=math.ceil(len(topreg)/5)
fig, axs = plt.subplots(nr, 5, figsize=(15, 10), dpi=100, sharey=False)
for i,ax in enumerate(axs): #i is the row of plot
  for j in range(0,len(ax)):
    if i*5+j >= len(topreg):break
    sns.distplot(auc_mtx[ topreg[i*5+j] ], ax=ax[j], norm_hist=True, bins=100)
    ax[j].plot( [ thresholds[ topreg[i*5+j] ] ]*2, ax[j].get_ylim(), 'r:')
    ax[j].title.set_text( topreg[i*5+j] )
    ax[j].set_xlabel('')
    
fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='large')
fig.text(0.5, -0.01, 'AUC', ha='center', va='center', rotation='horizontal', size='large')

fig.tight_layout()
fig.savefig(sys.argv[1]+'/'+sys.argv[3]+'.cellType-binaryPlot2.pdf', dpi=600, bbox_inches='tight')



