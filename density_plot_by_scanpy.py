#!/usr/bin/env python
#coding=gb18030
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import os
import argparse
import getpass
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description="using scanpy to show density plot",formatter_class=RawTextHelpFormatter)
parser.add_argument('--infile',help="h5ad file",required=True)
parser.add_argument('--outdir',help="output path(default: ./)",default=os.path.abspath('./'))
parser.add_argument('--groupby',help="A column name in the metadata, and the cells will be grouped by it.",required=True)
parser.add_argument('--outpre',help="the prefix for output files",default='out.')
parser.add_argument('--reduction_name',help="A reduction name in the object, such as umap, umapbyharmony15.",default='umap')
#parser.add_argument('--seqstrag',help="WGS, WES_illu, WES_agi or TS",default='WGS',choices=['WGS','WES_agi','WES_illu','TS'])
argv = vars(parser.parse_args())


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

os.chdir(argv['outdir'])

adata = sc.read(argv['infile'])

print(adata)

celltype_counts = adata.obs[argv['groupby']].value_counts()
print(celltype_counts)
celltypes_gt_100 = celltype_counts[celltype_counts > 100]
print(celltypes_gt_100)

adata=adata[adata.obs[argv['groupby']].isin(celltypes_gt_100.index)]
print(adata)

print(argv['groupby'])
print(argv['reduction_name'])
sc.tl.embedding_density(adata, groupby=argv['groupby'], basis=argv['reduction_name']) #computing density on 'umap'
with rc_context({'figure.figsize': (20, 15)}):
    sc.pl.embedding_density(adata, key=argv['reduction_name']+'_density_'+argv['groupby'], basis=argv['reduction_name'], bg_dotsize=160, fg_dotsize=360)
plt.savefig(argv['outpre']+".DensityPlot."+argv['groupby']+".sc."+argv['reduction_name']+".pdf")


