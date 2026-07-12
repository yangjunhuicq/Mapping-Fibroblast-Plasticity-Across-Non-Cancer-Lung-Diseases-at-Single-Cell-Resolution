#/home/lijia/yangjunhui/softwares/miniconda3/envs/scCODA/bin/python
import importlib
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import warnings
import pickle as pkl
from matplotlib.pyplot import rc_context

from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
from sccoda.util import comp_ana as mod

import sccoda.datasets as scd

warnings.filterwarnings("ignore")
import numpy as np

adata = ad.read_h5ad("/home/lijia/yangjunhui/other/Revision/sccoda/Figure1F/COVID19_Donor_PF_COPD.fib.ana.for_sccoda.h5ad")
print(adata)
print(adata.obs['Normal'].value_counts())
pd.set_option('display.max_rows', None)
print(adata.obs['orig.ident'].value_counts())
print(np.max(adata.X))
print(np.min(adata.X))
print(np.max(adata.raw.X))
print(np.min(adata.raw.X))

celltype_counts = adata.obs['orig.ident'].value_counts()
celltypes_gt_100 = celltype_counts[celltype_counts >= 60]
print(celltypes_gt_100)
adata=adata[adata.obs['orig.ident'].isin(celltypes_gt_100.index)]
print(adata)
print(adata.obs['orig.ident'].value_counts())

sample_covariate = adata.obs.groupby('orig.ident')['Normal'].first()
# 转换为 DataFrame
covariate_df = pd.DataFrame(sample_covariate)
print(covariate_df)

data_scanpy_1 = dat.from_scanpy(
    adata,
    cell_type_identifier="celltype",
    sample_identifier="orig.ident",
    covariate_df=covariate_df
)
print(data_scanpy_1)
print(data_scanpy_1.obs)
print(data_scanpy_1.var)
print(data_scanpy_1.X)

# Stacked barplot for each sample
data_scanpy_1.obs['origidents'] = data_scanpy_1.obs.index
print(data_scanpy_1.obs)

with rc_context({'figure.figsize': (20, 5)}):
    viz.stacked_barplot(data_scanpy_1, feature_name="origidents")
plt.savefig("stacked_barplot.orig.ident.pdf")

data_scanpy_1.obs['Normal'] = pd.Categorical(
    data_scanpy_1.obs['Normal'],
    categories=['Donor', 'COPD', 'PF', 'COVID19'],
    ordered=True
)
# Stacked barplot for the levels of "Normal"
with rc_context({'figure.figsize': (10, 7)}):
    viz.stacked_barplot(data_scanpy_1, feature_name="Normal")
plt.savefig("stacked_barplot.Disease.pdf")


# Grouped boxplots. No facets, relative abundance, no dots.
with rc_context({'figure.figsize': (10, 7)}):
    viz.boxplots(
        data_scanpy_1,
        feature_name="Normal",
        plot_facets=False,
        y_scale="relative",
        add_dots=False,
    )
plt.savefig("boxplots.Disease.pdf")

# Grouped boxplots. Facets, log scale, added dots and custom color palette.
with rc_context({'figure.figsize': (10, 7)}):
    viz.boxplots(
        data_scanpy_1,
        feature_name="Normal",
        plot_facets=True,
        y_scale="log",
        add_dots=True,
        cmap="Reds",
    )
plt.savefig("boxplots.Disease.facets.pdf")



#Finding a reference cell type
with rc_context({'figure.figsize': (10, 7)}):
    viz.rel_abundance_dispersion_plot(
        data=data_scanpy_1,
        abundant_threshold=0.7
    )
plt.savefig("rel_abundance0.7_dispersion_plot.pdf")



model_salm = mod.CompositionalAnalysis(data_scanpy_1, formula='C(Normal, Treatment("Donor"))', reference_cell_type="Adventitial Fibroblast")#formula='C(Normal, Treatment("Donor"))',  # 关键：指定 Donor 为参考

# Run MCMC
sim_results = model_salm.sample_hmc()

pd.set_option('display.max_columns', None)
print(sim_results.summary())
effects_df = sim_results.effect_df
effects_df.to_csv("effect_results.fdr0.05.csv")
intercepts_df = sim_results.intercept_df
intercepts_df.to_csv("intercept_results.fdr0.05.csv")
print(sim_results.credible_effects())
cred = sim_results.credible_effects()
cred.to_csv("credible_effects.fdr0.05.csv")

sim_results.set_fdr(est_fdr=0.1)
print(sim_results.summary())
print(sim_results.credible_effects())
effects_df = sim_results.effect_df
effects_df.to_csv("effect_results.fdr0.1.csv")
intercepts_df = sim_results.intercept_df
intercepts_df.to_csv("intercept_results.fdr0.1.csv")
cred = sim_results.credible_effects()
cred.to_csv("credible_effects.fdr0.1.csv")

sim_results.set_fdr(est_fdr=0.2)
print(sim_results.summary())
print(sim_results.credible_effects())
effects_df = sim_results.effect_df
effects_df.to_csv("effect_results.fdr0.2.csv")
intercepts_df = sim_results.intercept_df
intercepts_df.to_csv("intercept_results.fdr0.2.csv")
cred = sim_results.credible_effects()
cred.to_csv("credible_effects.fdr0.2.csv")
# saving
path = "all_orig.ident.result"
sim_results.save(path)
'''
# loading
with open(path, "rb") as f:
    sim_results_2 = pkl.load(f)

sim_results_2.summary()
'''
