
# Understanding adata_combined — The AnnData Object


import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from pathlib import Path
import matplotlib.pyplot as plt
import scanpy as sc
import scvi
from scipy import sparse
# Set random seed for reproducibility
scvi.settings.seed = 42


# Load data: (preprocessing,HVG selection, scVI, latent representation, UMAP, Leiden clusters)
integrated_path = Path("data/processed/combined_integrated.h5ad")
adata_combined = sc.read_h5ad("data/processed/combined_integrated.h5ad")
print(adata_combined.shape)
print(adata_combined.X)









# Gene expression matrix — .X:
print(adata_combined.shape)
print(adata_combined.X)


# To view a small portion as a regular DataFrame:
# First 5 nuclei and first 10 genes
df = pd.DataFrame(adata_combined.X[:5, :10].toarray(), index=adata_combined.obs_names[:5], columns=adata_combined.var_names[:10])
print(df)

# Nucleus metadata — .obs
print(adata_combined.obs.head())
print(adata_combined.obs.columns.tolist())
print(adata_combined.obs[["biosample_id", "donor_id", "group", "leiden"]].head(10))


# Gene metadata — .var
print(adata_combined.var.head())
print(adata_combined.var.columns.tolist())
print(adata_combined.var[adata_combined.var["highly_variable"]].head(10))

# scVI and UMAP representations — .obsm
print(adata_combined.obsm.keys())
# scVI latent representation
print(adata_combined.obsm["X_scVI"].shape)
# UMAP coordinates
print(adata_combined.obsm["X_umap"].shape)
print(adata_combined.obsm["X_scVI"][:5])

#  Quick reference
# What you want to inspect	Python command
# Dataset dimensions	adata_combined.shape
# Gene expression	adata_combined.X
# Nucleus IDs	adata_combined.obs_names
# Gene names	adata_combined.var_names
# Nucleus metadata	adata_combined.obs
# Gene metadata	adata_combined.var
# Nuclei per sample	adata_combined.obs["biosample_id"].value_counts()
# Nuclei per Leiden cluster	adata_combined.obs["leiden"].value_counts()
# scVI latent representation	adata_combined.obsm["X_scVI"]
# UMAP coordinates	adata_combined.obsm["X_umap"]




# The most important concept
# All these components are synchronized by nucleus.
# For example, row 100 in each of these structures refers to the same nucleus:
adata_combined.X[100, :]             # Gene expression
adata_combined.obs.iloc[100]         # Nucleus metadata
adata_combined.obsm["X_scVI"][100]   # 10 scVI features
adata_combined.obsm["X_umap"][100]   # 2 UMAP coordinates