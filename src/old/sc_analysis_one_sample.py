import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

import scanpy as sc
import scvi



# ==================================================================================================
# ==================================================================================================
# 1.SINGLE-CELL RNA-SEQ ANALYSIS — WORKFLOW SUMMARY (ONE SAMPLE)
# 2.METADATA INSPECTION AND VALIDATION — SINGLE SAMPLE (C51ctr)
# ==================================================================================================
# ==================================================================================================







# ==================================================================================================
# SINGLE-CELL RNA-SEQ ANALYSIS — WORKFLOW SUMMARY (ONE SAMPLE)
# ==================================================================================================

# 1. LOAD RAW DATA
# Load the raw gene-count matrix into AnnData.
# Rows (obs) = nuclei/cells; columns (var) = genes; X = raw UMI counts.


# 2. DOUBLET DETECTION (scVI + SOLO)
# A doublet occurs when two cells/nuclei are captured together and appear as one.
#
# - Temporarily filter low-detection genes and select 2,000 highly variable genes.
# - Train scVI to learn a low-dimensional representation of the expression data.
# - Train SOLO to distinguish likely singlets from doublets.
# - Inspect predicted doublet probabilities.
# - Remove high-confidence doublets (probability > 0.70).
# - Reload the complete raw count matrix before removing doublets,
#   preserving all original genes for downstream analysis.


# 3. QUALITY CONTROL (QC)
# Identify low-quality nuclei and genes using expression-based QC metrics.
#
# - Calculate the number of detected genes and total UMI counts per nucleus.
# - Calculate mitochondrial and ribosomal RNA percentages.
# - Remove genes detected in fewer than 3 nuclei.
# - Remove nuclei above the 98th percentile of detected genes.
# - Remove nuclei with mitochondrial RNA > 20%.
# - Remove nuclei with ribosomal RNA >= 2%.
#
# Purpose: Remove potential technical artifacts and low-quality nuclei.


# 4. NORMALIZATION
# Make expression values more comparable across nuclei.
#
# - Normalize each nucleus to 10,000 total counts.
# - Apply log1p transformation: log(1 + normalized expression).
# - Store normalized, log-transformed expression in adata.raw.
#
# Purpose: Reduce differences caused by sequencing depth and compress
# the range of expression values.
#
# NOTE: After normalization, values are no longer raw UMI counts.


# 5. FEATURE SELECTION AND SCALING
# Prepare the expression matrix for dimensionality reduction.
#
# - Select the 2,000 most highly variable genes (HVGs).
# - Regress out total counts, mitochondrial %, and ribosomal %.
# - Standardize gene expression (mean ~0, variance ~1).
#
# Purpose: Focus on informative genes and reduce selected technical variation.


# 6. DIMENSIONALITY REDUCTION (PCA + UMAP)
# Represent thousands of gene-expression measurements in fewer dimensions.
#
# - PCA: Reduce expression data to principal components.
# - Neighbors: Build a similarity graph using the first 30 PCs.
# - UMAP: Project nuclei into 2D for visualization.
#
# Purpose: Explore




# 1. LOAD FULL RAW DATA
# ==================================================================================================
# Open the first file (first sample)
adata = sc.read_csv(r'data\raw_counts\GSM5226574_C51ctr_raw_counts.csv.gz').T

# The adata object has
adata.obs  # the cells (here these nuclei) metadata
adata.var  # The genes
adata.X    # The counts







# 2. DOUBLET DETECTION
# ==================================================================================================
# These filtration are only for the doublets detection step
sc.pp.filter_genes(adata, min_cells=10)  # Keep only genes that are detected in at least 10 cells/nuclei.
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, flavor='seurat_v3')  # Keep the most variable genes

# Train scVI to learn a lower-dimensional representation of the nuclei
scvi.model.SCVI.setup_anndata(adata)
vae = scvi.model.SCVI(adata)  # model for learning representation of cells
vae.train()

# Train SOLO doublet classifier using the trained scVI model
solo = scvi.external.SOLO.from_scvi_model(vae)  # SOLO is the actual doublet classifier.
solo.train()

# Get doublet probabilities
df = solo.predict()
df['prediction'] = solo.predict(soft=False)

# INSPECT DOUBLET PROBABILITY
doublet_candidates = df[df["prediction"] == "doublet"]

plt.figure()
sns.histplot( doublet_candidates["doublet"], bins=30)
plt.axvline( 0.70, linestyle="--")
plt.xlabel("SOLO doublet probability")
plt.ylabel("Number of nuclei")
plt.title("SOLO predicted doublets")
plt.show()

# DEFINE HIGH-CONFIDENCE DOUBLETS
doublet_threshold = 0.70
doublets = df[df["doublet"] > doublet_threshold]
print("Total nuclei:", adata.n_obs)
print("SOLO predicted doublets:",(df["prediction"] == "doublet").sum())
print("High-confidence doublets:",len(doublets))
print("Fraction removed:", len(doublets) / adata.n_obs)

# Reload the complete sample
adata = sc.read_csv(r'data\raw_counts\GSM5226574_C51ctr_raw_counts.csv.gz').T

# REMOVE HIGH-CONFIDENCE DOUBLETS FROM THE FULL DATA
adata = adata[~adata.obs_names.isin(doublets.index)].copy()
print("Remaining nuclei:", adata.n_obs)



# 3. PREPROCESSING / QUALITY CONTROL
# ==================================================================================================

# Mark mitochondrial genes (human mitochondrial genes start with "MT-")
adata.var['mt'] = adata.var.index.str.startswith('MT-')

# Mark ribosomal genes
ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
ribo_genes = pd.read_table(ribo_url, skiprows=2, header=None)
adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)

# Calculate QC metrics for each nucleus:
# n_genes_by_counts, total_counts, pct_counts_mt, pct_counts_ribo, etc.
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)

# Inspect gene-level and nucleus-level QC metrics
adata.var.sort_values('n_cells_by_counts')
adata.var
adata.obs

# Remove genes detected in fewer than 3 nuclei
sc.pp.filter_genes(adata, min_cells=3)

# Inspect number of detected genes and total UMI counts per nucleus
adata.obs[["n_genes_by_counts", "total_counts"]].sort_values("n_genes_by_counts")

# Plot the main QC metrics before filtering nuclei
sc.pl.violin(adata,['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'],jitter=0.4,multi_panel=True)

# Remove nuclei in the highest 2% for number of detected genes
# (possible outliers / remaining doublets)
upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, 0.98)
adata = adata[adata.obs.n_genes_by_counts < upper_lim].copy()

# Remove nuclei with high mitochondrial RNA percentage
# (potentially damaged / low-quality nuclei)
adata = adata[adata.obs.pct_counts_mt < 20].copy()

# Remove nuclei with high ribosomal RNA percentage
adata = adata[adata.obs.pct_counts_ribo < 2].copy()




# 4. NORMALIZATION
# ==================================================================================================
# It is generally to normalize after integration

# Check the total number of UMI counts in each nucleus before normalization
# Different nuclei have different sequencing/capture depths.
adata.X.sum(axis=1)

# Normalize each nucleus to the same total count (10,000).
# This makes expression values more comparable between nuclei with different total UMI counts.
sc.pp.normalize_total(adata, target_sum=1e4)

# Check: each nucleus should now have approximately 10,000 total normalized counts.
# NOTE: these are now normalized values, NOT actual/raw UMI counts.
adata.X.sum(axis=1)

# Log-transform the normalized expression values: log(1 + x).
# This compresses large values and reduces the effect of highly expressed genes.
sc.pp.log1p(adata)

# Store the current normalized + log-transformed data as .raw.
# This preserves these expression values for later use (e.g. marker-gene analysis),
# even if adata.X is later scaled or subset to highly variable genes.
adata.raw = adata



# 5. CLUSTERING
# ==================================================================================================
# Select the 2,000 most highly variable genes (HVGs).
# These genes contain most of the useful variation for distinguishing cell populations.
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.var
sc.pl.highly_variable_genes(adata)

# Keep only the 2,000 HVGs for dimensionality reduction and clustering.
adata = adata[:, adata.var.highly_variable].copy()

# Regress out technical/QC effects:
# remove variation associated with sequencing depth, mitochondrial %, and ribosomal %.
sc.pp.regress_out(adata,['total_counts', 'pct_counts_mt', 'pct_counts_ribo'])

# Standardize each gene (approximately mean=0, variance=1).
# Clip extreme scaled values at 10.
sc.pp.scale(adata, max_value=10)

# PCA: reduce the 2,000 genes to a smaller number of principal components.
sc.tl.pca(adata, svd_solver='arpack')

# Plot how much variance is explained by each PC.
# Used to help decide how many PCs to retain.
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)

# Build the cell-neighbor graph using the first 30 PCs.
# Cells with similar expression profiles become neighbors.
sc.pp.neighbors(adata, n_pcs=30)

# Calculate UMAP coordinates from the neighbor graph.
sc.tl.umap(adata)

# Visualize the cells in 2 dimensions.
sc.pl.umap(adata)

# Leiden clustering: divide the neighbor graph into communities/clusters.
# Resolution controls how fine the clustering is:
# higher resolution -> more/smaller clusters.
sc.tl.leiden(adata, resolution=0.5)

# Inspect cluster labels stored in cell metadata.
adata.obs

# Plot the UMAP colored according to Leiden cluster.
sc.pl.umap(adata, color=['leiden'])












# ======================================================================================================================
# METADATA INSPECTION AND VALIDATION — SINGLE SAMPLE (C51ctr)
# ======================================================================================================================
# Goal:
# 1. Load the metadata table containing cell-level information for all samples.
# 2. Load one raw count matrix to examine the structure of the cell identifiers.
# 3. Verify that cell IDs in the raw count matrix match those in the metadata.
# 4. Identify cells that are missing from either dataset.
# 5. Remove cells without matching metadata.
# 6. Attach the corresponding metadata to the AnnData object.
#
# This is a validation step performed on one sample before automating the process for all samples.
# ======================================================================================================================


# ----------------------------------------------------------------------------------------------------------------------
# 1. LOAD AND INSPECT THE METADATA
# ----------------------------------------------------------------------------------------------------------------------

# Load the complete metadata table (all samples)
metadata = pd.read_csv(r"data\GSE171524_lung_metaData.txt", sep="\t")

# Remove the first row (excluded during initial metadata inspection)
metadata = metadata.iloc[1:].copy()

# Reset the DataFrame index after removing the first row
metadata.reset_index(drop=True, inplace=True)

# Inspect the dimensions and first rows
print(metadata.shape)
print(metadata.head())

# Count the number of cells associated with each biological sample
print(metadata["biosample_id"].value_counts())

# Examine the main metadata columns:
# NAME         → Unique cell identifier
# biosample_id → Sample identifier (e.g., C51ctr)
# donor_id     → Donor identifier
# group        → Experimental group (Control or COVID-19)
print(metadata[["NAME", "biosample_id", "donor_id", "group"]].head(20))


# ----------------------------------------------------------------------------------------------------------------------
# 2. LOAD ONE RAW COUNT MATRIX
# ----------------------------------------------------------------------------------------------------------------------

# Select sample C51ctr to validate the metadata-matching procedure
file_path = r'data\raw_counts\GSM5226574_C51ctr_raw_counts.csv.gz'

# Load the raw counts and transpose:
# Rows = cells, columns = genes
adata_test = sc.read_csv(file_path).T

# Inspect the sample name and dimensions (cells × genes)
print("Sample:", Path(file_path).name)
print("Shape:", adata_test.shape)

# Compare the format of cell identifiers in the two datasets
print("\nFirst 5 cell IDs from AnnData:")
print(adata_test.obs_names[:5].tolist())

print("\nFirst 5 cell IDs from metadata:")
print(metadata["NAME"].head().tolist())


# ----------------------------------------------------------------------------------------------------------------------
# 3. VALIDATE CELL IDENTIFIER MATCHING
# ----------------------------------------------------------------------------------------------------------------------

# Select only metadata rows belonging to sample C51ctr
sample_meta = metadata[metadata["biosample_id"] == "C51ctr"].copy()

# Convert cell identifiers into sets to compare unique cell IDs
raw_cells = set(adata_test.obs_names)
meta_cells = set(sample_meta["NAME"])

# Compare the number of unique cells in each dataset
print("Raw cells:", len(raw_cells))
print("Metadata cells:", len(meta_cells))

# Identify cells present in the raw count matrix but missing from the metadata
print("Only in raw data:", raw_cells - meta_cells)

# Identify cells present in the metadata but missing from the raw count matrix
print("Only in metadata:", meta_cells - raw_cells)

# Observed result for C51ctr:
# Raw count matrix: 6099 cells
# Metadata:         6098 cells
# One cell (TAGGTACCATGGCCAC-1_1) is present only in the raw data.
# All metadata cell IDs are present in the raw count matrix.


# ----------------------------------------------------------------------------------------------------------------------
# 4. REMOVE CELLS WITHOUT MATCHING METADATA
# ----------------------------------------------------------------------------------------------------------------------

# Keep only cells whose identifiers appear in the sample metadata
# .isin() creates a Boolean mask indicating which cells have matching metadata
adata_test = adata_test[adata_test.obs_names.isin(sample_meta["NAME"])].copy()

# Verify the dimensions after filtering
print(adata_test.shape)

# Expected: (6098, 34546)
# One cell was removed; the number of genes remains unchanged.


# ----------------------------------------------------------------------------------------------------------------------
# 5. ATTACH METADATA TO THE ANNDATA OBJECT
# ----------------------------------------------------------------------------------------------------------------------

# Set the cell identifier (NAME) as the metadata DataFrame index
sample_meta = sample_meta.set_index("NAME")

# Match metadata rows to the exact order of cells in AnnData
# This ensures that each cell receives its correct metadata information
adata_test.obs = sample_meta.loc[adata_test.obs_names].copy()

# Verify that the metadata was correctly attached
print(adata_test.obs[["biosample_id", "donor_id", "group"]].head())

# Final result:
# adata_test.X         → Raw gene expression counts (cells × genes)
# adata_test.obs       → Cell-level metadata (sample, donor, experimental group, etc.)
# adata_test.var       → Gene information
# adata_test.obs_names → Cell identifiers
#
# The AnnData object is now ready for sample-level preprocessing (doublet detection and QC).
# ======================================================================================================================