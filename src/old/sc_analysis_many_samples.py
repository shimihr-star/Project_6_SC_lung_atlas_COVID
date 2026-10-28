"""GSE171524 lung snRNA-seq: preprocessing, integration, and exploratory annotation.

Run sections 1–4 once for all samples; they include computationally expensive
per-sample scVI/SOLO training. The combined raw-count AnnData is saved to disk.
In a later session, load combined.h5ad and start at section 5.
"""

# ========================================================================
# 1. SETUP
# ========================================================================
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from pathlib import Path
import matplotlib.pyplot as plt
import scanpy as sc
import scvi
import seaborn as sns
from scipy import sparse
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
# Set random seed for reproducibility
scvi.settings.seed = 42

# Define the ribosomal list once
ribo_url = ("http://software.broadinstitute.org/gsea/msigdb/" "download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt")
ribo_genes = pd.read_table(ribo_url, skiprows=2, header=None)


# ______________________________________________________________________________________________________________________
# COVID-19 snRNA-seq project — Summary
# Our main goal: Understand how COVID-19 affects different cell populations in human lungs by comparing COVID-19 patients with healthy controls.
#
# We're using dataset GSE171524, containing 26 biological donors (19 COVID-19 and 7 controls) across 27 sample files.

# Load and clean the data (QC)
# We loaded raw RNA counts for individual nuclei, identified likely doublets using SOLO, and filtered low-quality nuclei.
#
# Attach metadata
# We linked nuclei to their sample, biological donor, and experimental group (COVID-19 or Control).
#
# Combine all samples
# We merged the processed samples into one dataset while preserving raw UMI counts.
#
# Integrate samples using scVI
# We used a deep-learning model to learn a representation of gene expression while accounting for differences between samples (batch effects).
#
# Cluster nuclei using Leiden
# We constructed a nearest-neighbor graph, identified 18 clusters, and visualized them using UMAP.
#
# Identify cell types
# We examined marker genes to assign biological identities to clusters, such as macrophages, AT1/AT2 epithelial cells, endothelial cells, and fibroblasts. These annotations are provisional.
#
# Compare cell-type proportions
# We calculated the proportion of each cell type per donor and set up statistical comparisons between COVID-19 and controls, including a sensitivity analysis for donor-dominated clusters.
# ______________________________________________________________________________________________________________________


# ========================================================================
# 2. PREPROCESSING FUNCTIONS
# ========================================================================
# attach_metadata: match nucleus barcodes to metadata for the current sample.
# Preserve QC columns and discard unmatched nuclei.

def attach_metadata(adata, metadata, sample_id):
    # Select metadata for the current sample
    sample_meta = metadata[metadata["biosample_id"] == sample_id].copy()

    # Use cell IDs as the index
    sample_meta = sample_meta.set_index("NAME")

    # Identify cells with matching metadata
    valid_cells = adata.obs_names.isin(sample_meta.index)
    print("Cells before metadata matching:", adata.n_obs)
    print("Cells without metadata:", (~valid_cells).sum())

    # Keep only cells with metadata
    adata = adata[valid_cells].copy()

    # Attach metadata while PRESERVING existing QC columns
    adata.obs = adata.obs.join(sample_meta, how="left")
    print("Cells after metadata matching:", adata.n_obs)
    return adata


# preprocess_sample: use scVI + SOLO to identify high-confidence doublets,
# then filter low-quality nuclei. Return unnormalized UMI counts.
# HVG filtering for SOLO occurs on a temporary copy, preserving all genes.

def preprocess_sample(file_path, doublet_threshold=0.70):
    # 2.1 Load raw counts

    adata = sc.read_csv(file_path).T

    # Convert the raw count matrix to sparse CSR format
    adata.X = sparse.csr_matrix(adata.X)
    # Verify matrix format
    print("Matrix format:", type(adata.X))
    print(f"\nProcessing: {Path(file_path).name}")
    print("Initial shape:", adata.shape)

    # 2.2 Detect doublets with scVI + SOLO

    # Keep a temporary object for SOLO because HVG filtering removes genes
    adata_solo = adata.copy()

    # These filters are only for doublet detection
    sc.pp.filter_genes(adata_solo, min_cells=10)
    sc.pp.highly_variable_genes(adata_solo, n_top_genes=2000, subset=True, flavor="seurat_v3")

    # Train scVI
    scvi.model.SCVI.setup_anndata(adata_solo)

    vae = scvi.model.SCVI(adata_solo)
    vae.train()

    # Train SOLO doublet classifier
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()

    # Get doublet probabilities
    df = solo.predict()
    df["prediction"] = solo.predict(soft=False)

    # Verify that SOLO cell IDs match the original AnnData cell IDs
    n_matched = df.index.isin(adata.obs_names).sum()
    print("SOLO cells:", len(df))
    print("Matching cell IDs:", n_matched)
    assert n_matched == len(df), "SOLO cell IDs do not match AnnData!"

    # Define high-confidence doublets
    doublets = df[df["doublet"] > doublet_threshold]

    print("SOLO predicted doublets:", (df["prediction"] == "doublet").sum())
    print("High-confidence doublets:", len(doublets))
    # Calculate the percentage of nuclei removed as doublets
    doublet_pct = 100 * len(doublets) / adata.n_obs
    print(f"Doublets removed: {doublet_pct:.2f}%")

    # Remove high-confidence doublets from the FULL raw-count object
    adata = adata[~adata.obs_names.isin(doublets.index)].copy()
    print("After doublet removal:", adata.shape)

    # 2.3 Filter low-quality nuclei
    # Number of nuclei before QC filtering
    n_before_qc = adata.n_obs

    # Mark mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # Check mitochondrial gene annotation
    print("Mitochondrial genes detected:", adata.var["mt"].sum())

    # Mark ribosomal genes
    adata.var["ribo"] = adata.var_names.isin(ribo_genes[0].values)
    # Check ribosomal gene annotation
    print("Ribosomal genes detected:", adata.var["ribo"].sum())

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)

    # Remove nuclei with fewer than 200 detected genes
    sc.pp.filter_cells(adata, min_genes=200)

    # Remove nuclei in the highest 2% for number of detected genes
    upper_lim = np.quantile(adata.obs["n_genes_by_counts"], 0.98)
    adata = adata[adata.obs["n_genes_by_counts"] < upper_lim].copy()

    # Remove high-mitochondrial nuclei
    adata = adata[adata.obs["pct_counts_mt"] < 20].copy()

    # Remove high-ribosomal nuclei
    adata = adata[adata.obs["pct_counts_ribo"] < 2].copy()
    print("Final shape after QC:", adata.shape)

    # Summarize QC filtering
    n_removed_qc = n_before_qc - adata.n_obs
    qc_pct = 100 * n_removed_qc / n_before_qc
    print(f"Nuclei removed during QC: {n_removed_qc} ({qc_pct:.2f}%)")

    # Summary of detected genes per nucleus after QC
    print("Median genes per nucleus:", adata.obs["n_genes_by_counts"].median())
    print("Minimum genes per nucleus:", adata.obs["n_genes_by_counts"].min())
    print("Maximum genes per nucleus:", adata.obs["n_genes_by_counts"].max())
    return adata


# ========================================================================
# 3. PROCESS AND VALIDATE SAMPLES
# ========================================================================
# 3.1 Load metadata and correct donor mapping
metadata = pd.read_csv(r"data\GSE171524_lung_metaData.txt", sep="\t", low_memory=False)
# Remove the extra annotation row
metadata = metadata.iloc[1:].copy()
metadata.reset_index(drop=True, inplace=True)

# L04 covaddon is an additional specimen from donor L04cov.
# Correct only the biological donor ID; preserve its distinct biosample_id.
metadata.loc[metadata["biosample_id"] == "L04covaddon", "donor_id"] = "L04cov"
print("Metadata shape:", metadata.shape)

# 3.2 Locate sample files
# Full dataset (27 sample files):
data_dir = Path(r"data\raw_counts")

sample_files = sorted(data_dir.glob("*_raw_counts.csv.gz"))
print("Number of sample files:", len(sample_files))
assert len(sample_files) > 0, "ERROR: No sample files found!"

# Verify that sample IDs extracted from filenames are unique
sample_ids = [file_path.name.split("_")[1] for file_path in sample_files]
assert len(sample_ids) == len(set(sample_ids)), ("ERROR: Duplicate sample IDs found in filenames!")

# 3.3 Process each sample (SOLO, QC, metadata)
processed_samples = {}
for file_path in sample_files:
    # Example: GSM5226574_C51ctr_raw_counts.csv.gz -> C51ctr
    sample_id = file_path.name.split("_")[1]
    print("\n" + "=" * 80)
    print("Processing sample:", sample_id)
    print("=" * 80)

    # STEP 1: Doublet detection and QC
    adata = preprocess_sample(file_path)

    # STEP 2: Attach metadata while preserving QC metrics
    adata = attach_metadata(adata, metadata, sample_id)

    # STEP 3: Store the processed sample
    processed_samples[sample_id] = adata
    print(f"Completed {sample_id}: " f"{adata.n_obs} nuclei x {adata.n_vars} genes")

# 3.4 Validate processed samples
print("\n" + "=" * 80)
print("PREPROCESSING VALIDATION")
print("=" * 80)
print("Samples found:", len(sample_files))
print("Samples processed:", len(processed_samples))
assert len(processed_samples) == len(sample_files), ("ERROR: Not all samples were processed!")

# Required metadata and QC columns
required_columns = [
    "biosample_id",
    "donor_id",
    "group",
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
    "pct_counts_ribo"]

for sample_id, adata in processed_samples.items():
    print(f"\nSample: {sample_id}")
    print("Shape:", adata.shape)

    # 4.1 Verify metadata and QC columns
    missing_columns = [
        col for col in required_columns
        if col not in adata.obs.columns]

    assert not missing_columns, (f"{sample_id}: Missing columns: {missing_columns}")

    # 4.2 Check missing metadata values
    for col in ["biosample_id", "donor_id", "group"]:
        n_missing = adata.obs[col].isna().sum()
        print(f"Missing {col}: {n_missing}")
        assert n_missing == 0, (f"{sample_id}: Missing values in {col}!")

    # Verify metadata belongs to the correct sample
    assert (adata.obs["biosample_id"] == sample_id).all(), (f"{sample_id}: Incorrect sample metadata!")

    # 4.3 Check duplicated cell and gene IDs
    duplicate_cells = adata.obs_names.duplicated().sum()
    duplicate_genes = adata.var_names.duplicated().sum()
    print("Duplicate cell IDs:", duplicate_cells)
    print("Duplicate gene IDs:", duplicate_genes)
    assert duplicate_cells == 0, (f"{sample_id}: Duplicate cell IDs!")
    assert duplicate_genes == 0, (f"{sample_id}: Duplicate gene IDs!")

    # 4.4 Verify minimum detected genes per nucleus
    n_low = (adata.obs["n_genes_by_counts"] < 200).sum()
    print("Nuclei with fewer than 200 genes:", n_low)
    assert n_low == 0, (f"{sample_id}: Low-complexity nuclei remain!")

    # 4.5 Verify raw-count matrix
    # For sparse matrices, check stored nonzero values.
    # For dense matrices, check all values.
    if hasattr(adata.X, "toarray"):
        X_values = adata.X.data
    else:
        X_values = np.asarray(adata.X).ravel()

    assert np.isfinite(X_values).all(), (f"{sample_id}: Non-finite count values!")
    assert (X_values >= 0).all(), (f"{sample_id}: Negative count values!")
    assert np.allclose(X_values, np.round(X_values)), (f"{sample_id}: Non-integer count values!")
    print("Raw count matrix: Valid")

    # 4.6 Display QC summary (including the upper tails, not only medians)

# 3.5 Check gene identifiers across samples
# Confirm that gene identifiers are consistent across all samples.

gene_sets = [set(adata.var_names)for adata in processed_samples.values()]

# Genes present in every processed sample
shared_genes = set.intersection(*gene_sets)

# Genes present in at least one processed sample
all_genes = set.union(*gene_sets)
print("\n" + "=" * 80)
print("GENE OVERLAP SUMMARY")
print("=" * 80)
print("Genes shared by all samples:", len(shared_genes))
print("Total unique genes:", len(all_genes))

# 3.6 Summarize QC across samples
# Create a compact table for comparing QC metrics between samples.
summary = pd.DataFrame([
    {
        "sample_id": sample_id,
        "group": ", ".join(map(str, adata.obs["group"].unique())),
        "n_nuclei": adata.n_obs,
        "n_genes": adata.n_vars,
        "median_genes": adata.obs["n_genes_by_counts"].median(),
        "median_UMIs": adata.obs["total_counts"].median(),
        "median_mt_pct": adata.obs["pct_counts_mt"].median(),
        "median_ribo_pct": adata.obs["pct_counts_ribo"].median(),
        "mt_95th_pct": adata.obs["pct_counts_mt"].quantile(0.95),
        "ribo_95th_pct": adata.obs["pct_counts_ribo"].quantile(0.95)
    }
    for sample_id, adata in processed_samples.items()])


print("\n" + "=" * 80)
print("FINAL SAMPLE SUMMARY")
print("=" * 80)
print(summary.to_string(index=False))
print("\nPreprocessing and validation completed successfully.")

# ========================================================================
# 4. COMBINE, QC, AND SAVE
# ========================================================================
# Stack nuclei without normalization; retain raw UMI counts and metadata.
# Batch correction happens later in scVI, not during concatenation.

# 4.1 Confirm compatible gene sets
reference_genes = processed_samples[next(iter(processed_samples))].var_names
for sample_id, adata in processed_samples.items():
    assert set(adata.var_names) == set(reference_genes), (f"{sample_id}: Gene sets do not match!")
print("All samples contain the same gene set.")

# 4.2 Concatenate and validate sample/donor identities
adata_combined = sc.concat(processed_samples, join="inner", label="sample_id", index_unique="-", merge="same")
print("Combined shape:", adata_combined.shape)

# Check donor-level replication: multiple sample files may belong to one donor.
print("Unique sample files:", adata_combined.obs["biosample_id"].nunique())
print("Unique biological donors:", adata_combined.obs["donor_id"].nunique())
assert adata_combined.obs["biosample_id"].nunique() == 27
assert adata_combined.obs["donor_id"].nunique() == 26
print("Donors by experimental group:")
print(adata_combined.obs.groupby("group", observed=True)["donor_id"].nunique())
assert (adata_combined.obs.groupby("donor_id", observed=True)[
            "group"].nunique() == 1).all(), "Donor assigned to multiple groups!"
print("Sample-to-donor mapping:")
print(adata_combined.obs[["biosample_id", "donor_id", "group"]].drop_duplicates().sort_values("biosample_id").to_string(
    index=False))

# Convert numeric metadata columns to float64 while preserving NaN values
numeric_cols = ["intubation_days", "interval_death_symptoms_onset_days", "pmi_h"]
for col in numeric_cols:
    adata_combined.obs[col] = pd.to_numeric(adata_combined.obs[col], errors="coerce")

# Verify column types
print(adata_combined.obs[numeric_cols].dtypes)

# 4.3 Verify concatenation

# Check number of nuclei
expected_nuclei = sum(adata.n_obs for adata in processed_samples.values())
assert adata_combined.n_obs == expected_nuclei, ("ERROR: Number of nuclei changed during concatenation!")

# Check unique nucleus IDs
assert adata_combined.obs_names.is_unique, ("ERROR: Duplicate nucleus IDs after concatenation!")

# Check sample identities
print("\nNuclei per sample:")
print(adata_combined.obs["sample_id"].value_counts())

# Check experimental groups
print("\nNuclei per group:")
print(adata_combined.obs["group"].value_counts())

# Check count matrix format
print("\nMatrix format:", type(adata_combined.X))
print("\nConcatenation completed successfully.")

# 4.4 Global gene filtering
# Remove genes detected in fewer than 3 nuclei
# across the entire combined dataset.
print("Before gene filtering:", adata_combined.shape)
sc.pp.filter_genes(adata_combined, min_cells=3)
print("After gene filtering:", adata_combined.shape)

# 4.5 Select highly variable genes
# Before training scVI for integration, we want to identify genes that show meaningful variation across nuclei.
# Since scVI expects raw counts, we'll use seurat_v3, which is designed to work with count data.
sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000, flavor="seurat_v3", subset=False)
print("Highly variable genes:", adata_combined.var["highly_variable"].sum())

# 4.6 Save combined raw counts and metadata
# Save the combined AnnData object, including:
# - Raw UMI counts (adata.X)
# - Cell metadata and QC metrics (adata.obs)
# - Gene annotations and HVG selection (adata.var)

output_path = Path("data/processed/combined.h5ad")
# Create output directory if it does not exist
output_path.parent.mkdir(parents=True, exist_ok=True)
# Save AnnData
# adata_combined.write_h5ad(output_path, compression="gzip")

# reload: (preprocessing,HVG selection)
# In a fresh session, load the saved object instead of rerunning SOLO:
adata_combined = sc.read_h5ad(output_path)
print("Dataset shape:", adata_combined.shape)
print("Unique samples:", adata_combined.obs["biosample_id"].nunique())
print("Unique donors:", adata_combined.obs["donor_id"].nunique())
print("HVGs:", adata_combined.var["highly_variable"].sum())

# ============================================================
# QC REVIEW — POST-FILTERING DATASET

#
# PURPOSE:
# Examine the quality-control metrics of nuclei retained
# after preprocessing, without rerunning scVI or SOLO.
#
# NOTE:
# These are surviving nuclei only. We cannot determine
# how many nuclei each original filter removed.

qc_columns = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"]

# Check which QC metrics are available.
print("Available QC columns:")
print([col for col in qc_columns if col in adata_combined.obs.columns])

# Display descriptive statistics.
print(adata_combined.obs[qc_columns].describe())

# Visualize the QC distributions.
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
for ax, col in zip(axes.flat, qc_columns):
    sns.histplot(adata_combined.obs[col], bins=50, ax=ax)
    ax.set_title(col)

plt.tight_layout()
plt.show()


qc_by_sample = (
    adata_combined.obs
    .groupby("biosample_id", observed=True)
    .agg(
        n_nuclei=("n_genes_by_counts", "size"),
        median_genes=("n_genes_by_counts", "median"),
        median_umis=("total_counts", "median"),
        median_mt=("pct_counts_mt", "median"),
        median_ribo=("pct_counts_ribo", "median"),
        max_mt=("pct_counts_mt", "max"),
    )
)

print(qc_by_sample.round(2).to_string())

# =============================================================================
# 5. scVI INTEGRATION AND CLUSTERING
# =============================================================================
# 5.1 Prepare HVG dataset
# Create a separate AnnData object containing only the selected HVGs.
# Keep the full dataset for downstream analysis.
adata_scvi = adata_combined[:, adata_combined.var["highly_variable"]].copy()
print("Full dataset:", adata_combined.shape)
print("scVI dataset:", adata_scvi.shape)

# 5.2 Configure scVI
# scVI expects raw UMI counts.
# Use sample ID as the batch key, not donor_id.
scvi.model.SCVI.setup_anndata(adata_scvi,batch_key="biosample_id")

# 5.3 Train scVI
# Learn a 10-dimensional latent representation of each nucleus.
model = scvi.model.SCVI(adata_scvi, n_latent=10)
print(model)
model.train()

# 5.4 Extract latent representation
# Store the learned features in the full AnnData object.
latent = model.get_latent_representation()
adata_combined.obsm["X_scVI"] = latent
print("Latent representation:", latent.shape)
print("Full dataset:", adata_combined.shape)

# 5.5 Construct nearest-neighbor graph
# Use the scVI latent representation instead of raw gene counts.
sc.pp.neighbors(adata_combined, use_rep="X_scVI", n_neighbors=15)
print("Neighbor graph created successfully.")

# 5.6 Compute UMAP
# Project the scVI latent representation into two dimensions.
sc.tl.umap(adata_combined)
print("UMAP coordinates:", adata_combined.obsm["X_umap"].shape)

# 5.7 Leiden clustering
# Identify communities using the nearest-neighbor graph.
# Resolution controls clustering granularity.
sc.tl.leiden(adata_combined, resolution=0.5, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
print("Number of clusters:", adata_combined.obs["leiden"].nunique())
print(adata_combined.obs["leiden"].value_counts().sort_index())

# VISUALIZE AND SAVE UMAP FIGURES
# Create the figures directory if it does not exist.
figures_dir = Path("figures")
figures_dir.mkdir(parents=True, exist_ok=True)

# Figure 1 — UMAP colored by sample
# Examine mixing between the 27 samples.
sc.pl.umap(adata_combined, color="biosample_id", frameon=False, title="UMAP — scVI Integration by Sample", show=False)
plt.savefig(figures_dir / "umap_by_sample.png", dpi=300, bbox_inches="tight")
plt.show()
plt.close()

# Figure 2 — UMAP colored by experimental group
# Compare the distribution of COVID and control nuclei.
sc.pl.umap(adata_combined, color="group", frameon=False, title="UMAP — COVID vs Control", show=False)
plt.savefig(figures_dir / "umap_covid_vs_control.png",dpi=300, bbox_inches="tight")
plt.show()
plt.close()

# Figure 3 — UMAP colored by Leiden cluster
# Display the graph-based clusters identified by Leiden.
sc.pl.umap(adata_combined, color="leiden", legend_loc="on data", frameon=False, title="UMAP — Leiden Clusters", show=False)
plt.savefig(figures_dir / "umap_leiden_clusters.png", dpi=300,bbox_inches="tight")
plt.show()
plt.close()

print("UMAP figures saved to:", figures_dir.resolve())



# 5.8 SAVE AND LOAD INTEGRATED DATA

# Define output paths
integrated_path = Path("data/processed/combined_integrated.h5ad")
model_path = Path("data/processed/scvi_model")

# Create directory if needed
# integrated_path.parent.mkdir(parents=True, exist_ok=True)

# Verify donor correction before saving
# assert adata_combined.obs["donor_id"].nunique() == 26, \
# "Donor metadata is not corrected!"

# Save integrated AnnData
# Includes raw counts, metadata, HVGs, scVI latent representation,
# neighbor graph, UMAP coordinates, and Leiden clusters.
# adata_combined.write_h5ad(integrated_path, compression="gzip")

# Save trained scVI model separately
# model.save(model_path, overwrite=True)
print("Integrated dataset saved:", integrated_path.resolve())
print("scVI model saved:", model_path.resolve())
print("Dataset shape:", adata_combined.shape)
print("Leiden clusters:", adata_combined.obs["leiden"].nunique())



# LOAD SAVED DATA (FOR FUTURE SESSIONS)
# Run this only in a new session.
# No need to repeat preprocessing, scVI training, or clustering.
adata_combined = sc.read_h5ad(integrated_path)
# print("Loaded dataset:", adata_combined.shape)

# ========================================================================
# 6. MARKER GENES AND CELL-TYPE ANNOTATION
# ========================================================================
# Create a copy to preserve the original raw UMI counts
adata_markers = adata_combined.copy()

# Normalize each nucleus to 10,000 total counts
sc.pp.normalize_total(adata_markers, target_sum=1e4)

# Apply log(1 + x) transformation
sc.pp.log1p(adata_markers)

# Verify that Leiden labels are available
print("Number of clusters:", adata_markers.obs["leiden"].nunique())
print("Marker analysis dataset:", adata_markers.shape)

# 6.2 Rank marker genes by cluster
sc.tl.rank_genes_groups(adata_markers, groupby="leiden", method="wilcoxon", reference="rest", use_raw=False)

# Display the top 10 marker genes for each cluster
sc.pl.rank_genes_groups(adata_markers, n_genes=10, sharey=False)

# 6.3 Print top-ranked markers
# Extract ranked marker-gene results
markers = sc.get.rank_genes_groups_df(adata_markers, group=None)

# Display the top 10 genes for each cluster
top_markers = (markers.groupby("group", sort=False).head(10))
print(top_markers[["group", "names", "scores", "logfoldchanges", "pvals_adj"]].to_string(index=False))

# 6.4 Plot known cell-type markers

marker_genes = [
    "SFTPC", "SFTPB", "ABCA3",  # Alveolar type II
    "AGER", "PDPN", "CAV1",  # Alveolar type I
    "FOXJ1", "PIFO",  # Ciliated epithelial
    "PECAM1", "VWF",  # Endothelial
    "MARCO", "CD68", "MRC1",  # Macrophages
    "CD3D", "CD3E", "CD247",  # T cells
    "NKG7", "GNLY",  # NK cells
    "KIT", "CPA3",  # Mast cells
    "ACTA2", "TAGLN",  # Smooth muscle
    "COL1A1", "DCN"  # Fibroblasts
]

# Keep only genes present in the dataset
marker_genes = [gene for gene in marker_genes if gene in adata_markers.var_names]
sc.pl.dotplot(adata_markers, var_names=marker_genes, groupby="leiden", standard_scale="var")


# 6.5 EXAMINE CLUSTER-SPECIFIC MARKER GENES
# Extract the top 15 marker genes for each Leiden cluster.
# These genes will help identify cell types that were not clearly
# recognized using the known-marker dot plot.
top15 = (markers.groupby("group", sort=False).head(15))

# Create a compact summary of marker genes for each cluster
marker_summary = (top15.groupby("group", sort=False)["names"].apply(lambda genes: ", ".join(genes)))

# Display the marker genes for all clusters
for cluster, genes in marker_summary.items():
    print(f"Cluster {cluster}: {genes}")



# 6.6 VALIDATE CELL-TYPE MARKERS
# Examine established marker genes to validate the biological identity
# of Leiden clusters, especially clusters with uncertain annotations.
#
# Dot size: percentage of nuclei expressing each gene.
# Dot color: average expression, scaled separately for each gene.

validation_markers = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "AT1": ["AGER", "AQP5", "PDPN", "CAV1"],
    "AT2": ["SFTPC", "SFTPB", "ABCA3", "SFTPA1"],
    "Secretory": ["SCGB1A1", "SCGB3A2", "MUC5B"],
    "Smooth muscle": ["ACTA2", "TAGLN", "MYH11"],
    "Pericyte": ["RGS5", "PDGFRB", "CSPG4"],
    "Fibroblast": ["DCN", "LUM", "COL1A1", "PDGFRA"],
    "Endothelial": ["PECAM1", "VWF", "EMCN", "CA4", "ACKR1"],
    "Plasma cell": ["JCHAIN", "MZB1", "XBP1", "SDC1"],
    "Erythroid": ["HBB", "HBA1", "ALAS2", "GYPA"]
}

# Keep only marker genes present in the dataset.
# Remove categories with no available genes.
validation_markers = {
    cell_type: available
    for cell_type, genes in validation_markers.items()
    if (available := [gene for gene in genes if gene in adata_markers.var_names])}

# Create output directory if needed.
figures_dir = Path("figures")
figures_dir.mkdir(parents=True, exist_ok=True)

# Plot marker expression across Leiden clusters.
sc.pl.dotplot(adata_markers, var_names=validation_markers, groupby="leiden", standard_scale="var", figsize=(17, 8), show=False)
plt.savefig(figures_dir / "cell_type_validation.png", dpi=300, bbox_inches="tight")
plt.close("all")


# 6.7 VALIDATE DONOR REPRESENTATION
# Check whether each cluster is represented across multiple donors.
# Also calculate the percentage contributed by the largest donor.
#
# This helps identify clusters dominated by one biological donor.

# Verify donor metadata correction.
assert adata_combined.obs["donor_id"].nunique() == 26, \
    "Expected 26 biological donors. Check donor metadata."

# Count nuclei from each donor in each Leiden cluster.
donor_counts = (adata_combined.obs .groupby(["leiden", "donor_id"], observed=True) .size() .unstack(fill_value=0))

# Summarize cluster sizes and donor representation.
cluster_summary = pd.DataFrame({
    "n_nuclei": donor_counts.sum(axis=1),
    "n_donors": (donor_counts > 0).sum(axis=1),
    "max_donor_pct": (
        donor_counts.max(axis=1)
        / donor_counts.sum(axis=1) * 100)})

print("\nCLUSTER DONOR SUMMARY")
print(cluster_summary.round(2).sort_index())


# 6.8 ASSIGN BIOLOGICAL CELL-TYPE LABELS
# Assign provisional biological identities based on:
# - Cluster-specific marker genes
# - Established lung cell-type markers
# - Marker validation dot plots
# Some cell types remain broadly classified until further validation.

cell_type_map = {
    "0":  "AT2",
    "1":  "Secretory-like epithelial",
    "2":  "Mural cell",
    "3":  "AT1",
    "4":  "Fibroblast",
    "5":  "Fibroblast-like stromal",
    "6":  "Alveolar macrophage",
    "7":  "Monocyte",
    "8":  "Macrophage",
    "9":  "Ciliated epithelial",
    "10": "Endothelial",
    "11": "Endothelial - unclassified",
    "12": "Plasma cell",
    "13": "B cell",
    "14": "Mast cell",
    "15": "T cell",
    "16": "Secretory epithelial - putative",
    "17": "Erythroid"
}

# Assign a biological label to every nucleus using its Leiden cluster.
adata_combined.obs["cell_type"] = (adata_combined.obs["leiden"].astype(str).map(cell_type_map).fillna("Unassigned"))

# Verify annotation completeness.
print("\nCELL-TYPE ANNOTATION SUMMARY")
print("Unassigned nuclei:",(adata_combined.obs["cell_type"] == "Unassigned").sum())
print("\nNuclei per cell type:")
print(adata_combined.obs["cell_type"].value_counts())

# 6.9 VISUALIZE ANNOTATED CELL TYPES
# Display the same UMAP coordinates, colored by biological cell type.
# This does not recalculate UMAP or Leiden clustering.

sc.pl.umap(adata_combined,color="cell_type",frameon=False, title="UMAP — Annotated Lung Cell Types", legend_loc="right margin", show=False)
plt.savefig(figures_dir / "umap_cell_types.png",dpi=300, bbox_inches="tight")
plt.close("all")
print("\nAnnotation figures saved to:", figures_dir.resolve())


# 6.10 SAVE ANNOTATED DATASET
# Save raw UMI counts, metadata, scVI latent representation,
# UMAP coordinates, Leiden clusters, and cell-type labels.
#
# Keep the original preprocessed and integrated files unchanged.
annotated_path = Path("data/processed/combined_annotated.h5ad")
annotated_path.parent.mkdir(parents=True, exist_ok=True)
adata_combined.write_h5ad(annotated_path,compression="gzip")
print("\nANNOTATED DATASET SAVED")
print("File:", annotated_path.resolve())
print("Dataset shape:", adata_combined.shape)
print("Biological donors:", adata_combined.obs["donor_id"].nunique())
print("Leiden clusters:", adata_combined.obs["leiden"].nunique())
print("Cell-type labels:", adata_combined.obs["cell_type"].nunique())


# 6.11 INVESTIGATE DONOR-DOMINATED CLUSTERS
# Identify the largest contributing donors for selected clusters
for cluster in ["16", "17", "13", "2"]:
    counts = donor_counts.loc[cluster].sort_values(ascending=False)

    print(f"\nCluster {cluster}")
    print(counts.head(5))



# 6.11.1 QC CHECK FOR DONOR-DOMINATED CLUSTERS
# Compare sequencing depth, gene complexity, and mitochondrial
# percentages for clusters 16 and 17.
qc_comparison = (
    adata_combined.obs
    .assign(
        cluster_group=lambda df: df["leiden"].astype(str).map({
            "16": "Cluster 16",
            "17": "Cluster 17"
        }).fillna("Other clusters")).groupby("cluster_group", observed=True).agg(
        n_nuclei=("total_counts", "size"),
        median_genes=("n_genes_by_counts", "median"),
        median_UMIs=("total_counts", "median"),
        median_mt_pct=("pct_counts_mt", "median"),
        mt_95th_pct=("pct_counts_mt", lambda x: x.quantile(0.95))
    ))

print(qc_comparison.round(2))




# ========================================================================
# 7. COVID VS. CONTROL — CELL-TYPE COMPOSITION
# ========================================================================
# 7.1 Calculate cell-type proportions for each biological donor
# Each donor is treated as one independent biological observation.
# L04cov and L04covaddon are already assigned to the same donor.
assert adata_combined.obs["donor_id"].nunique() == 26
assert adata_combined.obs["cell_type"].notna().all()

# Count nuclei for each donor and cell type
cell_counts = (adata_combined.obs.groupby(["donor_id", "group", "cell_type"], observed=True).size().unstack(fill_value=0))

# Convert counts to percentages within each donor
cell_proportions = cell_counts.div(cell_counts.sum(axis=1),axis=0) * 100
print("Proportion table shape:", cell_proportions.shape)
print(cell_proportions.round(2).head())

# Verify percentages sum to 100 for each donor
assert np.allclose(cell_proportions.sum(axis=1), 100)


# 7.2 VISUALIZE CELL-TYPE COMPOSITION PER DONOR
# Arrange donors by experimental group
cell_proportions = cell_proportions.sort_index(level=["group", "donor_id"])

# Create stacked bar chart
ax = cell_proportions.plot(kind="bar", stacked=True,figsize=(17, 7),width=0.85)
ax.set_ylabel("Percentage of nuclei (%)")
ax.set_xlabel("Biological donor")
ax.set_title("Lung Cell-Type Composition by Donor")
ax.set_ylim(0, 100)

# Display donor IDs without the group label
ax.set_xticklabels(cell_proportions.index.get_level_values("donor_id"), rotation=90)
ax.legend(title="Cell type", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
figures_dir = Path("figures")
figures_dir.mkdir(parents=True, exist_ok=True)
plt.savefig(figures_dir / "cell_type_composition_by_donor.png", dpi=300, bbox_inches="tight")
plt.close()


# 7.3 COMPARE COVID AND CONTROL CELL-TYPE PROPORTIONS
# Calculate mean proportion across donors within each group
group_means = (cell_proportions.groupby(level="group").mean())
# Transpose for easier comparison
comparison = group_means.T
# Calculate the difference in percentage points
comparison["COVID_minus_Control"] = ( comparison["COVID-19"] - comparison["Control"])
# Sort by the largest absolute difference
comparison = comparison.reindex(comparison["COVID_minus_Control"].abs().sort_values(ascending=False).index)
print(comparison.round(2).to_string())


# 7.4 VISUALIZE DONOR-LEVEL CELL-TYPE PROPORTIONS
# Compare the distributions across biological donors.
# Each point represents one donor, not one nucleus.

selected_types = ["Macrophage", "AT2", "AT1"]
plot_data = (cell_proportions[selected_types].reset_index().melt(id_vars=["donor_id", "group"], var_name="cell_type", value_name="proportion"))
fig, axes = plt.subplots(1, 3, figsize=(13, 5))
for ax, cell_type in zip(axes, selected_types):
    subset = plot_data[plot_data["cell_type"] == cell_type]
    sns.boxplot(data=subset, x="group",  y="proportion", ax=ax, showfliers=False)
    sns.stripplot(data=subset, x="group", y="proportion", ax=ax, color="black", alpha=0.7, size=5)
    ax.set_title(cell_type)
    ax.set_xlabel("")
    ax.set_ylabel("Percentage of nuclei (%)")

plt.tight_layout()
plt.savefig(figures_dir / "cell_type_proportions_by_group.png",dpi=300, bbox_inches="tight")
plt.close()


# 7.5 STATISTICAL COMPARISON OF CELL-TYPE PROPORTIONS
# Compare donor-level cell-type proportions between COVID and Control.
# Each donor is one biological observation.
# Mann-Whitney U test is used for exploratory comparisons.
# Benjamini-Hochberg correction controls the false discovery rate.


results = []
for cell_type in cell_proportions.columns:
    covid = cell_proportions.xs("COVID-19", level="group")[cell_type]
    control = cell_proportions.xs("Control", level="group")[cell_type]
    statistic, p_value = mannwhitneyu(covid,control, alternative="two-sided")
    results.append({"cell_type": cell_type,"COVID_mean": covid.mean(), "Control_mean": control.mean(),"difference_pp": covid.mean() - control.mean(), "p_value": p_value})

results_df = pd.DataFrame(results)
# Correct for multiple testing across all cell types
results_df["p_adj"] = multipletests(results_df["p_value"], method="fdr_bh")[1]

# Sort by adjusted p-value
results_df = results_df.sort_values("p_adj")
print(results_df.round(4).to_string(index=False))



# 7.6 SENSITIVITY ANALYSIS — EXCLUDE CLUSTERS 16 AND 17
# Test whether the COVID vs Control cell-type composition results
# remain stable after excluding two donor-dominated clusters.
#
# Cluster 16: Putative secretory epithelial population, dominated by C56ctr.
# Cluster 17: Small erythroid population, dominated by L17cov.
#
# IMPORTANT:
# - The original AnnData object is not modified.
# - Proportions are recalculated within each donor.
# - Each donor remains one biological observation.


# 1. Exclude the selected clusters from the sensitivity analysis
obs_sensitivity = adata_combined.obs.loc[~adata_combined.obs["leiden"].astype(str).isin(["16", "17"])].copy()

# 2. Count nuclei per donor and cell type
counts_sensitivity = (obs_sensitivity.groupby(["donor_id", "group", "cell_type"], observed=True).size().unstack(fill_value=0))

# 3. Recalculate cell-type proportions within each donor
proportions_sensitivity = counts_sensitivity.div(counts_sensitivity.sum(axis=1), axis=0) * 100

# Verify that all donors remain represented
assert proportions_sensitivity.shape[0] == 26

# Verify that proportions sum to 100%
assert np.allclose(proportions_sensitivity.sum(axis=1), 100)

# 4. Compare COVID and Control donors
sensitivity_results = []
for cell_type in proportions_sensitivity.columns:
    covid = proportions_sensitivity.xs("COVID-19", level="group")[cell_type]
    control = proportions_sensitivity.xs("Control", level="group")[cell_type]
    statistic, p_value = mannwhitneyu(covid,control,alternative="two-sided")
    sensitivity_results.append({ "cell_type": cell_type, "COVID_mean": covid.mean(), "Control_mean": control.mean(), "difference_pp": covid.mean() - control.mean(), "p_value": p_value})

sensitivity_df = pd.DataFrame(sensitivity_results)

# 5. Correct for multiple comparisons
sensitivity_df["p_adj"] = multipletests(sensitivity_df["p_value"], method="fdr_bh")[1]

# 6. Compare the original results with the sensitivity analysis
comparison_sensitivity = results_df[["cell_type", "difference_pp", "p_adj"]].merge(sensitivity_df[["cell_type", "difference_pp", "p_adj"]],on="cell_type",suffixes=("_original", "_sensitivity"))
comparison_sensitivity = comparison_sensitivity.sort_values("p_adj_original")
print("\nSENSITIVITY ANALYSIS RESULTS")
print(comparison_sensitivity.round(4).to_string(index=False))





# ========================================================================
# 8. PSEUDOBULK DIFFERENTIAL GENE EXPRESSION
# ========================================================================
# Our goal is to identify which genes behave differently in COVID-19 patients compared with controls, within a specific cell type.
# Which genes are upregulated or downregulated in lung macrophages from COVID-19 patients?

# pseudobulk
# our dataset contains thousands of individual nuclei, but only 26 biological donors.
# We cannot treat every nucleus as an independent biological replicate because nuclei from the same donor are related.
# Instead, we combine the raw gene counts of all nuclei belonging to the same cell type and donor.
# We treat each donor within the same experimental group as an independent biological replicate

# 8.1 Check donor representation for selected cell types
# Each donor will become one biological replicate in pseudobulk DGE.
selected_cell_types = ["Macrophage", "AT2"]  # We selected macrophages and AT2 cells as a starting point

# Count nuclei per donor and cell type
dge_counts = (adata_combined.obs[adata_combined.obs["cell_type"].isin(selected_cell_types)].groupby(["donor_id", "group", "cell_type"],observed=True).size().unstack(fill_value=0))
print("\nNUCLEI PER DONOR AND CELL TYPE")
print(dge_counts.to_string())
# Summarize the number of donors represented in each cell type
print("\nBIOLOGICAL REPLICATES")
for cell_type in selected_cell_types:
    for group in ["COVID-19", "Control"]:
        counts = dge_counts.xs(group, level="group")[cell_type]
        print(
            f"{cell_type} | {group}: "
            f"{(counts > 0).sum()} donors, "
            f"median {counts.median():.0f} nuclei per donor, "
            f"minimum {counts.min()} nuclei")


# 8.2 CREATE PSEUDOBULK COUNT MATRICES
# Sum raw UMI counts across nuclei of the same cell type and donor.
# Each row represents one biological donor.
# Each column represents one gene.
# Keep macrophage and AT2 pseudobulk matrices separate.
selected_cell_types = ["Macrophage", "AT2"]
pseudobulk = {}
for cell_type in selected_cell_types:

    # Select nuclei belonging to the current cell type
    adata_subset = adata_combined[adata_combined.obs["cell_type"] == cell_type]

    # Aggregate raw UMI counts separately for each donor
    donor_counts = []
    donor_ids = []
    for donor_id in adata_subset.obs["donor_id"].unique():
        donor_mask = (adata_subset.obs["donor_id"].to_numpy() == donor_id)

        # Sum gene counts across nuclei from the same donor
        summed_counts = np.asarray(adata_subset.X[donor_mask].sum(axis=0)).ravel()
        donor_counts.append(summed_counts)
        donor_ids.append(donor_id)

    # Create donor-by-gene count matrix
    counts_df = pd.DataFrame(donor_counts, index=donor_ids, columns=adata_subset.var_names)
    pseudobulk[cell_type] = counts_df
    print(f"\n{cell_type}")
    print("Pseudobulk shape:", counts_df.shape)
    print("Total counts:", counts_df.to_numpy().sum())


# 8.3 PSEUDOBULK QUALITY CONTROL
# Examine library sizes and gene detection for each donor.
# Identify samples with unusually low counts before DGE.
pseudobulk_qc = {}
for cell_type, counts_df in pseudobulk.items():

    # Calculate QC metrics for each donor
    qc_df = pd.DataFrame({
        "donor_id": counts_df.index,
        "total_UMIs": counts_df.sum(axis=1).values,
        "detected_genes": (counts_df > 0).sum(axis=1).values})

    # Add experimental group and number of nuclei
    qc_df["group"] = (qc_df["donor_id"].map(adata_combined.obs.drop_duplicates("donor_id") .set_index("donor_id")["group"]))
    qc_df["n_nuclei"] = (qc_df["donor_id"].map(dge_counts[cell_type].groupby(level="donor_id").sum()))
    pseudobulk_qc[cell_type] = qc_df
    print(f"\n{cell_type.upper()} PSEUDOBULK QC")
    print(qc_df.sort_values("total_UMIs").to_string(index=False))


# 8.4 PREPARE PSEUDOBULK MATRICES FOR DIFFERENTIAL EXPRESSION
# Filter very low-expression genes separately for each cell type.
# Preserve raw integer counts for DESeq2.
pseudobulk_filtered = {}
for cell_type, counts_df in pseudobulk.items():
    # Verify raw count values
    assert np.isfinite(counts_df.to_numpy()).all()
    assert (counts_df.to_numpy() >= 0).all()
    assert np.allclose(counts_df.to_numpy(),np.round(counts_df.to_numpy()))

    # Require at least 10 total counts across donors
    # and expression in at least 3 donors.
    keep_genes = ((counts_df.sum(axis=0) >= 10) &((counts_df > 0).sum(axis=0) >= 3))
    filtered = counts_df.loc[:, keep_genes].astype(np.int64)
    pseudobulk_filtered[cell_type] = filtered
    print(f"\n{cell_type}")
    print("Genes before filtering:", counts_df.shape[1])
    print("Genes after filtering:", filtered.shape[1])
    print("Donors retained:", filtered.shape[0])


# Check pseudobulk sample quality
for cell_type, qc_df in pseudobulk_qc.items():
    print(f"\n{'=' * 60}")
    print(f"{cell_type.upper()} — SAMPLE QUALITY")
    print(f"{'=' * 60}")
    print(qc_df[["donor_id", "group", "n_nuclei","total_UMIs", "detected_genes"]].sort_values("n_nuclei").to_string(index=False))
    print("\nSummary by group:")
    print(qc_df.groupby("group")[["n_nuclei", "total_UMIs", "detected_genes"]].agg(["count", "median", "min"]))



# 8.5 PREPARE SAMPLE METADATA FOR DESEQ2
# Now we need a matching metadata table for each cell type, identifying which donor belongs to which experimental group.
pseudobulk_metadata = {}
for cell_type, counts_df in pseudobulk_filtered.items():

    # Retrieve the experimental group for each donor
    sample_metadata = (
        adata_combined.obs[["donor_id", "group"]]
        .drop_duplicates()
        .set_index("donor_id")
        .loc[counts_df.index]
        .copy()
    )

    # Verify metadata matches the count matrix
    assert sample_metadata.index.equals(counts_df.index)
    assert sample_metadata["group"].notna().all()
    pseudobulk_metadata[cell_type] = sample_metadata
    print(f"\n{cell_type.upper()}")
    print("Count matrix:", counts_df.shape)
    print("Metadata:", sample_metadata.shape)
    print(sample_metadata["group"].value_counts())