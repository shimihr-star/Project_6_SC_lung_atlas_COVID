
# ============================================================
# PREPROCESSING FUNCTIONS
# ============================================================
#
# PURPOSE:
# Prepare individual lung snRNA-seq samples for analysis.
#
# This module handles:
# 1. Matching nucleus barcodes to sample metadata.
# 2. Doublet detection using scVI + SOLO.
# 3. Quality-control filtering.
#
# Raw UMI counts must be preserved for downstream
# pseudobulk differential expression.
# ============================================================



from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from scipy import sparse
from src.config import (SOLO_DOUBLET_THRESHOLD, MIN_GENES_PER_NUCLEUS, MAX_GENES_QUANTILE, MAX_MT_PERCENT,
                        MAX_RIBO_PERCENT, N_HIGHLY_VARIABLE_GENES,)





# ======================================================================================================================
def attach_metadata(adata, metadata, sample_id):
    """
    Attach donor and experimental-group metadata to nuclei.
    Only nuclei with matching metadata are retained.
    Existing QC metrics in adata.obs are preserved.
    """

    # Select metadata belonging to the current sample.
    sample_meta = metadata.loc[metadata["biosample_id"] == sample_id].copy()

    # Match nuclei using their barcode identifiers.
    sample_meta = sample_meta.set_index("NAME")

    # Check that metadata identifiers are unique.
    if not sample_meta.index.is_unique:
        raise ValueError(f"Duplicate metadata barcodes in {sample_id}")

    # Identify nuclei with matching metadata.
    valid_cells = adata.obs_names.isin(sample_meta.index)
    print(f"\nSample: {sample_id}")
    print("Nuclei before matching:", adata.n_obs)
    print("Unmatched nuclei:", (~valid_cells).sum())

    # Retain only nuclei with matching metadata.
    adata = adata[valid_cells].copy()

    # Join metadata while preserving existing QC columns.
    adata.obs = adata.obs.join(sample_meta, how="left")

    # Verify that required metadata is present.
    required_columns = ["biosample_id", "donor_id", "group"]
    if adata.n_obs == 0:
        raise ValueError(f"No matching nuclei in {sample_id}")

    if adata.obs[required_columns].isna().any().any():
        raise ValueError(f"Missing required metadata in {sample_id}")

    if not (adata.obs["biosample_id"] == sample_id).all():
        raise ValueError(f"Incorrect sample assignment in {sample_id}")

    print("Nuclei after matching:", adata.n_obs)
    return adata
# ======================================================================================================================


# ======================================================================================================================
def preprocess_sample(file_path, ribo_genes, doublet_threshold=SOLO_DOUBLET_THRESHOLD,):
    """
    Process one snRNA-seq sample.
    Steps:
    1. Load raw UMI counts.
    2. Detect doublets using scVI and SOLO.
    3. Calculate nucleus-level QC metrics.
    4. Filter low-quality nuclei.

    Returns an AnnData object containing raw UMI counts.
    """

    # --------------------------------------------------------
    # 1. Load raw counts
    # --------------------------------------------------------
    adata = sc.read_csv(file_path).T
    adata.X = sparse.csr_matrix(adata.X)
    print(f"\nProcessing: {Path(file_path).name}")
    print("Initial shape:", adata.shape)

    # --------------------------------------------------------
    # 2. Detect doublets with scVI + SOLO
    # --------------------------------------------------------
    # Use a temporary copy so HVG selection does not
    # remove genes from the original raw-count dataset.
    adata_solo = adata.copy()
    sc.pp.filter_genes(adata_solo, min_cells=10)
    sc.pp.highly_variable_genes(adata_solo,n_top_genes=N_HIGHLY_VARIABLE_GENES, subset=True, flavor="seurat_v3",)
    scvi.model.SCVI.setup_anndata(adata_solo)
    vae = scvi.model.SCVI(adata_solo)
    vae.train()
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()

    # Obtain doublet classification scores.
    predictions = solo.predict()
    predictions["prediction"] = solo.predict(soft=False)
    if not predictions.index.isin(adata.obs_names).all():
        raise ValueError("SOLO nucleus IDs do not match raw data.")
    doublets = predictions[predictions["doublet"] > doublet_threshold]
    print("SOLO predicted doublets:",(predictions["prediction"] == "doublet").sum())
    print("High-confidence doublets:", len(doublets))

    # Remove doublets from the full raw-count dataset.
    adata = adata[~adata.obs_names.isin(doublets.index)].copy()
    print("After doublet removal:", adata.shape)

    # --------------------------------------------------------
    # 3. Calculate QC metrics
    # --------------------------------------------------------
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.isin(ribo_genes)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True,)

    # --------------------------------------------------------
    # 4. Filter low-quality nuclei
    # --------------------------------------------------------
    n_before_qc = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=MIN_GENES_PER_NUCLEUS,)
    upper_lim = np.quantile(adata.obs["n_genes_by_counts"], MAX_GENES_QUANTILE,)
    adata = adata[adata.obs["n_genes_by_counts"] < upper_lim].copy()
    adata = adata[adata.obs["pct_counts_mt"] < MAX_MT_PERCENT].copy()
    adata = adata[adata.obs["pct_counts_ribo"] < MAX_RIBO_PERCENT].copy()

    # --------------------------------------------------------
    # 5. Report filtering results
    # --------------------------------------------------------
    n_removed = n_before_qc - adata.n_obs
    print("Final shape:", adata.shape)
    print("Nuclei removed during QC:", n_removed)
    print("Median detected genes:", adata.obs["n_genes_by_counts"].median(),)

    return adata
# ======================================================================================================================