
# ============================================================
# 01 — PREPROCESS ALL SAMPLES
# ============================================================
#
# PURPOSE:
# Process all lung snRNA-seq samples, attach donor metadata,
# combine raw counts, and save the processed dataset.
#
# INPUT:
# - Raw count files
# - GSE171524 metadata
#
# OUTPUT:
# - data/processed/combined.h5ad
#
# WARNING:
# This script trains scVI and SOLO for every sample.
# Run only when preprocessing needs to be repeated.
# ============================================================

import sys
from pathlib import Path

# Allow imports from src when running this script directly.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd
import scanpy as sc
import scvi
from src.config import (RAW_COUNTS_DIR, PROCESSED_DIR, METADATA_PATH, RANDOM_SEED, N_HIGHLY_VARIABLE_GENES)
from src.preprocessing import (preprocess_sample, attach_metadata,)

scvi.settings.seed = RANDOM_SEED


# ======================================================================================================================
def load_metadata():
    """Load metadata and correct the additional L04cov sample."""

    metadata = pd.read_csv(METADATA_PATH,sep="\t", low_memory=False,)

    # The original file contains an extra annotation row.
    metadata = metadata.iloc[1:].copy()
    metadata.reset_index(drop=True, inplace=True)

    # Two sample files belong to the same biological donor.
    metadata.loc[metadata["biosample_id"] == "L04covaddon","donor_id",] = "L04cov"

    return metadata
# ======================================================================================================================


# ======================================================================================================================
def process_all_samples(metadata, ribo_genes):
    """Preprocess every sample and attach its metadata."""

    sample_files = sorted(RAW_COUNTS_DIR.glob("*_raw_counts.csv.gz"))
    print("Sample files found:", len(sample_files))
    if len(sample_files) != 27:
        raise ValueError(f"Expected 27 sample files, found {len(sample_files)}")
    processed_samples = {}

    for file_path in sample_files:
        # Example filename:
        # GSM5226574_C51ctr_raw_counts.csv.gz
        sample_id = file_path.name.split("_")[1]
        print(f"\nProcessing sample: {sample_id}")

        # Doublet detection and QC.
        adata = preprocess_sample(file_path,ribo_genes=ribo_genes,)

        # Attach donor and experimental-group information.
        adata = attach_metadata(adata,metadata,sample_id,)
        processed_samples[sample_id] = adata
        print(f"Completed {sample_id}: "f"{adata.n_obs} nuclei")

    return processed_samples
# ======================================================================================================================


# ======================================================================================================================
def combine_and_save(processed_samples):
    """Combine processed samples and preserve raw counts."""

    # Confirm all samples contain the same gene identifiers.
    reference_genes = set(next(iter(processed_samples.values())).var_names)

    for sample_id, adata in processed_samples.items():
        if set(adata.var_names) != reference_genes:
            raise ValueError( f"Gene sets differ in {sample_id}")

    # Combine nuclei from all samples.
    adata_combined = sc.concat(processed_samples, join="inner", label="sample_id", index_unique="-", merge="same",)

    # Confirm biological sample and donor counts.
    assert adata_combined.obs["biosample_id"].nunique() == 27
    assert adata_combined.obs["donor_id"].nunique() == 26

    # Each donor must belong to only one experimental group.
    assert (adata_combined.obs.groupby("donor_id", observed=True)["group"].nunique().eq(1).all())

    # Preserve numeric metadata types.
    numeric_cols = ["intubation_days", "interval_death_symptoms_onset_days", "pmi_h",]
    for col in numeric_cols:
        adata_combined.obs[col] = pd.to_numeric(adata_combined.obs[col],errors="coerce",)

    # Remove genes detected in fewer than 3 nuclei.
    sc.pp.filter_genes(adata_combined, min_cells=3)

    # Mark highly variable genes without removing other genes.
    sc.pp.highly_variable_genes(adata_combined, n_top_genes=N_HIGHLY_VARIABLE_GENES, flavor="seurat_v3", subset=False,)

    # Save raw counts, metadata, QC metrics, and HVG annotations.
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    output_path = PROCESSED_DIR / "combined.h5ad"
    if output_path.exists():
        raise FileExistsError(f"{output_path} already exists. ""Refusing to overwrite the saved dataset.")

    adata_combined.write_h5ad(output_path,compression="gzip",)
    print("Saved:", output_path)
    print("Dataset shape:", adata_combined.shape)
    return adata_combined
# ======================================================================================================================


# ======================================================================================================================
def load_ribosomal_genes():
    """Load the ribosomal gene set used in the original pipeline."""

    ribo_url = ("http://software.broadinstitute.org/gsea/msigdb/""download_geneset.jsp?""geneSetName=KEGG_RIBOSOME&fileType=txt" )
    ribo_genes = pd.read_table(ribo_url, skiprows=2,header=None,)
    return ribo_genes[0].values
# ======================================================================================================================


# ======================================================================================================================
def main():
    metadata = load_metadata()
    ribo_genes = load_ribosomal_genes()
    processed_samples = process_all_samples(metadata,  ribo_genes,)
    combine_and_save(processed_samples)

if __name__ == "__main__":
    main()

# ======================================================================================================================