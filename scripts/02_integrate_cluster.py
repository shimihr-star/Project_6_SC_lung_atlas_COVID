
# ============================================================
# 02 — scVI INTEGRATION AND LEIDEN CLUSTERING
# ============================================================
#
# PURPOSE:
# Integrate all samples using scVI, construct a neighbor graph,
# perform UMAP and Leiden clustering, and save the results.
#
# INPUT:
# data/processed/combined.h5ad
#
# OUTPUT:
# data/processed/combined_integrated.h5ad
# data/processed/scvi_model/
# figures/umap_by_sample.png
# figures/umap_covid_vs_control.png
# figures/umap_leiden_clusters.png
#
# IMPORTANT:
# This script retrains scVI. Do not run unless you want
# to regenerate the integration results.
# ============================================================


import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt
import scanpy as sc
import scvi

from src.config import (PROCESSED_DIR, FIGURES_DIR, RANDOM_SEED,)
from src.integration import (train_scvi, cluster_nuclei, validate_integration,)




# ======================================================================================================================
def save_umap_figures(adata):
    """
    Generate UMAP plots colored by sample,
    experimental group, and Leiden cluster.
    """

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    plots = [("biosample_id", "umap_by_sample.png", "UMAP — scVI Integration by Sample"),
             ("group", "umap_covid_vs_control.png","UMAP — COVID vs Control"),
             ("leiden", "umap_leiden_clusters.png","UMAP — Leiden Clusters"),]

    for column, filename, title in plots:
        sc.pl.umap(adata, color=column, frameon=False, title=title, show=False,)
        plt.savefig(FIGURES_DIR / filename, dpi=300, bbox_inches="tight",)
        plt.close("all")
    print("UMAP figures saved to:", FIGURES_DIR)
# ======================================================================================================================


# ======================================================================================================================
def main():

    # --------------------------------------------------------
    # 1. Define input and output paths
    # --------------------------------------------------------
    input_path = PROCESSED_DIR / "combined.h5ad"
    integrated_path = PROCESSED_DIR / "combined_integrated.h5ad"
    model_path = PROCESSED_DIR / "scvi_model"

    # Prevent accidental overwriting of existing results.
    if integrated_path.exists() or model_path.exists():
        raise FileExistsError("Existing integration results detected. " "Refusing to retrain or overwrite them.")

    # --------------------------------------------------------
    # 2. Load preprocessed raw-count dataset
    # --------------------------------------------------------
    adata = sc.read_h5ad(input_path)
    print("Loaded dataset:", adata.shape)

    # --------------------------------------------------------
    # 3. Train scVI
    # --------------------------------------------------------
    model, latent = train_scvi(adata, n_latent=10,)

    # --------------------------------------------------------
    # 4. Construct graph, UMAP, and Leiden clusters
    # --------------------------------------------------------
    adata = cluster_nuclei(adata, latent, n_neighbors=15,leiden_resolution=0.5,)

    # --------------------------------------------------------
    # 5. Validate integration output
    # --------------------------------------------------------
    validate_integration(adata)

    # --------------------------------------------------------
    # 6. Save integrated dataset and trained model
    # --------------------------------------------------------
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(integrated_path, compression="gzip",)
    model.save(model_path, overwrite=False)
    print("Integrated dataset saved:", integrated_path)
    print("scVI model saved:", model_path)

    # --------------------------------------------------------
    # 7. Generate UMAP figures
    # --------------------------------------------------------
    save_umap_figures(adata)
    print("Integration pipeline completed.")


if __name__ == "__main__":
    main()
# ======================================================================================================================