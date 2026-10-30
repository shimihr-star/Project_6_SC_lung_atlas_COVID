
# ============================================================
# 05 — PSEUDOBULK DIFFERENTIAL GENE EXPRESSION ANALYSIS
# ============================================================
#
# PURPOSE:
# 1. Load the annotated single-nucleus RNA-seq dataset.
# 2. Verify the integrity of raw UMI counts.
# 3. Examine nucleus representation per donor and cell type.
# 4. Aggregate raw counts by donor and cell type (pseudobulk).
# 5. Calculate sample-level pseudobulk QC metrics.
# 6. Generate QC figures and save results.
#
# FUTURE STEPS:
# - Evaluate and filter low-quality samples.
# - Filter low-expression genes.
# - Normalize pseudobulk gene-expression data.
# - Perform differential gene expression analysis.
#
# IMPORTANT:
# Each donor represents an independent biological replicate.
# Nuclei from the same donor and cell type are aggregated.
# Raw UMI counts must be used for pseudobulk aggregation.
# The original annotated dataset remains unchanged.
# ============================================================

import sys
from pathlib import Path
import anndata as ad
# Allow imports from src when running this script directly.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from src.pseudobulk import (validate_raw_counts, count_nuclei_per_donor, aggregate_pseudobulk,
                            calculate_pseudobulk_qc,plot_pseudobulk_qc, plot_pseudobulk_pca)

# --------------------------------------------------------
# Define input and output directories
# --------------------------------------------------------
DATA_PATH = (PROJECT_ROOT / "data" / "processed" / "combined_annotated_v2.h5ad")
RESULTS_DIR = PROJECT_ROOT / "data" / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"


# ======================================================================================================================
def main():
    # --------------------------------------------------------
    # 1. Load annotated dataset
    # --------------------------------------------------------
    adata = ad.read_h5ad(DATA_PATH)
    print("Dataset shape:", adata.shape)
    print("Available layers:", list(adata.layers.keys()))

    # --------------------------------------------------------
    # 2. Validate raw counts
    # --------------------------------------------------------
    validate_raw_counts(adata)

    # --------------------------------------------------------
    # 3. Examine donor representation
    # --------------------------------------------------------
    cell_types = ["Macrophage", "AT2"]
    counts = count_nuclei_per_donor(adata, cell_types,)

    for cell_type in cell_types:
        result = counts[counts["cell_type"] == cell_type].sort_values("n_nuclei")
        print(f"\n{cell_type} nuclei per donor:")
        print(result.to_string(index=False))
        print("\nSummary by group:")
        summary = (result.groupby("group", observed=True)["n_nuclei"].agg(["count", "min", "median", "max"]))
        print(summary)

    # --------------------------------------------------------
    # 4. Save donor nucleus-count QC
    # --------------------------------------------------------
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    output_path = (RESULTS_DIR / "pseudobulk_donor_nuclei_qc.csv")
    counts.to_csv(output_path, index=False)
    print("\nSaved:", output_path)

    # --------------------------------------------------------
    # 5. Aggregate raw counts into pseudobulk samples
    # --------------------------------------------------------
    pseudobulk = aggregate_pseudobulk(adata,cell_types,)

    # --------------------------------------------------------
    # 6. Calculate and visualize sample-level QC
    # --------------------------------------------------------
    for cell_type, matrix in pseudobulk.items():
        print(f"\n{cell_type} pseudobulk matrix:")
        print("Shape:", matrix.shape)
        # Calculate QC metrics for each donor.
        qc = calculate_pseudobulk_qc(matrix=matrix, nuclei_counts=counts, cell_type=cell_type,)
        print("\nPseudobulk sample QC:")
        print(qc.to_string(index=False))
        # Generate a safe filename for the cell type.
        safe_name = cell_type.lower().replace(" ", "_")

        # ----------------------------------------------------
        # 7. Save sample-level QC table
        # ----------------------------------------------------
        qc_path = (RESULTS_DIR / f"{safe_name}_pseudobulk_qc.csv")
        qc.to_csv(qc_path, index=False)
        print("\nSaved QC table:", qc_path)

        # ----------------------------------------------------
        # 8. Generate QC diagnostic plots
        # ----------------------------------------------------
        plot_pseudobulk_qc(qc=qc, cell_type=cell_type, output_dir=FIGURES_DIR,)

        # ----------------------------------------------------
        # 9. Exploratory PCA of pseudobulk samples
        # ----------------------------------------------------
        plot_pseudobulk_pca(matrix=matrix, qc=qc, cell_type=cell_type, output_dir=FIGURES_DIR,)

        # ----------------------------------------------------
        # 10. Save aggregated raw count matrix
        # ----------------------------------------------------
        counts_path = (RESULTS_DIR / f"{safe_name}_pseudobulk_counts.csv")
        matrix.to_csv(counts_path)
        print("Saved pseudobulk counts:", counts_path)

    print("\nPseudobulk aggregation and sample QC completed.")


if __name__ == "__main__":
    main()

# ======================================================================================================================
