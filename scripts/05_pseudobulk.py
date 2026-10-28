
# ============================================================
# 05 — PSEUDOBULK DIFFERENTIAL GENE EXPRESSION ANALYSIS
# ============================================================
# PURPOSE:
# 1. Load the annotated single-nucleus RNA-seq dataset.
# 2. Verify the availability and integrity of raw UMI counts.
# 3. Examine nucleus representation per donor and cell type.
# 4. Aggregate raw counts by donor and cell type (pseudobulk).
# 5. Perform quality control and filter low-quality samples.
# 6. Normalize pseudobulk gene-expression data.
# 7. Compare COVID-19 and control donors using differential
#    gene expression analysis.
# 8. Save statistical results and visualization figures.
#
# IMPORTANT:
# Each donor represents an independent biological replicate.
# Nuclei from the same donor and cell type are aggregated.
# Raw UMI counts must be used for pseudobulk aggregation.
# The original annotated dataset must remain unchanged.
# ============================================================

import pandas as pd
from pathlib import Path
import anndata as ad
from src.pseudobulk import (
    validate_raw_counts,
    count_nuclei_per_donor,aggregate_pseudobulk)


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_PATH = (PROJECT_ROOT / "data" / "processed" / "combined_annotated_v2.h5ad")
RESULTS_DIR = PROJECT_ROOT / "data" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)



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
    counts = count_nuclei_per_donor(adata,cell_types,)

    for cell_type in cell_types:
        result = counts[counts["cell_type"] == cell_type].sort_values("n_nuclei")
        print(f"\n{cell_type} nuclei per donor:")
        print(result.to_string(index=False))
        print("\nSummary by group:")
        print(result.groupby("group")["n_nuclei"].agg(["count", "min", "median", "max"]))

    # --------------------------------------------------------
    # 4. Save QC results
    # --------------------------------------------------------
    output_path = RESULTS_DIR / "pseudobulk_donor_nuclei_qc.csv"
    counts.to_csv(output_path, index=False)
    print("\nSaved:", output_path)

    # --------------------------------------------------------
    # 5. Aggregate raw counts into pseudobulk samples
    # --------------------------------------------------------
    pseudobulk = aggregate_pseudobulk(adata, cell_types,)

    for cell_type, matrix in pseudobulk.items():

        print(f"\n{cell_type} pseudobulk matrix:")
        print("Shape:", matrix.shape)

        # Calculate QC metrics per donor.
        library_size = matrix.sum(axis=1)
        detected_genes = (matrix > 0).sum(axis=1)
        qc = pd.DataFrame({"donor_id": matrix.index,"library_size": library_size.to_numpy(),"detected_genes": detected_genes.to_numpy(),})
        print("\nPseudobulk QC:")
        print(qc.to_string(index=False))

        # Save pseudobulk QC metrics.
        safe_name = cell_type.lower().replace(" ", "_")
        qc.to_csv(RESULTS_DIR / f"{safe_name}_pseudobulk_qc.csv", index=False,)

        # Save aggregated raw count matrix.
        matrix.to_csv(RESULTS_DIR / f"{safe_name}_pseudobulk_counts.csv",)

if __name__ == "__main__":
    main()

# ======================================================================================================================