# ============================================================
# 04 — CELL-TYPE COMPOSITION ANALYSIS
# ============================================================
# PURPOSE:
# 1. Load the annotated single-nucleus RNA-seq dataset.
# 2. Calculate cell-type counts and percentages per donor.
# 3. Validate donor-level cell-type composition.
# 4. Extract donor-level clinical metadata.
# 5. Compare cell-type proportions between COVID-19 and
#    control donors.
# 6. Generate cell-type composition visualizations.
# 7. Perform Mann-Whitney U tests for each cell type.
# 8. Apply Benjamini-Hochberg FDR correction.
# 9. Examine potential confounding factors, including
#    age, sex, and postmortem interval.
# 10. Perform age- and sex-adjusted regression analysis
#     for macrophage and AT2 proportions.
# 11. Save statistical results, clinical summaries,
#     and visualization figures.
#
# IMPORTANT:
# Each donor represents an independent biological replicate.
# Cell-type proportions are calculated from recovered nuclei.
# Differences in proportions do not necessarily indicate
# absolute changes in cell abundance.
# Age and sex are treated as potential confounders.
# Missing postmortem interval data prevent PMI adjustment.
# Cell-type annotations remain provisional.
# The original annotated dataset remains unchanged.
# ============================================================

import sys
import pandas as pd
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
import scanpy as sc
from src.config import PROCESSED_DIR, RESULTS_DIR, FIGURES_DIR
from src.composition import (calculate_cell_type_composition, get_donor_metadata, summarize_composition_by_group,
                             plot_composition_by_group, plot_group_composition_barplot, test_cell_type_composition,
                             summarize_donor_characteristics, adjusted_composition_analysis)


def main():

    # --------------------------------------------------------
    # 1. Load annotated dataset
    # --------------------------------------------------------
    input_path = PROCESSED_DIR / "combined_annotated_v2.h5ad"
    adata = sc.read_h5ad(input_path)
    print("Loaded dataset:", adata.shape)

    print("\nAvailable metadata columns:")
    print(adata.obs.columns.tolist())
    # --------------------------------------------------------
    # 2. Calculate cell-type composition per donor
    # --------------------------------------------------------
    counts, percentages = calculate_cell_type_composition(adata)

    # --------------------------------------------------------
    # 3. Validate results
    # --------------------------------------------------------
    assert counts.values.sum() == adata.n_obs
    assert (abs(percentages.sum(axis=1) - 100) < 1e-6).all()
    print("Number of donors:", len(counts))
    print("Number of cell types:", len(counts.columns))


    # 3.1. Inspect donor coverage for additional cell types
    print("\nAvailable cell types:")
    print(counts.columns.tolist())
    print("\nDonor-level nuclei counts:")
    print(counts.to_string())

    # --------------------------------------------------------
    # 4. Extract donor-level clinical metadata
    # --------------------------------------------------------
    donor_metadata = get_donor_metadata(adata)

    # Ensure the same donors appear in both tables.
    assert set(donor_metadata.index) == set(percentages.index)

    # Display disease-group information.
    print("\nDonor group distribution:")
    print(donor_metadata["group"].value_counts())
    print("\nGroup and disease combinations:")
    print(donor_metadata.groupby(["group", "disease"],observed=True).size())

    # --------------------------------------------------------
    # 5. Save results
    # --------------------------------------------------------
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    counts_path = RESULTS_DIR / "donor_cell_type_counts.csv"
    percentages_path = RESULTS_DIR / "donor_cell_type_percentages.csv"
    counts.to_csv(counts_path)
    percentages.to_csv(percentages_path)
    print("\nSaved:", counts_path)
    print("Saved:", percentages_path)

    # --------------------------------------------------------
    # 6. Display a preview
    # --------------------------------------------------------
    print("\nCell-type percentages (first 5 donors):")
    print(percentages.head().round(2).to_string())

    # --------------------------------------------------------
    # 7. Compare COVID-19 and control composition
    # --------------------------------------------------------
    summary = summarize_composition_by_group(percentages,donor_metadata,)
    print("\nMedian cell-type percentages by group:")
    print(summary.round(2).to_string())

    # Save group-level summary.
    summary_path = RESULTS_DIR / "cell_type_group_summary.csv"
    summary.to_csv(summary_path)
    print("\nSaved:", summary_path)

    # --------------------------------------------------------
    # 8. Visualize donor-level cell-type composition
    # --------------------------------------------------------
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    figure_path = FIGURES_DIR / "cell_composition_by_group.png"
    plot_composition_by_group(percentages, donor_metadata, figure_path,)
    print("\nSaved composition figure:", figure_path)

    # --------------------------------------------------------
    # 9. Plot median cell-type composition by group
    # --------------------------------------------------------
    barplot_path = FIGURES_DIR / "cell_composition_group_barplot.png"
    plot_group_composition_barplot(summary,barplot_path,)

    # --------------------------------------------------------
    # 10. Statistical comparison of cell-type composition
    # --------------------------------------------------------
    stats_results = test_cell_type_composition( percentages, donor_metadata,)

    print("\nStatistical comparison: COVID-19 vs Control")
    print(stats_results.round(4).to_string(index=False))

    # Save statistical results.
    stats_path = RESULTS_DIR / "cell_type_composition_statistics.csv"
    stats_results.to_csv(stats_path, index=False,)
    print("\nSaved statistical results:", stats_path)

    # Display FDR-significant cell types.
    significant = stats_results[stats_results["fdr"] < 0.05]
    print("\nCell types with FDR < 0.05:")
    print(significant[["cell_type", "median_difference_pp", "p_value", "fdr"]].round(4).to_string(index=False))

    # --------------------------------------------------------
    # 11. Examine potential confounding factors
    # --------------------------------------------------------
    clinical = summarize_donor_characteristics(adata,counts,)
    print("\nDonor-level clinical characteristics:")
    print(clinical.to_string())
    print("\nSex distribution by group:")
    print(pd.crosstab(clinical["group"], clinical["sex"],))
    print("\nMissing values:")
    print(clinical.isna().sum())

    # Save clinical metadata.
    clinical_path = RESULTS_DIR / "donor_clinical_characteristics.csv"
    clinical.to_csv(clinical_path)
    print("\nSaved:", clinical_path)

    # --------------------------------------------------------
    # 12. Age- and sex-adjusted composition analysis
    # --------------------------------------------------------
    adjusted_results = adjusted_composition_analysis(percentages, clinical,)
    print("\nAge- and sex-adjusted results:")
    print(adjusted_results.round(4).to_string(index=False))
    adjusted_path = (RESULTS_DIR / "cell_type_composition_adjusted.csv")
    adjusted_results.to_csv(adjusted_path, index=False,)
    print("\nSaved:", adjusted_path)


if __name__ == "__main__":
    main()
