
# ============================================================
# 03 — CELL-TYPE ANNOTATION
# ============================================================
# PURPOSE:
# 1. Load the integrated dataset.
# 2. Identify marker genes for each Leiden cluster.
# 3. Compare markers with known lung cell-type markers.
# 4. Save marker analysis results.
# 5. Assign provisional biological cell-type labels.
# 6. Save and validate the annotated dataset.
#
# IMPORTANT:
# Preserves the original integrated dataset and raw UMI counts.
# Cell-type assignments remain provisional.
# ============================================================

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import scanpy as sc

from src.config import PROCESSED_DIR, RESULTS_DIR
from src.annotation import (prepare_marker_data, identify_cluster_markers, compare_markers_with_reference, assign_cell_types,)

# ======================================================================================================================
def main():

    output_path = PROCESSED_DIR / "combined_annotated_v2.h5ad"
    if output_path.exists():
        raise FileExistsError(f"Annotated dataset already exists: {output_path}")

    # --------------------------------------------------------
    # 1. Load integrated dataset
    # --------------------------------------------------------
    input_path = PROCESSED_DIR / "combined_integrated.h5ad"
    adata = sc.read_h5ad(input_path)
    print("Loaded dataset:", adata.shape)
    print("Leiden clusters:", adata.obs["leiden"].nunique())

    # --------------------------------------------------------
    # 2. Normalize a copy for marker analysis
    # --------------------------------------------------------
    adata_markers = prepare_marker_data(adata)

    # --------------------------------------------------------
    # 3. Identify cluster-specific marker genes
    # --------------------------------------------------------
    markers = identify_cluster_markers(adata_markers)

    # --------------------------------------------------------
    # 4. Compare markers with biological references
    # --------------------------------------------------------
    matches = compare_markers_with_reference(markers, top_n=50,)

    # --------------------------------------------------------
    # 5. Save results
    # --------------------------------------------------------
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    marker_path = RESULTS_DIR / "cluster_markers.csv"
    matches_path = RESULTS_DIR / "cell_type_candidates.csv"
    markers.to_csv(marker_path, index=False)
    matches.to_csv(matches_path, index=False)
    print("\nMarker genes saved:", marker_path)
    print("Candidate cell types saved:", matches_path)

    # --------------------------------------------------------
    # 6. Display the best candidate for each cluster
    # --------------------------------------------------------
    best_matches = (matches.sort_values(["cluster", "n_matching_markers"],ascending=[True, False],).drop_duplicates(subset="cluster"))
    print("\nCandidate cell types:")
    print(best_matches[["cluster", "candidate_cell_type", "n_matching_markers", "matching_genes",]].to_string(index=False))

    # --------------------------------------------------------
    # 7. Display top 10 marker genes for every cluster
    # --------------------------------------------------------

    print("\nTop 10 marker genes per Leiden cluster:")

    for cluster in sorted(markers["group"].unique(),key=int):
        cluster_markers = markers[(markers["group"] == cluster) & (markers["pvals_adj"] < 0.05) & (markers["logfoldchanges"] > 0)]
        top_genes = cluster_markers.head(10)["names"].tolist()
        print(f"Cluster {cluster}: {', '.join(top_genes)}")

    # --------------------------------------------------------
    # 8. Assign biological cell-type labels
    # --------------------------------------------------------
    adata_annotated = assign_cell_types(adata)

    # --------------------------------------------------------
    # 9. Validate annotations
    # --------------------------------------------------------
    assert adata_annotated.shape == adata.shape
    assert adata_annotated.obs_names.equals(adata.obs_names)
    assert adata_annotated.var_names.equals(adata.var_names)
    assert adata_annotated.obs["cell_type"].notna().all()

    # --------------------------------------------------------
    # 10. Save annotated dataset
    # --------------------------------------------------------
    if output_path.exists():
        raise FileExistsError(f"Annotated dataset already exists: {output_path}")

    adata_annotated.write_h5ad(output_path,compression="gzip",)
    print("\nAnnotated dataset saved:", output_path)


# This must remain OUTSIDE main().
if __name__ == "__main__":
    main()
