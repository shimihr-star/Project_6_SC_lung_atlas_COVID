
# ============================================================
# CHECK CELL-TYPE ANNOTATION
# ============================================================

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import scanpy as sc

from src.config import PROCESSED_DIR


def main():

    # 1. Load original and annotated datasets.
    original = sc.read_h5ad(
        PROCESSED_DIR / "combined_integrated.h5ad"
    )

    annotated = sc.read_h5ad(
        PROCESSED_DIR / "combined_annotated_v2.h5ad"
    )

    # 2. Check that nuclei and genes are preserved.
    assert original.shape == annotated.shape
    assert original.obs_names.equals(annotated.obs_names)
    assert original.var_names.equals(annotated.var_names)

    # 3. Check that Leiden labels are unchanged.
    assert original.obs["leiden"].astype(str).equals(
        annotated.obs["leiden"].astype(str)
    )

    # 4. Check that every nucleus has an annotation.
    assert "cell_type" in annotated.obs.columns
    assert annotated.obs["cell_type"].notna().all()

    # 5. Check that scVI and UMAP coordinates are retained.
    assert "X_scVI" in annotated.obsm
    assert "X_umap" in annotated.obsm

    # 6. Display annotation summary.
    print("Dataset shape:", annotated.shape)
    print("Leiden clusters:", annotated.obs["leiden"].nunique())
    print("Cell-type labels:", annotated.obs["cell_type"].nunique())

    print("\nCell-type distribution:")
    print(
        annotated.obs["cell_type"]
        .value_counts()
        .to_string()
    )

    print("\nStructural annotation checks passed.")


if __name__ == "__main__":
    main()
