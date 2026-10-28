
# ============================================================
# VISUALIZE LUNG CELL-TYPE MARKERS
# ============================================================

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import scanpy as sc
from src.config import PROCESSED_DIR, FIGURES_DIR
from src.annotation import (prepare_marker_data, plot_cell_type_markers,)


def main():

    # Load the existing integrated dataset.
    input_path = PROCESSED_DIR / "combined_integrated.h5ad"
    adata = sc.read_h5ad(input_path)
    print("Loaded dataset:", adata.shape)

    # Normalize a copy for visualization.
    adata_markers = prepare_marker_data(adata)

    # Generate marker dot plot.
    output_path = FIGURES_DIR / "cell_type_marker_dotplot.png"
    plot_cell_type_markers(adata_markers,output_path,)


if __name__ == "__main__":
    main()
