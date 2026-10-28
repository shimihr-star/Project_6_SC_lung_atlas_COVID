
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import scanpy as sc

from src.config import PROCESSED_DIR
from src.integration import validate_integration


# Load the previously saved integrated dataset.
integrated_path = PROCESSED_DIR / "combined_integrated.h5ad"

adata = sc.read_h5ad(integrated_path)

print("Loaded:", integrated_path)
print("Dataset shape:", adata.shape)

# Check the latent representation, UMAP and clusters.
validate_integration(adata)

# Check donor and sample identities.
assert adata.obs["biosample_id"].nunique() == 27
assert adata.obs["donor_id"].nunique() == 26

# Examine the number of nuclei in each cluster.
print("\nNuclei per Leiden cluster:")
print(adata.obs["leiden"].value_counts().sort_index())
