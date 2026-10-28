from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import scanpy as sc
from src.config import PROCESSED_DIR

# Load the existing checkpoint.
adata = sc.read_h5ad(PROCESSED_DIR / "combined.h5ad")

# Verify dataset dimensions and metadata.
print("Dataset shape:", adata.shape)
print("Samples:", adata.obs["biosample_id"].nunique())
print("Donors:", adata.obs["donor_id"].nunique())
print("HVGs:", adata.var["highly_variable"].sum())

assert adata.n_obs == 96915
assert adata.obs["biosample_id"].nunique() == 27
assert adata.obs["donor_id"].nunique() == 26
assert adata.obs_names.is_unique
assert adata.var_names.is_unique

print("Checkpoint validation passed.")