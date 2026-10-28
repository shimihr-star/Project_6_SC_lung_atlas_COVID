
from pathlib import Path

# ============================================================
# PROJECT CONFIGURATION
# ============================================================
#
# PURPOSE:
# Store shared project paths and reproducibility settings.
#
# All paths are defined relative to the project directory,
# allowing the project to run on Windows, macOS, or Linux.
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = PROJECT_ROOT / "data"
RAW_COUNTS_DIR = DATA_DIR / "raw_counts"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = DATA_DIR / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"

METADATA_PATH = DATA_DIR / "GSE171524_lung_metaData.txt"

RANDOM_SEED = 42

# Current preprocessing parameters.
# These preserve our existing methodology.
SOLO_DOUBLET_THRESHOLD = 0.70
MIN_GENES_PER_NUCLEUS = 200
MAX_GENES_QUANTILE = 0.98
MAX_MT_PERCENT = 20
MAX_RIBO_PERCENT = 2

N_HIGHLY_VARIABLE_GENES = 2000
