
import sys
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
import scanpy as sc
from src.config import PROCESSED_DIR, RESULTS_DIR
from src.annotation import check_donor_representation


def main():

    # Load existing integration results.
    adata = sc.read_h5ad(
        PROCESSED_DIR / "combined_integrated.h5ad"
    )

    # Calculate donor representation.
    counts, percentages, summary = check_donor_representation(adata)

    # Save results.
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    counts.to_csv(RESULTS_DIR / "cluster_donor_counts.csv")
    percentages.to_csv(RESULTS_DIR / "cluster_donor_percentages.csv")
    summary.to_csv(
        RESULTS_DIR / "cluster_donor_summary.csv",
        index=False,
    )

    print("\nDonor representation by cluster:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
