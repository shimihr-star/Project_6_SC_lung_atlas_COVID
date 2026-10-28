
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import scanpy as sc

from src.config import PROCESSED_DIR, RESULTS_DIR
from src.annotation import calculate_donor_enrichment


def main():

    adata = sc.read_h5ad(
        PROCESSED_DIR / "combined_integrated.h5ad"
    )

    enrichment = calculate_donor_enrichment(adata)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    enrichment.to_csv(
        RESULTS_DIR / "cluster_donor_enrichment.csv"
    )

    # Focus on the clusters with unusual donor representation.
    for cluster, donor in [
        ("13", "C53ctr"),
        ("16", "C56ctr"),
        ("17", "L17cov"),
    ]:
        ratio = enrichment.loc[cluster, donor]

        print(
            f"Cluster {cluster} | Donor {donor} | "
            f"Enrichment: {ratio:.2f}x"
        )


if __name__ == "__main__":
    main()
