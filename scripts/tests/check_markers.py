import pandas as pd

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

markers = pd.read_csv(
    PROJECT_ROOT / "data" / "results" / "cluster_markers.csv"
)

for cluster in sorted(markers["group"].unique(), key=int):

    cluster_markers = markers[
        (markers["group"] == cluster)
        & (markers["pvals_adj"] < 0.05)
        & (markers["logfoldchanges"] > 0)
    ]

    top_genes = cluster_markers.head(10)["names"].tolist()

    print(f"Cluster {cluster}: {', '.join(top_genes)}")