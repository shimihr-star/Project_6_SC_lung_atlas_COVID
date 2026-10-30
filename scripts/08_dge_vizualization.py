# =============================================================================
# 08_dge_vizualization.py
# Differential Gene Expression Visualization
# Macrophages and AT2 | COVID-19 vs Control
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load libraries
# -----------------------------------------------------------------------------
import pandas as pd
from pathlib import Path


# -----------------------------------------------------------------------------
# 2. Define project paths
# -----------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "data" / "results"
FIGURES_DIR = PROJECT_ROOT / "figures" / "dge_visualization"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


# -----------------------------------------------------------------------------
# 3. Load edgeR differential expression results
# -----------------------------------------------------------------------------
macrophage_dge = pd.read_csv(RESULTS_DIR / "macrophage_dge_edger.csv", index_col=0)
at2_dge = pd.read_csv(RESULTS_DIR / "at2_dge_edger.csv",index_col=0)


# -----------------------------------------------------------------------------
# 4. Validate loaded data
# -----------------------------------------------------------------------------

required_columns = {"logFC", "logCPM", "PValue", "FDR"}

for cell_type, df in [("Macrophage", macrophage_dge), ("AT2", at2_dge)]:
    assert required_columns.issubset(df.columns), (f"{cell_type}: Missing required columns")
    print(f"\n{cell_type}")
    print(f"Genes: {len(df)}")
    print(df.head())



# ----------------------------------------------------------------------------------------------------------------------
# 5. Prepare data for volcano plots
# ----------------------------------------------------------------------------------------------------------------------
import numpy as np
FDR_THRESHOLD = 0.05
LOGFC_THRESHOLD = 1.0


def prepare_volcano_data(df):
    # Copy the dataframe to preserve the original edgeR results.
    volcano_df = df.copy()

    # Store gene names as a regular column.
    volcano_df["gene"] = volcano_df.index

    # Calculate -log10(P-value).
    # Clipping prevents log10(0) if extremely small P-values occur.
    volcano_df["neg_log10_p"] = -np.log10(
        volcano_df["PValue"].clip(lower=np.finfo(float).tiny)
    )

    # Classify genes based on FDR and log2 fold change.
    volcano_df["status"] = "Not significant"

    volcano_df.loc[
        (volcano_df["FDR"] < FDR_THRESHOLD) &
        (volcano_df["logFC"] > LOGFC_THRESHOLD),
        "status"
    ] = "Upregulated"

    volcano_df.loc[
        (volcano_df["FDR"] < FDR_THRESHOLD) &
        (volcano_df["logFC"] < -LOGFC_THRESHOLD),
        "status"
    ] = "Downregulated"

    return volcano_df


macrophage_volcano = prepare_volcano_data(macrophage_dge)
at2_volcano = prepare_volcano_data(at2_dge)


# ----------------------------------------------------------------------------------------------------------------------
# 6. Check gene classifications
# ----------------------------------------------------------------------------------------------------------------------
for cell_type, df in [
    ("Macrophage", macrophage_volcano),
    ("AT2", at2_volcano)]:
    print(f"\n{cell_type} volcano classification:")
    print(df["status"].value_counts())


# ----------------------------------------------------------------------------------------------------------------------
# 7. Generate volcano plots
# ----------------------------------------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid")

COLORS = {
    "Not significant": "lightgray",
    "Downregulated": "steelblue",
    "Upregulated": "firebrick"
}


def plot_volcano(df, cell_type):
    fig, ax = plt.subplots(figsize=(9, 7))

    # Plot non-significant genes first so significant genes appear on top.
    for status in ["Not significant", "Downregulated", "Upregulated"]:
        subset = df[df["status"] == status]

        ax.scatter(
            subset["logFC"],
            subset["neg_log10_p"],
            color=COLORS[status],
            label=f"{status} (n={len(subset)})",
            alpha=0.6,
            s=12,
            edgecolors="none"
        )

    # Fold-change thresholds.
    ax.axvline(-LOGFC_THRESHOLD, color="gray", linestyle="--")
    ax.axvline(LOGFC_THRESHOLD, color="gray", linestyle="--")

    ax.set(
        title=f"{cell_type} — Differential Gene Expression",
        xlabel="Log2 Fold Change (COVID-19 vs Control)",
        ylabel="-Log10(P-value)"
    )

    ax.legend(loc="upper right", frameon=True)
    fig.tight_layout()

    return fig


# Generate plots.
fig_macrophage = plot_volcano(macrophage_volcano, "Macrophage")
fig_at2 = plot_volcano(at2_volcano, "AT2")

# Save plots.
fig_macrophage.savefig(
    FIGURES_DIR / "macrophage_volcano.png",
    dpi=300,
    bbox_inches="tight"
)

fig_at2.savefig(
    FIGURES_DIR / "at2_volcano.png",
    dpi=300,
    bbox_inches="tight"
)




# ----------------------------------------------------------------------------------------------------------------------
# 8. Label the top 10 most significant genes
# ----------------------------------------------------------------------------------------------------------------------

from adjustText import adjust_text


def add_gene_labels(fig, df, n_genes=10):
    ax = fig.axes[0]

    # Select significant genes, ordered by FDR and then P-value.
    top_genes = (
        df[df["status"] != "Not significant"]
        .sort_values(["FDR", "PValue"])
        .head(n_genes)
    )

    labels = []

    for _, gene in top_genes.iterrows():
        label = ax.text(
            gene["logFC"],
            gene["neg_log10_p"],
            gene["gene"],
            fontsize=9
        )
        labels.append(label)

    # Adjust labels to reduce overlap.
    adjust_text(
        labels,
        ax=ax,
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.6)
    )


add_gene_labels(fig_macrophage, macrophage_volcano)
add_gene_labels(fig_at2, at2_volcano)

print("Macrophage labels:", len(fig_macrophage.axes[0].texts))
print("AT2 labels:", len(fig_at2.axes[0].texts))

print("Macrophage figure:", fig_macrophage.number)
print("AT2 figure:", fig_at2.number)


# Save updated figures.
fig_macrophage.savefig(
    FIGURES_DIR / "macrophage_volcano_labeled.png",
    dpi=300,
    bbox_inches="tight"
)

fig_at2.savefig(
    FIGURES_DIR / "at2_volcano_labeled.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()



# ----------------------------------------------------------------------------------------------------------------------
# 9. Load TMM-normalized log2 CPM expression matrices
# ----------------------------------------------------------------------------------------------------------------------

macrophage_expr = pd.read_csv(
    RESULTS_DIR / "macrophage_log2cpm.csv",
    index_col=0
)

at2_expr = pd.read_csv(
    RESULTS_DIR / "at2_log2cpm.csv",
    index_col=0
)


# ----------------------------------------------------------------------------------------------------------------------
# 10. Validate normalized expression matrices
# ----------------------------------------------------------------------------------------------------------------------

for cell_type, expr, dge_results in [
    ("Macrophage", macrophage_expr, macrophage_dge),
    ("AT2", at2_expr, at2_dge)
]:
    print(f"\n{cell_type} normalized expression:")
    print("Donors:", expr.shape[0])
    print("Genes:", expr.shape[1])
    print("All DGE genes present:", set(dge_results.index).issubset(expr.columns))
    print("Missing values:", expr.isna().sum().sum())
    print("Donor IDs:", expr.index.tolist())



# ----------------------------------------------------------------------------------------------------------------------
# 11. Prepare heatmap data
# ----------------------------------------------------------------------------------------------------------------------

def prepare_heatmap_data(expr, dge, n_per_direction=10):

    # Select genes significantly upregulated in COVID-19.
    up_genes = (
        dge[
            (dge["FDR"] < FDR_THRESHOLD) &
            (dge["logFC"] > LOGFC_THRESHOLD)
        ]
        .sort_values("FDR")
        .head(n_per_direction)
        .index.tolist()
    )

    # Select genes significantly downregulated in COVID-19.
    down_genes = (
        dge[
            (dge["FDR"] < FDR_THRESHOLD) &
            (dge["logFC"] < -LOGFC_THRESHOLD)
        ]
        .sort_values("FDR")
        .head(n_per_direction)
        .index.tolist()
    )

    selected_genes = up_genes + down_genes

    # Heatmap format: genes as rows, donors as columns.
    heatmap_data = expr[selected_genes].T

    # Calculate a Z-score for each gene across donors.
    gene_means = heatmap_data.mean(axis=1)
    gene_stds = heatmap_data.std(axis=1)

    heatmap_z = heatmap_data.sub(gene_means, axis=0).div(
        gene_stds.replace(0, np.nan),
        axis=0
    ).fillna(0)

    return heatmap_z


macrophage_heatmap = prepare_heatmap_data(
    macrophage_expr, macrophage_dge
)

at2_heatmap = prepare_heatmap_data(
    at2_expr, at2_dge
)

print("\nMacrophage heatmap:", macrophage_heatmap.shape)
print("AT2 heatmap:", at2_heatmap.shape)


# ----------------------------------------------------------------------------------------------------------------------
# 12. Plot heatmaps
# ----------------------------------------------------------------------------------------------------------------------

def plot_dge_heatmap(data, cell_type):

    # Donor IDs beginning with C are controls;
    # donor IDs beginning with L are COVID-19 samples.
    donor_groups = [
        "Control" if donor.startswith("C") else "COVID-19"
        for donor in data.columns
    ]

    group_colors = pd.Series(
        donor_groups,
        index=data.columns
    ).map({
        "Control": "steelblue",
        "COVID-19": "firebrick"
    })

    grid = sns.clustermap(
        data,
        cmap="vlag",
        center=0,
        vmin=-2.5,
        vmax=2.5,
        col_colors=group_colors,
        row_cluster=False,
        col_cluster=False,
        figsize=(13, 9),
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "Gene-wise Z-score"}
    )

    grid.fig.suptitle(
        f"{cell_type} — Top Differentially Expressed Genes",
        y=1.02,
        fontsize=15
    )

    grid.ax_heatmap.set_xlabel("Donors")
    grid.ax_heatmap.set_ylabel("Genes")

    grid.ax_heatmap.tick_params(
        axis="x",
        labelrotation=90,
        labelsize=8
    )

    output_file = FIGURES_DIR / f"{cell_type.lower()}_heatmap.png"

    grid.fig.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(grid.fig)

    print(f"Saved: {output_file}")


plot_dge_heatmap(macrophage_heatmap, "Macrophage")
plot_dge_heatmap(at2_heatmap, "AT2")
