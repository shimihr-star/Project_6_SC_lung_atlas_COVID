# ----------------------------------------------------------------------------------------------------------------------
# 1. Imports and project paths
# ----------------------------------------------------------------------------------------------------------------------
from pathlib import Path
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "data" / "results"
FIGURES_DIR = PROJECT_ROOT / "figures" / "pathway_enrichment"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

print("GSEApy version:", gp.__version__)


# ----------------------------------------------------------------------------------------------------------------------
# 2. Load differential expression results
# ----------------------------------------------------------------------------------------------------------------------
macrophage_dge = pd.read_csv(RESULTS_DIR / "macrophage_dge_edger.csv", index_col=0)

at2_dge = pd.read_csv(RESULTS_DIR / "at2_dge_edger.csv",index_col=0)
print("\nMacrophage DGE shape:", macrophage_dge.shape)
print("AT2 DGE shape:", at2_dge.shape)
print("\nMacrophage columns:")
print(macrophage_dge.columns.tolist())
print("\nFirst five macrophage genes:")
print(macrophage_dge.head())

# ----------------------------------------------------------------------------------------------------------------------
# 3. Prepare ranked gene lists for GSEA
# ----------------------------------------------------------------------------------------------------------------------
def prepare_gsea_ranking(dge):
    df = dge.copy()

    # Remove genes with missing or invalid statistics
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["logFC", "F"])

    # Calculate signed ranking statistic
    df["rank_score"] = np.sign(df["logFC"]) * np.sqrt(df["F"].clip(lower=0))

    # Remove duplicated gene names and sort by ranking score
    df = df[~df.index.duplicated(keep="first")]
    df = df.sort_values("rank_score", ascending=False)
    return df[["rank_score"]]


macrophage_ranking = prepare_gsea_ranking(macrophage_dge)
at2_ranking = prepare_gsea_ranking(at2_dge)

# Save ranked gene lists for reproducibility
macrophage_ranking.to_csv(RESULTS_DIR / "macrophage_gsea_ranking.csv")
at2_ranking.to_csv(RESULTS_DIR / "at2_gsea_ranking.csv")

print("\nMacrophage ranking:", macrophage_ranking.shape)
print("AT2 ranking:", at2_ranking.shape)
print("\nTop 5 macrophage upregulated genes:")
print(macrophage_ranking.head())
print("\nTop 5 macrophage downregulated genes:")
print(macrophage_ranking.tail())
print("\nMissing values:")
print(macrophage_ranking.isna().sum())
print("\nDuplicated gene names:")
print(macrophage_ranking.index.duplicated().sum())

# ----------------------------------------------------------------------------------------------------------------------
# 4. Check available pathway libraries
# ----------------------------------------------------------------------------------------------------------------------
libraries = gp.get_library_name(organism="Human")
hallmark_libraries = [name for name in libraries if "hallmark" in name.lower()]

print("\nAvailable Hallmark libraries:")
for name in hallmark_libraries:
    print(name)


# ----------------------------------------------------------------------------------------------------------------------
# 5. Run preranked GSEA — Macrophages
# ----------------------------------------------------------------------------------------------------------------------
HALLMARK_LIBRARY = "MSigDB_Hallmark_2020"
macrophage_gsea = gp.prerank(rnk=macrophage_ranking,gene_sets=HALLMARK_LIBRARY, outdir=str(RESULTS_DIR
                             / "gsea_macrophage"),min_size=10, max_size=500, permutation_num=1000, seed=42, threads=4, verbose=True)

print("\nMacrophage GSEA completed.")
print("\nGSEA result columns:")
print(macrophage_gsea.res2d.columns.tolist())
print("\nTop 10 enriched pathways:")
print(macrophage_gsea.res2d.sort_values("FDR q-val").head(10).to_string(index=False))


# 5.1. Plot top enriched pathways — Macrophages
def plot_gsea_top_pathways(gsea_results, cell_type, top_n=10):
    df = gsea_results.copy()

    # Convert GSEA statistics to numeric values
    df["NES"] = pd.to_numeric(df["NES"])
    df["FDR q-val"] = pd.to_numeric(df["FDR q-val"])

    # Select the most significant pathways
    top = (df.sort_values("FDR q-val").head(top_n).sort_values("NES"))

    # Positive enrichment = red, negative enrichment = blue
    colors = ["firebrick" if nes > 0 else "steelblue"for nes in top["NES"]]
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.barh(top["Term"],top["NES"],color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Normalized Enrichment Score (NES)")
    ax.set_ylabel("")
    ax.set_title(f"{cell_type} — Top {top_n} Hallmark Pathways")

    # Add NES values beside the bars
    for i, nes in enumerate(top["NES"]):
        ax.text(nes + (0.04 if nes > 0 else -0.04),i,f"{nes:.2f}", va="center", ha="left" if nes > 0 else "right")

    ax.margins(x=0.15)
    plt.tight_layout()
    output_path = FIGURES_DIR / f"{cell_type.lower()}_gsea_top_pathways.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved GSEA plot: {output_path}")

plot_gsea_top_pathways(macrophage_gsea.res2d,cell_type="Macrophage",top_n=10)
# ----------------------------------------------------------------------------------------------------------------------
# 6. Run preranked GSEA — AT2 cells
# ----------------------------------------------------------------------------------------------------------------------
at2_gsea = gp.prerank(
    rnk=at2_ranking,
    gene_sets=HALLMARK_LIBRARY,
    outdir=str(RESULTS_DIR / "gsea_at2"),
    min_size=10,
    max_size=500,
    permutation_num=1000,
    seed=42,
    threads=4, verbose=True)

print("\nAT2 GSEA completed.")
print("\nTop 10 enriched AT2 pathways:")


# 6.1. Plot top enriched pathways — AT2
plot_gsea_top_pathways(at2_gsea.res2d,cell_type="AT2",top_n=10)


# ----------------------------------------------------------------------------------------------------------------------
# 7. Prepare macrophage vs AT2 GSEA comparison
# ----------------------------------------------------------------------------------------------------------------------
def prepare_gsea_comparison(gsea_results, cell_type):
    df = gsea_results[["Term", "NES", "FDR q-val"]].copy()
    df["NES"] = pd.to_numeric(df["NES"])
    df["FDR q-val"] = pd.to_numeric(df["FDR q-val"])

    df = df.rename(columns={
        "NES": f"{cell_type}_NES",
        "FDR q-val": f"{cell_type}_FDR"})

    return df.set_index("Term")


macrophage_summary = prepare_gsea_comparison(macrophage_gsea.res2d, "Macrophage")
at2_summary = prepare_gsea_comparison(at2_gsea.res2d, "AT2")
comparison = macrophage_summary.join(at2_summary, how="outer")
print("\nGSEA comparison shape:", comparison.shape)
print("\nFirst five pathways:")
print(comparison.head())
print("\nMissing values:")
print(comparison.isna().sum())


# ----------------------------------------------------------------------------------------------------------------------
# 8. Plot comparative Hallmark GSEA heatmap
# ----------------------------------------------------------------------------------------------------------------------

# Extract NES and FDR matrices
nes_matrix = comparison[["Macrophage_NES", "AT2_NES"]].copy()
fdr_matrix = comparison[["Macrophage_FDR", "AT2_FDR"]].copy()

nes_matrix.columns = ["Macrophage", "AT2"]
fdr_matrix.columns = ["Macrophage", "AT2"]

# Order pathways by strongest absolute enrichment
pathway_order = (nes_matrix.abs().max(axis=1).sort_values(ascending=False).index)
nes_matrix = nes_matrix.loc[pathway_order]
fdr_matrix = fdr_matrix.loc[pathway_order]

# Create annotations with significance markers
annotations = pd.DataFrame("",index=nes_matrix.index, columns=nes_matrix.columns)

for pathway in nes_matrix.index:
    for cell_type in nes_matrix.columns:

        nes = nes_matrix.loc[pathway, cell_type]
        fdr = fdr_matrix.loc[pathway, cell_type]

        if pd.notna(nes):
            marker = "*" if pd.notna(fdr) and fdr < 0.05 else ""
            annotations.loc[pathway, cell_type] = f"{nes:.2f}{marker}"

# Create heatmap
fig, ax = plt.subplots(figsize=(10, 17))

sns.heatmap(
    nes_matrix,
    annot=annotations,
    fmt="",
    cmap="RdBu_r",
    center=0,
    vmin=-3,
    vmax=3,
    linewidths=0.5,
    linecolor="white",
    mask=nes_matrix.isna(),
    cbar_kws={"label": "Normalized Enrichment Score (NES)"},
    ax=ax
)

# Show missing values in gray
ax.set_facecolor("lightgray")
ax.set_title("Hallmark GSEA: Macrophages vs AT2 Cells\n""COVID-19 vs Control",fontsize=14,pad=15)
ax.set_xlabel("")
ax.set_ylabel("")
ax.tick_params(axis="x", rotation=0)
ax.tick_params(axis="y", labelsize=9)
fig.text(0.5, 0.01,"* FDR < 0.05 | Gray = pathway not tested",ha="center",fontsize=10)

plt.tight_layout(rect=[0, 0.025, 1, 1])
output_path = FIGURES_DIR / "macrophage_at2_gsea_heatmap.png"
fig.savefig(output_path, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"\nSaved comparative GSEA heatmap: {output_path}")


# ----------------------------------------------------------------------------------------------------------------------
# 9. Select pathways significant in at least one cell type
# ----------------------------------------------------------------------------------------------------------------------

significant_mask = ((comparison["Macrophage_FDR"] < 0.05) |(comparison["AT2_FDR"] < 0.05))
significant_comparison = comparison.loc[significant_mask].copy()
print("\nPathways significant in at least one cell type:",len(significant_comparison))
print("\nSignificant in both:", ((significant_comparison["Macrophage_FDR"] < 0.05) &(significant_comparison["AT2_FDR"] < 0.05)).sum())
print("\nSignificant only in macrophages:", ((significant_comparison["Macrophage_FDR"] < 0.05) & ~(significant_comparison["AT2_FDR"] < 0.05)).sum())
print("\nSignificant only in AT2:",(~(significant_comparison["Macrophage_FDR"] < 0.05) & (significant_comparison["AT2_FDR"] < 0.05)).sum())


# ----------------------------------------------------------------------------------------------------------------------
# 10. Focused GSEA heatmap — significant pathways
# ----------------------------------------------------------------------------------------------------------------------
# Extract NES and FDR values
sig_nes = significant_comparison[["Macrophage_NES", "AT2_NES"]].copy()
sig_fdr = significant_comparison[["Macrophage_FDR", "AT2_FDR"]].copy()
sig_nes.columns = ["Macrophage", "AT2"]
sig_fdr.columns = ["Macrophage", "AT2"]

# Sort pathways by maximum absolute NES
pathway_order = ( sig_nes.abs().max(axis=1).sort_values(ascending=False).index)
sig_nes = sig_nes.loc[pathway_order]
sig_fdr = sig_fdr.loc[pathway_order]

# Add NES values and significance markers
annotations = pd.DataFrame("",index=sig_nes.index,columns=sig_nes.columns)

for pathway in sig_nes.index:
    for cell_type in sig_nes.columns:

        nes = sig_nes.loc[pathway, cell_type]
        fdr = sig_fdr.loc[pathway, cell_type]

        if pd.notna(nes):
            marker = "*" if pd.notna(fdr) and fdr < 0.05 else ""
            annotations.loc[pathway, cell_type] = f"{nes:.2f}{marker}"

# Plot
fig, ax = plt.subplots(figsize=(10, 9))

sns.heatmap(sig_nes, annot=annotations, fmt="", cmap="RdBu_r", center=0,vmin=-3,vmax=3,  linewidths=0.5,
            linecolor="white", mask=sig_nes.isna(), cbar_kws={"label": "Normalized Enrichment Score (NES)"}, ax=ax)

ax.set_facecolor("lightgray")
ax.set_title("Significant Hallmark Pathways\n""Macrophages vs AT2 — COVID-19 vs Control",fontsize=14, pad=15)
ax.set_xlabel("")
ax.set_ylabel("")
ax.tick_params(axis="x", rotation=0)
ax.tick_params(axis="y", labelsize=10)
fig.text(0.5, 0.01,"* FDR < 0.05 | Positive NES: COVID-19 | Negative NES: Control",ha="center",fontsize=10)
plt.tight_layout(rect=[0, 0.03, 1, 1])
output_path = (FIGURES_DIR / "macrophage_at2_gsea_significant_heatmap.png")
fig.savefig(output_path, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved focused GSEA heatmap: {output_path}")


# ----------------------------------------------------------------------------------------------------------------------
# 11. Export GSEA results
# ----------------------------------------------------------------------------------------------------------------------

def export_gsea_results(gsea_results, cell_type):
    df = gsea_results.copy()

    # Convert statistical columns to numeric
    numeric_columns = ["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"]

    for col in numeric_columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Sort by statistical significance
    df = df.sort_values("FDR q-val")

    # Save complete GSEA results
    output_path = (RESULTS_DIR / f"{cell_type.lower()}_gsea_hallmark.csv")
    df.to_csv(output_path, index=False)
    print(f"Saved {cell_type} GSEA results: {output_path}")
    return df

macrophage_gsea_df = export_gsea_results(macrophage_gsea.res2d,"Macrophage")
at2_gsea_df = export_gsea_results(at2_gsea.res2d,"AT2")

# Save comparison table
comparison.to_csv(RESULTS_DIR / "macrophage_at2_gsea_comparison.csv")
print("\nGSEA results exported successfully.")

