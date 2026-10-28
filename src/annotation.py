
# ============================================================
# CELL-TYPE ANNOTATION — GSE171524
# ============================================================
#
# PURPOSE:
# Identify marker genes for Leiden clusters to help
# determine their biological cell-type identities.
#
# INPUT:
# Integrated AnnData with raw UMI counts and Leiden labels.
#
# OUTPUT:
# Normalized AnnData for marker analysis and a marker table.
#
# IMPORTANT:
# Preserve the original raw-count AnnData.
# ============================================================

import scanpy as sc



# ======================================================================================================================
def prepare_marker_data(adata):
    """
    Create a normalized copy for marker-gene analysis.
    The original AnnData retains raw UMI counts.
    """

    # Work on a copy to preserve the original counts.
    adata_markers = adata.copy()

    # Normalize each nucleus to 10,000 total counts.
    sc.pp.normalize_total(adata_markers,target_sum=1e4,)

    # Apply log(1 + normalized expression).
    sc.pp.log1p(adata_markers)

    return adata_markers
# ======================================================================================================================


# ======================================================================================================================
def identify_cluster_markers(adata_markers):
    """
    Identify genes enriched in each Leiden cluster
    compared with all other nuclei.
    Uses the Wilcoxon rank-sum test.
    """

    sc.tl.rank_genes_groups(adata_markers, groupby="leiden", method="wilcoxon", reference="rest", use_raw=False,)

    # Convert Scanpy results into a pandas DataFrame.
    markers = sc.get.rank_genes_groups_df(adata_markers,group=None,)

    return markers
# ======================================================================================================================




# ======================================================================================================================
# Known lung cell-type markers from established biological literature.
# These are reference markers, NOT genes discovered by our analysis.
# ======================================================================================================================

KNOWN_CELL_MARKERS = {"AT2": ["SFTPC", "SFTPB", "ABCA3", "SFTPA1"], "AT1": ["AGER", "AQP5", "PDPN", "CAV1"],
                      "Macrophage": ["MARCO", "CD68", "MRC1"], "Monocyte": ["FCN1", "S100A8", "S100A9"],
                      "T cell": ["CD3D", "CD3E", "CD247"], "B cell": ["MS4A1", "CD79A", "CD79B"],
                      "Endothelial": ["PECAM1", "VWF", "EMCN"], "Fibroblast": ["DCN", "LUM", "COL1A1"],
                      "Ciliated epithelial": ["FOXJ1", "PIFO", "DNAH5"], "Mast cell": ["KIT", "TPSAB1", "CPA3"],
                      "Plasma cell": ["JCHAIN", "MZB1", "SDC1"], "Erythroid": ["HBB", "HBA1", "HBA2"],
                      "Mural / Pericyte": ["RGS5", "PDGFRB", "CSPG4", "NOTCH3"],
                      "Smooth muscle": ["ACTA2", "TAGLN", "MYH11", "CNN1"],
                      "Secretory epithelial": ["SCGB1A1", "SCGB3A1", "SCGB3A2", "BPIFB1"],
                      "Monocyte - inflammatory": ["FCN1", "VCAN", "LYZ", "S100A8", "S100A9"],
                      "Alveolar macrophage": ["PPARG", "FABP4", "MARCO", "MRC1"],
                      "Capillary endothelial": ["CA4", "RGCC", "EMCN", "EDNRB"],
                      "Venous endothelial": ["ACKR1", "SELE", "VCAM1"]}

# ======================================================================================================================


# ======================================================================================================================
def compare_markers_with_reference(markers, top_n=50):
    """
    Compare the top Wilcoxon marker genes for each Leiden cluster
    with known lung cell-type marker genes.

    Returns a DataFrame showing the number of matching markers
    for each cluster and candidate cell type.

    IMPORTANT:
    This is an exploratory matching method, not a definitive
    automatic cell-type annotation.
    """

    import pandas as pd

    results = []

    for cluster in sorted(markers["group"].unique(), key=lambda x: int(x)):

        # Select significant, positively enriched genes.
        cluster_markers = markers[(markers["group"] == cluster) & (markers["pvals_adj"] < 0.05) & (markers["logfoldchanges"] > 0)]

        # Take the highest-ranked genes.
        top_genes = set(cluster_markers.head(top_n)["names"])

        # Compare them against each known cell-type signature.
        for cell_type, reference_genes in KNOWN_CELL_MARKERS.items():
            matched_genes = sorted(top_genes.intersection(reference_genes))
            results.append({"cluster": cluster, "candidate_cell_type": cell_type,
                            "n_matching_markers": len(matched_genes), "matching_genes": ", ".join(matched_genes),})

    results_df = pd.DataFrame(results)

    return results_df.sort_values(["cluster", "n_matching_markers"], ascending=[True, False],
                                  key=lambda col: col.astype(int) if col.name == "cluster" else col,)
# ======================================================================================================================


# ======================================================================================================================
def plot_cell_type_markers(adata_markers, output_path):
    """
    Generate a dot plot of known lung cell-type markers
    across all Leiden clusters.

    Dot size: Percentage of nuclei expressing each gene.
    Dot color: Mean normalized expression, scaled per gene.

    Saves the figure without modifying the dataset.
    """

    from pathlib import Path
    import matplotlib.pyplot as plt

    # Include only marker genes present in our dataset.
    available_markers = {cell_type: [gene for gene in genes if gene in adata_markers.var_names]for cell_type, genes in KNOWN_CELL_MARKERS.items()}

    # Remove cell types without available marker genes.
    available_markers = {cell_type: genes for cell_type, genes in available_markers.items() if genes}

    if not available_markers:
        raise ValueError("No reference marker genes found in dataset.")

    # Generate the dot plot.
    sc.pl.dotplot(adata_markers, var_names=available_markers, groupby="leiden", standard_scale="var", show=False, save=False, figsize=(18, 9),)

    # Save the figure.
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight",)
    plt.close("all")
    print("Marker dot plot saved:", output_path)
# ======================================================================================================================



# ======================================================================================================================
def check_donor_representation(adata):
    """
    Calculate the contribution of each donor to every Leiden cluster.

    Helps identify clusters dominated by individual donors.

    Returns:
        donor_counts: Number of nuclei per donor and cluster.
        donor_percentages: Percentage contribution of each donor.
        summary: Largest donor contribution for each cluster.
    """

    import pandas as pd

    # Count nuclei for every cluster-donor combination.
    donor_counts = pd.crosstab( adata.obs["leiden"], adata.obs["donor_id"],)

    # Calculate each donor's percentage within a cluster.
    donor_percentages = donor_counts.div(donor_counts.sum(axis=1),axis=0,) * 100

    # Identify the donor contributing the most nuclei.
    summary = pd.DataFrame({
        "cluster": donor_counts.index.astype(str),
        "n_nuclei": donor_counts.sum(axis=1).values,
        "n_donors": (donor_counts > 0).sum(axis=1).values,
        "top_donor": donor_counts.idxmax(axis=1).values,
        "top_donor_percentage": donor_percentages.max(axis=1).values,})

    summary = summary.sort_values("cluster", key=lambda col: col.astype(int),)

    return donor_counts, donor_percentages, summary
# ======================================================================================================================



# ======================================================================================================================
def calculate_donor_enrichment(adata):
    """
    Calculate donor enrichment in each Leiden cluster.

    Enrichment = donor fraction within cluster /
                 donor fraction in the entire dataset.

    Interpretation:
        1.0 = expected representation
        >1  = donor is overrepresented
        <1  = donor is underrepresented
    """

    import pandas as pd

    # Count nuclei per cluster and donor.
    counts = pd.crosstab(adata.obs["leiden"], adata.obs["donor_id"],)

    # Fraction of each cluster contributed by each donor.
    cluster_fractions = counts.div(counts.sum(axis=1), axis=0,)

    # Overall fraction of nuclei contributed by each donor.
    overall_fractions = (adata.obs["donor_id"].value_counts(normalize=True).reindex(counts.columns))

    # Calculate enrichment relative to overall donor abundance.
    enrichment = cluster_fractions.div(overall_fractions, axis=1,)
    return enrichment
# ======================================================================================================================



# ======================================================================================================================
def assign_cell_types(adata):
    """
    Assign provisional biological cell-type labels to Leiden clusters.

    Labels are based on marker-gene analysis, reference markers,
    dot plots, and donor-representation checks.

    The original Leiden labels and raw UMI counts are preserved.
    """

    cell_type_map = {
        "0": "AT2",
        "1": "Secretory epithelial",
        "2": "Mural / smooth muscle",
        "3": "AT1 - putative",
        "4": "Fibroblast",
        "5": "Fibroblast-like stromal",
        "6": "Alveolar macrophage",
        "7": "Monocyte-like",
        "8": "Macrophage",
        "9": "Ciliated epithelial",
        "10": "Endothelial",
        "11": "Endothelial-like - unclassified",
        "12": "Plasma cell",
        "13": "B cell",
        "14": "Mast cell",
        "15": "T cell",
        "16": "Unclassified - donor-enriched",
        "17": "Erythroid",
    }

    # Check that every Leiden cluster has a defined annotation.
    clusters = set(adata.obs["leiden"].astype(str).unique())
    missing_clusters = clusters - set(cell_type_map)

    if missing_clusters:
        raise ValueError(f"Missing annotations for clusters: {missing_clusters}")

    # Work on a copy to preserve the original dataset.
    adata_annotated = adata.copy()

    # Map Leiden cluster numbers to biological names.
    adata_annotated.obs["cell_type"] = (adata_annotated.obs["leiden"].astype(str).map(cell_type_map))

    if adata_annotated.obs["cell_type"].isna().any():
        raise ValueError("Some nuclei have no cell-type annotation.")

    return adata_annotated
# ======================================================================================================================
