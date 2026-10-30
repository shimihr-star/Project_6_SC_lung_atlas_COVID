import numpy as np
import pandas as pd
from scipy import sparse
from pathlib import Path

# ======================================================================================================================
def validate_raw_counts(adata):
    """
    Check whether adata.X contains nonnegative integer counts.
    """

    X = adata.X
    values = X.data if sparse.issparse(X) else np.asarray(X).ravel()
    is_nonnegative = np.all(values >= 0)
    is_integer = np.all(np.isclose(values, np.round(values)))
    print("Nonnegative counts:", is_nonnegative)
    print("Integer-valued counts:", is_integer)
    print("Maximum count:", values.max())

    if not is_nonnegative or not is_integer:
        raise ValueError("adata.X does not appear to contain raw counts.")

# ======================================================================================================================



# ======================================================================================================================
def count_nuclei_per_donor(adata, cell_types):
    """
    Count nuclei for each donor and selected cell type.
    """

    subset = adata.obs[adata.obs["cell_type"].isin(cell_types)]
    counts = (subset.groupby(["cell_type", "donor_id", "group"],observed=True,).size().reset_index(name="n_nuclei"))

    return counts

# ======================================================================================================================



# ======================================================================================================================

def aggregate_pseudobulk(adata, cell_types):
    """
    Aggregate raw UMI counts by donor and cell type.

    Returns:
        pseudobulk: dictionary of DataFrames
            Keys: cell-type names
            Rows: donors
            Columns: genes
            Values: summed raw counts
    """

    import pandas as pd
    import numpy as np
    from scipy import sparse

    pseudobulk = {}

    for cell_type in cell_types:
        # Select nuclei belonging to this cell type.
        cell_adata = adata[adata.obs["cell_type"] == cell_type]
        donors = cell_adata.obs["donor_id"].astype(str)
        donor_ids = sorted(donors.unique())
        aggregated_counts = []

        for donor_id in donor_ids:

            # Select nuclei from one donor.
            donor_mask = (donors == donor_id).to_numpy()
            donor_X = cell_adata.X[donor_mask]

            # Sum raw counts across nuclei for each gene.
            gene_counts = donor_X.sum(axis=0)

            # Convert sparse matrix result to 1D array.
            gene_counts = np.asarray(gene_counts).ravel()
            aggregated_counts.append(gene_counts)

        # Construct donor-by-gene count matrix.
        pseudobulk[cell_type] = pd.DataFrame(np.vstack(aggregated_counts).astype(np.int64), index=donor_ids,
                                             columns=adata.var_names,)

    return pseudobulk
# ======================================================================================================================

# ======================================================================================================================
def calculate_pseudobulk_qc(matrix, nuclei_counts, cell_type):
    """
    Calculate sample-level QC metrics for one cell type.

    Parameters
    ----------
    matrix : pd.DataFrame
        Pseudobulk count matrix (donors x genes).

    nuclei_counts : pd.DataFrame
        Output of count_nuclei_per_donor().

    cell_type : str
        Cell type being analyzed.

    Returns
    -------
    pd.DataFrame
        One row per donor with:
        - donor_id
        - group
        - n_nuclei
        - library_size
        - detected_genes
        - counts_per_nucleus
    """

    # Select nucleus counts for the current cell type.
    donor_info = nuclei_counts.loc[nuclei_counts["cell_type"] == cell_type,["donor_id", "group", "n_nuclei"]].copy()

    # Calculate expression-based QC metrics.
    qc = pd.DataFrame({"donor_id": matrix.index.astype(str),"library_size": matrix.sum(axis=1).to_numpy(), "detected_genes": (matrix > 0).sum(axis=1).to_numpy()})

    # Add donor metadata and nucleus counts.
    qc = qc.merge(donor_info, on="donor_id", how="left", validate="one_to_one")

    # Ensure every pseudobulk sample has donor metadata.
    if qc[["group", "n_nuclei"]].isna().any().any():
        raise ValueError(f"Missing donor metadata for {cell_type}")

    # Calculate average UMI counts per nucleus.
    qc["counts_per_nucleus"] = (qc["library_size"] / qc["n_nuclei"])
    return qc
# ======================================================================================================================


# ======================================================================================================================
def plot_pseudobulk_qc(qc, cell_type, output_dir):
    """
    Generate diagnostic plots for pseudobulk samples.

    Plots:
    1. Number of nuclei versus library size.
    2. Number of nuclei versus detected genes.
    3. Library size versus detected genes.

    Samples are colored by experimental group.
    """

    import matplotlib.pyplot as plt
    import seaborn as sns

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    metrics = [ ("n_nuclei", "library_size"),  ("n_nuclei", "detected_genes"), ("library_size", "detected_genes")]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for ax, (x, y) in zip(axes, metrics):

        sns.scatterplot(data=qc, x=x, y=y, hue="group", ax=ax, s=70)

        # Log scales accommodate large differences in sample size.
        if x in ["n_nuclei", "library_size"]:
            ax.set_xscale("log")

        if y == "library_size":
            ax.set_yscale("log")

        ax.set_title(f"{x} vs {y}")

    fig.suptitle(f"{cell_type} pseudobulk sample QC")
    fig.tight_layout()
    filename = cell_type.lower().replace(" ", "_")
    output_path = output_dir / f"{filename}_pseudobulk_qc.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved QC figure:", output_path)
# ======================================================================================================================



# ======================================================================================================================
def plot_pseudobulk_pca(matrix, qc, cell_type, output_dir):
    """
    Perform exploratory PCA on pseudobulk gene expression.

    Steps:
    1. Filter genes with very low expression.
    2. Calculate CPM and log2-transform.
    3. Select the 2,000 most variable genes.
    4. Center genes and perform PCA.
    5. Plot donors colored by experimental group.

    This normalization is for exploratory PCA only,
    not for final differential expression analysis.
    """

    from pathlib import Path
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.decomposition import PCA

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # 1. Filter very low-expression genes
    # --------------------------------------------------------
    # Retain genes with at least 10 counts in at least
    # 3 donor samples. This is preliminary PCA filtering.
    keep_genes = (matrix >= 10).sum(axis=0) >= 3
    filtered = matrix.loc[:, keep_genes]

    if filtered.shape[1] < 2:
        raise ValueError(f"Too few genes for PCA: {cell_type}")

    print(f"\n{cell_type} PCA:")
    print("Genes before filtering:", matrix.shape[1])
    print("Genes after filtering:", filtered.shape[1])

    # --------------------------------------------------------
    # 2. Calculate log2(CPM + 1)
    # --------------------------------------------------------
    library_sizes = matrix.sum(axis=1)

    if (library_sizes <= 0).any():
        raise ValueError("Cannot normalize samples with zero counts.")

    cpm = filtered.div(library_sizes, axis=0) * 1e6
    log_cpm = np.log2(cpm + 1)

    # --------------------------------------------------------
    # 3. Select most variable genes
    # --------------------------------------------------------
    gene_variances = log_cpm.var(axis=0)
    top_genes = gene_variances.nlargest(min(2000, len(gene_variances))).index
    pca_input = log_cpm[top_genes]

    # --------------------------------------------------------
    # 4. Perform PCA
    # --------------------------------------------------------
    pca = PCA(n_components=2)
    coordinates = pca.fit_transform(pca_input)
    pca_df = pd.DataFrame(coordinates,columns=["PC1", "PC2"],index=matrix.index)
    pca_df["donor_id"] = pca_df.index

    # Add donor group and sample-level QC metrics.
    pca_df = pca_df.merge(qc[["donor_id", "group", "n_nuclei", "library_size"]],on="donor_id",how="left",validate="one_to_one")
    if pca_df["group"].isna().any():
        raise ValueError("Missing group labels in PCA results.")

    # --------------------------------------------------------
    # 5. Generate PCA plot
    # --------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 7))

    sns.scatterplot(data=pca_df, x="PC1", y="PC2",hue="group",s=100,ax=ax)

    # Label each donor so potential outliers are identifiable.
    for _, row in pca_df.iterrows():
        ax.annotate(row["donor_id"], (row["PC1"], row["PC2"]), xytext=(4, 4), textcoords="offset points", fontsize=8)

    explained = pca.explained_variance_ratio_ * 100
    ax.set_xlabel(f"PC1 ({explained[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({explained[1]:.1f}% variance)")
    ax.set_title(f"{cell_type} pseudobulk PCA")

    fig.tight_layout()
    safe_name = cell_type.lower().replace(" ", "_")
    figure_path = output_dir / f"{safe_name}_pseudobulk_pca.png"
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Save PCA coordinates for later investigation.
    coordinates_path = ( output_dir.parent / "data" / "results" / f"{safe_name}_pseudobulk_pca.csv")
    coordinates_path.parent.mkdir(parents=True, exist_ok=True)
    pca_df.to_csv(coordinates_path, index=False)
    print("Explained variance:", explained)
    print("Saved PCA figure:", figure_path)
    print("Saved PCA coordinates:", coordinates_path)

    return pca_df
# ======================================================================================================================
