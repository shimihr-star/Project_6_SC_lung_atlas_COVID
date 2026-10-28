import numpy as np
import pandas as pd
from scipy import sparse


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