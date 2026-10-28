
# ============================================================
# scVI INTEGRATION AND CLUSTERING
# ============================================================
#
# PURPOSE:
# Learn a low-dimensional representation of snRNA-seq data,
# construct a nearest-neighbor graph, and identify clusters.
#
# INPUT:
# Preprocessed AnnData with raw UMI counts and HVG annotations.
#
# OUTPUT:
# AnnData with scVI latent coordinates, UMAP and Leiden labels.
#
# NOTE:
# Keep raw counts in the original AnnData object.
# ============================================================



import scanpy as sc
import scvi

# ======================================================================================================================
def train_scvi(adata, n_latent=10):
    """
    Train scVI using highly variable genes.

    Returns:
        Trained scVI model and latent representation.
    """

    # Select HVGs without modifying the full dataset.
    adata_scvi = adata[:, adata.var["highly_variable"]].copy()
    print("Full dataset:", adata.shape)
    print("scVI dataset:", adata_scvi.shape)

    # Use sample identity as the batch covariate.
    scvi.model.SCVI.setup_anndata(adata_scvi, bbatch_key="biosample_id",)

    # Train the model on raw UMI counts.
    model = scvi.model.SCVI(adata_scvi,n_latent=n_latent,)
    model.train()

    # Extract the learned representation of each nucleus.
    latent = model.get_latent_representation()
    print("Latent representation:", latent.shape)
    return model, latent
# ======================================================================================================================


# ======================================================================================================================
def cluster_nuclei(adata, latent, n_neighbors=15, leiden_resolution=0.5,):
    """
    Construct a neighbor graph, calculate UMAP,
    and perform Leiden clustering.
    """

    # Store scVI latent coordinates in the full dataset.
    adata.obsm["X_scVI"] = latent

    # Build a graph using the scVI representation.
    sc.pp.neighbors(adata,use_rep="X_scVI", n_neighbors=n_neighbors,)

    # Calculate two-dimensional UMAP coordinates.
    sc.tl.umap(adata)

    # Identify clusters from the neighbor graph.
    sc.tl.leiden(adata, resolution=leiden_resolution, key_added="leiden", flavor="igraph", n_iterations=2, directed=False,)
    print("Number of Leiden clusters:",adata.obs["leiden"].nunique(),)
    return adata
# ======================================================================================================================


# ======================================================================================================================
def validate_integration(adata):
    """
    Check that integration and clustering produced
    the expected data structures.
    """

    assert "X_scVI" in adata.obsm
    assert "X_umap" in adata.obsm
    assert "leiden" in adata.obs
    assert adata.obsm["X_scVI"].shape == (adata.n_obs, 10)
    assert adata.obsm["X_umap"].shape == (adata.n_obs, 2)
    assert adata.obs["leiden"].notna().all()
    print("Integration validation passed.")
    print("Nuclei:", adata.n_obs)
    print("Clusters:", adata.obs["leiden"].nunique())
# ======================================================================================================================