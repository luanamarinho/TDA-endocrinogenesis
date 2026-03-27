from anndata import AnnData
import scanpy as sc
import umap
 
def qc_filter(
        adata: AnnData, 
        min_genes: int = 1200, 
        min_cells: int = 20,
        max_mt_perc: float = 20) -> AnnData:
    """
    Apply sequential QC filters to remove low-quality cells and lowly-expressed genes
    from a raw count matrix.

    Parameters
    ----------
    adata : AnnData
        The annotaded cell x gene matrix with UMI counts. `adata.X` and `adata.layers` will be modified in place.
    min_genes : int
        Minimum number of genes detected per cell (default: 1200).
    min_cells : int
        Minimum number of cells a gene must be detected in (default: 20).
    max_mt_perc : float
        Maximum percentage of mitochondrial counts per cell. Cells at or above
        this threshold are considered stressed or dying (default: 20).

    Returns
    -------
    AnnData
        Filtered AnnData object.
    """
    
    sc.pp.filter_cells(adata, min_genes=min_genes, inplace = True)
    sc.pp.filter_genes(adata, min_cells=min_cells, inplace = True)

    mt_genes = adata.var_names.str.upper().str.startswith("MT-")
    mt_frac = 100 * (np.asarray(adata[:, mt_genes].X.sum(axis=1)) / adata.X.sum(axis=1))
    adata = adata[mt_frac < max_mt_perc]
    
    return adata
 

def normalize_counts(adata: AnnData, **kwargs) -> AnnData:
    """
    Wrapper for sequencing depth normalization using `scanpy.pp.normalize_total`.
    Additional keyword arguments are passed directly to `scanpy.pp.normalize_total`.

    Parameters
    ----------
    adata : AnnData
        The annotaded cell x gene matrix with UMI counts. `adata.X` and `adata.layers` will be modified in place.
    **kwargs: optional
        Additional arguments forwarded to `sc.pp.normalize_total` (e.g. `target_sum`, `exclude_highly_expressed`, `max_fraction`).

    Returns
    -------
    AnnData
        Normalized AnnData object.
    """
    # TODO: add more normalization methods
    sc.pp.normalize_total(adata, **kwargs)
    return adata

def dimensionality_reduction(
    adata: AnnData,
    method: str = "umap",
    metric: str = "euclidean",
    random_state: int = 42,
    **kwargs) -> umap.UMAP:
    """
    Apply dimensionality reduction to an AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix. `adata.X` is used as input.
    method : str
        Dimensionality reduction method (default: "umap").
    metric : str
        Distance metric used by the reducer (default: "euclidean").
    random_state : int
        Random seed for reproducibility (default: 42).
    **kwargs
        Additional arguments forwarded to the dimensionality reduction method.

    Returns
    -------
    umap.UMAP
        Fitted UMAP mapper object.
    """
    # TODO: add more DR methods
    if method.lower() != "umap":
        raise ValueError(f"Dimensionality reduction method: {method!r} is currently not supported")

    mapper = umap.UMAP(
        metric=metric,
        random_state=random_state,
        **kwargs
    )
    return mapper.fit(adata.X)
