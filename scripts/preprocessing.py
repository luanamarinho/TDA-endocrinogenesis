from anndata import AnnData
import scanpy as sc
import umap
import numpy as np
 

def qc_filter(
    adata: AnnData,
    min_genes: int = 1200,
    min_cells: int = 20,
    max_mt_perc: float = 20
) -> AnnData:
    
    """
    Apply sequential QC filters to remove low-quality cells and lowly-expressed genes
    from raw count data.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix of shape (n_obs, n_vars), where rows correspond to
        cells and columns to genes. Expected to contain raw counts in ``adata.X``.
    min_genes : int
        Minimum number of genes that must be detected (count > 0) in a cell for it
        to be retained. Cells below this threshold are typically empty droplets or
        low-quality cells. Default is 1200.
    min_cells : int
        Minimum number of cells in which a gene must be detected for it to be
        retained. Genes below this threshold are considered lowly expressed and
        uninformative. Default is 20.
    max_mt_perc : float
        Maximum allowable mitochondrial gene count as a percentage of total counts
        per cell. Cells exceeding this threshold are likely damaged or dying.
        Default is 20 (%).

    Returns
    -------
    AnnData
        Filtered subset of the annotated data matrix containing only cells and genes that passed
        all QC thresholds.
    """
    
    sc.pp.filter_cells(adata, min_genes=min_genes, inplace = True, **kwargs)
    sc.pp.filter_genes(adata, min_cells=min_cells, inplace = True, **kwargs)

    mt_genes = adata.var_names.str.upper().str.startswith("MT-")
    mt_frac = 100 * (np.asarray(adata[:, mt_genes].X.sum(axis=1)) / adata.X.sum(axis=1))

    adata = adata[mt_frac < max_mt_perc]
    
    return adata
 


def normalize_counts(adata: AnnData, **kwargs) -> AnnData:
    """
    Wrapper for sequencing depth normalization using ``scanpy.pp.normalize_total``.
    Additional keyword arguments are passed directly to the scanpy function.

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
    ## TODO: add more normalization methods
    sc.pp.normalize_total(adata, **kwargs)
    return adata

def dimensionality_reduction(
    data: np.ndarray,
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
    if method != "umap":
        raise ValueError(f"Dimensionality reduction method: {method!r} is currently not supported")

    mapper = umap.UMAP(
        metric=metric,
        random_state=random_state,
        **kwargs
    )
    return mapper.fit(data)
