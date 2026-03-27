from anndata import AnnData
import scanpy as sc
import umap
import numpy as np
from typing import Dict, Any
 

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
    
    sc.pp.filter_cells(adata, min_genes=min_genes, inplace = True)
    sc.pp.filter_genes(adata, min_cells=min_cells, inplace = True)

    mt_genes = adata.var_names.str.upper().str.startswith("MT-")
    mt_frac = 100 * (np.asarray(adata[:, mt_genes].X.sum(axis=1)) / adata.X.sum(axis=1))

    adata = adata[mt_frac < max_mt_perc]
    
    return adata
 


def normalize_counts(adata: AnnData, layer:str = None,**kwargs) -> AnnData:
    """
    Wrapper for sequencing depth normalization using ``scanpy.pp.normalize_total``.
    Additional keyword arguments are passed directly to the scanpy function.

    Each cell's counts are scaled so that the total count per cell equals a target
    sum (default: median total count across cells), making cells comparable
    while correcting for differences in sequencing depth. Normalization is applied in-place.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix of shape (n_obs, n_vars), where rows correspond to
        cells and columns to genes. ``adata.X`` must contain raw UMI counts prior
        to calling this function; normalized counts will overwrite ``adata.X``
        in-place.
    
    **kwargs
        Additional keyword arguments forwarded directly to
        ``sc.pp.normalize_total``. Commonly used options include:

        - ``target_sum`` (*float*, default ``None``): total count to normalize
          each cell to. If ``None``, the median total count across cells is used.
        - ``exclude_highly_expressed`` (*bool*, default ``False``): exclude highly
          expressed genes from the normalization factor computation, to avoid
          a single dominant gene skewing the scaling. See (einreb et al., 2017.
        - ``max_fraction`` (*float*, default ``0.05``): if
          ``exclude_highly_expressed=True``, genes accounting for more than this
          fraction of total counts in any cell are excluded from the scaling
          factor.
        

        See the `scanpy docs
        <https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html>`_
        for the full list of options.

    Returns
    -------
    AnnData
        The input ``adata`` object with ``adata.X`` (or the specified layer)
        overwritten with depth-normalized counts. Modification is in-place;
        the returned object is the same as the input.
    """

    ## TODO: add more normalization methods
    sc.pp.normalize_total(adata, **kwargs)
    return adata

def dimensionality_reduction(
    data: np.ndarray,
    method: str = "umap",
    metric: str = "euclidean",
    random_state: int = 42,
    **kwargs) -> Dict[str, Any]:
    """
    Apply dimensionality reduction to an AnnData object. Currently only ``umap`` from ``umap-learn`` is supported.

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
        Additional arguments passed to the dimensionality reduction method.

    Returns
    -------
    Dict[str, Any]
        Keys:

        - ``"embedding"`` (*np.ndarray*, shape (n_samples, n_components)):
          the low-dimensional embedding coordinates.
        - ``"params"`` (*dict*): hyperparameters of the fitted reducer, as
          returned by ``get_params()``.
    """
    # TODO: add more DR methods
    if method != "umap":
        raise ValueError(f"Dimensionality reduction method: {method!r} is currently not supported")

    mapper = umap.UMAP(
        metric=metric,
        random_state=random_state,
        **kwargs
    )

    return {"embedding": mapper.fit_transform(data), "params": mapper.get_params()}


def select_hvg(adata: AnnData, flavor: str = "cell_ranger", n_top_genes: int = 4000) -> AnnData:
    if flavor in ("seurat_v3", "seurat_v3_paper"):
        sc.pp.highly_variable_genes(adata, layer = "raw", n_top_genes=n_top_genes)
    else:
        sc.pp.highly_variable_genes(adata, layer = "logNormal", n_top_genes=n_top_genes)
    return adata
