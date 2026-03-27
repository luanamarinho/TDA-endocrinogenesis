import argparse
import logging
from pathlib import Path
from anndata import AnnData
import scanpy as sc
import scvelo as scv
from time import perf_counter
import umap
import numpy as np
import utils

 
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
        Raw count matrix. Will be modified in place.
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
 

def dimensionality_reduction(
    data,
    perform_DR = True,
    method = "umap",
    metric = "cosine",
    random_state = 42,
    low_memory = True,
):
    if data is None:
        raise ValueError("`data` must not be None.")
    if not perform_DR:
        return None
    if method.lower() != "umap":
        raise ValueError(f"Currently unsupported dimensionality reduction method: {method!r}")

    logger.info("Starting dimensionality reduction. Method=%s", method)
    t0 = perf_counter()

    mapper = umap.UMAP(
        metric=metric,
        random_state=random_state,
        low_memory=low_memory
    )
    embedding = mapper.fit_transform(data)
    dt = perf_counter() - t0

    logger.info(
        "Finished dimensionality reduction in %.2fs. Final embedding shape: %s",
        dt,
        embedding.shape,
    )
    return {
        "map": embedding,
        "umap_param": mapper.get_params()
    }
