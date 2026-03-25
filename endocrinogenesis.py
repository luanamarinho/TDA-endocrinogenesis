from logging import exception, log
import scanpy as sc
import scvelo as scv
from time import perf_counter, time


"""
Preprocess the pancreas dataset following the scVelo workflow.

- Filter genes by a minimum number of shared counts, `min_shared_counts`, then retain the top `n_top_genes`
  most variable genes.
- “Shared counts” refer to signal present in both `spliced` and `unspliced` layers.
  The shared-count matrix is defined entry-wise as:
    X_ij = spliced_ij + unspliced_ij  if spliced_ij > 0 and unspliced_ij > 0,
    X_ij = 0                         otherwise.
- Cells are normalized to the median library size.
- Normalized and filtered data is `log1p` transformed.
"""

def pre_treatment(adata, min_shared_counts=20, n_top_genes=2000):
    try:
        if adata:
            log.info("Initiating pre-treatment of the raw counts matrix of shape:", adata.shape)
            t0 = perf_counter()
            scv.pp.filter_genes(adata, min_shared_counts=min_shared_counts)
            scv.pp.normalize_per_cell(adata)
            sc.pp.filter_genes_dispersion(adata, n_top_genes=n_top_genes)
            sc.pp.log1p(adata)
            dt = perf_counter() - t0
            log.info(f"Data was succesfully pre-treated after {dt} sec. Final shape: {}")
            log.info("Dimensions of the pancreas dataset after filtering and normalizing: ", adata.shape)
    except:
        exception





        


def dimensionality_reduction (adata, perform_DR = False, method = "umap"):
    if adataif perform_DR:
        if method.lower() == "umap":



adata = scv.datasets.pancreas()

pre_treatment(adata)

