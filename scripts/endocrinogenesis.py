import logging
import scanpy as sc
import scvelo as scv
from time import perf_counter
import umap
import numpy as np


"""
Preprocess the pancreas dataset following the scVelo workflow.

- Filter genes by a minimum number of shared counts (`min_shared_counts`), then retain
  the top `n_top_genes` most variable genes.
- “Shared counts” refer to signal present in both `spliced` and `unspliced` layers.
  The shared-count matrix is defined entry-wise as:
    X_ij = spliced_ij + unspliced_ij  if spliced_ij > 0 and unspliced_ij > 0,
    X_ij = 0                         otherwise.
- Normalize each cell by library size (total-count scaling) to the median pre-normalization
  total count per cell.
- Apply a `log1p` transform to the normalized data.
"""

logger = logging.getLogger(__name__)


def pre_treatment(adata, min_shared_counts = 20, n_top_genes = 2000):
    if adata is None:
        raise ValueError("`adata` must not be None.")

    logger.info("Starting preprocessing. Initial shape: %s", adata.shape)
    t0 = perf_counter()

    scv.pp.filter_genes(adata, min_shared_counts=min_shared_counts)
    scv.pp.normalize_per_cell(adata)
    sc.pp.filter_genes_dispersion(adata, n_top_genes=n_top_genes)
    sc.pp.log1p(adata)

    dt = perf_counter() - t0
    logger.info("Finished preprocessing in %.2fs. Final shape: %s", dt, adata.shape)
    return adata


def dimensionality_reduction(data, perform_DR = True, method = "umap"):
  if data is None:
    raise ValueError("`adata` must not be None.")
  if perform_DR:
    if method.lower() == "umap":
      logger.info(f"Starting dimensionality reduction. Method:")
      t0 = perf_counter()
      
      mapper = umap.UMAP(metric="cosine", random_state=42, low_memory=True)
      map = mapper.fit_transform(data)
      dt = perf_counter() - t0
      logger.info("Finished dimensionality reduction in %.2fs. Final shape: %s", dt, map.shape)
      return {"map":map, "n_neighbors":mapper.n_neighbors, "metric":mapper.metric, "learning_rate":mapper.learning_rate, "local_connectivity":mapper.local_connectivity}


def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")
    adata = scv.datasets.pancreas()
    pre_treatment(adata)
    map = dimensionality_reduction(adata.X)
    np.savez(
      "map_collection.npz",
      map=map["map"],
      params={"n_neighbors": map["n_neighbors"], "metric": map["metric"], "learning_rate": map["learning_rate"], "local_connectivity":map["local_connectivity"]},
    )    



if __name__ == "__main__":
    main()


