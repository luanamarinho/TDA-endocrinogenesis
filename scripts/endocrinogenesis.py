import argparse
import logging
from pathlib import Path
import scanpy as sc
import scvelo as scv
from time import perf_counter
import umap
import numpy as np
import utils


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


def pre_treatment(adata, min_shared_counts = 20, n_top_genes = 2000, keep_raw = False):
    # TODO: include QC
    
    if adata is None:
        raise ValueError("`adata` must not be None.")

    if keep_raw:
      adata.raw = adata.X.copy()

    logger.info("Starting preprocessing. Initial shape: %s", adata.shape)
    t0 = perf_counter()

    scv.pp.filter_genes(adata, min_shared_counts=min_shared_counts)
    scv.pp.normalize_per_cell(adata)
    sc.pp.filter_genes_dispersion(adata, n_top_genes=n_top_genes)
    sc.pp.log1p(adata)

    dt = perf_counter() - t0
    logger.info("Finished preprocessing in %.2fs. Final shape: %s", dt, adata.shape)

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


def main():
    parser = argparse.ArgumentParser(description="Preprocess scVelo pancreas dataset and compute UMAP.")
    parser.add_argument("--min-shared-counts", type=int, default=20)
    parser.add_argument("--n-top-genes", type=int, default=2000)

    parser.add_argument("--method", type=str, default="umap")
    parser.add_argument("--umap-metric", type=str, default="euclidean")
    parser.add_argument("--umap-random-state", type=int, default=42)


    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")

    root = Path(__file__).resolve().parents[1]
    cache_dir = utils._cache_file_folder(root, "data/Pancreas/cache")

    umap_params = {
        "metric": args.umap_metric,
        "random_state": args.umap_random_state
    }

    umap_cache_path = cache_dir / (
        f"umap_metric{args.umap_metric}_rs{args.umap_random_state}.npz"
    )

    preprocess_params = {
        "min_shared_counts": args.min_shared_counts,
        "n_top_genes": args.n_top_genes,
    }

    preprocessed_cache_path = cache_dir / (
        f"pancreas_preprocessed_minshared{args.min_shared_counts}_top{args.n_top_genes}.h5ad"
    )

    if preprocessed_cache_path.exists():
        logger.info("Loading preprocessed cache: %s", preprocessed_cache_path)
        adata = utils._load_preprocessed_cache(preprocessed_cache_path)
        if not utils._preprocess_params_match(adata, preprocess_params):
            logger.info("Preprocessed cache params mismatch. Recomputing preprocessing...")
            adata = scv.datasets.pancreas()
            adata = pre_treatment(adata, min_shared_counts=args.min_shared_counts, n_top_genes=args.n_top_genes)
            logger.info("Overwriting preprocessed cache (gzip): %s", preprocessed_cache_path)
            utils._save_preprocessed_cache(adata, preprocessed_cache_path, preprocess_params=preprocess_params)
    else:
      logger.info("No preprocessed cached filed was found. Loading data and initializing pretreatment.")
      adata = scv.datasets.pancreas()
      adata = pre_treatment(adata, min_shared_counts=args.min_shared_counts, n_top_genes=args.n_top_genes)
      logger.info("Saving new preprocessed cache (gzip): %s", preprocessed_cache_path)
      utils._save_preprocessed_cache(adata, preprocessed_cache_path, preprocess_params=preprocess_params)


    if umap_cache_path.exists():
        logger.info("Loading cached UMAP collection: %s", umap_cache_path)
        umap_collection = utils._load_umap_cache(umap_cache_path)
        if not utils._umap_params_match(umap_collection["params"], umap_params):
            logger.info("UMAP cache params mismatch. Recomputing UMAP...")
            umap_collection = dimensionality_reduction(
                adata.X,
                method=args.method,
                metric=args.umap_metric,
                random_state=args.umap_random_state,
                low_memory=args.low_memory
            )
            logger.info("Overwriting UMAP collection cache: %s", umap_cache_path)
            utils._save_umap_cache(umap_cache_path, umap_collection)
    else:
        logger.info("No compatible UMAP cached file was found. Computing UMAP...")
        umap_collection = dimensionality_reduction(
            adata.X,
            method=args.method,
            metric=args.umap_metric,
            random_state=args.umap_random_state
        )
        logger.info("Caching new UMAP collection: %s", umap_cache_path)
        utils._save_umap_cache(umap_cache_path, umap_collection)

    


if __name__ == "__main__":
    main()


