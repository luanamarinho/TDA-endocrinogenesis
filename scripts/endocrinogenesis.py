import argparse
import logging
from pathlib import Path
import scanpy as sc
import scvelo as scv
from time import perf_counter
import umap
import numpy as np
import utils


def main():
    parser = argparse.ArgumentParser(description="Preprocess scVelo pancreas dataset and compute UMAP.")
    parser.add_argument("--min_genes", type=int, default=1200)
    parser.add_argument("--min_cells", type=int, default=20)
    parser.add_argument("--max_mt_perc", type=float, default=20)

    parser.add_argument("--exclude_highly_expressed", type=bool, default=True)
    parser.add_argument("--max_fraction", type=float, default=0.05)
    parser.add_argument("--target_sum", type=float, default=None)
    parser.add_argument("--key_added", type=str, default="norm")

    parser.add_argument("--n_top_genes", type=int, default=4000)
    parser.add_argument("--flavor", type=str, default="cell_ranger")


    parser.add_argument("--method", type=str, default="umap")
    parser.add_argument("--umap-metric", type=str, default="euclidean")
    parser.add_argument("--umap-random-state", type=int, default=42)
    parser.add_argument("--force_recompute", action="store_true", default=False)


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

    if not args.force_recompute:
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


