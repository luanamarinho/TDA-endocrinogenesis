import argparse
import logging
from pathlib import Path
import scanpy as sc
import scvelo as scv
from time import perf_counter
import utils
import preprocessing as prep
from anndata import AnnData
import logging
from typing import Dict, Any


logger = logging.getLogger(__name__)

def bastidas_pontes_pipeline(
    adata: AnnData,
    min_genes: int = 1200,
    min_cells: int = 20,
    max_mt_perc: float = 20,
    exclude_highly_expressed: bool = True,
    max_fraction: float = 0.05,
    target_sum: float = None,
    key_added: str = "norm_factor",
    n_top_genes: int = 4000,
    flavor: str = "cell_ranger",
    dr_method: str = "umap",
    umap_metric: str = "euclidean",
    umap_random_state: int = 42,
    **kwargs) -> Dict[str, Any]:

    logger.info("Starting cell and gene QC")
    prep.qc_filter(adata, min_genes=min_genes, min_cells=min_cells, max_mt_perc=max_mt_perc)

    if flavor in ("seurat_v3", "seurat_v3_paper"):
        logger.info(f"Identifying HVG using raw count data. Flavor:{flavor}")
        sc.pp.highly_variable_genes(adata, flavor = flavor, n_top_genes=n_top_genes)

    logger.info("Starting cell normalization")
    prep.normalize_counts(
        adata,
        exclude_highly_expressed = exclude_highly_expressed,
        max_fraction = max_fraction,
        target_sum = target_sum)
    
    logger.info("Starting log transformation of normalized data")
    sc.pp.log1p(adata)
    
    if flavor not in ("seurat_v3", "seurat_v3_paper"):
        logger.info(f"Identifying HVG using log-normalized data. Flavor:{flavor}")
        sc.pp.highly_variable_genes(adata, flavor = flavor, n_top_genes=n_top_genes)

    hvg = adata.var["highly_variable"]
    adata = adata[:, hvg].copy()

    embedding_collection = prep.dimensionality_reduction(data = adata.X, method = dr_method, metric = umap_metric, random_state = umap_random_state)

    adata_collection = {"logNormal":adata.X , "X_umap_original":adata.obsm["X_umap"]}
    
    adata_metadata = {"clusters_coarse": adata.obs["clusters_coarse"], "clusters": adata.obs["clusters"], "highly_variable_genes": hvg}

    
    return {"embedding_collection": embedding_collection, "adata_collection": adata_collection, "adata_metadata":adata_metadata}


def main():
    parser = argparse.ArgumentParser(description="Preprocess scVelo pancreas dataset and compute UMAP.")
    parser.add_argument("--min_genes", type=int, default=1200)
    parser.add_argument("--min_cells", type=int, default=20)
    parser.add_argument("--max_mt_perc", type=float, default=20)

    parser.add_argument("--exclude_highly_expressed", type=eval, default=True)
    parser.add_argument("--max_fraction", type=float, default=0.05)
    parser.add_argument("--target_sum", type=float, default=None)
    parser.add_argument("--key_added", type=str, default="norm")

    parser.add_argument("--n_top_genes", type=int, default=4000)
    parser.add_argument("--flavor", type=str, default="cell_ranger")


    parser.add_argument("--dr_method", type=str, default="umap")
    parser.add_argument("--umap_metric", type=str, default="euclidean")
    parser.add_argument("--umap_random_state", type=int, default=42)

    args = parser.parse_args()

    if args.dr_method != "umap":
        raise ValueError("Currently only UMAP is supported")    

    pipeline_params = vars(args)

    log_format="%(levelname)s:%(name)s:%(message)s"

    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        handlers=[
            logging.FileHandler("pipeline.log"),  # Saves to file
            logging.StreamHandler()              # Displays in terminal
        ]
    )

    root = Path(__file__).resolve().parents[1]
    cache_dir = utils._cache_file_folder(root, "data/Pancreas/cache")

    pipeline_cache_fname = utils.get_cache_id(pipeline_params)

    pipeline_cache_path = cache_dir / (
        f"pipeline_bundle_{pipeline_cache_fname}.joblib"
    )
    if pipeline_cache_path.exists():
        logger.info("Skipping pipeline and loading preprocessed cache: %s", pipeline_cache_path)
        pipeline = utils._load_pipeline(pipeline_cache_path)
    else:
      logger.info("No valid cached pipeline file was found. Loading pancreas data and running the bastisdas pipeline.")
      adata = scv.datasets.pancreas()
      adata = bastidas_pontes_pipeline(adata, **pipeline_params)
      logger.info("Caching new pipeline: %s", pipeline_cache_path)
      utils._save_pipeline(adata, pipeline_cache_path)


if __name__ == "__main__":
    main()


