import scanpy as sc
from pathlib import Path
import numpy as np
import json


def _cache_file_folder(root: Path, cache_subdir: str):
    cache_dir = root / cache_subdir
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def _save_preprocessed_cache(adata, cache_path: Path, preprocess_params: dict):
    adata.uns["preprocess_params"] = preprocess_params
    adata.write_h5ad(cache_path, compression="gzip")


def _load_adata_cache(cache_path: Path):
    adata = sc.read_h5ad(cache_path)
    return adata


def _preprocess_params_match(adata, passed_params: dict):
    current = adata.uns["preprocess_params"]
    return current == passed_params


def _save_umap_cache(cache_path: Path, umap_collection: dict):
    np.savez_compressed(
        cache_path,
        map=umap_collection["map"],
        params_json=json.dumps(umap_collection["umap_param"], sort_keys=True)
    )


def _load_umap_cache(cache_path: Path):
    collection = np.load(cache_path)
    params_json_raw = collection["params_json"]
    params_json = params_json_raw.item() if hasattr(params_json_raw, "item") else str(params_json_raw)
    params = json.loads(params_json) if params_json else {}
    return {
        "map": collection["map"],
        "params": params,
    }

def _umap_params_match(cached_params: dict, passed_params: dict):
    if not (cached_params and passed_params):
        return False
    return all(cached_params.get(k) == v for k, v in passed_params.items())

def build_cache_filename(args) -> str:
    parts = [
        f"mg{args.min_genes}",
        f"mc{args.min_cells}",
        f"mt{args.max_mt_perc}",
        f"ehe{int(args.exclude_highly_expressed)}",
        f"mf{args.max_fraction}",
        f"ts{args.target_sum}",
        f"ntg{args.n_top_genes}",
        f"fl{args.flavor}",
    ]
    return "pancreas_" + "_".join(parts) + ".h5ad"