import scanpy as sc
from pathlib import Path
import numpy as np
import json
import hashlib
from joblib import load, dump


def _cache_file_folder(root: Path, cache_subdir: str):
    cache_dir = root / cache_subdir
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def _save_pipeline(pipeline_output, cache_path: Path):
    dump(pipeline_output, cache_path, compress=3)


def _load_pipeline(cache_path: Path):
    pipeline_bundle = load(cache_path)
    return pipeline_bundle


def get_cache_id(params: dict) -> str:
    """Creates a stable 12-character hash from a dictionary of parameters."""
    serializable = {
        k: v for k, v in params.items() 
        if isinstance(v, (str, int, float, bool, type(None)))
    }
    
    param_string = json.dumps(serializable, sort_keys=True, separators=(',', ':'))
    return hashlib.sha256(param_string.encode()).hexdigest()[:12]
