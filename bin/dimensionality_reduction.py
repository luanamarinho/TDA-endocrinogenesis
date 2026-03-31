#!/usr/bin/env python3
"""
Apply UMAP dimensionality reduction to an HVG-filtered AnnData object.
Outputs:
  - embedding.npy  : the low-dimensional embedding array (n_cells × n_components)
  - umap_params.json : the fitted UMAP hyper-parameters
"""
import argparse
import json
import logging
import sys
import numpy as np
import anndata
import umap as umap_lib

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="UMAP dimensionality reduction on an AnnData h5ad file."
    )
    parser.add_argument("--input",        required=True,  help="HVG-filtered .h5ad file")
    parser.add_argument("--embedding",    required=True,  help="Output embedding .npy file")
    parser.add_argument("--params_json",  required=True,  help="Output UMAP params JSON")
    parser.add_argument("--metric",       default="euclidean")
    parser.add_argument("--random_state", type=int, default=42)
    parser.add_argument("--n_components", type=int, default=40)
    args = parser.parse_args()

    logger.info("Reading %s", args.input)
    adata = anndata.read_h5ad(args.input)

    logger.info(
        "Fitting UMAP (n_components=%d, metric=%s, random_state=%d) on matrix %s",
        args.n_components, args.metric, args.random_state, adata.X.shape,
    )
    mapper = umap_lib.UMAP(
        metric=args.metric,
        random_state=args.random_state,
        n_components=args.n_components,
    )
    embedding = mapper.fit_transform(adata.X)
    logger.info("Embedding shape: %s", embedding.shape)

    logger.info("Saving embedding to %s", args.embedding)
    np.save(args.embedding, embedding)

    params = {k: (v if isinstance(v, (str, int, float, bool, type(None))) else str(v))
              for k, v in mapper.get_params().items()}
    logger.info("Saving UMAP params to %s", args.params_json)
    with open(args.params_json, "w") as fh:
        json.dump(params, fh, indent=2)

    logger.info("Done")


if __name__ == "__main__":
    main()
