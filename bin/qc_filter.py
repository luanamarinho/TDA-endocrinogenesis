#!/usr/bin/env python3
"""
Apply QC filters (cell count, gene count, mitochondrial fraction) to an AnnData object.
"""
import argparse
import logging
import sys
import numpy as np
import scanpy as sc
import anndata

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def qc_filter(
    adata: anndata.AnnData,
    min_genes: int,
    min_cells: int,
    max_mt_perc: float,
) -> anndata.AnnData:
    logger.info(
        "Before QC: %d cells × %d genes", adata.n_obs, adata.n_vars
    )

    sc.pp.filter_cells(adata, min_genes=min_genes, inplace=True)
    logger.info("After min_genes=%d filter: %d cells", min_genes, adata.n_obs)

    sc.pp.filter_genes(adata, min_cells=min_cells, inplace=True)
    logger.info("After min_cells=%d filter: %d genes", min_cells, adata.n_vars)

    mt_genes = adata.var_names.str.upper().str.startswith("MT-")
    mt_frac = 100 * (
        np.asarray(adata[:, mt_genes].X.sum(axis=1)).flatten()
        / np.asarray(adata.X.sum(axis=1)).flatten()
    )
    adata = adata[mt_frac < max_mt_perc].copy()
    logger.info(
        "After max_mt_perc=%.1f%% filter: %d cells", max_mt_perc, adata.n_obs
    )

    logger.info(
        "Final QC shape: %d cells × %d genes", adata.n_obs, adata.n_vars
    )
    return adata


def main():
    parser = argparse.ArgumentParser(
        description="Apply QC filters to an AnnData h5ad file."
    )
    parser.add_argument("--input",       required=True,  help="Input .h5ad file")
    parser.add_argument("--output",      required=True,  help="Output .h5ad file")
    parser.add_argument("--min_genes",   type=int,   default=1200)
    parser.add_argument("--min_cells",   type=int,   default=20)
    parser.add_argument("--max_mt_perc", type=float, default=20.0)
    args = parser.parse_args()

    logger.info("Reading %s", args.input)
    adata = anndata.read_h5ad(args.input)

    adata = qc_filter(
        adata,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        max_mt_perc=args.max_mt_perc,
    )

    logger.info("Writing filtered data to %s", args.output)
    adata.write_h5ad(args.output)
    logger.info("Done")


if __name__ == "__main__":
    main()
