#!/usr/bin/env python3
"""
Normalize counts, log-transform, and select highly variable genes (HVG).
Outputs the HVG-filtered AnnData and a TSV listing the selected genes.
"""
import argparse
import logging
import sys
import scanpy as sc
import anndata
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Normalize, log-transform and select HVGs."
    )
    parser.add_argument("--input",                    required=True)
    parser.add_argument("--output",                   required=True,  help="HVG-filtered .h5ad")
    parser.add_argument("--hvg_table",                required=True,  help="HVG summary TSV")
    parser.add_argument("--exclude_highly_expressed", type=lambda x: x.lower() == "true",
                                                      default=True)
    parser.add_argument("--max_fraction",             type=float, default=0.05)
    parser.add_argument("--target_sum",               type=float, default=None)
    parser.add_argument("--n_top_genes",              type=int,   default=4000)
    parser.add_argument("--flavor",                   type=str,   default="cell_ranger",
                        choices=["cell_ranger", "seurat", "seurat_v3", "seurat_v3_paper"])
    args = parser.parse_args()

    logger.info("Reading %s", args.input)
    adata = anndata.read_h5ad(args.input)

    # For seurat_v3 / seurat_v3_paper, HVG must be called on raw counts *before* normalisation
    if args.flavor in ("seurat_v3", "seurat_v3_paper"):
        logger.info("Selecting HVGs on raw counts (flavor=%s)", args.flavor)
        sc.pp.highly_variable_genes(
            adata, flavor=args.flavor, n_top_genes=args.n_top_genes
        )

    logger.info(
        "Normalizing counts (exclude_highly_expressed=%s, max_fraction=%.2f, target_sum=%s)",
        args.exclude_highly_expressed, args.max_fraction, args.target_sum,
    )
    sc.pp.normalize_total(
        adata,
        exclude_highly_expressed=args.exclude_highly_expressed,
        max_fraction=args.max_fraction,
        target_sum=args.target_sum,
    )

    logger.info("Log1p transforming")
    sc.pp.log1p(adata)

    if args.flavor not in ("seurat_v3", "seurat_v3_paper"):
        logger.info(
            "Selecting HVGs on log-normalized data (flavor=%s, n_top=%d)",
            args.flavor, args.n_top_genes,
        )
        sc.pp.highly_variable_genes(
            adata, flavor=args.flavor, n_top_genes=args.n_top_genes
        )

    n_hvg = adata.var["highly_variable"].sum()
    logger.info("Selected %d highly variable genes", n_hvg)

    # Subset to HVGs
    adata = adata[:, adata.var["highly_variable"]].copy()
    logger.info("Shape after HVG subsetting: %d cells × %d genes", adata.n_obs, adata.n_vars)

    logger.info("Writing HVG-filtered AnnData to %s", args.output)
    adata.write_h5ad(args.output)

    logger.info("Writing HVG table to %s", args.hvg_table)
    adata.var.to_csv(args.hvg_table, sep="\t")

    logger.info("Done")


if __name__ == "__main__":
    main()
