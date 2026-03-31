#!/usr/bin/env python3
"""
Load the scVelo pancreas dataset and save as an AnnData h5ad file.
"""
import argparse
import logging
import sys
import scvelo as scv

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Load scVelo pancreas dataset.")
    parser.add_argument("--output", required=True, help="Output .h5ad file path")
    args = parser.parse_args()

    logger.info("Loading scVelo pancreas dataset")
    adata = scv.datasets.pancreas()

    logger.info("Saving to %s  (shape: %s)", args.output, adata.shape)
    adata.write_h5ad(args.output)
    logger.info("Done")


if __name__ == "__main__":
    main()
