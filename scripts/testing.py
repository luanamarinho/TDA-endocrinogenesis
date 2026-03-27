import argparse
import logging
from pathlib import Path
import scanpy as sc
import scvelo as scv
from time import perf_counter
import umap
import numpy as np
from scripts import utils
from scripts import preprocessing as prep
from anndata import AnnData
import importlib


importlib.reload(prep) 



adata = scv.datasets.pancreas()

prep.qc_filter(adata)

rmatrix = adata.X.copy()

prep.normalize_counts(adata)

adata.layers["rmatrix"] = rmatrix

adata.layers["logNormal"] = sc.pp.log1p(adata.X)

sc.pp.highly_variable_genes(adata,layer="rmatrix", n_top_genes=2000, flavor = "seurat_v3")
hgv_seurat = adata.var["highly_variable_genes"]

sc.pp.highly_variable_genes(adata,layer="logNormal", n_top_genes=2000, flavor = "cell_ranger")
hvg_ranger = adata.var["highly_variable_genes"]

np.mean(adata.layers["rmatrix"] - adata.layers["logNormal"])

sRaw = adata.layers["rmatrix"].sum(axis=1).flatten()
sNorm = adata.layers["logNormal"].sum(axis=1).flatten()

