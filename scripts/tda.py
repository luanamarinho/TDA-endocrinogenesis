import numpy as np
import pandas as pd  
from gtda.plotting import plot_point_cloud
import scvelo as scv
from pathlib import Path
import utils
import scanpy as sc

from gtda.mapper import (
    CubicalCover,
    make_mapper_pipeline,
    plot_static_mapper_graph,
    plot_interactive_mapper_graph,
    MapperInteractivePlotter
)

from sklearn.cluster import DBSCAN


root = Path.cwd().resolve()
cache_dir = utils._cache_file_folder(root, "data/Pancreas")
pancreasData_cache_path = cache_dir / (
        "endocrinogenesis_day15.h5ad"
    )
if pancreasData_cache_path.exists():
    adata = sc.read_h5ad(pancreasData_cache_path)
else:
    adata = scv.datasets.pancreas()
data = adata.obsm["X_umap"]
print("Shape of pancreas data UMAP: ", data.shape)



# Creating mapper
cover = CubicalCover(n_intervals=10, overlap_frac=0.3)
clusterer = DBSCAN()

pipe = make_mapper_pipeline(
    cover=cover,
    clusterer=clusterer,
    verbose=False,
    n_jobs=1,
)


# TDA
MIP = MapperInteractivePlotter(pipe, data)

MIP.plot(color_data=data)
