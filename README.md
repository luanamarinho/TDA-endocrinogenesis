# TDA-endocrinogenesis

Applying giotto-tda Mapper to resolve mature and transient pancreatic cellular states.

We analyze single-cell gene expression to visualize unipotent populations (Ngn3+, Fev+) alongside differentiated islet cells (Alpha, Beta, Epsilon). This unified topological framework provides high-resolution insights into developmental differentiation trajectories.

We aimed to follow the preprocessing described in [Bastidas-Ponce et al. (2019)](https://doi.org/10.1242/dev.173849), though certain parameters were interpreted where the original methodology was ambiguous.


# Running the Pipeline
Navigate to the project root and execute the script via the terminal.

# 1. Running with Default Parameters
To run the pipeline using all default settings (e.g., 4000 HVGs, UMAP embedding): ``python endocrinogenesis.py``

# 2. Custom Parameters
Override specific settings directly via the CLI:
``python endocrinogenesis.py --n_top_genes 2000 --exclude_highly_expressed False --flavor seurat_v3``

Note: If the pipeline has not been run previously with the specified arguments, a new bundle file will be cached in `data/Pancreas/cache` within the project folder.

