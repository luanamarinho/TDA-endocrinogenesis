process DIMENSIONALITY_REDUCTION {
    tag "umap | metric=${params.umap_metric}"
    label 'process_high'

    publishDir "${params.outdir}/embedding", mode: 'copy', pattern: "*.{npy,json}"

    conda "conda-forge::umap-learn=0.5.11 conda-forge::numpy=2.4.3 conda-forge::joblib=1.5.3 conda-forge::anndata=0.12.10"

    input:
    path adata

    output:
    path "embedding.npy",      emit: embedding
    path "umap_params.json",   emit: params_json

    script:
    """
    dimensionality_reduction.py \\
        --input        ${adata} \\
        --embedding    embedding.npy \\
        --params_json  umap_params.json \\
        --metric       ${params.umap_metric} \\
        --random_state ${params.umap_random_state} \\
        --n_components ${params.n_components}
    """

    stub:
    """
    touch embedding.npy
    touch umap_params.json
    """
}
