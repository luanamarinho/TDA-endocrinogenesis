process NORMALIZE_HVG {
    tag "normalize | hvg"
    label 'process_medium'

    publishDir "${params.outdir}/normalize", mode: 'copy', pattern: "*.{h5ad,tsv}"

    conda "conda-forge::scvelo=0.3.4 conda-forge::anndata=0.10.9 conda-forge::scanpy=1.10.4"

    input:
    path adata

    output:
    path "pancreas_hvg.h5ad",  emit: adata
    path "hvg_table.tsv",      emit: hvg_table

    script:
    def target_sum_arg = params.target_sum != null ? "--target_sum ${params.target_sum}" : ""
    """
    normalize_hvg.py \\
        --input                    ${adata} \\
        --output                   pancreas_hvg.h5ad \\
        --hvg_table                hvg_table.tsv \\
        --exclude_highly_expressed ${params.exclude_highly_expressed} \\
        --max_fraction             ${params.max_fraction} \\
        ${target_sum_arg} \\
        --n_top_genes              ${params.n_top_genes} \\
        --flavor                   ${params.flavor}
    """

    stub:
    """
    touch pancreas_hvg.h5ad
    touch hvg_table.tsv
    """
}
