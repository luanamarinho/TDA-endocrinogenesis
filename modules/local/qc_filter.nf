process QC_FILTER {
    tag "qc"
    label 'process_low'

    publishDir "${params.outdir}/qc", mode: 'copy', pattern: "*.h5ad"

    conda "conda-forge::scvelo=0.3.4 conda-forge::anndata=0.10.9 conda-forge::scanpy=1.10.4"

    input:
    path adata

    output:
    path "pancreas_qc.h5ad", emit: adata

    script:
    """
    qc_filter.py \\
        --input          ${adata} \\
        --output         pancreas_qc.h5ad \\
        --min_genes      ${params.min_genes} \\
        --min_cells      ${params.min_cells} \\
        --max_mt_perc    ${params.max_mt_perc}
    """

    stub:
    """
    touch pancreas_qc.h5ad
    """
}
