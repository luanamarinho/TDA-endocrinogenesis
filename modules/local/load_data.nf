process LOAD_DATA {
    tag "pancreas"
    label 'process_low'

    conda "conda-forge::scvelo=0.3.4 conda-forge::anndata=0.10.9 conda-forge::scanpy=1.10.4"

    output:
    path "pancreas_raw.h5ad", emit: adata

    script:
    """
    load_data.py --output pancreas_raw.h5ad
    """

    stub:
    """
    touch pancreas_raw.h5ad
    """
}
