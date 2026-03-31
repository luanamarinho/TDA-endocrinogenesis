process LOAD_DATA {
    tag "pancreas"
    label 'process_low'

    conda "conda-forge::scvelo=0.3.4 conda-forge::anndata=0.10.9 conda-forge::scanpy=1.10.4"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer'
        ? 'oras://community.wave.seqera.io/library/python_scanpy_scvelo_anndata:8880c8fc704f18c2'
        : 'community.wave.seqera.io/library/python_scanpy_scvelo_anndata:8880c8fc704f18c2' }"

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
