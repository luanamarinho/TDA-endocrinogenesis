#!/usr/bin/env nextflow

// ============================================================
// Giotto-TDA — Endocrinogenesis scRNA-seq pipeline
//
// Steps:
//   1. LOAD_DATA              — Download the scVelo pancreas dataset
//   2. QC_FILTER              — Cell / gene / MT-fraction filtering
//   3. NORMALIZE_HVG          — Depth normalisation + log1p + HVG selection
//   4. DIMENSIONALITY_REDUCTION — UMAP embedding
// ============================================================

nextflow.enable.dsl = 2

include { LOAD_DATA              } from './modules/local/load_data'
include { QC_FILTER              } from './modules/local/qc_filter'
include { NORMALIZE_HVG          } from './modules/local/normalize_hvg'
include { DIMENSIONALITY_REDUCTION } from './modules/local/dimensionality_reduction'

// ── Parameter validation ─────────────────────────────────────
def validate_params() {
    def valid_flavors = ['cell_ranger', 'seurat', 'seurat_v3', 'seurat_v3_paper']
    if (!valid_flavors.contains(params.flavor)) {
        error "Invalid --flavor '${params.flavor}'. Choose from: ${valid_flavors.join(', ')}"
    }
    if (params.dr_method != 'umap') {
        error "Only --dr_method 'umap' is currently supported"
    }
    if (params.n_components < 2) {
        error "--n_components must be >= 2"
    }
}

// ── Workflow ─────────────────────────────────────────────────
workflow {
    validate_params()

    // 1 — Load raw pancreas data from scvelo
    LOAD_DATA()

    // 2 — QC filter cells and genes
    QC_FILTER(LOAD_DATA.out.adata)

    // 3 — Normalize, log-transform, and select HVGs
    NORMALIZE_HVG(QC_FILTER.out.adata)

    // 4 — UMAP dimensionality reduction
    DIMENSIONALITY_REDUCTION(NORMALIZE_HVG.out.adata)

    // ── Summary emit ──────────────────────────────────────────
    NORMALIZE_HVG.out.hvg_table
        | view { file -> "HVG table   : ${params.outdir}/normalize/${file.name}" }

    DIMENSIONALITY_REDUCTION.out.embedding
        | view { file -> "Embedding   : ${params.outdir}/embedding/${file.name}" }

    DIMENSIONALITY_REDUCTION.out.params_json
        | view { file -> "UMAP params : ${params.outdir}/embedding/${file.name}" }
}


