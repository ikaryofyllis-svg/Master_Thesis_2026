# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 14:13:21 2026

@author: jkary
"""

"""Βήμα 00G — Δημιουργία configs για raw inputs, RNA preprocessing και ATAC preprocessing"""

from pathlib import Path


# ============================================================
# 0. Refactored project root
# ============================================================

refactor_dir = Path(
    r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project/Refactored"
)

configs_dir = refactor_dir / "configs"
configs_dir.mkdir(parents=True, exist_ok=True)

print("Refactored directory:")
print(refactor_dir)

print("\nConfigs directory:")
print(configs_dir)


# ============================================================
# 1. Config για raw input files
# ============================================================

raw_inputs_yaml = configs_dir / "raw_inputs.yaml"

raw_inputs_content = """raw_inputs:
  rna_10x_h5: "../data/raw/pbmc_rna_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5"

  atac_10x_h5: "../data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_filtered_peak_bc_matrix.h5"
  atac_fragments: "../data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_fragments.tsv.gz"
  atac_fragments_index: "../data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_fragments.tsv.gz.tbi"
  atac_peaks_bed: "../data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_peaks.bed"
  atac_peak_annotation: "../data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_peak_annotation.tsv"
  atac_peak_motif_mapping: "../data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_peak_motif_mapping.bed"

outputs:
  raw_rna_h5ad: "data/processed/rna/pbmc_rna_raw.h5ad"
  raw_atac_h5ad: "data/processed/atac/pbmc_atac_raw.h5ad"
  file_check_table: "reports/tables/raw_input_file_check.csv"
"""


if raw_inputs_yaml.exists():
    print("SKIP existing:", raw_inputs_yaml)
else:
    raw_inputs_yaml.write_text(raw_inputs_content, encoding="utf-8")
    print("CREATED:", raw_inputs_yaml)


# ============================================================
# 2. Config για RNA preprocessing
# ============================================================

rna_yaml = configs_dir / "preprocessing_rna.yaml"

rna_content = """rna_preprocessing:
  min_genes_per_cell: 200
  min_cells_per_gene: 3
  max_pct_mito: 15
  max_total_counts: 40000

normalization:
  target_sum: 10000
  log1p: true

highly_variable_genes:
  n_top_genes: 3000
  flavor: "seurat"

dimensionality_reduction:
  scale_max_value: 10
  pca_solver: "arpack"
  n_neighbors: 15
  n_pcs: 30
  leiden_resolution: 0.6
  leiden_key: "leiden_0_6"

outputs:
  rna_qc_table: "reports/tables/rna_qc_summary.csv"
  rna_filtered_h5ad: "data/processed/rna/pbmc_rna_filtered.h5ad"
  rna_hvg_h5ad: "data/processed/rna/pbmc_rna_hvg.h5ad"
  rna_clustered_h5ad: "data/processed/rna/pbmc_rna_clustered.h5ad"
  rna_umap_clusters: "reports/figures/rna/pbmc_rna_umap_leiden_0_6.png"
"""


if rna_yaml.exists():
    print("SKIP existing:", rna_yaml)
else:
    rna_yaml.write_text(rna_content, encoding="utf-8")
    print("CREATED:", rna_yaml)


# ============================================================
# 3. Config για ATAC preprocessing
# ============================================================

atac_yaml = configs_dir / "preprocessing_atac.yaml"

atac_content = """atac_preprocessing:
  use_gex_only_false: true
  keep_feature_type: "Peaks"

qc:
  min_peaks_per_cell: 1000
  max_counts_per_cell: 60000
  min_cells_per_peak: 20

outputs:
  atac_qc_table_before_filtering: "reports/tables/atac_qc_summary_before_filtering.csv"
  atac_qc_quantiles: "reports/tables/atac_qc_quantiles.csv"
  atac_filtering_summary: "reports/tables/atac_filtering_summary.csv"
  atac_filtered_h5ad: "data/processed/atac/pbmc_atac_filtered.h5ad"
"""


if atac_yaml.exists():
    print("SKIP existing:", atac_yaml)
else:
    atac_yaml.write_text(atac_content, encoding="utf-8")
    print("CREATED:", atac_yaml)


# ============================================================
# 4. Τελικό summary
# ============================================================

print("\nDone.")
print("Δημιουργήθηκαν/ελέγχθηκαν τα preprocessing configs:")
print(raw_inputs_yaml)
print(rna_yaml)
print(atac_yaml)