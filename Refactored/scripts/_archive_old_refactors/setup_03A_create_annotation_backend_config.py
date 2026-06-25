"""Βήμα 03A.0 — Δημιουργία generic config για reference-based RNA annotation backends"""

from pathlib import Path
import yaml


# ============================================================
# 0. Refactored project root
# ============================================================

refactor_dir = Path(
    r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project/Refactored"
)

configs_dir = refactor_dir / "configs"
configs_dir.mkdir(parents=True, exist_ok=True)

config_out = configs_dir / "annotation_backend.yaml"

print("Writing annotation backend config to:")
print(config_out)


# ============================================================
# 1. Generic annotation backend config
# ============================================================

config = {
    "annotation": {
        "dataset_name": "pbmc_rna",
        "organism": "human",
        "biological_domain": "PBMC / immune",
        "method": "celltypist",

        "inputs": {
            "full_gene_filtered_h5ad": "data/processed/rna/pbmc_rna_filtered.h5ad",
            "clustered_h5ad": "data/processed/rna/pbmc_rna_clustered.h5ad",
            "cluster_key": "leiden_0_6",
        },

        "celltypist": {
            "model_name": "Immune_All_Low.pkl",
            "majority_voting": True,
            "label_column_preference": "majority_voting",
            "normalize_total_target_sum": 10000,
            "log1p": True,
        },

        "external_predictions_csv": {
            "csv_path": "",
            "barcode_col": "cell_barcode",
            "label_col": "predicted_label",
            "confidence_col": "",
            "notes": (
                "Use this mode for annotations from Azimuth, SingleR, scANVI, "
                "manual atlas mapping, or any external model."
            ),
        },

        "existing_obs_column": {
            "h5ad_path": "",
            "label_col": "",
            "confidence_col": "",
            "notes": (
                "Use this mode if annotation labels are already stored inside "
                "an AnnData obs column."
            ),
        },

        "cluster_aggregation": {
            "high_consensus_threshold": 0.80,
            "medium_consensus_threshold": 0.60,
        },

        "outputs": {
            "annotated_h5ad": "data/processed/rna/pbmc_rna_reference_annotated.h5ad",
            "cell_level_predictions": "reports/tables/rna_reference_cell_level_predictions.csv",
            "candidate_cluster_annotations": "reports/tables/rna_reference_candidate_cluster_annotations.csv",
            "label_count_table": "reports/tables/rna_reference_label_counts_by_cluster.csv",
            "label_fraction_table": "reports/tables/rna_reference_label_fractions_by_cluster.csv",
            "umap_reference_labels": "reports/figures/rna/pbmc_rna_umap_reference_labels.png",
            "umap_leiden_review": "reports/figures/rna/pbmc_rna_umap_leiden_reference_review.png",
        },
    }
}


# ============================================================
# 2. Save YAML
# ============================================================

with open(config_out, "w", encoding="utf-8") as handle:
    yaml.safe_dump(
        config,
        handle,
        sort_keys=False,
        allow_unicode=True,
        width=120,
    )

print("\nDone.")
print("Δημιουργήθηκε generic annotation backend config.")
print("Για PBMC τώρα χρησιμοποιούμε method: celltypist.")
print("Για άλλα datasets μπορείς να αλλάξεις method σε external_predictions_csv ή existing_obs_column.")