"""Βήμα 03A — Generic reference-based RNA annotation με configurable backend"""

from pathlib import Path
import sys

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


# ============================================================
# 0. Refactored project root και Python import path
# ============================================================

refactor_dir = Path(
    r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project/Refactored"
)

src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))


# ============================================================
# 1. Imports από refactored package
# ============================================================

from pbmcgrn.config import load_yaml
from pbmcgrn.io import write_csv, write_json
from pbmcgrn.paths import ProjectPaths


paths = ProjectPaths(refactor_root=refactor_dir)
paths.ensure_output_dirs()

print("Refactored project root:")
print(paths.refactor_root)


# ============================================================
# 2. Helper για resolving paths από config
# ============================================================

def resolve_from_refactor(path_string):
    """
    Resolve relative paths from the Refactored project root.
    """

    path = Path(path_string)

    if path.is_absolute():
        return path

    return (paths.refactor_root / path).resolve()


# ============================================================
# 3. Load annotation backend config
# ============================================================

config_path = paths.configs / "annotation_backend.yaml"

if not config_path.exists():
    raise FileNotFoundError(
        f"Missing config:\n{config_path}\n"
        "Τρέξε πρώτα setup_03A_create_annotation_backend_config.py"
    )

config = load_yaml(config_path)
annotation_cfg = config["annotation"]

method = annotation_cfg["method"]
inputs = annotation_cfg["inputs"]
outputs = annotation_cfg["outputs"]

cluster_key = inputs["cluster_key"]

print("\nLoaded annotation backend config:")
print(config_path)

print("\nSelected annotation method:")
print(method)


# ============================================================
# 4. Load clustered RNA object
# ============================================================

clustered_h5ad = resolve_from_refactor(inputs["clustered_h5ad"])

print("\nLoading clustered RNA object:")
print(clustered_h5ad)

if not clustered_h5ad.exists():
    raise FileNotFoundError(f"Missing clustered RNA object:\n{clustered_h5ad}")

rna_clustered = sc.read_h5ad(clustered_h5ad)

if cluster_key not in rna_clustered.obs.columns:
    raise ValueError(f"Cluster key not found in clustered object: {cluster_key}")

print("\nClustered RNA object:")
print(rna_clustered)


# ============================================================
# 5. Backend 1: CellTypist
# ============================================================

def run_celltypist_backend(annotation_cfg, rna_clustered):
    """
    Run CellTypist as one possible annotation backend.

    Important:
    CellTypist is not hardcoded as the pipeline logic.
    It is only one backend selected from annotation_backend.yaml.
    """

    try:
        import celltypist
        from celltypist import models
    except ImportError as error:
        raise ImportError(
            "CellTypist is not installed.\n"
            "Install with:\n"
            "python -m pip install celltypist"
        ) from error

    ct_cfg = annotation_cfg["celltypist"]
    inputs = annotation_cfg["inputs"]

    full_gene_h5ad = resolve_from_refactor(inputs["full_gene_filtered_h5ad"])

    print("\nRunning CellTypist backend.")
    print("Full-gene filtered RNA input:")
    print(full_gene_h5ad)

    if not full_gene_h5ad.exists():
        raise FileNotFoundError(f"Missing full-gene RNA object:\n{full_gene_h5ad}")

    adata = sc.read_h5ad(full_gene_h5ad)

    print("\nCellTypist full-gene input before normalization:")
    print(adata)

    sc.pp.normalize_total(
        adata,
        target_sum=ct_cfg["normalize_total_target_sum"],
    )

    if ct_cfg["log1p"]:
        sc.pp.log1p(adata)

    common_cells = adata.obs_names.intersection(rna_clustered.obs_names)

    if len(common_cells) != adata.n_obs:
        raise ValueError(
            "Filtered full-gene RNA and clustered RNA do not contain identical cells.\n"
            f"Full-gene cells: {adata.n_obs}\n"
            f"Common cells: {len(common_cells)}"
        )

    rna_clustered_ordered = rna_clustered[adata.obs_names].copy()

    adata.obs[cluster_key] = (
        rna_clustered_ordered.obs[cluster_key]
        .astype(str)
        .values
    )

    if "X_umap" in rna_clustered_ordered.obsm:
        adata.obsm["X_umap"] = rna_clustered_ordered.obsm["X_umap"].copy()

    model_name = ct_cfg["model_name"]

    print("\nCellTypist model:")
    print(model_name)

    try:
        models.download_models(
            force_update=False,
            model=[model_name],
        )
    except TypeError:
        models.download_models(force_update=False)

    model = models.Model.load(model_name)

    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=ct_cfg["majority_voting"],
    )

    adata_pred = predictions.to_adata()

    adata_pred.obs[cluster_key] = adata.obs[cluster_key].values

    if "X_umap" in adata.obsm:
        adata_pred.obsm["X_umap"] = adata.obsm["X_umap"].copy()

    preferred_col = ct_cfg["label_column_preference"]

    if preferred_col in adata_pred.obs.columns:
        label_col = preferred_col
    elif "majority_voting" in adata_pred.obs.columns:
        label_col = "majority_voting"
    elif "predicted_labels" in adata_pred.obs.columns:
        label_col = "predicted_labels"
    else:
        raise ValueError(
            "No usable CellTypist label column found. "
            f"Available columns: {adata_pred.obs.columns.tolist()}"
        )

    adata_pred.obs["reference_label"] = adata_pred.obs[label_col].astype(str)
    adata_pred.obs["reference_method"] = "celltypist"
    adata_pred.obs["reference_model"] = model_name
    adata_pred.obs["reference_label_source_column"] = label_col

    return adata_pred, {
        "backend": "celltypist",
        "model_name": model_name,
        "label_col": label_col,
        "majority_voting": ct_cfg["majority_voting"],
    }


# ============================================================
# 6. Backend 2: External predictions CSV
# ============================================================

def run_external_predictions_csv_backend(annotation_cfg, rna_clustered):
    """
    Use per-cell predictions from an external annotation model.

    This supports:
    - Azimuth exports
    - SingleR exports
    - scANVI predictions
    - custom atlas mapping
    - any future neural/brain annotation model
    """

    ext_cfg = annotation_cfg["external_predictions_csv"]

    csv_path = resolve_from_refactor(ext_cfg["csv_path"])

    if not csv_path.exists():
        raise FileNotFoundError(
            f"External predictions CSV not found:\n{csv_path}"
        )

    barcode_col = ext_cfg["barcode_col"]
    label_col = ext_cfg["label_col"]
    confidence_col = ext_cfg["confidence_col"]

    predictions = pd.read_csv(csv_path)

    required = [barcode_col, label_col]
    missing = [col for col in required if col not in predictions.columns]

    if missing:
        raise ValueError(
            f"Missing required external prediction columns: {missing}\n"
            f"Available columns: {predictions.columns.tolist()}"
        )

    adata_pred = rna_clustered.copy()

    pred_map = predictions.set_index(barcode_col)[label_col].astype(str)

    adata_pred.obs["reference_label"] = (
        adata_pred.obs_names
        .map(pred_map)
    )

    missing_labels = int(adata_pred.obs["reference_label"].isna().sum())

    if missing_labels > 0:
        raise ValueError(
            f"Some cells did not receive external labels: {missing_labels}"
        )

    if confidence_col and confidence_col in predictions.columns:
        conf_map = predictions.set_index(barcode_col)[confidence_col]
        adata_pred.obs["reference_confidence"] = adata_pred.obs_names.map(conf_map)

    adata_pred.obs["reference_method"] = "external_predictions_csv"
    adata_pred.obs["reference_model"] = str(csv_path.name)
    adata_pred.obs["reference_label_source_column"] = label_col

    return adata_pred, {
        "backend": "external_predictions_csv",
        "csv_path": str(csv_path),
        "label_col": label_col,
        "confidence_col": confidence_col,
    }


# ============================================================
# 7. Backend 3: Existing obs column
# ============================================================

def run_existing_obs_column_backend(annotation_cfg, rna_clustered):
    """
    Use labels already stored in an AnnData obs column.
    """

    obs_cfg = annotation_cfg["existing_obs_column"]

    h5ad_path_string = obs_cfg["h5ad_path"]
    label_col = obs_cfg["label_col"]

    if h5ad_path_string:
        h5ad_path = resolve_from_refactor(h5ad_path_string)
        adata_pred = sc.read_h5ad(h5ad_path)
    else:
        adata_pred = rna_clustered.copy()

    if label_col not in adata_pred.obs.columns:
        raise ValueError(
            f"Label column not found in AnnData obs: {label_col}\n"
            f"Available columns: {adata_pred.obs.columns.tolist()}"
        )

    if cluster_key not in adata_pred.obs.columns:
        adata_pred.obs[cluster_key] = rna_clustered.obs[cluster_key].astype(str).values

    adata_pred.obs["reference_label"] = adata_pred.obs[label_col].astype(str)
    adata_pred.obs["reference_method"] = "existing_obs_column"
    adata_pred.obs["reference_model"] = "existing_obs_column"
    adata_pred.obs["reference_label_source_column"] = label_col

    return adata_pred, {
        "backend": "existing_obs_column",
        "label_col": label_col,
    }


# ============================================================
# 8. Run selected backend
# ============================================================

if method == "celltypist":
    adata_pred, backend_info = run_celltypist_backend(
        annotation_cfg=annotation_cfg,
        rna_clustered=rna_clustered,
    )

elif method == "external_predictions_csv":
    adata_pred, backend_info = run_external_predictions_csv_backend(
        annotation_cfg=annotation_cfg,
        rna_clustered=rna_clustered,
    )

elif method == "existing_obs_column":
    adata_pred, backend_info = run_existing_obs_column_backend(
        annotation_cfg=annotation_cfg,
        rna_clustered=rna_clustered,
    )

else:
    raise ValueError(
        f"Unsupported annotation method: {method}\n"
        "Supported methods: celltypist, external_predictions_csv, existing_obs_column"
    )

print("\nReference-annotated object:")
print(adata_pred)

print("\nBackend info:")
print(backend_info)


# ============================================================
# 9. Save cell-level predictions
# ============================================================

cell_predictions = adata_pred.obs.copy()
cell_predictions.insert(0, "cell_barcode", cell_predictions.index)

cell_predictions_out = resolve_from_refactor(outputs["cell_level_predictions"])

write_csv(
    cell_predictions,
    cell_predictions_out,
    index=False,
)

print("\nSaved cell-level reference predictions to:")
print(cell_predictions_out)


# ============================================================
# 10. Aggregate reference labels per Leiden cluster
# ============================================================

high_threshold = annotation_cfg["cluster_aggregation"]["high_consensus_threshold"]
medium_threshold = annotation_cfg["cluster_aggregation"]["medium_consensus_threshold"]

count_table = pd.crosstab(
    adata_pred.obs[cluster_key].astype(str),
    adata_pred.obs["reference_label"].astype(str),
)

count_table.index.name = "cluster"

fraction_table = count_table.div(
    count_table.sum(axis=1),
    axis=0,
)

count_table_out = resolve_from_refactor(outputs["label_count_table"])
fraction_table_out = resolve_from_refactor(outputs["label_fraction_table"])

count_table.to_csv(count_table_out)
fraction_table.to_csv(fraction_table_out)

rows = []

for cluster, row in count_table.iterrows():
    n_cells = int(row.sum())

    sorted_counts = row.sort_values(ascending=False)
    sorted_fractions = sorted_counts / n_cells

    labels = sorted_counts.index.tolist()
    counts = sorted_counts.values.tolist()
    fractions = sorted_fractions.values.tolist()

    top_label = labels[0]
    top_count = int(counts[0])
    top_fraction = float(fractions[0])

    second_label = labels[1] if len(labels) > 1 else ""
    second_count = int(counts[1]) if len(labels) > 1 else 0
    second_fraction = float(fractions[1]) if len(labels) > 1 else 0.0

    if top_fraction >= high_threshold:
        consensus = "high_consensus"
    elif top_fraction >= medium_threshold:
        consensus = "medium_consensus"
    else:
        consensus = "mixed_or_low_consensus"

    rows.append(
        {
            "cluster": str(cluster),
            "n_cells": n_cells,
            "reference_method": method,
            "reference_model": backend_info.get("model_name", backend_info.get("csv_path", "")),
            "top_reference_label": str(top_label),
            "top_reference_count": top_count,
            "top_reference_fraction": top_fraction,
            "second_reference_label": str(second_label),
            "second_reference_count": second_count,
            "second_reference_fraction": second_fraction,
            "top_5_reference_labels": "; ".join(map(str, labels[:5])),
            "top_5_reference_fractions": "; ".join([f"{x:.3f}" for x in fractions[:5]]),
            "candidate_label": str(top_label),
            "candidate_label_source": f"{method}:{backend_info}",
            "candidate_confidence": consensus,
            "manual_review_label": "",
            "manual_review_status": "pending",
            "review_notes": "",
            "final_label": "",
        }
    )

candidate_table = pd.DataFrame(rows)

candidate_table = candidate_table.sort_values(
    by="cluster",
    key=lambda x: x.astype(int),
).reset_index(drop=True)

candidate_table_out = resolve_from_refactor(outputs["candidate_cluster_annotations"])

write_csv(
    candidate_table,
    candidate_table_out,
    index=False,
)

print("\nCandidate cluster annotation table:")
print(candidate_table.to_string(index=False))

print("\nSaved candidate cluster annotations to:")
print(candidate_table_out)


# ============================================================
# 11. Save reference-annotated h5ad
# ============================================================

annotated_h5ad_out = resolve_from_refactor(outputs["annotated_h5ad"])
annotated_h5ad_out.parent.mkdir(parents=True, exist_ok=True)

adata_pred.write_h5ad(annotated_h5ad_out)

print("\nSaved reference-annotated h5ad to:")
print(annotated_h5ad_out)


# ============================================================
# 12. Save UMAP plots
# ============================================================

umap_reference_out = resolve_from_refactor(outputs["umap_reference_labels"])
umap_leiden_out = resolve_from_refactor(outputs["umap_leiden_review"])

umap_reference_out.parent.mkdir(parents=True, exist_ok=True)
umap_leiden_out.parent.mkdir(parents=True, exist_ok=True)

if "X_umap" in adata_pred.obsm:
    sc.pl.umap(
        adata_pred,
        color="reference_label",
        legend_loc="right margin",
        frameon=False,
        show=False,
    )

    plt.savefig(
        umap_reference_out,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    print("\nSaved reference-label UMAP to:")
    print(umap_reference_out)

    sc.pl.umap(
        adata_pred,
        color=cluster_key,
        legend_loc="on data",
        legend_fontsize=6,
        legend_fontweight="normal",
        frameon=False,
        show=False,
    )

    plt.savefig(
        umap_leiden_out,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    print("\nSaved Leiden-review UMAP to:")
    print(umap_leiden_out)

else:
    print("\nNo X_umap found. UMAP plots were skipped.")


# ============================================================
# 13. Save manifest
# ============================================================

manifest = {
    "step": "03A",
    "description": (
        "Generic reference-based RNA annotation using configurable backend. "
        "The selected backend provides candidate cell-type labels; final labels "
        "must be reviewed at cluster level before downstream GRN inference."
    ),
    "config": str(config_path),
    "method": method,
    "backend_info": backend_info,
    "biological_interpretation": {
        "annotation_role": (
            "Reference labels are used to define cell populations for stratified "
            "downstream analysis."
        ),
        "not_used_as": (
            "Reference labels are not used as GRN ground truth and are not used "
            "as regulatory edge labels."
        ),
        "risk": (
            "Incorrect annotation can change cell-type composition and therefore "
            "pseudo-bulk GRN results."
        ),
        "mitigation": (
            "Use cluster-level consensus, DE/marker review, and optional sensitivity "
            "analysis across annotation strategies."
        ),
    },
    "inputs": {
        "clustered_h5ad": str(clustered_h5ad),
    },
    "outputs": {
        "annotated_h5ad": str(annotated_h5ad_out),
        "cell_level_predictions": str(cell_predictions_out),
        "candidate_cluster_annotations": str(candidate_table_out),
        "label_count_table": str(count_table_out),
        "label_fraction_table": str(fraction_table_out),
        "umap_reference_labels": str(umap_reference_out),
        "umap_leiden_review": str(umap_leiden_out),
    },
    "summary": {
        "n_cells": int(adata_pred.n_obs),
        "n_genes": int(adata_pred.n_vars),
        "n_clusters": int(candidate_table.shape[0]),
        "n_reference_labels": int(count_table.shape[1]),
    },
}

manifest_out = paths.manifests / "step03A_reference_based_rna_annotation_manifest.json"

write_json(
    manifest,
    manifest_out,
)

print("\nSaved manifest to:")
print(manifest_out)


# ============================================================
# 14. Τελικό summary
# ============================================================

print("\nDone.")
print("Το generic Step 03A ολοκληρώθηκε.")
print("Το annotation backend είναι πλέον configurable και όχι hardcoded για CellTypist.")
print("Επόμενο βήμα: cluster-level review και harmonization σε thesis cell types.")