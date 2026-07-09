# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:44:30 2026

@author: jkary
"""
# -*- coding: utf-8 -*-
"""
Step 03C — RNA CellTypist confirmation.

Purpose:
- Load reviewed annotated RNA object from Step 03B.
- Run CellTypist as complementary confirmation.
- Compare reviewed/manual annotation with CellTypist predictions on the same UMAP.
- Save per-cell predictions, cluster-level summaries, comparison table, figures, and manifest.

Important:
- This script does NOT replace reviewed annotation.
- This script does NOT create final labels.
- CellTypist is used only as confirmation evidence.
"""

from pathlib import Path
import sys
from datetime import datetime

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

import celltypist
from celltypist import models


try:
    refactor_dir = Path(__file__).resolve().parents[1]  # Refactored root
except NameError:
    refactor_dir = Path.cwd().resolve()  # interactive fallback

src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))  # allow pbmcgrn imports

from pbmcgrn.config import load_yaml
from pbmcgrn.io import write_csv, write_json
from pbmcgrn.paths import ProjectPaths


def assert_path_inside_project(path: Path, project_root: Path, label: str) -> None:
    """Fail if a path is outside active project."""

    resolved_path = Path(path).resolve()
    resolved_project_root = Path(project_root).resolve()

    try:
        resolved_path.relative_to(resolved_project_root)
    except ValueError as exc:
        raise ValueError(
            f"{label} is outside active project.\n"
            f"Path: {resolved_path}\n"
            f"Project root: {resolved_project_root}"
        ) from exc


def save_show_close(path: Path, dpi: int = 300, show: bool = True) -> None:
    """Save current matplotlib figure and optionally show it."""

    plt.savefig(path, dpi=dpi, bbox_inches="tight")  # save figure

    if show:
        plt.show()  # show interactively

    plt.close()  # close figure to avoid memory buildup


paths = ProjectPaths(refactor_root=refactor_dir)  # project-aware paths
paths.ensure_output_dirs()  # ensure output folders exist

print("\nActive project:")
print(paths.project_id)

print("\nProject root:")
print(paths.project_root)


# ============================================================
# 1. Read selected cluster key from preprocessing_rna.yaml
# ============================================================

preprocessing_config_path = paths.configs / "preprocessing_rna.yaml"  # Step 02B source of truth

assert_path_inside_project(
    preprocessing_config_path,
    paths.project_root,
    "preprocessing_rna.yaml",
)

preprocessing_config = load_yaml(preprocessing_config_path)  # load selected preprocessing config

dim_cfg = preprocessing_config["dimensionality_reduction"]  # PCA/neighbors/Leiden settings

cluster_key = str(dim_cfg["leiden_key"])  # selected cluster key, e.g. leiden_0_4

print("\nUsing selected cluster key:")
print(cluster_key)


# ============================================================
# 2. Resolve paths
# ============================================================

rna_in = paths.processed_rna / "pbmc_rna_annotated_reviewed.h5ad"  # reviewed object from Step 03B

celltypist_annotated_out = paths.prepare_output_path(
    "data/processed/rna/pbmc_rna_annotated_reviewed_celltypist.h5ad",
    label="pbmc_rna_annotated_reviewed_celltypist",
)  # reviewed object plus CellTypist columns

per_cell_out = paths.prepare_output_path(
    "reports/tables/rna/celltypist_per_cell_predictions.csv",
    label="celltypist_per_cell_predictions",
)  # per-cell prediction table

cluster_summary_out = paths.prepare_output_path(
    "reports/tables/rna/celltypist_cluster_label_summary.csv",
    label="celltypist_cluster_label_summary",
)  # label distribution per cluster

majority_out = paths.prepare_output_path(
    "reports/tables/rna/celltypist_cluster_majority_labels.csv",
    label="celltypist_cluster_majority_labels",
)  # majority CellTypist label per cluster

comparison_out = paths.prepare_output_path(
    "reports/tables/rna/celltypist_reviewed_annotation_comparison.csv",
    label="celltypist_reviewed_annotation_comparison",
)  # reviewed annotation vs CellTypist

figures_dir = paths.prepare_output_path(
    "reports/figures/rna/annotation",
    label="annotation_figures_dir",
)  # annotation figure folder

umap_celltypist_out = figures_dir / f"rna_umap_celltypist_majority_by_{cluster_key}.png"

umap_reviewed_vs_celltypist_out = (
    figures_dir / f"rna_umap_reviewed_vs_celltypist_by_{cluster_key}.png"
)

umap_cluster_reviewed_celltypist_out = (
    figures_dir / f"rna_umap_cluster_reviewed_celltypist_by_{cluster_key}.png"
)

manifest_out = paths.prepare_output_path(
    "reports/manifests/step03C_rna_celltypist_confirmation_manifest.json",
    label="step03C_rna_celltypist_confirmation_manifest",
)

for label, path in {
    "rna_in": rna_in,
    "celltypist_annotated_out": celltypist_annotated_out,
    "per_cell_out": per_cell_out,
    "cluster_summary_out": cluster_summary_out,
    "majority_out": majority_out,
    "comparison_out": comparison_out,
    "figures_dir": figures_dir,
    "umap_celltypist_out": umap_celltypist_out,
    "umap_reviewed_vs_celltypist_out": umap_reviewed_vs_celltypist_out,
    "umap_cluster_reviewed_celltypist_out": umap_cluster_reviewed_celltypist_out,
    "manifest_out": manifest_out,
}.items():
    assert_path_inside_project(path, paths.project_root, label)

figures_dir.mkdir(parents=True, exist_ok=True)
per_cell_out.parent.mkdir(parents=True, exist_ok=True)
celltypist_annotated_out.parent.mkdir(parents=True, exist_ok=True)
manifest_out.parent.mkdir(parents=True, exist_ok=True)


# ============================================================
# 3. Load reviewed RNA object
# ============================================================

if not rna_in.exists():
    raise FileNotFoundError(
        f"Reviewed annotated RNA object not found:\n{rna_in}\n\n"
        "Run Step 03B first to create pbmc_rna_annotated_reviewed.h5ad."
    )

rna = sc.read_h5ad(rna_in)  # load reviewed annotated object

if cluster_key not in rna.obs.columns:
    raise ValueError(
        f"Cluster key not found in rna.obs: {cluster_key}\n"
        f"Available Leiden columns: {[col for col in rna.obs.columns if 'leiden' in col.lower()]}"
    )

required_reviewed_columns = [
    "cell_type_broad",
    "cell_type_detailed",
    "annotation_confidence",
    "include_in_downstream",
]  # columns expected from Step 03B

missing_reviewed_columns = [
    col for col in required_reviewed_columns if col not in rna.obs.columns
]

if missing_reviewed_columns:
    raise ValueError(
        f"Reviewed annotation columns missing from rna.obs: {missing_reviewed_columns}\n"
        "This script expects the Step 03B reviewed annotated object."
    )

if "X_umap" not in rna.obsm:
    raise ValueError("X_umap not found in rna.obsm. Step 02B must be completed first.")

print("\nLoaded reviewed RNA object:")
print(rna)

print("\nReviewed annotation counts:")
print(rna.obs["cell_type_broad"].value_counts())

print("\nCluster sizes:")
print(rna.obs[cluster_key].value_counts().sort_index())


# ============================================================
# 4. Run CellTypist
# ============================================================

celltypist_model = "Immune_All_Low.pkl"  # broad immune model, good first confirmation layer

print("\nCellTypist version:")
print(celltypist.__version__)

print("\nChecking/downloading CellTypist models...")
models.download_models(force_update=False)  # download model if needed

print("\nRunning CellTypist with model:")
print(celltypist_model)

predictions = celltypist.annotate(
    rna,
    model=celltypist_model,
    majority_voting=True,
)  # run CellTypist; reviewed annotation remains primary

pred_adata = predictions.to_adata()  # convert prediction result to AnnData

celltypist_obs = pred_adata.obs.copy()  # prediction obs table

print("\nCellTypist output columns:")
print(celltypist_obs.columns.tolist())

if "predicted_labels" not in celltypist_obs.columns:
    raise ValueError(
        "CellTypist output does not contain 'predicted_labels'. "
        f"Available columns: {celltypist_obs.columns.tolist()}"
    )

rna.obs["celltypist_predicted_label"] = (
    celltypist_obs["predicted_labels"].astype(str).values
)  # raw per-cell CellTypist label

if "majority_voting" in celltypist_obs.columns:
    rna.obs["celltypist_majority_voting"] = (
        celltypist_obs["majority_voting"].astype(str).values
    )  # majority-voted CellTypist label
else:
    rna.obs["celltypist_majority_voting"] = rna.obs[
        "celltypist_predicted_label"
    ]  # fallback

score_candidates = ["conf_score", "probability", "score"]  # possible confidence-like columns

score_col = next(
    (col for col in score_candidates if col in celltypist_obs.columns),
    None,
)  # find confidence score if available

if score_col is not None:
    rna.obs["celltypist_confidence"] = celltypist_obs[score_col].values  # save score
else:
    rna.obs["celltypist_confidence"] = pd.NA  # keep column for consistency

print("\nCellTypist predicted label counts:")
print(rna.obs["celltypist_predicted_label"].value_counts())

print("\nCellTypist majority-voting label counts:")
print(rna.obs["celltypist_majority_voting"].value_counts())


# ============================================================
# 5. Save per-cell predictions
# ============================================================

per_cell = rna.obs[
    [
        cluster_key,
        "cell_type_broad",
        "cell_type_detailed",
        "annotation_confidence",
        "include_in_downstream",
        "celltypist_predicted_label",
        "celltypist_majority_voting",
        "celltypist_confidence",
    ]
].copy()  # combine reviewed labels and CellTypist labels per cell

per_cell.insert(
    0,
    "cell_barcode",
    per_cell.index.astype(str),
)  # preserve cell barcode/index

write_csv(per_cell, per_cell_out, index=False)

print("\nSaved CellTypist per-cell predictions to:")
print(per_cell_out)


# ============================================================
# 6. Cluster-level CellTypist summaries
# ============================================================

cluster_summary = (
    rna.obs
    .groupby([cluster_key, "celltypist_majority_voting"], dropna=False)
    .size()
    .reset_index(name="n_cells")
)  # count CellTypist labels inside each Leiden cluster

cluster_totals = (
    rna.obs
    .groupby(cluster_key, dropna=False)
    .size()
    .rename("cluster_n_cells")
    .reset_index()
)  # total cells per cluster

cluster_summary = cluster_summary.merge(
    cluster_totals,
    on=cluster_key,
    how="left",
)  # add total cluster size

cluster_summary["fraction_in_cluster"] = (
    cluster_summary["n_cells"] / cluster_summary["cluster_n_cells"]
)  # fraction of each CellTypist label inside cluster

cluster_summary = cluster_summary.sort_values(
    [cluster_key, "n_cells"],
    ascending=[True, False],
)  # majority label first

write_csv(cluster_summary, cluster_summary_out, index=False)

print("\nSaved CellTypist cluster label summary to:")
print(cluster_summary_out)

majority = (
    cluster_summary
    .sort_values([cluster_key, "n_cells"], ascending=[True, False])
    .groupby(cluster_key, as_index=False)
    .head(1)
    .copy()
)  # one majority label per cluster

majority = majority.rename(
    columns={
        "celltypist_majority_voting": "celltypist_majority_label",
        "n_cells": "celltypist_majority_n_cells",
        "fraction_in_cluster": "celltypist_majority_fraction",
    }
)  # clearer names

write_csv(majority, majority_out, index=False)

print("\nSaved CellTypist majority labels to:")
print(majority_out)

print("\nCellTypist majority label per cluster:")
print(
    majority[
        [
            cluster_key,
            "celltypist_majority_label",
            "celltypist_majority_n_cells",
            "cluster_n_cells",
            "celltypist_majority_fraction",
        ]
    ]
)


# ============================================================
# 7. Compare reviewed annotation vs CellTypist majority labels
# ============================================================

reviewed_cluster_summary = (
    rna.obs
    .groupby(
        [
            cluster_key,
            "cell_type_broad",
            "cell_type_detailed",
            "annotation_confidence",
            "include_in_downstream",
        ],
        dropna=False,
    )
    .size()
    .reset_index(name="reviewed_n_cells")
)  # one row per reviewed cluster annotation

reviewed_cluster_summary[cluster_key] = reviewed_cluster_summary[cluster_key].astype(str)
majority[cluster_key] = majority[cluster_key].astype(str)

comparison = reviewed_cluster_summary.merge(
    majority[
        [
            cluster_key,
            "celltypist_majority_label",
            "celltypist_majority_n_cells",
            "cluster_n_cells",
            "celltypist_majority_fraction",
        ]
    ],
    on=cluster_key,
    how="left",
)  # reviewed annotation + CellTypist majority annotation

comparison["celltypist_role"] = "confirmation_evidence_not_final_truth"  # thesis audit note

comparison["manual_review_needed"] = (
    comparison["celltypist_majority_fraction"] < 0.75
)  # practical flag for mixed/uncertain CellTypist clusters

write_csv(comparison, comparison_out, index=False)

print("\nSaved reviewed annotation vs CellTypist comparison to:")
print(comparison_out)


# ============================================================
# 8. Save annotated object with CellTypist columns
# ============================================================

rna.write_h5ad(celltypist_annotated_out)  # save reviewed object plus CellTypist evidence columns

print("\nSaved reviewed RNA object with CellTypist columns to:")
print(celltypist_annotated_out)


# ============================================================
# 9. UMAP figures
# ============================================================

sc.pl.umap(
    rna,
    color="celltypist_majority_voting",
    legend_loc="right margin",
    frameon=False,
    title="CellTypist majority voting",
    show=False,
)
save_show_close(umap_celltypist_out, dpi=300, show=True)

sc.pl.umap(
    rna,
    color=["cell_type_broad", "celltypist_majority_voting"],
    legend_loc="right margin",
    frameon=False,
    title=[
        "Reviewed broad annotation",
        "CellTypist majority voting",
    ],
    show=False,
)
save_show_close(umap_reviewed_vs_celltypist_out, dpi=300, show=True)

sc.pl.umap(
    rna,
    color=[cluster_key, "cell_type_broad", "celltypist_majority_voting"],
    legend_loc="right margin",
    frameon=False,
    title=[
        "Leiden clusters",
        "Reviewed broad annotation",
        "CellTypist majority voting",
    ],
    show=False,
)
save_show_close(umap_cluster_reviewed_celltypist_out, dpi=300, show=True)


# ============================================================
# 10. Manifest
# ============================================================

manifest = {
    "step": "03C",
    "description": "CellTypist confirmation of reviewed RNA annotation",
    "timestamp": datetime.now().isoformat(timespec="seconds"),
    "project": paths.as_dict(),
    "inputs": {
        "reviewed_annotated_h5ad": str(rna_in),
        "preprocessing_config": str(preprocessing_config_path),
    },
    "celltypist": {
        "version": str(celltypist.__version__),
        "model": celltypist_model,
        "majority_voting": True,
        "role": "confirmation_evidence_not_final_truth",
    },
    "clustering_source": {
        "cluster_key": cluster_key,
        "n_pcs": int(dim_cfg["n_pcs"]),
        "n_neighbors": int(dim_cfg["n_neighbors"]),
        "leiden_resolution": float(dim_cfg["leiden_resolution"]),
    },
    "outputs": {
        "annotated_reviewed_celltypist_h5ad": str(celltypist_annotated_out),
        "per_cell_predictions": str(per_cell_out),
        "cluster_label_summary": str(cluster_summary_out),
        "cluster_majority_labels": str(majority_out),
        "reviewed_vs_celltypist_comparison": str(comparison_out),
        "umap_celltypist_majority": str(umap_celltypist_out),
        "umap_reviewed_vs_celltypist": str(umap_reviewed_vs_celltypist_out),
        "umap_cluster_reviewed_celltypist": str(umap_cluster_reviewed_celltypist_out),
    },
    "summary": {
        "n_cells": int(rna.n_obs),
        "n_genes": int(rna.n_vars),
        "n_clusters": int(rna.obs[cluster_key].nunique()),
        "n_reviewed_broad_labels": int(rna.obs["cell_type_broad"].nunique()),
        "n_celltypist_majority_labels": int(rna.obs["celltypist_majority_voting"].nunique()),
        "n_clusters_flagged_for_review": int(comparison["manual_review_needed"].sum()),
    },
    "rules": {
        "reviewed_annotation_remains_primary": True,
        "celltypist_used_as_confirmation_only": True,
        "does_not_replace_final_labels": True,
        "does_not_modify_step03B_outputs": True,
    },
}

write_json(manifest, manifest_out)

print("\nSaved Step 03C manifest to:")
print(manifest_out)

print("\nDone.")
print("Step 03C CellTypist confirmation completed.")
print("Open the reviewed-vs-CellTypist UMAP and comparison table next.")