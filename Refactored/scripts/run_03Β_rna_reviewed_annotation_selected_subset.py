# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 15:32:06 2026

@author: jkary
"""

"""
Step 03B — RNA reviewed cluster annotation.

Purpose:
- Read selected cluster key dynamically from preprocessing_rna.yaml.
- Read manual reviewed cluster annotation table.
- Add reviewed annotation columns to RNA object.
- Keep ambiguous/small/mixed clusters visible in the full annotated object.
- Create a selected downstream subset only from include_in_downstream == yes.
- Save audit tables, annotated object, selected object, and manifest.
"""

from pathlib import Path
import sys
from datetime import datetime

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

try:
    refactor_dir = Path(__file__).resolve().parents[1]  # Refactored root
except NameError:
    refactor_dir = Path.cwd().resolve()

src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

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

    plt.savefig(path, dpi=dpi, bbox_inches="tight")  # save reproducible figure

    if show:
        plt.show()  # display interactively

    plt.close()  # avoid memory buildup


paths = ProjectPaths(refactor_root=refactor_dir)  # project-aware paths
paths.ensure_output_dirs()  # create project output dirs

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

cluster_key = str(dim_cfg["leiden_key"])  # dynamic selected cluster key

print("\nUsing selected cluster key:")
print(cluster_key)


# ============================================================
# 2. Resolve paths
# ============================================================

rna_in = paths.processed_rna / "pbmc_rna_clustered.h5ad"  # selected clustered RNA object

reviewed_annotation_path = paths.prepare_output_path(
    "reports/tables/rna/rna_reviewed_cluster_annotation.csv",
    label="rna_reviewed_cluster_annotation",
)  # manual reviewed annotation table

annotated_out = paths.prepare_output_path(
    "data/processed/rna/pbmc_rna_annotated_reviewed.h5ad",
    label="pbmc_rna_annotated_reviewed",
)  # full annotated object, includes ambiguous clusters

selected_out = paths.prepare_output_path(
    "data/processed/rna/pbmc_rna_selected_high_confidence.h5ad",
    label="pbmc_rna_selected_high_confidence",
)  # downstream subset, only include_in_downstream == yes

annotation_summary_out = paths.prepare_output_path(
    "reports/tables/rna/rna_reviewed_annotation_summary.csv",
    label="rna_reviewed_annotation_summary",
)  # cluster-level annotation summary

cell_counts_out = paths.prepare_output_path(
    "reports/tables/rna/rna_reviewed_cell_type_counts.csv",
    label="rna_reviewed_cell_type_counts",
)  # cell counts by reviewed label

selected_counts_out = paths.prepare_output_path(
    "reports/tables/rna/rna_selected_high_confidence_cell_counts.csv",
    label="rna_selected_high_confidence_cell_counts",
)  # selected subset counts

figures_dir = paths.prepare_output_path(
    "reports/figures/rna/annotation",
    label="annotation_figures_dir",
)  # project-local annotation figure directory

umap_broad_out = figures_dir / f"rna_umap_cell_type_broad_by_{cluster_key}.png"
umap_detailed_out = figures_dir / f"rna_umap_cell_type_detailed_by_{cluster_key}.png"
umap_include_out = figures_dir / f"rna_umap_include_downstream_by_{cluster_key}.png"
umap_cluster_and_annotation_out = figures_dir / f"rna_umap_cluster_and_annotation_by_{cluster_key}.png"

manifest_out = paths.prepare_output_path(
    "reports/manifests/step03B_rna_reviewed_annotation_manifest.json",
    label="step03B_manifest",
)

for label, path in {
    "rna_in": rna_in,
    "reviewed_annotation_path": reviewed_annotation_path,
    "annotated_out": annotated_out,
    "selected_out": selected_out,
    "annotation_summary_out": annotation_summary_out,
    "cell_counts_out": cell_counts_out,
    "selected_counts_out": selected_counts_out,
    "figures_dir": figures_dir,
    "manifest_out": manifest_out,
}.items():
    assert_path_inside_project(path, paths.project_root, label)

figures_dir.mkdir(parents=True, exist_ok=True)
annotated_out.parent.mkdir(parents=True, exist_ok=True)
selected_out.parent.mkdir(parents=True, exist_ok=True)


# ============================================================
# 3. Load RNA object and reviewed annotation table
# ============================================================

if not rna_in.exists():
    raise FileNotFoundError(f"RNA clustered object not found:\n{rna_in}")

if not reviewed_annotation_path.exists():
    raise FileNotFoundError(
        f"Reviewed annotation table not found:\n{reviewed_annotation_path}\n\n"
        "Create this table from rna_informed_annotation_template.csv, then fill manual labels."
    )

rna = sc.read_h5ad(rna_in)  # load selected clustered RNA object

if cluster_key not in rna.obs.columns:
    raise ValueError(
        f"Cluster key from preprocessing_rna.yaml not found in rna.obs: {cluster_key}"
    )

reviewed = pd.read_csv(reviewed_annotation_path, dtype={"cluster": str})  # manual reviewed mapping

required_columns = [
    "cluster",
    "cell_type_broad",
    "cell_type_detailed",
    "annotation_confidence",
    "include_in_downstream",
    "exclusion_reason",
    "needs_subclustering",
    "notes",
]

missing_columns = [col for col in required_columns if col not in reviewed.columns]

if missing_columns:
    raise ValueError(
        f"Reviewed annotation table is missing required columns: {missing_columns}"
    )

reviewed["cluster"] = reviewed["cluster"].astype(str)  # stable mapping key
rna.obs[cluster_key] = rna.obs[cluster_key].astype(str)  # stable cluster labels


# ============================================================
# 4. Validate every cluster has a reviewed decision
# ============================================================

clusters_in_object = set(rna.obs[cluster_key].astype(str).unique())
clusters_in_table = set(reviewed["cluster"].astype(str).unique())

missing_from_table = sorted(clusters_in_object - clusters_in_table)
extra_in_table = sorted(clusters_in_table - clusters_in_object)

if missing_from_table:
    raise ValueError(
        f"These clusters exist in the RNA object but are missing from reviewed table:\n"
        f"{missing_from_table}"
    )

if extra_in_table:
    raise ValueError(
        f"These clusters are in reviewed table but not in RNA object:\n"
        f"{extra_in_table}"
    )


# ============================================================
# 5. Add reviewed annotations to rna.obs
# ============================================================

mapping_broad = reviewed.set_index("cluster")["cell_type_broad"].to_dict()
mapping_detailed = reviewed.set_index("cluster")["cell_type_detailed"].to_dict()
mapping_confidence = reviewed.set_index("cluster")["annotation_confidence"].to_dict()
mapping_include = reviewed.set_index("cluster")["include_in_downstream"].to_dict()
mapping_exclusion = reviewed.set_index("cluster")["exclusion_reason"].fillna("").to_dict()
mapping_subcluster = reviewed.set_index("cluster")["needs_subclustering"].to_dict()
mapping_notes = reviewed.set_index("cluster")["notes"].fillna("").to_dict()

rna.obs["cell_type_broad"] = rna.obs[cluster_key].map(mapping_broad)  # broad reviewed annotation
rna.obs["cell_type_detailed"] = rna.obs[cluster_key].map(mapping_detailed)  # cluster-aware annotation
rna.obs["annotation_confidence"] = rna.obs[cluster_key].map(mapping_confidence)  # high/medium/low
rna.obs["include_in_downstream"] = rna.obs[cluster_key].map(mapping_include)  # yes/no
rna.obs["exclusion_reason"] = rna.obs[cluster_key].map(mapping_exclusion)  # why excluded
rna.obs["needs_subclustering"] = rna.obs[cluster_key].map(mapping_subcluster)  # yes/no/maybe
rna.obs["annotation_notes"] = rna.obs[cluster_key].map(mapping_notes)  # manual evidence notes

missing_annotation_cells = rna.obs["cell_type_broad"].isna().sum()

if missing_annotation_cells > 0:
    raise ValueError(
        f"{missing_annotation_cells} cells have missing reviewed annotation."
    )


# ============================================================
# 6. Save full reviewed annotated object
# ============================================================

rna.write_h5ad(annotated_out)  # full object retains all clusters, including excluded/ambiguous

print("\nSaved reviewed annotated full RNA object to:")
print(annotated_out)


# ============================================================
# 7. Build selected high-confidence subset
# ============================================================

include_mask = (
    rna.obs["include_in_downstream"]
    .astype(str)
    .str.lower()
    .isin(["yes", "true", "1"])
)  # keep only reviewed good clusters

rna_selected = rna[include_mask].copy()  # selected downstream object

rna_selected.write_h5ad(selected_out)  # save selected subset

print("\nSelected high-confidence RNA subset:")
print(rna_selected)

print("\nSaved selected high-confidence RNA object to:")
print(selected_out)


# ============================================================
# 8. Save audit tables
# ============================================================

annotation_summary = (
    rna.obs
    .groupby(
        [
            cluster_key,
            "cell_type_broad",
            "cell_type_detailed",
            "annotation_confidence",
            "include_in_downstream",
            "exclusion_reason",
            "needs_subclustering",
        ],
        dropna=False,
    )
    .size()
    .reset_index(name="n_cells")
    .sort_values(cluster_key)
)  # one row per cluster annotation decision

write_csv(annotation_summary, annotation_summary_out, index=False)

cell_counts = (
    rna.obs["cell_type_broad"]
    .value_counts()
    .rename_axis("cell_type_broad")
    .reset_index(name="n_cells")
)  # full object counts

write_csv(cell_counts, cell_counts_out, index=False)

selected_counts = (
    rna_selected.obs["cell_type_broad"]
    .value_counts()
    .rename_axis("cell_type_broad")
    .reset_index(name="n_cells")
)  # selected object counts

write_csv(selected_counts, selected_counts_out, index=False)

print("\nSaved reviewed annotation summary to:")
print(annotation_summary_out)

print("\nSaved full cell type counts to:")
print(cell_counts_out)

print("\nSaved selected high-confidence counts to:")
print(selected_counts_out)


# ============================================================
# 9. UMAP visualizations with reviewed annotations
# ============================================================

sc.pl.umap(
    rna,
    color="cell_type_broad",
    legend_loc="on data",
    legend_fontsize=7,
    legend_fontweight="normal",
    frameon=False,
    title="RNA reviewed broad annotation",
    show=False,
)
save_show_close(umap_broad_out, dpi=300, show=True)

sc.pl.umap(
    rna,
    color="cell_type_detailed",
    legend_loc="right margin",
    frameon=False,
    title="RNA reviewed detailed annotation",
    show=False,
)
save_show_close(umap_detailed_out, dpi=300, show=True)

sc.pl.umap(
    rna,
    color="include_in_downstream",
    legend_loc="on data",
    legend_fontsize=7,
    legend_fontweight="normal",
    frameon=False,
    title="RNA clusters included in downstream analysis",
    show=False,
)
save_show_close(umap_include_out, dpi=300, show=True)

sc.pl.umap(
    rna,
    color=[cluster_key, "cell_type_broad"],
    legend_loc="on data",
    legend_fontsize=6,
    legend_fontweight="normal",
    frameon=False,
    title=["Selected Leiden clusters", "Reviewed broad annotation"],
    show=False,
)
save_show_close(umap_cluster_and_annotation_out, dpi=300, show=True)


# ============================================================
# 10. Manifest
# ============================================================

manifest = {
    "step": "03B",
    "description": "Reviewed cluster annotation and downstream inclusion/exclusion decision",
    "timestamp": datetime.now().isoformat(timespec="seconds"),
    "project": paths.as_dict(),
    "inputs": {
        "rna_clustered_h5ad": str(rna_in),
        "reviewed_annotation_table": str(reviewed_annotation_path),
        "preprocessing_config": str(preprocessing_config_path),
    },
    "clustering_source": {
        "cluster_key": cluster_key,
        "n_pcs": int(dim_cfg["n_pcs"]),
        "n_neighbors": int(dim_cfg["n_neighbors"]),
        "leiden_resolution": float(dim_cfg["leiden_resolution"]),
    },
    "outputs": {
        "annotated_reviewed_h5ad": str(annotated_out),
        "selected_high_confidence_h5ad": str(selected_out),
        "annotation_summary": str(annotation_summary_out),
        "cell_counts": str(cell_counts_out),
        "selected_counts": str(selected_counts_out),
        "umap_broad": str(umap_broad_out),
        "umap_detailed": str(umap_detailed_out),
        "umap_include": str(umap_include_out),
        "umap_cluster_and_annotation": str(umap_cluster_and_annotation_out),
    },
    "summary": {
        "n_cells_full": int(rna.n_obs),
        "n_cells_selected": int(rna_selected.n_obs),
        "selected_fraction": float(rna_selected.n_obs / rna.n_obs),
        "n_clusters_full": int(rna.obs[cluster_key].nunique()),
        "n_clusters_selected": int(rna_selected.obs[cluster_key].nunique()),
        "excluded_cells": int(rna.n_obs - rna_selected.n_obs),
    },
    "rules": {
        "full_annotated_object_keeps_all_clusters": True,
        "selected_object_keeps_only_include_in_downstream_yes": True,
        "ambiguous_or_small_clusters_are_not_deleted_from_full_object": True,
    },
}

write_json(manifest, manifest_out)

print("\nSaved Step 03B manifest to:")
print(manifest_out)

print("\nDone.")
print("Step 03B reviewed annotation completed.")
print("Full object keeps all clusters.")
print("Selected subset contains only include_in_downstream == yes.")

from pathlib import Path  # paths

print(Path.cwd())  # δείχνει τον current working directory

import pandas as pd  # χειρισμός CSV table

reviewed = pd.read_csv(
    "reports/tables/rna/rna_reviewed_cluster_annotation.csv",
    dtype={"cluster": str}
)  # διάβασε το manual reviewed annotation table

kept = reviewed[
    reviewed["include_in_downstream"].astype(str).str.lower().isin(["yes", "true", "1"])
]  # κράτα μόνο clusters για downstream

excluded = reviewed[
    ~reviewed["include_in_downstream"].astype(str).str.lower().isin(["yes", "true", "1"])
]  # βρες clusters που αποκλείονται

print("Kept downstream clusters:")
print(kept[["cluster", "cell_type_broad", "cell_type_detailed", "annotation_confidence"]])

print("\nExcluded clusters:")
print(excluded[["cluster", "cell_type_broad", "exclusion_reason", "annotation_confidence"]])
