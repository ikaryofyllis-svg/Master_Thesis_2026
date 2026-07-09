# -*- coding: utf-8 -*-
"""
Step 04D2 — focused scATAC selected-cluster UMAP using LSI.

Purpose:
- Load original clustered ATAC object from Step 04B.
- Select cells from chosen original ATAC Leiden clusters.
- Recompute neighbors/UMAP on the selected subset using ATAC LSI representation.
- Compare multiple n_neighbors settings, e.g. 15 vs 30.
- Optionally copy marker score columns from Step 04D gene activity evidence object.
- Save subset AnnData, audit tables, comparison figures, and manifest.

Important:
- This is visualization / diagnostic only.
- It does not replace the global ATAC clustering.
- It does not create final cell type annotation.
- UMAP is recomputed from scATAC LSI, not from gene activity X.
- Gene activity marker scores are copied only for interpretation.
- All outputs must stay inside the active project.
"""

from pathlib import Path
import sys
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


# ============================================================
# 0. Project-aware setup
# ============================================================

try:
    refactor_dir = Path(__file__).resolve().parents[1]
except NameError:
    refactor_dir = Path.cwd().resolve()

src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

from pbmcgrn.config import load_yaml
from pbmcgrn.io import write_csv, write_json
from pbmcgrn.paths import ProjectPaths


def assert_path_inside_project(path: Path, project_root: Path, label: str) -> None:
    """Fail if path is outside active project."""

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


def resolve_project_path(paths: ProjectPaths, relative_path: str, label: str) -> Path:
    """Resolve a YAML relative path inside the active project."""

    out_path = paths.prepare_output_path(relative_path, label=label)

    assert_path_inside_project(out_path, paths.project_root, label)

    return out_path


def save_show_close(path: Path, dpi: int = 300, show: bool = True) -> None:
    """Save current matplotlib figure and optionally show it."""

    plt.savefig(path, dpi=dpi, bbox_inches="tight")

    print("\nSaved figure:")
    print(path)

    if show:
        plt.show(block=True)

    plt.close()


def sort_cluster_values(values) -> list[str]:
    """Sort cluster labels numerically when possible."""

    unique_values = pd.Series(values).astype(str).unique().tolist()

    def key_fn(x):
        try:
            return int(x)
        except ValueError:
            return x

    return sorted(unique_values, key=key_fn)


def copy_gene_activity_evidence_obs(
    atac_subset,
    gene_activity_evidence,
    prefix: str,
    max_cols: int,
):
    """Copy compatible marker-score obs columns from gene activity evidence object."""

    common_cells = atac_subset.obs_names.intersection(gene_activity_evidence.obs_names)

    if len(common_cells) == 0:
        print("\nWarning: no common cells between ATAC subset and gene activity evidence object.")
        return atac_subset, []

    score_cols = [
        col for col in gene_activity_evidence.obs.columns
        if col.startswith(prefix)
    ]

    score_cols = score_cols[:max_cols]

    for col in score_cols:
        atac_subset.obs[col] = np.nan
        atac_subset.obs.loc[common_cells, col] = gene_activity_evidence.obs.loc[common_cells, col].values

    print("\nCopied marker score columns from gene activity evidence:")
    print(score_cols)

    return atac_subset, score_cols


def plot_umap_with_cluster_centers(
    adata,
    umap_key: str,
    color_key: str,
    title: str,
    output_path: Path,
    dpi: int,
    show: bool,
    point_size: float,
    alpha: float,
    label_font_size: int,
) -> None:
    """Plot a UMAP embedding colored by original cluster with cluster labels at medians."""

    if umap_key not in adata.obsm.keys():
        raise KeyError(f"{umap_key} not found in adata.obsm.")

    if color_key not in adata.obs.columns:
        raise KeyError(f"{color_key} not found in adata.obs.")

    plot_df = pd.DataFrame(
        adata.obsm[umap_key],
        index=adata.obs_names,
        columns=["UMAP1", "UMAP2"],
    )

    plot_df[color_key] = adata.obs[color_key].astype(str).values

    cluster_codes = pd.Categorical(plot_df[color_key]).codes

    centers = (
        plot_df
        .groupby(color_key, observed=True)[["UMAP1", "UMAP2"]]
        .median()
        .reset_index()
    )

    plt.figure(figsize=(8, 7))

    plt.scatter(
        plot_df["UMAP1"],
        plot_df["UMAP2"],
        c=cluster_codes,
        s=point_size,
        alpha=alpha,
        cmap="tab20",
        linewidths=0,
    )

    for _, row in centers.iterrows():
        plt.text(
            row["UMAP1"],
            row["UMAP2"],
            str(row[color_key]),
            fontsize=label_font_size,
            fontweight="bold",
            ha="center",
            va="center",
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "edgecolor": "black",
                "alpha": 0.8,
                "linewidth": 0.5,
            },
        )

    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title(title, fontsize=12, fontweight="bold")
    plt.tight_layout()

    save_show_close(output_path, dpi=dpi, show=show)


def plot_umap_comparison(
    adata,
    umap_keys: list[str],
    titles: list[str],
    color_key: str,
    output_path: Path,
    dpi: int,
    show: bool,
    point_size: float,
    alpha: float,
    label_font_size: int,
) -> None:
    """Plot multiple UMAP embeddings side-by-side."""

    n = len(umap_keys)

    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5))

    if n == 1:
        axes = [axes]

    color_values = adata.obs[color_key].astype(str)
    cluster_codes = pd.Categorical(color_values).codes

    for ax, umap_key, title in zip(axes, umap_keys, titles):
        if umap_key not in adata.obsm.keys():
            raise KeyError(f"{umap_key} not found in adata.obsm.")

        plot_df = pd.DataFrame(
            adata.obsm[umap_key],
            index=adata.obs_names,
            columns=["UMAP1", "UMAP2"],
        )

        plot_df[color_key] = color_values.values

        ax.scatter(
            plot_df["UMAP1"],
            plot_df["UMAP2"],
            c=cluster_codes,
            s=point_size,
            alpha=alpha,
            cmap="tab20",
            linewidths=0,
        )

        centers = (
            plot_df
            .groupby(color_key, observed=True)[["UMAP1", "UMAP2"]]
            .median()
            .reset_index()
        )

        for _, row in centers.iterrows():
            ax.text(
                row["UMAP1"],
                row["UMAP2"],
                str(row[color_key]),
                fontsize=label_font_size,
                fontweight="bold",
                ha="center",
                va="center",
                bbox={
                    "boxstyle": "round,pad=0.25",
                    "facecolor": "white",
                    "edgecolor": "black",
                    "alpha": 0.8,
                    "linewidth": 0.5,
                },
            )

        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")

    plt.tight_layout()

    save_show_close(output_path, dpi=dpi, show=show)


# ============================================================
# 1. Initialize project
# ============================================================

paths = ProjectPaths(refactor_root=refactor_dir)
paths.ensure_output_dirs()

print("\nActive project:")
print(paths.project_id)

print("\nProject root:")
print(paths.project_root)


# ============================================================
# 2. Load config
# ============================================================

config_path = paths.configs / "atac_focused_umap_04D2.yaml"

assert_path_inside_project(
    config_path,
    paths.project_root,
    "atac_focused_umap_04D2.yaml",
)

if not config_path.exists():
    raise FileNotFoundError(f"Config not found:\n{config_path}")

config = load_yaml(config_path)

input_cfg = config["input"]
cluster_cfg = config["cluster"]
representation_cfg = config["representation"]
neighbors_cfg = config["neighbors"]
neighbors_grid = config["neighbors_grid"]
umap_cfg = config["umap"]
clustering_cfg = config["clustering"]
plots_cfg = config["plots"]
outputs_cfg = config["outputs"]

print("\nLoaded config:")
print(config_path)


# ============================================================
# 3. Resolve paths
# ============================================================

atac_clustered_in = paths.project_root / input_cfg["atac_clustered_h5ad"]

gene_activity_evidence_in = paths.project_root / input_cfg["gene_activity_evidence_h5ad"]

subset_h5ad_out = resolve_project_path(
    paths,
    outputs_cfg["subset_h5ad"],
    "subset_h5ad",
)

subset_cell_table_out = resolve_project_path(
    paths,
    outputs_cfg["subset_cell_table"],
    "subset_cell_table",
)

subset_summary_out = resolve_project_path(
    paths,
    outputs_cfg["subset_summary"],
    "subset_summary",
)

neighbor_run_summary_out = resolve_project_path(
    paths,
    outputs_cfg["neighbor_run_summary"],
    "neighbor_run_summary",
)

original_vs_recomputed_umap_out = resolve_project_path(
    paths,
    outputs_cfg["original_vs_recomputed_umap"],
    "original_vs_recomputed_umap",
)

recomputed_umap_neighbors_15_out = resolve_project_path(
    paths,
    outputs_cfg["recomputed_umap_neighbors_15"],
    "recomputed_umap_neighbors_15",
)

recomputed_umap_neighbors_30_out = resolve_project_path(
    paths,
    outputs_cfg["recomputed_umap_neighbors_30"],
    "recomputed_umap_neighbors_30",
)

recomputed_umap_neighbors_comparison_out = resolve_project_path(
    paths,
    outputs_cfg["recomputed_umap_neighbors_comparison"],
    "recomputed_umap_neighbors_comparison",
)

marker_score_umap_neighbors_15_out = resolve_project_path(
    paths,
    outputs_cfg["marker_score_umap_neighbors_15"],
    "marker_score_umap_neighbors_15",
)

marker_score_umap_neighbors_30_out = resolve_project_path(
    paths,
    outputs_cfg["marker_score_umap_neighbors_30"],
    "marker_score_umap_neighbors_30",
)

manifest_out = resolve_project_path(
    paths,
    outputs_cfg["manifest"],
    "step04D2_manifest",
)

for label, path in {
    "atac_clustered_in": atac_clustered_in,
    "gene_activity_evidence_in": gene_activity_evidence_in,
    "subset_h5ad_out": subset_h5ad_out,
    "subset_cell_table_out": subset_cell_table_out,
    "subset_summary_out": subset_summary_out,
    "neighbor_run_summary_out": neighbor_run_summary_out,
    "original_vs_recomputed_umap_out": original_vs_recomputed_umap_out,
    "recomputed_umap_neighbors_15_out": recomputed_umap_neighbors_15_out,
    "recomputed_umap_neighbors_30_out": recomputed_umap_neighbors_30_out,
    "recomputed_umap_neighbors_comparison_out": recomputed_umap_neighbors_comparison_out,
    "marker_score_umap_neighbors_15_out": marker_score_umap_neighbors_15_out,
    "marker_score_umap_neighbors_30_out": marker_score_umap_neighbors_30_out,
    "manifest_out": manifest_out,
}.items():
    assert_path_inside_project(path, paths.project_root, label)

for path in [
    subset_h5ad_out,
    subset_cell_table_out,
    subset_summary_out,
    neighbor_run_summary_out,
    original_vs_recomputed_umap_out,
    recomputed_umap_neighbors_15_out,
    recomputed_umap_neighbors_30_out,
    recomputed_umap_neighbors_comparison_out,
    marker_score_umap_neighbors_15_out,
    marker_score_umap_neighbors_30_out,
    manifest_out,
]:
    path.parent.mkdir(parents=True, exist_ok=True)

if not atac_clustered_in.exists():
    raise FileNotFoundError(f"ATAC clustered object not found:\n{atac_clustered_in}")

if not gene_activity_evidence_in.exists():
    print(
        "\nWarning: gene activity evidence object not found. "
        "Marker score transfer will be skipped."
    )
    gene_activity_evidence_in = None


# ============================================================
# 4. Load ATAC clustered object
# ============================================================

atac = sc.read_h5ad(atac_clustered_in)

print("\nLoaded ATAC clustered object:")
print(atac)

cluster_key = cluster_cfg.get("key", "leiden_0_6")

if cluster_key not in atac.obs.columns:
    raise KeyError(
        f"Cluster key '{cluster_key}' not found in ATAC obs.\n"
        f"Available obs columns: {atac.obs.columns.tolist()}"
    )

use_rep = representation_cfg.get("use_rep", "X_lsi_2_50")

if use_rep not in atac.obsm.keys():
    raise KeyError(
        f"Requested representation '{use_rep}' not found in atac.obsm.\n"
        f"Available obsm keys: {list(atac.obsm.keys())}"
    )

atac.obs[cluster_key] = atac.obs[cluster_key].astype(str)

selected_clusters = [str(x) for x in cluster_cfg["selected_clusters"]]

print("\nSelected original ATAC clusters:")
print(selected_clusters)


# ============================================================
# 5. Subset selected cells
# ============================================================

cell_mask = atac.obs[cluster_key].isin(selected_clusters)

n_cells_before = int(atac.n_obs)
n_cells_after = int(cell_mask.sum())

if n_cells_after == 0:
    raise ValueError(
        "No cells matched selected clusters.\n"
        f"Selected clusters: {selected_clusters}\n"
        f"Available clusters: {sort_cluster_values(atac.obs[cluster_key])}"
    )

atac_subset = atac[cell_mask].copy()

atac_subset.obs["original_atac_cluster"] = atac_subset.obs[cluster_key].astype(str)

print("\nSubset ATAC object:")
print(atac_subset)

print("\nCells before:", n_cells_before)
print("Cells after:", n_cells_after)
print("Retention fraction:", n_cells_after / n_cells_before)


# ============================================================
# 6. Copy original global UMAP
# ============================================================

has_original_umap = "X_umap" in atac_subset.obsm.keys()

if has_original_umap:
    atac_subset.obsm["X_umap_global_04B"] = atac_subset.obsm["X_umap"].copy()
    print("\nCopied original/global UMAP to obsm['X_umap_global_04B'].")

else:
    print("\nWarning: original X_umap not found. Original-vs-new comparison will skip global panel.")


# ============================================================
# 7. Optional: copy gene activity marker score columns
# ============================================================

marker_score_columns = []

if gene_activity_evidence_in is not None:
    gene_activity_evidence = sc.read_h5ad(gene_activity_evidence_in)

    score_prefix = plots_cfg.get("marker_score_columns_prefix", "ga_score_")
    max_score_cols = int(plots_cfg.get("max_marker_score_columns", 12))

    atac_subset, marker_score_columns = copy_gene_activity_evidence_obs(
        atac_subset=atac_subset,
        gene_activity_evidence=gene_activity_evidence,
        prefix=score_prefix,
        max_cols=max_score_cols,
    )


# ============================================================
# 8. Save audit tables
# ============================================================

subset_cell_table = atac_subset.obs.copy()
subset_cell_table.insert(0, "cell_id", atac_subset.obs_names)

write_csv(subset_cell_table, subset_cell_table_out, index=False)

subset_summary = (
    atac_subset.obs
    .groupby("original_atac_cluster", observed=True)
    .size()
    .reset_index(name="n_cells")
)

subset_summary["fraction_of_subset"] = subset_summary["n_cells"] / atac_subset.n_obs

subset_summary = subset_summary.sort_values(
    "original_atac_cluster",
    key=lambda s: s.astype(int),
)

write_csv(subset_summary, subset_summary_out, index=False)

print("\nSubset summary:")
print(subset_summary)


# ============================================================
# 9. Recompute neighbors/UMAP for each n_neighbors
# ============================================================

metric = neighbors_cfg.get("metric", "cosine")
random_state = int(umap_cfg.get("random_state", 42))
run_leiden = bool(clustering_cfg.get("run_leiden", True))
leiden_resolution = float(clustering_cfg.get("leiden_resolution", 0.5))
leiden_key_template = clustering_cfg.get(
    "leiden_key_template",
    "leiden_04D2_lsi_neighbors_{n_neighbors}_res_0_5",
)

neighbor_run_records = []

for n_neighbors in neighbors_grid:
    n_neighbors = int(n_neighbors)

    print("\nRecomputing neighbors/UMAP with:")
    print("use_rep:", use_rep)
    print("n_neighbors:", n_neighbors)
    print("metric:", metric)

    sc.pp.neighbors(
        atac_subset,
        n_neighbors=n_neighbors,
        use_rep=use_rep,
        metric=metric,
        key_added=f"neighbors_lsi_{n_neighbors}",
    )

    sc.tl.umap(
        atac_subset,
        neighbors_key=f"neighbors_lsi_{n_neighbors}",
        random_state=random_state,
    )

    umap_key = f"X_umap_lsi_neighbors_{n_neighbors}"

    atac_subset.obsm[umap_key] = atac_subset.obsm["X_umap"].copy()

    leiden_key = None
    n_leiden_clusters = None

    if run_leiden:
        leiden_key = leiden_key_template.format(n_neighbors=n_neighbors)

        sc.tl.leiden(
            atac_subset,
            resolution=leiden_resolution,
            key_added=leiden_key,
            neighbors_key=f"neighbors_lsi_{n_neighbors}",
        )

        n_leiden_clusters = int(atac_subset.obs[leiden_key].nunique())

        print(f"\nLeiden counts for {leiden_key}:")
        print(atac_subset.obs[leiden_key].value_counts().sort_index())

    neighbor_run_records.append(
        {
            "n_neighbors": n_neighbors,
            "use_rep": use_rep,
            "metric": metric,
            "umap_key": umap_key,
            "neighbors_key": f"neighbors_lsi_{n_neighbors}",
            "run_leiden": run_leiden,
            "leiden_key": leiden_key,
            "leiden_resolution": leiden_resolution if run_leiden else None,
            "n_leiden_clusters": n_leiden_clusters,
        }
    )

neighbor_run_summary = pd.DataFrame(neighbor_run_records)

write_csv(neighbor_run_summary, neighbor_run_summary_out, index=False)

print("\nNeighbor run summary:")
print(neighbor_run_summary)


# ============================================================
# 10. Plots
# ============================================================

plot_dpi = int(plots_cfg.get("dpi", 300))
plot_show = bool(plots_cfg.get("show", True))
point_size = float(plots_cfg.get("point_size", 6))
alpha = float(plots_cfg.get("alpha", 0.75))
label_font_size = int(plots_cfg.get("label_font_size", 10))

# Individual UMAPs
plot_umap_with_cluster_centers(
    adata=atac_subset,
    umap_key="X_umap_lsi_neighbors_15",
    color_key="original_atac_cluster",
    title="Selected ATAC clusters — LSI UMAP, n_neighbors=15",
    output_path=recomputed_umap_neighbors_15_out,
    dpi=plot_dpi,
    show=plot_show,
    point_size=point_size,
    alpha=alpha,
    label_font_size=label_font_size,
)

plot_umap_with_cluster_centers(
    adata=atac_subset,
    umap_key="X_umap_lsi_neighbors_30",
    color_key="original_atac_cluster",
    title="Selected ATAC clusters — LSI UMAP, n_neighbors=30",
    output_path=recomputed_umap_neighbors_30_out,
    dpi=plot_dpi,
    show=plot_show,
    point_size=point_size,
    alpha=alpha,
    label_font_size=label_font_size,
)

# Side-by-side comparison: 15 vs 30
plot_umap_comparison(
    adata=atac_subset,
    umap_keys=["X_umap_lsi_neighbors_15", "X_umap_lsi_neighbors_30"],
    titles=["LSI UMAP n_neighbors=15", "LSI UMAP n_neighbors=30"],
    color_key="original_atac_cluster",
    output_path=recomputed_umap_neighbors_comparison_out,
    dpi=plot_dpi,
    show=plot_show,
    point_size=point_size,
    alpha=alpha,
    label_font_size=label_font_size,
)

# Original global vs recomputed comparison
if has_original_umap:
    plot_umap_comparison(
        adata=atac_subset,
        umap_keys=[
            "X_umap_global_04B",
            "X_umap_lsi_neighbors_15",
            "X_umap_lsi_neighbors_30",
        ],
        titles=[
            "Original global ATAC UMAP",
            "Subset LSI UMAP n_neighbors=15",
            "Subset LSI UMAP n_neighbors=30",
        ],
        color_key="original_atac_cluster",
        output_path=original_vs_recomputed_umap_out,
        dpi=plot_dpi,
        show=plot_show,
        point_size=point_size,
        alpha=alpha,
        label_font_size=label_font_size,
    )

# Marker score plots if available
if len(marker_score_columns) > 0:
    sc.pl.embedding(
        atac_subset,
        basis="umap_lsi_neighbors_15",
        color=marker_score_columns,
        ncols=3,
        cmap="viridis",
        frameon=False,
        show=False,
    )

    save_show_close(
        marker_score_umap_neighbors_15_out,
        dpi=plot_dpi,
        show=plot_show,
    )

    sc.pl.embedding(
        atac_subset,
        basis="umap_lsi_neighbors_30",
        color=marker_score_columns,
        ncols=3,
        cmap="viridis",
        frameon=False,
        show=False,
    )

    save_show_close(
        marker_score_umap_neighbors_30_out,
        dpi=plot_dpi,
        show=plot_show,
    )

else:
    print("\nNo marker score columns available. Marker-score UMAPs skipped.")


# ============================================================
# 11. Save subset object
# ============================================================

atac_subset.uns["step04D2_lsi_selected_cluster_umap_comparison"] = {
    "description": "Selected ATAC clusters with recomputed UMAP using ATAC LSI representation",
    "source_atac_clustered_h5ad": str(atac_clustered_in),
    "source_gene_activity_evidence_h5ad": str(gene_activity_evidence_in) if gene_activity_evidence_in is not None else None,
    "cluster_key": cluster_key,
    "selected_clusters": selected_clusters,
    "n_cells_before": n_cells_before,
    "n_cells_after": n_cells_after,
    "retention_fraction": n_cells_after / n_cells_before,
    "use_rep": use_rep,
    "neighbors_grid": [int(x) for x in neighbors_grid],
    "metric": metric,
    "random_state": random_state,
    "run_leiden": run_leiden,
    "leiden_resolution": leiden_resolution if run_leiden else None,
    "marker_score_columns": marker_score_columns,
    "note": (
        "Visualization-only focused UMAP. "
        "UMAPs are based on scATAC LSI representation, not gene activity X. "
        "Do not use as final annotation without marker evidence and RNA-ATAC matching."
    ),
}

atac_subset.write_h5ad(subset_h5ad_out)

print("\nSaved focused LSI subset object:")
print(subset_h5ad_out)


# ============================================================
# 12. Manifest
# ============================================================

manifest = {
    "step": "04D2",
    "description": "Focused scATAC selected-cluster UMAP comparison using LSI representation",
    "timestamp": datetime.now().isoformat(timespec="seconds"),
    "project": paths.as_dict(),
    "inputs": {
        "atac_clustered_h5ad": str(atac_clustered_in),
        "gene_activity_evidence_h5ad": str(gene_activity_evidence_in) if gene_activity_evidence_in is not None else None,
        "config": str(config_path),
    },
    "settings": {
        "cluster_key": cluster_key,
        "selected_clusters": selected_clusters,
        "n_cells_before": n_cells_before,
        "n_cells_after": n_cells_after,
        "retention_fraction": n_cells_after / n_cells_before,
        "use_rep": use_rep,
        "neighbors_grid": [int(x) for x in neighbors_grid],
        "metric": metric,
        "random_state": random_state,
        "run_leiden": run_leiden,
        "leiden_resolution": leiden_resolution if run_leiden else None,
        "has_original_umap": has_original_umap,
        "n_marker_score_columns_copied": len(marker_score_columns),
        "marker_score_columns": marker_score_columns,
    },
    "outputs": {
        "subset_h5ad": str(subset_h5ad_out),
        "subset_cell_table": str(subset_cell_table_out),
        "subset_summary": str(subset_summary_out),
        "neighbor_run_summary": str(neighbor_run_summary_out),
        "original_vs_recomputed_umap": str(original_vs_recomputed_umap_out) if has_original_umap else None,
        "recomputed_umap_neighbors_15": str(recomputed_umap_neighbors_15_out),
        "recomputed_umap_neighbors_30": str(recomputed_umap_neighbors_30_out),
        "recomputed_umap_neighbors_comparison": str(recomputed_umap_neighbors_comparison_out),
        "marker_score_umap_neighbors_15": str(marker_score_umap_neighbors_15_out) if len(marker_score_columns) > 0 else None,
        "marker_score_umap_neighbors_30": str(marker_score_umap_neighbors_30_out) if len(marker_score_columns) > 0 else None,
        "manifest": str(manifest_out),
    },
}

write_json(manifest, manifest_out)

print("\nSaved Step 04D2 manifest:")
print(manifest_out)

print("\nStep 04D2 complete.")