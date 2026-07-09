# -*- coding: utf-8 -*-
"""
Step 04D — scATAC marker-based annotation evidence.

Purpose:
- Load scATAC gene activity object from Step 04C.
- Evaluate known PBMC marker genes across ATAC clusters.
- Compute marker availability.
- Compute mean gene activity per cluster.
- Compute percent-positive raw gene activity per cluster.
- Compute marker-set/module scores.
- Suggest likely ATAC cluster identities from gene activity evidence.
- Save tables, figures, marker-evidence AnnData, and manifest.

Important:
- This is an evidence/review step, not final annotation.
- ATAC marker evidence should support RNA–ATAC matching.
- It does not replace RNA-based annotation.
- All outputs must stay inside the active project.
"""

from pathlib import Path
import sys
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from scipy import sparse


# ============================================================
# 0. Project-aware setup
# ============================================================

try:
    refactor_dir = Path(__file__).resolve().parents[1]  # Refactored root
except NameError:
    refactor_dir = Path.cwd().resolve()  # interactive fallback

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

    if show:
        plt.show()

    plt.close()


def flatten_marker_sets(marker_sets: dict) -> list[str]:
    """Flatten marker sets into unique marker list while preserving order."""

    seen = set()
    markers = []

    for _, genes in marker_sets.items():
        for gene in genes:
            gene = str(gene).strip()

            if gene and gene not in seen:
                seen.add(gene)
                markers.append(gene)

    return markers


def to_dense_1d(x) -> np.ndarray:
    """Convert dense/sparse vector-like matrix to 1D array."""

    if sparse.issparse(x):
        return np.asarray(x.toarray()).ravel()

    return np.asarray(x).ravel()


def sort_cluster_values(values: pd.Series) -> list[str]:
    """Sort cluster labels numerically when possible."""

    unique_values = values.astype(str).unique().tolist()

    def sort_key(x):
        try:
            return int(x)
        except ValueError:
            return x

    return sorted(unique_values, key=sort_key)


def compute_marker_availability(marker_sets: dict, var_names: pd.Index) -> pd.DataFrame:
    """Create marker availability table."""

    rows = []

    for marker_set, genes in marker_sets.items():
        for gene in genes:
            rows.append(
                {
                    "marker_set": marker_set,
                    "marker_gene": gene,
                    "available_in_gene_activity": gene in var_names,
                }
            )

    return pd.DataFrame(rows)


def compute_marker_mean_by_cluster(
    adata,
    cluster_key: str,
    markers: list[str],
) -> pd.DataFrame:
    """Compute mean normalized/log gene activity per cluster."""

    marker_df = sc.get.obs_df(
        adata,
        keys=[cluster_key] + markers,
    )

    marker_df[cluster_key] = marker_df[cluster_key].astype(str)

    mean_df = (
        marker_df
        .groupby(cluster_key, observed=True)[markers]
        .mean()
        .reset_index()
    )

    return mean_df


def compute_marker_pct_positive_by_cluster(
    adata,
    cluster_key: str,
    markers: list[str],
    layer: str | None,
    threshold: float,
) -> pd.DataFrame:
    """Compute fraction of cells with raw marker activity > threshold."""

    rows = []

    cluster_values = adata.obs[cluster_key].astype(str)

    for cluster in sort_cluster_values(cluster_values):
        cell_mask = cluster_values == cluster
        adata_cluster = adata[cell_mask.values, :]

        for gene in markers:
            gene_idx = adata_cluster.var_names.get_loc(gene)

            if layer is not None and layer in adata_cluster.layers.keys():
                values = adata_cluster.layers[layer][:, gene_idx]
            else:
                values = adata_cluster.X[:, gene_idx]

            values = to_dense_1d(values)

            rows.append(
                {
                    "cluster": cluster,
                    "marker_gene": gene,
                    "pct_positive": float(np.mean(values > threshold)),
                    "mean_value_used": float(np.mean(values)),
                    "n_cells": int(cell_mask.sum()),
                    "layer_used": layer if layer is not None else "X",
                    "threshold": threshold,
                }
            )

    long_df = pd.DataFrame(rows)

    wide_df = (
        long_df
        .pivot(index="cluster", columns="marker_gene", values="pct_positive")
        .reset_index()
    )

    return wide_df


def add_marker_set_scores(
    adata,
    marker_sets: dict,
    score_prefix: str,
) -> tuple[object, pd.DataFrame]:
    """Add marker-set scores to adata.obs."""

    rows = []

    for marker_set, genes in marker_sets.items():
        available = [gene for gene in genes if gene in adata.var_names]
        missing = [gene for gene in genes if gene not in adata.var_names]
        score_col = f"{score_prefix}{marker_set}"

        if len(available) >= 2:
            sc.tl.score_genes(
                adata,
                gene_list=available,
                score_name=score_col,
                use_raw=False,
            )
            status = "score_genes"

        elif len(available) == 1:
            gene = available[0]
            adata.obs[score_col] = to_dense_1d(adata[:, gene].X)
            status = "single_marker_copied"

        else:
            adata.obs[score_col] = np.nan
            status = "no_available_markers"

        rows.append(
            {
                "marker_set": marker_set,
                "score_column": score_col,
                "n_requested_markers": len(genes),
                "n_available_markers": len(available),
                "available_markers": ";".join(available),
                "missing_markers": ";".join(missing),
                "status": status,
            }
        )

    return adata, pd.DataFrame(rows)


def compute_score_means_by_cluster(
    adata,
    cluster_key: str,
    score_columns: list[str],
) -> pd.DataFrame:
    """Compute mean marker-set score per cluster."""

    score_df = adata.obs[[cluster_key] + score_columns].copy()

    score_df[cluster_key] = score_df[cluster_key].astype(str)

    summary = (
        score_df
        .groupby(cluster_key, observed=True)[score_columns]
        .mean()
        .reset_index()
    )

    return summary


def make_cluster_suggestions(
    score_summary: pd.DataFrame,
    cluster_key: str,
    score_prefix: str,
    n_top: int,
) -> pd.DataFrame:
    """Suggest likely labels from marker-set score rankings."""

    score_cols = [
        col for col in score_summary.columns
        if col.startswith(score_prefix)
    ]

    rows = []

    for _, row in score_summary.iterrows():
        cluster = str(row[cluster_key])

        scores = row[score_cols].astype(float).sort_values(ascending=False)

        top = scores.head(n_top)

        top_labels = [
            col.replace(score_prefix, "")
            for col in top.index
        ]

        rows.append(
            {
                "cluster": cluster,
                "suggested_label": top_labels[0] if top_labels else "Unknown",
                "top_labels": ";".join(top_labels),
                "top_scores": ";".join([str(round(float(v), 4)) for v in top.values]),
                "n_top_labels": len(top_labels),
                "note": "Suggestion from ATAC gene activity marker-set scores; manual review required.",
            }
        )

    return pd.DataFrame(rows)

def plot_umap_with_cluster_suggestions(
    adata,
    cluster_key: str,
    suggestions: pd.DataFrame,
    output_path: Path,
    dpi: int = 300,
    show: bool = True,
    font_size: int = 10,
) -> None:
    """Plot UMAP with each ATAC cluster labeled by suggested cell type."""

    if "X_umap" not in adata.obsm.keys():
        print("\nWarning: X_umap not found. Cannot make cluster suggestion UMAP.")
        return

    suggestion_map = (
        suggestions
        .assign(cluster=suggestions["cluster"].astype(str))
        .set_index("cluster")["suggested_label"]
        .to_dict()
    )

    obs = adata.obs.copy()
    obs[cluster_key] = obs[cluster_key].astype(str)

    umap_df = pd.DataFrame(
        adata.obsm["X_umap"],
        index=adata.obs_names,
        columns=["UMAP1", "UMAP2"],
    )

    plot_df = pd.concat([umap_df, obs[[cluster_key]]], axis=1)

    cluster_centers = (
        plot_df
        .groupby(cluster_key, observed=True)[["UMAP1", "UMAP2"]]
        .median()
        .reset_index()
    )

    cluster_centers["suggested_label"] = (
        cluster_centers[cluster_key]
        .map(suggestion_map)
        .fillna("Unknown")
    )

    plt.figure(figsize=(10, 8))

    cluster_codes = pd.Categorical(plot_df[cluster_key]).codes

    plt.scatter(
        plot_df["UMAP1"],
        plot_df["UMAP2"],
        c=cluster_codes,
        s=5,
        alpha=0.7,
        cmap="tab20",
        linewidths=0,
    )

    for _, row in cluster_centers.iterrows():
        label = f"{row[cluster_key]}: {row['suggested_label']}"

        plt.text(
            row["UMAP1"],
            row["UMAP2"],
            label,
            fontsize=font_size,
            fontweight="bold",
            ha="center",
            va="center",
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "edgecolor": "black",
                "alpha": 0.85,
                "linewidth": 0.6,
            },
        )

    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title(
        "ATAC clusters with suggested cell types",
        fontsize=12,
        fontweight="bold",
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")

    print("\nSaved cluster suggestion UMAP:")
    print(output_path)

    if show:
        plt.show(block=True)

    plt.close()
# ============================================================
# 1. Initialize paths
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

config_path = paths.configs / "atac_marker_annotation.yaml"

assert_path_inside_project(
    config_path,
    paths.project_root,
    "atac_marker_annotation.yaml",
)

if not config_path.exists():
    raise FileNotFoundError(f"Config not found:\n{config_path}")

config = load_yaml(config_path)

input_cfg = config["input"]
cluster_cfg = config["cluster"]
markers_cfg = config["markers"]
scoring_cfg = config["scoring"]
plots_cfg = config["plots"]
outputs_cfg = config["outputs"]

print("\nLoaded config:")
print(config_path)


# ============================================================
# 3. Resolve paths
# ============================================================

gene_activity_in = paths.project_root / input_cfg["gene_activity_h5ad"]

marker_availability_out = resolve_project_path(
    paths,
    outputs_cfg["marker_availability"],
    "marker_availability",
)

marker_mean_by_cluster_out = resolve_project_path(
    paths,
    outputs_cfg["marker_mean_by_cluster"],
    "marker_mean_by_cluster",
)

marker_pct_positive_by_cluster_out = resolve_project_path(
    paths,
    outputs_cfg["marker_pct_positive_by_cluster"],
    "marker_pct_positive_by_cluster",
)

marker_set_availability_out = resolve_project_path(
    paths,
    outputs_cfg["marker_set_availability"],
    "marker_set_availability",
)

marker_module_score_by_cluster_out = resolve_project_path(
    paths,
    outputs_cfg["marker_module_score_by_cluster"],
    "marker_module_score_by_cluster",
)

cluster_annotation_suggestions_out = resolve_project_path(
    paths,
    outputs_cfg["cluster_annotation_suggestions"],
    "cluster_annotation_suggestions",
)

marker_dotplot_out = resolve_project_path(
    paths,
    outputs_cfg["marker_dotplot"],
    "marker_dotplot",
)

marker_umap_out = resolve_project_path(
    paths,
    outputs_cfg["marker_umap"],
    "marker_umap",
)

module_score_umap_out = resolve_project_path(
    paths,
    outputs_cfg["module_score_umap"],
    "module_score_umap",
)

module_score_heatmap_out = resolve_project_path(
    paths,
    outputs_cfg["module_score_heatmap"],
    "module_score_heatmap",
)
cluster_suggestion_umap_out = resolve_project_path(
    paths,
    outputs_cfg.get(
        "cluster_suggestion_umap",
        "reports/figures/atac/annotation/atac_umap_cluster_suggested_celltypes.png",
    ),
    "cluster_suggestion_umap",
)

evidence_h5ad_out = resolve_project_path(
    paths,
    outputs_cfg["evidence_h5ad"],
    "evidence_h5ad",
)

manifest_out = resolve_project_path(
    paths,
    outputs_cfg["manifest"],
    "step04D_manifest",
)

for label, path in {
    "gene_activity_in": gene_activity_in,
    "marker_availability_out": marker_availability_out,
    "marker_mean_by_cluster_out": marker_mean_by_cluster_out,
    "marker_pct_positive_by_cluster_out": marker_pct_positive_by_cluster_out,
    "marker_set_availability_out": marker_set_availability_out,
    "marker_module_score_by_cluster_out": marker_module_score_by_cluster_out,
    "cluster_annotation_suggestions_out": cluster_annotation_suggestions_out,
    "marker_dotplot_out": marker_dotplot_out,
    "marker_umap_out": marker_umap_out,
    "module_score_umap_out": module_score_umap_out,
    "module_score_heatmap_out": module_score_heatmap_out,
    "evidence_h5ad_out": evidence_h5ad_out,
    "manifest_out": manifest_out,
    "cluster_suggestion_umap_out": cluster_suggestion_umap_out,
}.items():
    assert_path_inside_project(path, paths.project_root, label)

for path in [
    marker_availability_out,
    marker_mean_by_cluster_out,
    marker_pct_positive_by_cluster_out,
    marker_set_availability_out,
    marker_module_score_by_cluster_out,
    cluster_annotation_suggestions_out,
    marker_dotplot_out,
    marker_umap_out,
    module_score_umap_out,
    module_score_heatmap_out,
    evidence_h5ad_out,
    manifest_out,
    cluster_suggestion_umap_out,
]:
    path.parent.mkdir(parents=True, exist_ok=True)

if not gene_activity_in.exists():
    raise FileNotFoundError(f"Gene activity object not found:\n{gene_activity_in}")


# ============================================================
# 4. Load gene activity object
# ============================================================

adata = sc.read_h5ad(gene_activity_in)

print("\nLoaded ATAC gene activity object:")
print(adata)

cluster_key = cluster_cfg.get("key", "leiden_0_6")

if cluster_key not in adata.obs.columns:
    raise KeyError(
        f"Cluster key not found in adata.obs: {cluster_key}\n"
        f"Available obs columns: {adata.obs.columns.tolist()}"
    )

adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)

has_umap = "X_umap" in adata.obsm.keys()

if not has_umap:
    print("\nWarning: X_umap not found. UMAP plots will be skipped.")


# ============================================================
# 5. Marker availability
# ============================================================

marker_sets = markers_cfg["marker_sets"]

all_markers = flatten_marker_sets(marker_sets)

marker_availability = compute_marker_availability(
    marker_sets=marker_sets,
    var_names=adata.var_names,
)

write_csv(marker_availability, marker_availability_out, index=False)

available_markers = [
    marker for marker in all_markers
    if marker in adata.var_names
]

missing_markers = [
    marker for marker in all_markers
    if marker not in adata.var_names
]

print("\nRequested unique markers:", len(all_markers))
print("Available unique markers:", len(available_markers))
print("Missing unique markers:", len(missing_markers))

print("\nAvailable markers:")
print(available_markers)

print("\nMissing markers:")
print(missing_markers)

if len(available_markers) == 0:
    raise ValueError(
        "No requested marker genes were found in gene_activity.var_names. "
        "Check gene symbols, genome annotation, and peak-to-gene annotation."
    )


# ============================================================
# 6. Marker mean activity by cluster
# ============================================================

marker_mean_by_cluster = compute_marker_mean_by_cluster(
    adata=adata,
    cluster_key=cluster_key,
    markers=available_markers,
)

write_csv(marker_mean_by_cluster, marker_mean_by_cluster_out, index=False)

print("\nMarker mean gene activity by cluster:")
print(marker_mean_by_cluster)


# ============================================================
# 7. Marker percent-positive by cluster
# ============================================================

raw_layer = scoring_cfg.get("raw_layer_for_pct_positive", "raw_gene_activity_counts")

if raw_layer is not None and raw_layer not in adata.layers.keys():
    print(
        f"\nWarning: layer '{raw_layer}' not found. "
        "Using adata.X for percent-positive calculations."
    )
    raw_layer = None

pct_positive_threshold = float(scoring_cfg.get("pct_positive_threshold", 0.0))

marker_pct_positive_by_cluster = compute_marker_pct_positive_by_cluster(
    adata=adata,
    cluster_key=cluster_key,
    markers=available_markers,
    layer=raw_layer,
    threshold=pct_positive_threshold,
)

write_csv(
    marker_pct_positive_by_cluster,
    marker_pct_positive_by_cluster_out,
    index=False,
)

print("\nMarker percent-positive by cluster:")
print(marker_pct_positive_by_cluster)


# ============================================================
# 8. Marker-set scores
# ============================================================

score_prefix = scoring_cfg.get("module_score_prefix", "ga_score_")

adata, marker_set_availability = add_marker_set_scores(
    adata=adata,
    marker_sets=marker_sets,
    score_prefix=score_prefix,
)

write_csv(marker_set_availability, marker_set_availability_out, index=False)

score_columns = marker_set_availability["score_column"].tolist()

marker_module_score_by_cluster = compute_score_means_by_cluster(
    adata=adata,
    cluster_key=cluster_key,
    score_columns=score_columns,
)

write_csv(
    marker_module_score_by_cluster,
    marker_module_score_by_cluster_out,
    index=False,
)

print("\nMarker-set score means by cluster:")
print(marker_module_score_by_cluster)


# ============================================================
# 9. Cluster annotation suggestions
# ============================================================

n_top_suggested_labels = int(scoring_cfg.get("n_top_suggested_labels", 3))

cluster_annotation_suggestions = make_cluster_suggestions(
    score_summary=marker_module_score_by_cluster,
    cluster_key=cluster_key,
    score_prefix=score_prefix,
    n_top=n_top_suggested_labels,
)

write_csv(
    cluster_annotation_suggestions,
    cluster_annotation_suggestions_out,
    index=False,
)

print("\nCluster annotation suggestions:")
print(cluster_annotation_suggestions)


# ============================================================
# 10. Plots
# ============================================================

plot_dpi = int(plots_cfg.get("dpi", 300))
plot_show = bool(plots_cfg.get("show", True))

cmap = plots_cfg.get("cmap", "viridis")
max_umap_markers = int(plots_cfg.get("max_umap_markers", 16))
print("\nDEBUG has_umap:", "X_umap" in adata.obsm.keys())
print("DEBUG suggestions columns:", cluster_annotation_suggestions.columns.tolist())
print("DEBUG suggestion UMAP output:", cluster_suggestion_umap_out)
print(cluster_annotation_suggestions.head())
if has_umap:
    plot_umap_with_cluster_suggestions(
        adata=adata,
        cluster_key=cluster_key,
        suggestions=cluster_annotation_suggestions,
        output_path=cluster_suggestion_umap_out,
        dpi=plot_dpi,
        show=plot_show,
        font_size=10,
    )

# Dotplot of marker genes
sc.pl.dotplot(
    adata,
    var_names=available_markers,
    groupby=cluster_key,
    standard_scale="var",
    show=False,
)

save_show_close(
    marker_dotplot_out,
    dpi=plot_dpi,
    show=plot_show,
)

# UMAP marker panels
if has_umap:
    umap_markers = available_markers[:max_umap_markers]

    sc.pl.umap(
        adata,
        color=[cluster_key] + umap_markers,
        ncols=4,
        cmap=cmap,
        frameon=False,
        show=False,
    )

    save_show_close(
        marker_umap_out,
        dpi=plot_dpi,
        show=plot_show,
    )

    sc.pl.umap(
        adata,
        color=[cluster_key] + score_columns,
        ncols=3,
        cmap=cmap,
        frameon=False,
        show=False,
    )

    save_show_close(
        module_score_umap_out,
        dpi=plot_dpi,
        show=plot_show,
    )

# Heatmap of marker-set scores
heatmap_df = marker_module_score_by_cluster.copy()
heatmap_df = heatmap_df.set_index(cluster_key)

heatmap_df.columns = [
    col.replace(score_prefix, "")
    for col in heatmap_df.columns
]

plt.figure(
    figsize=(
        max(8, 0.6 * heatmap_df.shape[1]),
        max(5, 0.35 * heatmap_df.shape[0]),
    )
)

plt.imshow(
    heatmap_df.values,
    aspect="auto",
    cmap="viridis",
)

plt.colorbar(label="Mean marker-set score")

plt.xticks(
    ticks=np.arange(heatmap_df.shape[1]),
    labels=heatmap_df.columns,
    rotation=45,
    ha="right",
)

plt.yticks(
    ticks=np.arange(heatmap_df.shape[0]),
    labels=heatmap_df.index,
)

plt.xlabel("Marker set")
plt.ylabel("ATAC cluster")
plt.title("ATAC gene activity marker-set scores by cluster")
plt.tight_layout()

save_show_close(
    module_score_heatmap_out,
    dpi=plot_dpi,
    show=plot_show,
)


# ============================================================
# 11. Save marker-evidence AnnData
# ============================================================

adata.uns["step04D_marker_annotation_evidence"] = {
    "description": "Marker-based ATAC cluster interpretation evidence from gene activity",
    "cluster_key": cluster_key,
    "n_requested_unique_markers": len(all_markers),
    "n_available_unique_markers": len(available_markers),
    "n_missing_unique_markers": len(missing_markers),
    "score_prefix": score_prefix,
    "raw_layer_for_pct_positive": raw_layer,
    "pct_positive_threshold": pct_positive_threshold,
    "note": (
        "This object contains marker-set scores as evidence only. "
        "It is not a final ATAC annotation object."
    ),
}

adata.write_h5ad(evidence_h5ad_out)

print("\nSaved marker-evidence gene activity object:")
print(evidence_h5ad_out)


# ============================================================
# 12. Manifest
# ============================================================

manifest = {
    "step": "04D",
    "description": "scATAC marker-based annotation evidence from gene activity",
    "timestamp": datetime.now().isoformat(timespec="seconds"),
    "project": paths.as_dict(),
    "inputs": {
        "gene_activity_h5ad": str(gene_activity_in),
        "config": str(config_path),
    },
    "settings": {
        "cluster_key": cluster_key,
        "n_marker_sets": len(marker_sets),
        "n_requested_unique_markers": len(all_markers),
        "n_available_unique_markers": len(available_markers),
        "n_missing_unique_markers": len(missing_markers),
        "raw_layer_for_pct_positive": raw_layer,
        "pct_positive_threshold": pct_positive_threshold,
        "score_prefix": score_prefix,
        "n_top_suggested_labels": n_top_suggested_labels,
        "has_umap": has_umap,
    },
    "outputs": {
        "marker_availability": str(marker_availability_out),
        "marker_mean_by_cluster": str(marker_mean_by_cluster_out),
        "marker_pct_positive_by_cluster": str(marker_pct_positive_by_cluster_out),
        "marker_set_availability": str(marker_set_availability_out),
        "marker_module_score_by_cluster": str(marker_module_score_by_cluster_out),
        "cluster_annotation_suggestions": str(cluster_annotation_suggestions_out),
        "marker_dotplot": str(marker_dotplot_out),
        "marker_umap": str(marker_umap_out) if has_umap else None,
        "module_score_umap": str(module_score_umap_out) if has_umap else None,
        "module_score_heatmap": str(module_score_heatmap_out),
        "evidence_h5ad": str(evidence_h5ad_out),
        "cluster_suggestion_umap": str(cluster_suggestion_umap_out) if has_umap else None,
    },
}

write_json(manifest, manifest_out)

print("\nSaved manifest:")
print(manifest_out)

print("\nStep 04D complete.")