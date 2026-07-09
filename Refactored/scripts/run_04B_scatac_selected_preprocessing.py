# -*- coding: utf-8 -*-
"""
Step 04B — selected scATAC preprocessing.

Purpose:
- Read selected ATAC preprocessing parameters from preprocessing_atac.yaml.
- Load raw ATAC object.
- Compute ATAC QC metrics.
- Apply selected filtering thresholds.
- Run TF-IDF normalization.
- Run LSI using TruncatedSVD.
- Build neighbors graph using selected LSI components.
- Run UMAP and Leiden clustering.
- Save filtered, LSI, clustered objects, audit tables, figures, and manifest.

Important:
- Step 04A is diagnostic only.
- Step 04B is the selected ATAC preprocessing run.
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
from sklearn.decomposition import TruncatedSVD


try:
    refactor_dir = Path(__file__).resolve().parents[1]  # Refactored root
except NameError:
    refactor_dir = Path.cwd().resolve()  # interactive fallback

src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))  # allow local pbmcgrn imports

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


def resolve_project_path(paths: ProjectPaths, relative_path: str, label: str) -> Path:
    """Resolve a YAML relative path inside the active project."""

    out_path = paths.prepare_output_path(relative_path, label=label)  # project-aware output path

    assert_path_inside_project(out_path, paths.project_root, label)  # safety check

    return out_path


def add_basic_atac_qc_metrics(atac):
    """Add basic ATAC QC metrics to obs and var."""

    X = atac.X  # cells x peaks matrix

    atac.obs["total_counts"] = np.asarray(X.sum(axis=1)).ravel()  # total fragments/counts per cell

    atac.obs["n_peaks_by_counts"] = np.asarray((X > 0).sum(axis=1)).ravel()  # nonzero peaks per cell

    atac.var["n_cells_by_counts"] = np.asarray((X > 0).sum(axis=0)).ravel()  # cells per peak

    atac.var["total_counts"] = np.asarray(X.sum(axis=0)).ravel()  # total counts per peak

    return atac


def validate_single_value(value, label: str) -> None:
    """Fail if a config value is a diagnostic list instead of a selected scalar."""

    if isinstance(value, list):
        raise ValueError(
            f"{label} must be a single selected value in Step 04B, not a list: {value}"
        )


def compute_tfidf_lsi(atac, n_components: int, scale_factor: float, binarize: bool, random_state: int):
    """Compute TF-IDF and LSI for scATAC peak matrix."""

    X = atac.X  # input count/accessibility matrix

    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)  # convert dense matrix to sparse CSR
    else:
        X = X.tocsr()  # ensure efficient row slicing

    if binarize:
        X_use = X.copy()  # copy count matrix
        X_use.data = np.ones_like(X_use.data)  # convert nonzero counts to 1
    else:
        X_use = X.copy()  # keep counts if requested

    cell_sums = np.asarray(X_use.sum(axis=1)).ravel()  # total accessible peaks/counts per cell
    cell_sums[cell_sums == 0] = 1  # avoid division by zero

    tf = X_use.multiply(1 / cell_sums[:, None])  # term frequency: normalize each cell

    n_cells = X_use.shape[0]  # number of cells

    peak_sums = np.asarray(X_use.sum(axis=0)).ravel()  # number/count of cells per peak

    idf = np.log(1 + n_cells / (1 + peak_sums))  # inverse document frequency

    tfidf = tf.multiply(idf)  # TF-IDF weighted matrix

    tfidf = tfidf * scale_factor  # scale for numerical stability

    svd = TruncatedSVD(
        n_components=n_components,
        random_state=random_state,
    )  # sparse-friendly SVD for LSI

    X_lsi = svd.fit_transform(tfidf)  # cells x LSI components

    atac.obsm["X_lsi"] = X_lsi  # store LSI embedding

    atac.uns["lsi"] = {
        "n_components": int(n_components),
        "scale_factor": float(scale_factor),
        "binarize": bool(binarize),
        "random_state": int(random_state),
        "explained_variance_ratio": svd.explained_variance_ratio_.tolist(),
        "singular_values": svd.singular_values_.tolist(),
    }  # store LSI metadata

    return atac, svd


paths = ProjectPaths(refactor_root=refactor_dir)  # project-aware paths
paths.ensure_output_dirs()  # create standard output dirs

print("\nActive project:")
print(paths.project_id)

print("\nProject root:")
print(paths.project_root)


# ============================================================
# 1. Load preprocessing_atac.yaml
# ============================================================

config_path = paths.configs / "preprocessing_atac.yaml"  # selected ATAC config

assert_path_inside_project(config_path, paths.project_root, "preprocessing_atac.yaml")

if not config_path.exists():
    raise FileNotFoundError(f"preprocessing_atac.yaml not found:\n{config_path}")

config = load_yaml(config_path)  # load selected ATAC config

input_cfg = config["input"]  # input paths
qc_cfg = config["qc"]  # selected QC thresholds
tfidf_cfg = config["tfidf"]  # TF-IDF settings
lsi_cfg = config["lsi"]  # LSI settings
neighbors_cfg = config["neighbors"]  # graph settings
umap_cfg = config["umap"]  # UMAP settings
clustering_cfg = config["clustering"]  # Leiden settings
outputs_cfg = config["outputs"]  # output paths

for label, value in {
    "qc.min_peaks_per_cell": qc_cfg["min_peaks_per_cell"],
    "qc.max_counts_per_cell": qc_cfg["max_counts_per_cell"],
    "qc.max_peaks_per_cell": qc_cfg.get("max_peaks_per_cell"),
    "qc.min_cells_per_peak": qc_cfg["min_cells_per_peak"],
    "lsi.n_components": lsi_cfg["n_components"],
    "neighbors.n_neighbors": neighbors_cfg["n_neighbors"],
    "clustering.leiden_resolution": clustering_cfg["leiden_resolution"],
}.items():
    validate_single_value(value, label)  # ensure Step 04B uses selected values, not diagnostic lists

print("\nLoaded ATAC preprocessing config:")
print(config_path)


# ============================================================
# 2. Resolve paths
# ============================================================

atac_in = paths.project_root / input_cfg["raw_atac_h5ad"]  # raw ATAC object from YAML

filtered_out = resolve_project_path(
    paths,
    outputs_cfg["atac_filtered_h5ad"],
    "atac_filtered_h5ad",
)

lsi_out = resolve_project_path(
    paths,
    outputs_cfg["atac_lsi_h5ad"],
    "atac_lsi_h5ad",
)

clustered_out = resolve_project_path(
    paths,
    outputs_cfg["atac_clustered_h5ad"],
    "atac_clustered_h5ad",
)

thresholds_used_out = resolve_project_path(
    paths,
    outputs_cfg["atac_qc_thresholds_used"],
    "atac_qc_thresholds_used",
)

filtering_summary_out = resolve_project_path(
    paths,
    outputs_cfg["atac_filtering_summary"],
    "atac_filtering_summary",
)

lsi_variance_table_out = resolve_project_path(
    paths,
    outputs_cfg["atac_lsi_variance_table"],
    "atac_lsi_variance_table",
)

cluster_summary_out = resolve_project_path(
    paths,
    outputs_cfg["atac_cluster_summary"],
    "atac_cluster_summary",
)

lsi_variance_plot_out = resolve_project_path(
    paths,
    outputs_cfg["atac_lsi_variance_plot"],
    "atac_lsi_variance_plot",
)

umap_clusters_out = resolve_project_path(
    paths,
    outputs_cfg["atac_umap_clusters"],
    "atac_umap_clusters",
)

manifest_out = resolve_project_path(
    paths,
    outputs_cfg["manifest"],
    "step04B_manifest",
)

for label, path in {
    "atac_in": atac_in,
    "filtered_out": filtered_out,
    "lsi_out": lsi_out,
    "clustered_out": clustered_out,
    "thresholds_used_out": thresholds_used_out,
    "filtering_summary_out": filtering_summary_out,
    "lsi_variance_table_out": lsi_variance_table_out,
    "cluster_summary_out": cluster_summary_out,
    "lsi_variance_plot_out": lsi_variance_plot_out,
    "umap_clusters_out": umap_clusters_out,
    "manifest_out": manifest_out,
}.items():
    assert_path_inside_project(path, paths.project_root, label)

filtered_out.parent.mkdir(parents=True, exist_ok=True)
thresholds_used_out.parent.mkdir(parents=True, exist_ok=True)
lsi_variance_plot_out.parent.mkdir(parents=True, exist_ok=True)
manifest_out.parent.mkdir(parents=True, exist_ok=True)


# ============================================================
# 3. Load raw ATAC object
# ============================================================

if not atac_in.exists():
    raise FileNotFoundError(f"Raw ATAC object not found:\n{atac_in}")

atac = sc.read_h5ad(atac_in)  # load raw ATAC object

print("\nRaw ATAC object:")
print(atac)

atac = add_basic_atac_qc_metrics(atac)  # compute QC metrics


# ============================================================
# 4. Apply selected ATAC filtering
# ============================================================

min_peaks_per_cell = int(qc_cfg["min_peaks_per_cell"])  # lower complexity threshold
max_counts_per_cell = float(qc_cfg["max_counts_per_cell"])  # high-count outlier threshold
max_peaks_per_cell = qc_cfg.get("max_peaks_per_cell")  # optional high-peak outlier threshold
min_cells_per_peak = int(qc_cfg["min_cells_per_peak"])  # rare peak threshold

n_cells_before = int(atac.n_obs)  # raw cells
n_peaks_before = int(atac.n_vars)  # raw peaks

cell_filter = (
    (atac.obs["n_peaks_by_counts"] >= min_peaks_per_cell)
    & (atac.obs["total_counts"] <= max_counts_per_cell)
)  # selected cell filtering rule

if max_peaks_per_cell is not None:
    cell_filter = cell_filter & (
        atac.obs["n_peaks_by_counts"] <= int(max_peaks_per_cell)
    )  # optional upper complexity threshold

peak_filter = atac.var["n_cells_by_counts"] >= min_cells_per_peak  # selected peak filtering rule

atac_filtered = atac[cell_filter, peak_filter].copy()  # filtered ATAC object

atac_filtered = add_basic_atac_qc_metrics(atac_filtered)  # recompute QC after filtering

n_cells_after = int(atac_filtered.n_obs)  # filtered cells
n_peaks_after = int(atac_filtered.n_vars)  # filtered peaks

print("\nFiltered ATAC object:")
print(atac_filtered)

thresholds_used = pd.DataFrame(
    {
        "parameter": [
            "min_peaks_per_cell",
            "max_counts_per_cell",
            "max_peaks_per_cell",
            "min_cells_per_peak",
        ],
        "value": [
            min_peaks_per_cell,
            max_counts_per_cell,
            max_peaks_per_cell,
            min_cells_per_peak,
        ],
    }
)  # selected thresholds audit table

write_csv(thresholds_used, thresholds_used_out, index=False)

filtering_summary = pd.DataFrame(
    {
        "metric": [
            "cells_before",
            "cells_after",
            "cells_removed",
            "cell_retention_fraction",
            "peaks_before",
            "peaks_after",
            "peaks_removed",
            "peak_retention_fraction",
            "min_peaks_per_cell",
            "max_counts_per_cell",
            "max_peaks_per_cell",
            "min_cells_per_peak",
            "mean_total_counts_after",
            "median_total_counts_after",
            "max_total_counts_after",
            "mean_n_peaks_by_counts_after",
            "median_n_peaks_by_counts_after",
            "max_n_peaks_by_counts_after",
        ],
        "value": [
            n_cells_before,
            n_cells_after,
            n_cells_before - n_cells_after,
            n_cells_after / n_cells_before,
            n_peaks_before,
            n_peaks_after,
            n_peaks_before - n_peaks_after,
            n_peaks_after / n_peaks_before,
            min_peaks_per_cell,
            max_counts_per_cell,
            max_peaks_per_cell,
            min_cells_per_peak,
            float(atac_filtered.obs["total_counts"].mean()),
            float(atac_filtered.obs["total_counts"].median()),
            float(atac_filtered.obs["total_counts"].max()),
            float(atac_filtered.obs["n_peaks_by_counts"].mean()),
            float(atac_filtered.obs["n_peaks_by_counts"].median()),
            float(atac_filtered.obs["n_peaks_by_counts"].max()),
        ],
    }
)  # filtering audit summary

write_csv(filtering_summary, filtering_summary_out, index=False)

print("\nFiltering summary:")
print(filtering_summary)

print("\nSaved filtering summary to:")
print(filtering_summary_out)

atac_filtered.layers["counts"] = atac_filtered.X.copy()  # preserve filtered raw counts

atac_filtered.write_h5ad(filtered_out)  # save filtered object

print("\nSaved filtered ATAC object to:")
print(filtered_out)


# ============================================================
# 5. TF-IDF normalization and LSI
# ============================================================

n_components = int(lsi_cfg["n_components"])  # number of LSI components
lsi_random_state = int(lsi_cfg["random_state"])  # reproducible SVD
tfidf_scale_factor = float(tfidf_cfg["scale_factor"])  # numerical scaling
tfidf_binarize = bool(tfidf_cfg["binarize"])  # use binary accessibility

atac_lsi, svd = compute_tfidf_lsi(
    atac_filtered,
    n_components=n_components,
    scale_factor=tfidf_scale_factor,
    binarize=tfidf_binarize,
    random_state=lsi_random_state,
)  # compute LSI representation

lsi_variance = pd.DataFrame(
    {
        "component": np.arange(1, n_components + 1),
        "explained_variance_ratio": svd.explained_variance_ratio_,
        "singular_value": svd.singular_values_,
    }
)  # LSI variance table

write_csv(lsi_variance, lsi_variance_table_out, index=False)

plt.figure(figsize=(7, 4))
plt.plot(
    lsi_variance["component"],
    lsi_variance["explained_variance_ratio"],
    marker="o",
    linewidth=1,
)  # LSI variance elbow-style plot
plt.xlabel("LSI component")
plt.ylabel("Explained variance ratio")
plt.title("ATAC LSI explained variance")
plt.tight_layout()
save_show_close(lsi_variance_plot_out, dpi=300, show=True)

atac_lsi.write_h5ad(lsi_out)  # save LSI object

print("\nSaved LSI ATAC object to:")
print(lsi_out)

print("\nSaved LSI variance table to:")
print(lsi_variance_table_out)


# ============================================================
# 6. Neighbors, UMAP, Leiden
# ============================================================

lsi_start_component = int(lsi_cfg["lsi_start_component"])  # human-readable start component, e.g. 2
lsi_end_component = int(lsi_cfg["lsi_end_component"])  # human-readable end component, e.g. 50

if lsi_start_component < 1:
    raise ValueError("lsi_start_component must be >= 1.")

if lsi_end_component > n_components:
    raise ValueError(
        f"lsi_end_component={lsi_end_component} cannot exceed n_components={n_components}."
    )

slice_start = lsi_start_component - 1  # convert to Python zero-based index
slice_end = lsi_end_component  # Python end is exclusive, so this keeps component 50

lsi_rep_key = f"X_lsi_{lsi_start_component}_{lsi_end_component}"  # e.g. X_lsi_2_50

atac_lsi.obsm[lsi_rep_key] = atac_lsi.obsm["X_lsi"][:, slice_start:slice_end]  # selected LSI components

use_rep = str(neighbors_cfg["use_rep"])  # representation from YAML

if use_rep != lsi_rep_key:
    print(
        f"\nWarning: YAML neighbors.use_rep={use_rep}, but selected LSI key is {lsi_rep_key}."
    )
    print("Using selected LSI key from lsi_start_component/lsi_end_component.")
    use_rep = lsi_rep_key

n_neighbors = int(neighbors_cfg["n_neighbors"])  # graph neighbors
metric = str(neighbors_cfg["metric"])  # cosine for LSI

sc.pp.neighbors(
    atac_lsi,
    n_neighbors=n_neighbors,
    use_rep=use_rep,
    metric=metric,
)  # build ATAC neighbor graph

sc.tl.umap(
    atac_lsi,
    random_state=int(umap_cfg["random_state"]),
)  # compute UMAP

leiden_resolution = float(clustering_cfg["leiden_resolution"])  # selected clustering resolution
leiden_key = str(clustering_cfg["leiden_key"])  # selected cluster key

sc.tl.leiden(
    atac_lsi,
    resolution=leiden_resolution,
    key_added=leiden_key,
)  # Leiden clustering

cluster_summary = (
    atac_lsi.obs[leiden_key]
    .value_counts()
    .rename_axis(leiden_key)
    .reset_index(name="n_cells")
    .sort_values(leiden_key)
)  # cluster size summary

write_csv(cluster_summary, cluster_summary_out, index=False)

print("\nATAC cluster summary:")
print(cluster_summary)

sc.pl.umap(
    atac_lsi,
    color=leiden_key,
    legend_loc="on data",
    legend_fontsize=6,
    legend_fontweight="normal",
    frameon=False,
    title=f"ATAC Leiden clusters ({leiden_key})",
    show=False,
)  # UMAP colored by selected ATAC clusters

save_show_close(umap_clusters_out, dpi=300, show=True)

atac_lsi.write_h5ad(clustered_out)  # save final clustered ATAC object

print("\nSaved clustered ATAC object to:")
print(clustered_out)


# ============================================================
# 7. Manifest
# ============================================================

manifest = {
    "step": "04B",
    "description": "Selected scATAC preprocessing: filtering, TF-IDF, LSI, neighbors, UMAP, Leiden",
    "timestamp": datetime.now().isoformat(timespec="seconds"),
    "project": paths.as_dict(),
    "inputs": {
        "raw_atac_h5ad": str(atac_in),
        "preprocessing_atac_config": str(config_path),
    },
    "selected_qc_thresholds": {
        "min_peaks_per_cell": int(min_peaks_per_cell),
        "max_counts_per_cell": float(max_counts_per_cell),
        "max_peaks_per_cell": None if max_peaks_per_cell is None else int(max_peaks_per_cell),
        "min_cells_per_peak": int(min_cells_per_peak),
    },
    "tfidf": {
        "binarize": bool(tfidf_binarize),
        "scale_factor": float(tfidf_scale_factor),
    },
    "lsi": {
        "n_components": int(n_components),
        "random_state": int(lsi_random_state),
        "selected_lsi_components": f"{lsi_start_component}-{lsi_end_component}",
        "selected_lsi_rep": use_rep,
    },
    "neighbors_umap_clustering": {
        "n_neighbors": int(n_neighbors),
        "metric": metric,
        "umap_random_state": int(umap_cfg["random_state"]),
        "leiden_resolution": float(leiden_resolution),
        "leiden_key": leiden_key,
    },
    "outputs": {
        "filtered_atac_h5ad": str(filtered_out),
        "lsi_atac_h5ad": str(lsi_out),
        "clustered_atac_h5ad": str(clustered_out),
        "qc_thresholds_used": str(thresholds_used_out),
        "filtering_summary": str(filtering_summary_out),
        "lsi_variance_table": str(lsi_variance_table_out),
        "cluster_summary": str(cluster_summary_out),
        "lsi_variance_plot": str(lsi_variance_plot_out),
        "umap_clusters": str(umap_clusters_out),
    },
    "summary": {
        "n_cells_before": int(n_cells_before),
        "n_cells_after": int(n_cells_after),
        "cell_retention_fraction": float(n_cells_after / n_cells_before),
        "n_peaks_before": int(n_peaks_before),
        "n_peaks_after": int(n_peaks_after),
        "peak_retention_fraction": float(n_peaks_after / n_peaks_before),
        "n_clusters": int(atac_lsi.obs[leiden_key].nunique()),
    },
    "rules": {
        "uses_selected_yaml_values_only": True,
        "does_not_use_diagnostic_lists": True,
        "outputs_project_local_only": True,
        "raw_atac_object_not_overwritten": True,
    },
}

write_json(manifest, manifest_out)

print("\nSaved Step 04B manifest to:")
print(manifest_out)

print("\nDone.")
print("Step 04B selected scATAC preprocessing completed.")


