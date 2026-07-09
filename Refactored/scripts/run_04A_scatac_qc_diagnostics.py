# -*- coding: utf-8 -*-
"""
Step 04A — scATAC QC diagnostics.

Purpose:
- Load raw ATAC AnnData object.
- Compute basic ATAC QC metrics.
- Summarize cell-level and peak-level QC.
- Save diagnostic tables, plots, and manifest.
- Do NOT perform final filtering in this step.

This is the ATAC equivalent of RNA Step 02A diagnostics.
"""

from pathlib import Path
import sys
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


try:
    refactor_dir = Path(__file__).resolve().parents[1]  # Refactored root
except NameError:
    refactor_dir = Path.cwd().resolve()  # interactive fallback

src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))  # allow local pbmcgrn imports

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

    plt.close()  # close figure to avoid memory buildup


def get_processed_atac_dir(paths: ProjectPaths) -> Path:
    """Return project-local processed ATAC directory."""

    if hasattr(paths, "processed_atac"):
        return paths.processed_atac  # use ProjectPaths attribute if it exists

    return paths.project_root / "data" / "processed" / "atac"  # safe fallback inside project


def add_basic_atac_qc_metrics(atac):
    """Add basic ATAC QC metrics to obs and var."""

    X = atac.X  # cell x peak matrix

    atac.obs["total_counts"] = np.asarray(X.sum(axis=1)).ravel()  # total ATAC counts/fragments per cell

    atac.obs["n_peaks_by_counts"] = np.asarray((X > 0).sum(axis=1)).ravel()  # accessible peaks per cell

    atac.var["n_cells_by_counts"] = np.asarray((X > 0).sum(axis=0)).ravel()  # cells with accessibility per peak

    atac.var["total_counts"] = np.asarray(X.sum(axis=0)).ravel()  # total counts per peak

    return atac


paths = ProjectPaths(refactor_root=refactor_dir)  # project-aware paths
paths.ensure_output_dirs()  # create standard output dirs

print("\nActive project:")
print(paths.project_id)

print("\nProject root:")
print(paths.project_root)


# ============================================================
# 1. Resolve project-aware paths
# ============================================================

processed_atac = get_processed_atac_dir(paths)  # project-local ATAC processed folder

atac_in = processed_atac / "pbmc_atac_raw.h5ad"  # raw ATAC object

qc_summary_out = paths.prepare_output_path(
    "reports/tables/atac/atac_qc_summary_before_filtering.csv",
    label="atac_qc_summary_before_filtering",
)  # summary of basic ATAC QC metrics

qc_quantiles_out = paths.prepare_output_path(
    "reports/tables/atac/atac_qc_quantiles_before_filtering.csv",
    label="atac_qc_quantiles_before_filtering",
)  # quantiles for threshold selection

figures_dir = paths.prepare_output_path(
    "reports/figures/atac/qc",
    label="atac_qc_figures_dir",
)  # ATAC QC figure folder

total_counts_fig = figures_dir / "atac_total_counts_per_cell_before_filtering.png"

n_peaks_fig = figures_dir / "atac_n_peaks_by_counts_per_cell_before_filtering.png"

cells_per_peak_fig = figures_dir / "atac_cells_per_peak_before_filtering.png"

manifest_out = paths.prepare_output_path(
    "reports/manifests/step04A_scatac_qc_diagnostics_manifest.json",
    label="step04A_scatac_qc_diagnostics_manifest",
)

for label, path in {
    "atac_in": atac_in,
    "qc_summary_out": qc_summary_out,
    "qc_quantiles_out": qc_quantiles_out,
    "figures_dir": figures_dir,
    "total_counts_fig": total_counts_fig,
    "n_peaks_fig": n_peaks_fig,
    "cells_per_peak_fig": cells_per_peak_fig,
    "manifest_out": manifest_out,
}.items():
    assert_path_inside_project(path, paths.project_root, label)

figures_dir.mkdir(parents=True, exist_ok=True)
qc_summary_out.parent.mkdir(parents=True, exist_ok=True)
manifest_out.parent.mkdir(parents=True, exist_ok=True)


# ============================================================
# 2. Load raw ATAC object
# ============================================================

if not atac_in.exists():
    raise FileNotFoundError(
        f"Raw ATAC object not found:\n{atac_in}\n\n"
        "Expected project-local path:\n"
        "data/processed/atac/pbmc_atac_raw.h5ad"
    )

atac = sc.read_h5ad(atac_in)  # load raw ATAC AnnData object

print("\nLoaded raw ATAC object:")
print(atac)

print("\nATAC obs columns:")
print(atac.obs.columns.tolist())

print("\nATAC var columns:")
print(atac.var.columns.tolist())


# ============================================================
# 3. Compute basic ATAC QC metrics
# ============================================================

atac = add_basic_atac_qc_metrics(atac)  # add total_counts, n_peaks_by_counts, n_cells_by_counts

print("\nAdded basic ATAC QC metrics:")
print(["total_counts", "n_peaks_by_counts", "n_cells_by_counts"])


# ============================================================
# 4. Save QC summary table
# ============================================================

qc_summary = pd.DataFrame(
    {
        "metric": [
            "n_cells",
            "n_peaks",
            "mean_total_counts_per_cell",
            "median_total_counts_per_cell",
            "min_total_counts_per_cell",
            "max_total_counts_per_cell",
            "mean_n_peaks_by_counts_per_cell",
            "median_n_peaks_by_counts_per_cell",
            "min_n_peaks_by_counts_per_cell",
            "max_n_peaks_by_counts_per_cell",
            "mean_cells_per_peak",
            "median_cells_per_peak",
            "min_cells_per_peak",
            "max_cells_per_peak",
        ],
        "value": [
            int(atac.n_obs),
            int(atac.n_vars),
            float(atac.obs["total_counts"].mean()),
            float(atac.obs["total_counts"].median()),
            float(atac.obs["total_counts"].min()),
            float(atac.obs["total_counts"].max()),
            float(atac.obs["n_peaks_by_counts"].mean()),
            float(atac.obs["n_peaks_by_counts"].median()),
            float(atac.obs["n_peaks_by_counts"].min()),
            float(atac.obs["n_peaks_by_counts"].max()),
            float(atac.var["n_cells_by_counts"].mean()),
            float(atac.var["n_cells_by_counts"].median()),
            float(atac.var["n_cells_by_counts"].min()),
            float(atac.var["n_cells_by_counts"].max()),
        ],
    }
)  # compact QC summary

write_csv(qc_summary, qc_summary_out, index=False)  # save summary table

print("\nATAC QC summary before filtering:")
print(qc_summary)

print("\nSaved ATAC QC summary to:")
print(qc_summary_out)


# ============================================================
# 5. Save QC quantiles table
# ============================================================

quantile_levels = [
    0.00,
    0.01,
    0.02,
    0.05,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.98,
    0.99,
    1.00,
]  # useful thresholds for manual filtering decisions

cell_quantiles = atac.obs[
    ["total_counts", "n_peaks_by_counts"]
].quantile(quantile_levels)  # cell-level QC quantiles

peak_quantiles = atac.var[
    ["n_cells_by_counts", "total_counts"]
].quantile(quantile_levels)  # peak-level QC quantiles

cell_quantiles = cell_quantiles.reset_index().rename(columns={"index": "quantile"})  # table format

peak_quantiles = peak_quantiles.reset_index().rename(columns={"index": "quantile"})  # table format

cell_quantiles["level"] = "cell"  # mark cell-level metrics

peak_quantiles["level"] = "peak"  # mark peak-level metrics

qc_quantiles = pd.concat(
    [cell_quantiles, peak_quantiles],
    axis=0,
    ignore_index=True,
    sort=False,
)  # combined quantile table

write_csv(qc_quantiles, qc_quantiles_out, index=False)  # save quantiles

print("\nATAC QC quantiles before filtering:")
print(qc_quantiles)

print("\nSaved ATAC QC quantiles to:")
print(qc_quantiles_out)


# ============================================================
# 6. Save QC diagnostic plots
# ============================================================

plt.figure(figsize=(6, 4))
plt.hist(atac.obs["total_counts"], bins=100)  # total counts distribution
plt.xlabel("Total ATAC counts per cell")
plt.ylabel("Number of cells")
plt.title("ATAC total counts per cell before filtering")
plt.tight_layout()
save_show_close(total_counts_fig, dpi=300, show=True)

plt.figure(figsize=(6, 4))
plt.hist(atac.obs["n_peaks_by_counts"], bins=100)  # accessible peaks distribution
plt.xlabel("Accessible peaks per cell")
plt.ylabel("Number of cells")
plt.title("ATAC accessible peaks per cell before filtering")
plt.tight_layout()
save_show_close(n_peaks_fig, dpi=300, show=True)

plt.figure(figsize=(6, 4))
plt.hist(atac.var["n_cells_by_counts"], bins=100)  # cells per peak distribution
plt.xlabel("Cells per peak")
plt.ylabel("Number of peaks")
plt.title("ATAC cells per peak before filtering")
plt.tight_layout()
save_show_close(cells_per_peak_fig, dpi=300, show=True)

print("\nSaved ATAC QC figures to:")
print(figures_dir)


# ============================================================
# 7. Manifest
# ============================================================

manifest = {
    "step": "04A",
    "description": "scATAC QC diagnostics before selected filtering",
    "timestamp": datetime.now().isoformat(timespec="seconds"),
    "project": paths.as_dict(),
    "inputs": {
        "raw_atac_h5ad": str(atac_in),
    },
    "outputs": {
        "qc_summary": str(qc_summary_out),
        "qc_quantiles": str(qc_quantiles_out),
        "total_counts_histogram": str(total_counts_fig),
        "n_peaks_by_counts_histogram": str(n_peaks_fig),
        "cells_per_peak_histogram": str(cells_per_peak_fig),
    },
    "summary": {
        "n_cells": int(atac.n_obs),
        "n_peaks": int(atac.n_vars),
        "median_total_counts_per_cell": float(atac.obs["total_counts"].median()),
        "median_n_peaks_by_counts_per_cell": float(atac.obs["n_peaks_by_counts"].median()),
        "median_cells_per_peak": float(atac.var["n_cells_by_counts"].median()),
    },
    "rules": {
        "diagnostic_only": True,
        "does_not_filter_cells_or_peaks": True,
        "does_not_overwrite_raw_atac_object": True,
        "outputs_project_local_only": True,
    },
}

write_json(manifest, manifest_out)  # save reproducibility manifest

print("\nSaved Step 04A manifest to:")
print(manifest_out)

print("\nDone.")
print("Step 04A scATAC QC diagnostics completed.")