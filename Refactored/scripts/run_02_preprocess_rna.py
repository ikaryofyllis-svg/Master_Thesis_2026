"""Βήμα 02B — Selected RNA preprocessing rerun / validation.

This script performs the selected RNA preprocessing run after Step 02A.

Main purpose:
- validate that the final RNA QC thresholds selected in Step 02A
  match the thresholds stored in preprocessing_rna.yaml
- rerun the RNA preprocessing pipeline using those selected thresholds
- save all outputs inside the active project workspace only
- write a manifest with the selected-threshold handoff metadata

Pipeline:
raw RNA
→ QC metrics
→ selected QC filtering
→ normalization/log1p
→ HVG selection
→ PCA/neighbors/UMAP/Leiden
→ tables, plots, h5ad objects, manifest
"""

from __future__ import annotations  # cleaner type hints

from pathlib import Path  # safe path handling
import sys  # lets us add Refactored/src to Python path
from typing import Any  # general type hints

import pandas as pd  # tables and summaries
import scanpy as sc  # single-cell RNA-seq analysis
import matplotlib.pyplot as plt  # save UMAP plots

from itertools import product  # create all parameter combinations from lists
from datetime import datetime  # timestamp for YAML backup
import shutil  # copy YAML before overwrite
import yaml  # write updated YAML safely

# ============================================================
# 0. Define Refactored project root and Python import path
# ============================================================


refactor_dir = Path(__file__).resolve().parents[1]  # Refactored root folder

src_dir = refactor_dir / "src"  # shared package source folder

if str(src_dir) not in sys.path:  # avoid duplicate path insertion
    sys.path.insert(0, str(src_dir))  # allows imports from pbmcgr


# ============================================================
# 1. Import reusable project helpers
# ============================================================

from pbmcgrn.config import load_yaml  # YAML config loader
from pbmcgrn.io import write_csv, write_json  # reusable output writers
from pbmcgrn.paths import ProjectPaths  # project-aware path manager

from pbmcgrn.preprocessing.rna import (
    calculate_rna_qc_metrics,  # calculates RNA QC metrics
    summarize_rna_qc,  # summarizes QC before/after filtering
    filter_rna_cells_and_genes,  # filters cells and genes
    normalize_log_hvg,  # normalization, log1p, HVG selection
    run_rna_dimensionality_reduction_and_clustering,  # old all-in-one PCA/UMAP/Leiden
    summarize_clusters,  # cluster size summary

    run_rna_pca_only,  # new: scale + PCA only for elbow diagnostics
    summarize_pca_variance,  # new: PCA variance/cumulative variance table
    run_rna_neighbors_umap_leiden,  # new: neighbors + UMAP + Leiden after PCA
)


# ============================================================
# 2. Step 02B helper functions: selected threshold validation
# ============================================================

def require_single_value(value: Any, parameter_name: str) -> Any:
    """Ensure a YAML parameter is a single final value, not a diagnostic list."""

    if isinstance(value, list):  # Step 02A allowed candidate lists, but Step 02B must not
        if len(value) != 1:  # multiple values means the selected run is ambiguous
            raise ValueError(
                f"{parameter_name} must be a single final value for Step 02B, "
                f"but got list: {value}. "
                "Run Step 02A final selection or edit preprocessing_rna.yaml."
            )  # fail early to protect reproducibility

        return value[0]  # allow a single-item list, but convert it to a scalar

    return value  # already scalar


def coerce_rna_qc_thresholds(rna_cfg: dict[str, Any]) -> dict[str, Any]:
    """Build a clean, typed threshold dictionary from preprocessing_rna.yaml."""

    selected_thresholds = {
        "mito_prefix": str(
            require_single_value(rna_cfg["mito_prefix"], "mito_prefix")
        ),  # mitochondrial gene prefix, usually MT- for human data

        "min_genes_per_cell": int(
            require_single_value(rna_cfg["min_genes_per_cell"], "min_genes_per_cell")
        ),  # minimum detected genes per cell

        "min_cells_per_gene": int(
            require_single_value(rna_cfg["min_cells_per_gene"], "min_cells_per_gene")
        ),  # minimum cells where a gene must be detected

        "max_pct_mito": float(
            require_single_value(rna_cfg["max_pct_mito"], "max_pct_mito")
        ),  # maximum mitochondrial percentage

        "max_total_counts": float(
            require_single_value(rna_cfg["max_total_counts"], "max_total_counts")
        ),  # maximum total counts per cell
    }

    return selected_thresholds  # final typed thresholds


def load_step02a_selected_thresholds(selected_thresholds_path: Path) -> dict[str, Any]:
    """Load the final selected threshold row written by Step 02A."""

    if not selected_thresholds_path.exists():  # Step 02B requires the 02A handoff file
        raise FileNotFoundError(
            "Step 02B requires the selected threshold file from Step 02A, "
            f"but it was not found:\n{selected_thresholds_path}\n\n"
            "Run Step 02A and update preprocessing_rna.yaml with the final selected candidate."
        )  # avoids silent use of unverified YAML values

    selected_df = pd.read_csv(selected_thresholds_path)  # read selected threshold table

    if selected_df.shape[0] != 1:  # selected run must have exactly one final threshold row
        raise ValueError(
            f"Expected exactly one selected threshold row in:\n{selected_thresholds_path}\n"
            f"Found {selected_df.shape[0]} rows instead."
        )  # prevents ambiguous selected preprocessing

    required_columns = [
        "mito_prefix",
        "min_genes_per_cell",
        "min_cells_per_gene",
        "max_pct_mito",
        "max_total_counts",
    ]  # columns expected from Step 02A selected-threshold table

    missing_columns = [
        col for col in required_columns if col not in selected_df.columns
    ]  # detect corrupted/incomplete 02A output

    if missing_columns:  # fail clearly if the table is not the expected file
        raise ValueError(
            f"Missing required columns in Step 02A selected threshold file: {missing_columns}\n"
            f"File: {selected_thresholds_path}"
        )  # protects the selected-run handoff

    selected_row = selected_df.iloc[0]  # take the one selected row

    selected_thresholds = {
        "mito_prefix": str(selected_row["mito_prefix"]),  # selected mitochondrial prefix
        "min_genes_per_cell": int(selected_row["min_genes_per_cell"]),  # selected min genes
        "min_cells_per_gene": int(selected_row["min_cells_per_gene"]),  # selected min cells per gene
        "max_pct_mito": float(selected_row["max_pct_mito"]),  # selected mito threshold
        "max_total_counts": float(selected_row["max_total_counts"]),  # selected total-count threshold
    }

    return selected_thresholds  # selected thresholds from Step 02A


def values_match(value_a: Any, value_b: Any, tolerance: float = 1e-9) -> bool:
    """Compare threshold values safely, including floats."""

    if isinstance(value_a, float) or isinstance(value_b, float):  # float comparison needs tolerance
        return abs(float(value_a) - float(value_b)) <= tolerance  # avoid false mismatch from float formatting

    return value_a == value_b  # exact comparison for strings and integers


def compare_threshold_dicts(
    yaml_thresholds: dict[str, Any],
    step02a_thresholds: dict[str, Any],
) -> pd.DataFrame:
    """Compare YAML thresholds with Step 02A selected thresholds."""

    rows = []  # mismatch/detail rows

    for parameter in yaml_thresholds:  # check all selected threshold parameters
        yaml_value = yaml_thresholds[parameter]  # value stored in preprocessing_rna.yaml
        step02a_value = step02a_thresholds[parameter]  # value selected by Step 02A

        rows.append(
            {
                "parameter": parameter,  # threshold name
                "yaml_value": yaml_value,  # value that Step 02B would use
                "step02a_selected_value": step02a_value,  # selected value from Step 02A
                "matches": values_match(yaml_value, step02a_value),  # pass/fail per parameter
            }
        )  # one row per threshold

    comparison_df = pd.DataFrame(rows)  # tidy validation table

    return comparison_df  # comparison report


def assert_path_inside_project(path: Path, project_root: Path, label: str) -> None:
    """Fail if an input/output path is outside the active project root."""

    resolved_path = Path(path).resolve()  # absolute normalized path
    resolved_project_root = Path(project_root).resolve()  # absolute project root

    try:
        resolved_path.relative_to(resolved_project_root)  # succeeds only if path is inside project
    except ValueError as exc:
        raise ValueError(
            f"{label} is outside the active project workspace.\n"
            f"Path: {resolved_path}\n"
            f"Active project root: {resolved_project_root}"
        ) from exc  # project-local output guard


def resolve_selected_rna_qc_thresholds(
    config: dict[str, Any],
    paths: ProjectPaths,
) -> tuple[dict[str, Any], str, Path, pd.DataFrame]:
    """Validate Step 02A → Step 02B selected RNA QC threshold handoff."""

    rna_cfg = config["rna_preprocessing"]  # YAML RNA preprocessing section

    yaml_thresholds = coerce_rna_qc_thresholds(rna_cfg)  # typed single-value thresholds from YAML

    selected_thresholds_path = (
        paths.project_root
        / "reports"
        / "tables"
        / "rna"
        / "rna_selected_qc_thresholds.csv"
    )  # expected Step 02A selected-threshold output

    assert_path_inside_project(
        selected_thresholds_path,
        paths.project_root,
        label="Step 02A selected threshold file",
    )  # ensure handoff file is project-local

    step02a_thresholds = load_step02a_selected_thresholds(
        selected_thresholds_path
    )  # load selected thresholds from Step 02A

    threshold_comparison = compare_threshold_dicts(
        yaml_thresholds=yaml_thresholds,
        step02a_thresholds=step02a_thresholds,
    )  # compare YAML with Step 02A output

    if not bool(threshold_comparison["matches"].all()):  # any mismatch means ambiguous run
        raise ValueError(
            "Step 02B cannot continue because preprocessing_rna.yaml does not match "
            "the selected thresholds from Step 02A.\n\n"
            f"{threshold_comparison}\n\n"
            "Fix preprocessing_rna.yaml or rerun Step 02A final selection."
        )  # protects thesis reproducibility

    threshold_source = (
        "preprocessing_rna.yaml validated against Step 02A selected thresholds"
    )  # source description for logs/manifest

    print("\nSelected RNA QC thresholds for Step 02B:")
    print(pd.DataFrame([yaml_thresholds]))  # log exactly what will be used

    print("\nStep 02A selected threshold validation:")
    print(threshold_comparison)  # log validation table

    print("\nValidated against Step 02A selected thresholds file:")
    print(selected_thresholds_path)  # show handoff file path

    return yaml_thresholds, threshold_source, selected_thresholds_path, threshold_comparison


def validate_log1p_setting(norm_cfg: dict[str, Any]) -> None:
    """Protect against mismatch between YAML and current reusable normalize_log_hvg implementation."""

    if bool(norm_cfg["log1p"]) is not True:  # current reusable function always applies sc.pp.log1p
        raise ValueError(
            "normalization.log1p is set to false in preprocessing_rna.yaml, "
            "but the current normalize_log_hvg() function applies log1p. "
            "Set normalization.log1p to true or refactor normalize_log_hvg()."
        )  # avoids silent config/function mismatch

def parse_int_or_int_list(user_text: str) -> int | list[int]:
    """Parse user input as one integer or comma-separated integer list."""

    values = [x.strip() for x in user_text.split(",")]  # split by comma

    parsed_values = [int(x) for x in values if x != ""]  # convert each value to int

    if len(parsed_values) == 1:
        return parsed_values[0]  # single final value

    return parsed_values  # diagnostic list


def parse_float_or_float_list(user_text: str) -> float | list[float]:
    """Parse user input as one float or comma-separated float list."""

    values = [x.strip() for x in user_text.split(",")]  # split by comma

    parsed_values = [float(x) for x in values if x != ""]  # convert each value to float

    if len(parsed_values) == 1:
        return parsed_values[0]  # single final value

    return parsed_values  # diagnostic list


def as_list(value: Any) -> list[Any]:
    """Convert a scalar or list into a list."""

    if isinstance(value, list):
        return value  # already list

    return [value]  # scalar becomes one-item list


def safe_resolution_label(value: float) -> str:
    """Make Leiden resolution safe for filenames and obs keys."""

    return str(value).replace(".", "_")  # 0.6 -> 0_6


def summarize_grid_result(
    adata: sc.AnnData,
    cluster_key: str,
    n_pcs: int,
    n_neighbors: int,
    leiden_resolution: float,
) -> dict[str, Any]:
    """Summarize one clustering parameter combination."""

    cluster_counts = adata.obs[cluster_key].value_counts()  # cells per cluster

    return {
        "n_pcs": n_pcs,  # PCA dimensions used
        "n_neighbors": n_neighbors,  # neighbor graph size
        "leiden_resolution": leiden_resolution,  # Leiden resolution used
        "cluster_key": cluster_key,  # cluster column in obs
        "n_clusters": int(cluster_counts.shape[0]),  # number of clusters
        "min_cluster_size": int(cluster_counts.min()),  # smallest cluster
        "max_cluster_size": int(cluster_counts.max()),  # largest cluster
        "median_cluster_size": float(cluster_counts.median()),  # median cluster size
        "small_clusters_lt_20_cells": int((cluster_counts < 20).sum()),  # tiny clusters
        "small_clusters_lt_50_cells": int((cluster_counts < 50).sum()),  # small clusters
        "small_clusters_lt_100_cells": int((cluster_counts < 100).sum()),  # small-ish clusters
    }


def summarize_leiden_clusters_for_grid(
    adata: sc.AnnData,
    cluster_key: str,
    n_pcs: int,
    n_neighbors: int,
    leiden_resolution: float,
) -> pd.DataFrame:
    """Create cluster size table for one grid combination."""

    counts = (
        adata.obs[cluster_key]
        .value_counts()
        .sort_index()
        .rename_axis("cluster")
        .reset_index(name="n_cells")
    )  # one row per cluster

    counts["fraction_cells"] = counts["n_cells"] / adata.n_obs  # cluster fraction
    counts["n_pcs"] = n_pcs  # PCA dimensions
    counts["n_neighbors"] = n_neighbors  # neighbors parameter
    counts["leiden_resolution"] = leiden_resolution  # Leiden resolution
    counts["cluster_key"] = cluster_key  # obs key

    return counts


def backup_and_update_preprocessing_yaml(
    config: dict[str, Any],
    config_path: Path,
    selected_n_pcs: int,
    selected_n_neighbors: int,
    selected_leiden_resolution: float,
) -> Path:
    """Backup preprocessing_rna.yaml and update final PCA/neighbors/Leiden values."""

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")  # backup timestamp

    backup_path = config_path.with_suffix(
        f".before_step02B_clustering_update_{timestamp}.yaml"
    )  # backup YAML path

    shutil.copy2(config_path, backup_path)  # save backup before overwrite

    leiden_key = f"leiden_{safe_resolution_label(selected_leiden_resolution)}"  # final Leiden key

    config["dimensionality_reduction"]["n_pcs"] = int(selected_n_pcs)  # update final n_pcs
    config["dimensionality_reduction"]["n_neighbors"] = int(selected_n_neighbors)  # update final neighbors
    config["dimensionality_reduction"]["leiden_resolution"] = float(selected_leiden_resolution)  # update final resolution
    config["dimensionality_reduction"]["leiden_key"] = leiden_key  # update final Leiden key

    with open(config_path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)  # write updated YAML

    return backup_path  # path to backup YAML
# ============================================================
# 3. Initialize project-aware paths
# ============================================================

paths = ProjectPaths(refactor_root=refactor_dir)  # reads active_project_id from project registry
paths.ensure_output_dirs()  # creates project-local output folders

print("\nRefactored project root:")
print(paths.refactor_root)  # should be .../Refactored

print("\nActive project id:")
print(paths.project_id)  # should be pbmc_10k_v1

print("\nActive project root:")
print(paths.project_root)  # should be .../Refactored/projects/pbmc_10k_v1


# ============================================================
# 4. Load RNA preprocessing config
# ============================================================

config_path = paths.configs / "preprocessing_rna.yaml"  # project-local config file

print("\nLoading RNA preprocessing config:")
print(config_path)  # exact YAML file used by Step 02B

if not config_path.exists():  # fail clearly if config is missing
    raise FileNotFoundError(
        f"RNA preprocessing config not found:\n{config_path}"
    )  # no config means no reproducible selected run

assert_path_inside_project(
    config_path,
    paths.project_root,
    label="RNA preprocessing config",
)  # config must belong to active project

config = load_yaml(config_path)  # load preprocessing_rna.yaml

rna_cfg = config["rna_preprocessing"]  # raw YAML RNA QC/filtering parameters
norm_cfg = config["normalization"]  # normalization parameters
hvg_cfg = config["highly_variable_genes"]  # HVG parameters
dim_cfg = config["dimensionality_reduction"]  # PCA/UMAP/Leiden parameters
outputs = config["outputs"]  # output paths from YAML

validate_log1p_setting(norm_cfg)  # prevents YAML/function mismatch

print("\nLoaded RNA preprocessing config successfully.")


# ============================================================
# 5. Validate Step 02A selected thresholds against YAML
# ============================================================

selected_qc_thresholds, threshold_source, selected_thresholds_file, threshold_comparison = (
    resolve_selected_rna_qc_thresholds(
        config=config,
        paths=paths,
    )
)  # validates Step 02A → Step 02B threshold handoff

selected_thresholds_used_out = paths.prepare_output_path(
    "reports/tables/rna/rna_step02B_selected_qc_thresholds_used.csv",
    label="rna_step02B_selected_qc_thresholds_used",
)  # audit table for thresholds actually used in Step 02B

threshold_validation_out = paths.prepare_output_path(
    "reports/tables/rna/rna_step02B_threshold_validation.csv",
    label="rna_step02B_threshold_validation",
)  # validation report comparing YAML vs Step 02A


# ============================================================
# 6. Resolve and validate project-aware output paths
# ============================================================

rna_filtered_out = paths.prepare_output_path(
    outputs["rna_filtered_h5ad"],
    label="rna_filtered_h5ad",
)  # project-local filtered RNA h5ad

rna_qc_table = paths.prepare_output_path(
    outputs["rna_qc_table"],
    label="rna_qc_table",
)  # project-local QC summary table

rna_hvg_out = paths.prepare_output_path(
    outputs["rna_hvg_h5ad"],
    label="rna_hvg_h5ad",
)  # project-local HVG RNA h5ad

rna_clustered_out = paths.prepare_output_path(
    outputs["rna_clustered_h5ad"],
    label="rna_clustered_h5ad",
)  # project-local clustered RNA h5ad

rna_umap_out = paths.prepare_output_path(
    outputs["rna_umap_clusters"],
    label="rna_umap_clusters",
)  # project-local UMAP figure

filtering_summary_out = paths.prepare_output_path(
    outputs["rna_filtering_summary"],
    label="rna_filtering_summary",
)  # project-local filtering summary

hvg_summary_out = paths.prepare_output_path(
    outputs["rna_hvg_summary"],
    label="rna_hvg_summary",
)  # project-local HVG summary

cluster_summary_out = paths.prepare_output_path(
    outputs["rna_cluster_summary"],
    label="rna_cluster_summary",
)  # project-local cluster summary

manifest_out = paths.prepare_output_path(
    outputs["manifest"],
    label="step02_manifest",
)  # project-local manifest

rna_pca_elbow_out = paths.prepare_output_path(
    outputs["rna_pca_elbow"],
    label="rna_pca_elbow",
)  # project-local PCA elbow plot

rna_pca_variance_out = paths.prepare_output_path(
    outputs["rna_pca_variance_table"],
    label="rna_pca_variance_table",
)  # project-local PCA variance table

rna_pca_qc_out = paths.prepare_output_path(
    outputs["rna_pca_qc_plot"],
    label="rna_pca_qc_plot",
)  # project-local PCA QC plot

rna_umap_qc_out = paths.prepare_output_path(
    outputs["rna_umap_qc_plot"],
    label="rna_umap_qc_plot",
)  # project-local UMAP QC plot

rna_leiden_diagnostics_out = paths.prepare_output_path(
    outputs["rna_leiden_diagnostics"],
    label="rna_leiden_diagnostics",
)  # project-local Leiden diagnostics table

rna_umap_diagnostics_dir = paths.prepare_output_path(
    outputs["rna_umap_diagnostics_dir"],
    label="rna_umap_diagnostics_dir",
)  # project-local UMAP diagnostics directory

clustering_grid_summary_out = paths.prepare_output_path(
    "reports/tables/rna/rna_step02B_clustering_parameter_grid_summary.csv",
    label="rna_step02B_clustering_parameter_grid_summary",
)  # diagnostic grid summary

clustering_grid_cluster_sizes_out = paths.prepare_output_path(
    "reports/tables/rna/rna_step02B_clustering_parameter_grid_cluster_sizes.csv",
    label="rna_step02B_clustering_parameter_grid_cluster_sizes",
)  # cluster sizes per grid combination

selected_clustering_params_out = paths.prepare_output_path(
    "reports/tables/rna/rna_step02B_selected_clustering_parameters.csv",
    label="rna_step02B_selected_clustering_parameters",
)  # final selected PCA/neighbors/Leiden values

clustering_grid_figures_dir = paths.prepare_output_path(
    "reports/figures/rna/clustering_grid",
    label="rna_step02B_clustering_grid_figures_dir",
)  # diagnostic UMAP folder

clustering_grid_figures_dir.mkdir(parents=True, exist_ok=True)  # create folder

project_local_outputs = {
    "rna_filtered_h5ad": rna_filtered_out,
    "rna_qc_table": rna_qc_table,
    "rna_hvg_h5ad": rna_hvg_out,
    "rna_clustered_h5ad": rna_clustered_out,
    "rna_umap_clusters": rna_umap_out,
    "rna_filtering_summary": filtering_summary_out,
    "rna_hvg_summary": hvg_summary_out,
    "rna_cluster_summary": cluster_summary_out,
    "manifest": manifest_out,
    "rna_step02B_selected_qc_thresholds_used": selected_thresholds_used_out,
    "rna_step02B_threshold_validation": threshold_validation_out,
    "rna_step02B_clustering_parameter_grid_summary": clustering_grid_summary_out,
    "rna_step02B_clustering_parameter_grid_cluster_sizes": clustering_grid_cluster_sizes_out,
    "rna_step02B_selected_clustering_parameters": selected_clustering_params_out,
    "rna_step02B_clustering_grid_figures_dir": clustering_grid_figures_dir,
}  # all Step 02B outputs

for label, output_path in project_local_outputs.items():  # audit all outputs
    assert_path_inside_project(output_path, paths.project_root, label=label)  # enforce project-local writes

print("\nResolved project-aware outputs:")
for label, output_path in project_local_outputs.items():  # print all resolved outputs
    print(f"{label}: {output_path}")  # clean path log



# ============================================================
# 7. Save selected-threshold validation audit tables
# ============================================================

selected_thresholds_used = pd.DataFrame(
    [
        {
            **selected_qc_thresholds,  # thresholds actually used by Step 02B
            "threshold_source": threshold_source,  # provenance description
            "step02a_selected_thresholds_file": str(selected_thresholds_file),  # handoff file
            "config_file": str(config_path),  # YAML used
            "script": "run_02_preprocess_rna.py",  # current script
            "step": "02B",  # selected preprocessing rerun
        }
    ]
)  # one-row audit table

write_csv(
    selected_thresholds_used,
    selected_thresholds_used_out,
    index=False,
)  # save thresholds used by Step 02B

write_csv(
    threshold_comparison,
    threshold_validation_out,
    index=False,
)  # save YAML vs Step 02A validation table

print("\nSaved Step 02B selected thresholds used table to:")
print(selected_thresholds_used_out)  # audit output path

print("\nSaved Step 02B threshold validation table to:")
print(threshold_validation_out)  # validation output path


# ============================================================
# 8. Load raw RNA h5ad
# ============================================================

rna_raw_path = paths.processed_rna / "pbmc_rna_raw.h5ad"  # project-local raw RNA object from Step 01

print("\nLoading raw RNA object:")
print(rna_raw_path)  # should be inside active project

if not rna_raw_path.exists():  # fail if Step 01 output is missing
    raise FileNotFoundError(
        f"Raw RNA h5ad not found:\n{rna_raw_path}\n"
        "Run Step 01 first or migrate the existing raw RNA output into the active project."
    )  # no raw input means no preprocessing

assert_path_inside_project(
    rna_raw_path,
    paths.project_root,
    label="raw RNA input",
)  # raw input must be project-local

rna = sc.read_h5ad(rna_raw_path)  # load raw RNA AnnData object

print("\nRaw RNA object:")
print(rna)  # shows n_cells × n_genes


# ============================================================
# 9. Calculate QC metrics before filtering
# ============================================================

print("\nCalculating RNA QC metrics before filtering...")

rna_qc = calculate_rna_qc_metrics(
    adata=rna,  # raw RNA AnnData from Step 01
    mito_prefix=selected_qc_thresholds["mito_prefix"],  # selected mitochondrial prefix from Step 02A/YAML
)  # adds total_counts, n_genes_by_counts, pct_counts_mt

qc_before = summarize_rna_qc(
    rna_qc,
    label="before_filtering",
)  # QC summary before filtering

print("\nRNA QC summary before filtering:")
print(qc_before)  # useful for reporting and sanity checks


# ============================================================
# 10. Filter RNA cells and genes using selected thresholds
# ============================================================

print("\nFiltering RNA cells and genes with selected Step 02A thresholds...")

rna_filtered, filtering_summary = filter_rna_cells_and_genes(
    adata=rna_qc,  # raw RNA object with QC metrics
    min_genes_per_cell=selected_qc_thresholds["min_genes_per_cell"],  # selected low-gene cutoff
    min_cells_per_gene=selected_qc_thresholds["min_cells_per_gene"],  # selected gene detection cutoff
    max_pct_mito=selected_qc_thresholds["max_pct_mito"],  # selected mito cutoff
    max_total_counts=selected_qc_thresholds["max_total_counts"],  # selected high-count cutoff
)  # selected QC filtering

qc_after = summarize_rna_qc(
    rna_filtered,
    label="after_filtering",
)  # QC summary after filtering

filtering_summary["threshold_source"] = threshold_source  # provenance in filtering table
filtering_summary["step02a_selected_thresholds_file"] = str(selected_thresholds_file)  # selected file path

print("\nRNA filtering summary:")
print(filtering_summary)  # cells/genes retained and thresholds used

print("\nRNA QC summary after filtering:")
print(qc_after)  # QC summary after filtering


# ============================================================
# 11. Save filtered RNA object and QC tables
# ============================================================

rna_filtered.write_h5ad(rna_filtered_out)  # save filtered RNA object inside active project

qc_combined = pd.concat(
    [qc_before, qc_after],
    ignore_index=True,
)  # combine before/after QC summaries

write_csv(qc_combined, rna_qc_table, index=False)  # save combined QC table
write_csv(filtering_summary, filtering_summary_out, index=False)  # save filtering summary

print("\nSaved filtered RNA object to:")
print(rna_filtered_out)  # output path

print("\nSaved RNA QC table to:")
print(rna_qc_table)  # output path

print("\nSaved RNA filtering summary to:")
print(filtering_summary_out)  # output path


# ============================================================
# 12. Normalize, log-transform and select HVGs
# ============================================================

print("\nNormalizing RNA, applying log1p, and selecting HVGs...")

rna_hvg, hvg_summary = normalize_log_hvg(
    adata=rna_filtered,  # selected filtered RNA AnnData
    target_sum=norm_cfg["target_sum"],  # normalization target, usually 10000
    apply_log1p=norm_cfg["log1p"],  # validated true above
    n_top_genes=hvg_cfg["n_top_genes"],  # number of HVGs, e.g. 3000
    hvg_flavor=hvg_cfg["flavor"],  # HVG method, e.g. seurat_v3
)  # normalization, log1p, HVG selection

hvg_summary["threshold_source"] = threshold_source  # tracks selected preprocessing provenance

print("\nRNA HVG object:")
print(rna_hvg)  # cells × HVGs

print("\nHVG summary:")
print(hvg_summary)  # HVG selection summary

rna_hvg.write_h5ad(rna_hvg_out)  # save HVG AnnData object
write_csv(hvg_summary, hvg_summary_out, index=False)  # save HVG summary table

print("\nSaved RNA HVG object to:")
print(rna_hvg_out)  # output path

print("\nSaved HVG summary to:")
print(hvg_summary_out)  # output path


# ============================================================
# 13. RNA PCA diagnostics
# ============================================================

print("\nRunning RNA PCA diagnostics...")

rna_pca = run_rna_pca_only(
    adata=rna_hvg,  # HVG-filtered RNA AnnData from Step 02B
    scale_max_value=dim_cfg["scale_max_value"],  # clipping after scaling
    pca_solver=dim_cfg["pca_solver"],  # PCA solver, e.g. arpack
    n_comps=dim_cfg["n_comps"],  # compute enough PCs for elbow plot
    random_state=dim_cfg["random_state"],  # reproducibility
)  # scale and PCA only, no clustering yet

print("\nRNA PCA object:")
print(rna_pca)  # confirms PCA was added


# ============================================================
# 14. Save PCA elbow plot and PCA variance table
# ============================================================

print("\nSaving RNA PCA elbow plot...")

sc.pl.pca_variance_ratio(
    rna_pca,
    n_pcs=dim_cfg["n_comps"],
    log=True,
    show=False,
)  # elbow plot for selecting n_pcs

plt.savefig(
    rna_pca_elbow_out,
    dpi=300,
    bbox_inches="tight",
)  # save PCA elbow plot
plt.show()  # show PCA elbow plot on screen for local inspection
plt.close()  # free memory

pca_variance_table = summarize_pca_variance(
    rna_pca,
)  # make PCA variance/cumulative variance table

write_csv(
    pca_variance_table,
    rna_pca_variance_out,
    index=False,
)  # save table for thesis/report
# ============================================================
# 15. Interactive PCA/neighbors/Leiden diagnostics and final selection
# ============================================================

print("\nPCA elbow plot has been saved and shown.")
print("Use the elbow plot to choose candidate n_pcs values.")
print("Example: 30")
print("Example list for diagnostics: 20,30,40")

n_pcs_input = input(
    "\nEnter n_pcs value or comma-separated list for diagnostics: "
).strip()  # user-guided PCA choice after elbow plot

n_neighbors_input = input(
    "Enter n_neighbors value or comma-separated list for diagnostics: "
).strip()  # user-guided neighbors choice

leiden_resolution_input = input(
    "Enter Leiden resolution value or comma-separated list for diagnostics: "
).strip()  # user-guided Leiden resolution choice

n_pcs_candidate = parse_int_or_int_list(n_pcs_input)  # scalar or list
n_neighbors_candidate = parse_int_or_int_list(n_neighbors_input)  # scalar or list
leiden_resolution_candidate = parse_float_or_float_list(leiden_resolution_input)  # scalar or list

n_pcs_grid = as_list(n_pcs_candidate)  # ensure list for grid
n_neighbors_grid = as_list(n_neighbors_candidate)  # ensure list for grid
leiden_resolution_grid = as_list(leiden_resolution_candidate)  # ensure list for grid

print("\nCandidate clustering grid:")
print(f"n_pcs_grid: {n_pcs_grid}")
print(f"n_neighbors_grid: {n_neighbors_grid}")
print(f"leiden_resolution_grid: {leiden_resolution_grid}")

grid_summary_rows = []  # one summary row per parameter combination
grid_cluster_size_tables = []  # detailed cluster sizes per combination

for n_pcs, n_neighbors, leiden_resolution in product(
    n_pcs_grid,
    n_neighbors_grid,
    leiden_resolution_grid,
):
    resolution_label = safe_resolution_label(leiden_resolution)  # filename-safe label

    cluster_key = (
        f"leiden_pcs{n_pcs}_neighbors{n_neighbors}_res{resolution_label}"
    )  # unique diagnostic cluster key

    print(
        f"\nRunning diagnostic combination: "
        f"n_pcs={n_pcs}, n_neighbors={n_neighbors}, "
        f"leiden_resolution={leiden_resolution}"
    )  # progress log

    adata_grid = rna_pca.copy()  # copy PCA object for this diagnostic run

    sc.pp.neighbors(
        adata_grid,
        n_neighbors=int(n_neighbors),
        n_pcs=int(n_pcs),
        random_state=dim_cfg["random_state"],
    )  # build neighbor graph

    sc.tl.umap(
        adata_grid,
        random_state=dim_cfg["random_state"],
    )  # compute UMAP

    sc.tl.leiden(
        adata_grid,
        resolution=float(leiden_resolution),
        key_added=cluster_key,
        random_state=dim_cfg["random_state"],
    )  # compute Leiden clusters

    grid_summary_rows.append(
        summarize_grid_result(
            adata=adata_grid,
            cluster_key=cluster_key,
            n_pcs=int(n_pcs),
            n_neighbors=int(n_neighbors),
            leiden_resolution=float(leiden_resolution),
        )
    )  # save compact diagnostic summary

    grid_cluster_size_tables.append(
        summarize_leiden_clusters_for_grid(
            adata=adata_grid,
            cluster_key=cluster_key,
            n_pcs=int(n_pcs),
            n_neighbors=int(n_neighbors),
            leiden_resolution=float(leiden_resolution),
        )
    )  # save detailed cluster sizes

    umap_grid_out = (
        clustering_grid_figures_dir
        / f"umap_pcs{n_pcs}_neighbors{n_neighbors}_res{resolution_label}.png"
    )  # diagnostic UMAP output path

    sc.pl.umap(
        adata_grid,
        color=cluster_key,
        legend_loc="on data",
        legend_fontsize=5,
        legend_fontweight="normal",
        frameon=False,
        title=f"PCs={n_pcs}, neighbors={n_neighbors}, Leiden={leiden_resolution}",
        show=False,
    )  # diagnostic UMAP plot

    plt.savefig(
        umap_grid_out,
        dpi=250,
        bbox_inches="tight",
    )  # save diagnostic UMAP

    plt.show()  # show diagnostic UMAP for local inspection

    plt.close()  # free memory

    print(f"Saved diagnostic UMAP to: {umap_grid_out}")  # log output path


grid_summary = pd.DataFrame(grid_summary_rows)  # table of all combinations

grid_cluster_sizes = pd.concat(
    grid_cluster_size_tables,
    ignore_index=True,
)  # detailed cluster sizes

write_csv(
    grid_summary,
    clustering_grid_summary_out,
    index=False,
)  # save diagnostic grid summary

write_csv(
    grid_cluster_sizes,
    clustering_grid_cluster_sizes_out,
    index=False,
)  # save diagnostic cluster sizes

print("\nClustering diagnostic grid summary:")
print(
    grid_summary[
        [
            "n_pcs",
            "n_neighbors",
            "leiden_resolution",
            "n_clusters",
            "min_cluster_size",
            "small_clusters_lt_20_cells",
            "small_clusters_lt_50_cells",
            "small_clusters_lt_100_cells",
        ]
    ]
)  # terminal decision table

print("\nSaved clustering grid summary to:")
print(clustering_grid_summary_out)

print("\nSaved clustering grid cluster sizes to:")
print(clustering_grid_cluster_sizes_out)


# ============================================================
# 15B. Final user selection after diagnostics
# ============================================================

print("\nNow select the final clustering parameters to lock into preprocessing_rna.yaml.")
print("Choose one final value for each parameter based on PCA elbow + UMAPs + cluster sizes.")

selected_n_pcs = int(
    input("\nFinal selected n_pcs: ").strip()
)  # final PCA dimensions

selected_n_neighbors = int(
    input("Final selected n_neighbors: ").strip()
)  # final neighbors

selected_leiden_resolution = float(
    input("Final selected Leiden resolution: ").strip()
)  # final Leiden resolution

selected_leiden_key = f"leiden_{safe_resolution_label(selected_leiden_resolution)}"  # final cluster key

selected_clustering_parameters = pd.DataFrame(
    [
        {
            "n_pcs": selected_n_pcs,
            "n_neighbors": selected_n_neighbors,
            "leiden_resolution": selected_leiden_resolution,
            "leiden_key": selected_leiden_key,
            "selection_source": "interactive Step 02B clustering diagnostics",
            "pca_elbow_plot": str(rna_pca_elbow_out),
            "pca_variance_table": str(rna_pca_variance_out),
            "grid_summary": str(clustering_grid_summary_out),
            "grid_cluster_sizes": str(clustering_grid_cluster_sizes_out),
        }
    ]
)  # final selected parameter audit table

write_csv(
    selected_clustering_parameters,
    selected_clustering_params_out,
    index=False,
)  # save final clustering parameter selection

print("\nSaved selected clustering parameters to:")
print(selected_clustering_params_out)

update_yaml_answer = input(
    "\nUpdate preprocessing_rna.yaml with these final clustering parameters? Type yes/no: "
).strip().lower()  # explicit YAML update decision

if update_yaml_answer == "yes":
    yaml_backup_path = backup_and_update_preprocessing_yaml(
        config=config,
        config_path=config_path,
        selected_n_pcs=selected_n_pcs,
        selected_n_neighbors=selected_n_neighbors,
        selected_leiden_resolution=selected_leiden_resolution,
    )  # backup and update YAML

    print("\nUpdated preprocessing_rna.yaml with final clustering parameters.")
    print("YAML backup saved to:")
    print(yaml_backup_path)

    dim_cfg["n_pcs"] = selected_n_pcs  # update in-memory config
    dim_cfg["n_neighbors"] = selected_n_neighbors  # update in-memory config
    dim_cfg["leiden_resolution"] = selected_leiden_resolution  # update in-memory config
    dim_cfg["leiden_key"] = selected_leiden_key  # update in-memory config

else:
    raise ValueError(
        "Step 02B stopped because final clustering parameters were not written to YAML. "
        "For thesis reproducibility, update YAML before producing the final selected clustered object."
    )


# ============================================================
# 15C. Run final selected neighbors, UMAP and Leiden
# ============================================================

print("\nRunning final selected RNA neighbors, UMAP and Leiden after YAML update...")

rna_clustered = run_rna_neighbors_umap_leiden(
    adata=rna_pca,
    n_neighbors=selected_n_neighbors,
    n_pcs=selected_n_pcs,
    leiden_resolution=selected_leiden_resolution,
    leiden_key=selected_leiden_key,
    random_state=dim_cfg["random_state"],
)  # final selected graph, UMAP and Leiden clustering

# ============================================================
# 16. Save final clustered RNA object and cluster summary
# ============================================================

print("\nSummarizing final RNA Leiden clusters...")

cluster_summary = summarize_clusters(
    rna_clustered,
    cluster_key=dim_cfg["leiden_key"],
)  # create table with number/fraction of cells per Leiden cluster

cluster_summary["threshold_source"] = threshold_source  # record where QC thresholds came from

print("\nRNA cluster summary:")
print(cluster_summary)  # inspect cluster sizes before saving

rna_clustered.write_h5ad(
    rna_clustered_out,
)  # save final RNA object with PCA, neighbors, UMAP and Leiden labels

write_csv(
    cluster_summary,
    cluster_summary_out,
    index=False,
)  # save final cluster summary table

print("\nSaved clustered RNA object to:")
print(rna_clustered_out)  # show clustered h5ad output path

print("\nSaved RNA cluster summary to:")
print(cluster_summary_out)  # show cluster table output path

# ============================================================
# 17. Save final RNA UMAP cluster plot
# ============================================================

print("\nSaving final RNA UMAP cluster plot...")

sc.pl.umap(
    rna_clustered,
    color=dim_cfg["leiden_key"],
    legend_loc="on data",
    legend_fontsize=6,
    legend_fontweight="normal",
    frameon=False,
    show=False,
)  # draw final UMAP colored by selected Leiden clusters

plt.savefig(
    rna_umap_out,
    dpi=300,
    bbox_inches="tight",
)  # save final UMAP cluster figure
plt.show()  # show final UMAP on screen for local inspection
plt.close()  # free matplotlib memory

print("\nSaved final RNA UMAP cluster plot to:")
print(rna_umap_out)  # show UMAP output path

# ============================================================
# 18. Final summary
# ============================================================

print("\nDone.")
print("Step 02B selected RNA preprocessing rerun completed.")
print("Step 02A selected thresholds were validated against preprocessing_rna.yaml.")
print("All outputs were written inside the active project workspace.")
print(f"Active project: {paths.project_id}")
print(f"Project root: {paths.project_root}")
print("\nNext step:")
print("Step 03C — RNA PCA / neighbors / Leiden diagnostics.")