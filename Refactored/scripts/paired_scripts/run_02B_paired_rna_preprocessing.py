# -*- coding: utf-8 -*-
"""
Run 02B — Paired RNA preprocessing, synchronized filtering,
PCA/UMAP/Leiden diagnostics and final clustering.

Inputs
------
1. Raw paired RNA:
   data/processed/rna/pbmc_multiome_rna_raw.h5ad

2. Raw paired ATAC:
   data/processed/atac/pbmc_multiome_atac_raw.h5ad

3. Raw full Multiome:
   data/processed/multiome/pbmc_multiome_raw.h5ad

4. Step 02A selected thresholds:
   qc_diagnostics.outputs.selected_thresholds

5. Step 02A Scrublet scores:
   qc_diagnostics.outputs.doublet_scores_table

6. configs/preprocessing_rna.yaml

Workflow
--------
raw paired RNA / ATAC / Multiome
→ validate identical barcodes and order
→ load final Step 02A QC and Scrublet decisions
→ apply one shared cell-retention mask
→ save synchronized RNA / ATAC / Multiome objects
→ filter RNA genes
→ calculate HVGs from raw counts
→ normalize total counts and log1p
→ preserve full log-normalized genes in rna.raw
→ subset to HVGs
→ scale and PCA
→ interactive neighbors / UMAP / Leiden grid
→ final clustering selection
→ update preprocessing_rna.yaml
→ save final clustered RNA h5ad

Important
---------
- Scrublet is NOT rerun.
- ATAC peaks are NOT filtered here.
- RNA, ATAC and Multiome always retain identical cells and cell order.
- Final RNA clustered object is HVG-based, but rna.raw contains the
  full log-normalized gene matrix for marker annotation in Step 03A.
"""

from __future__ import annotations

from datetime import datetime
from itertools import product
from pathlib import Path
from typing import Any
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml


# ============================================================
# 0. Locate Refactored root robustly
# ============================================================

def find_refactored_root(start_path: Path) -> Path:
    """Find the parent directory containing both src/ and projects/."""

    start_path = start_path.resolve()

    for candidate in [start_path, *start_path.parents]:
        if (
            (candidate / "src").is_dir()
            and (candidate / "projects").is_dir()
        ):
            return candidate

    raise FileNotFoundError(
        "Could not locate the Refactored root.\n"
        f"Search started from:\n{start_path}\n\n"
        "Expected a parent directory containing both src/ and projects/."
    )


try:
    script_path = Path(__file__).resolve()
    refactor_dir = find_refactored_root(script_path.parent)
except NameError:
    script_path = Path("run_02B_paired_rna_preprocessing.py")
    refactor_dir = find_refactored_root(Path.cwd())

src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))


# ============================================================
# 1. Project helpers
# ============================================================

from pbmcgrn.config import load_yaml
from pbmcgrn.io import write_csv, write_json
from pbmcgrn.paths import ProjectPaths
from pbmcgrn.preprocessing.rna import calculate_rna_qc_metrics


# ============================================================
# 2. Generic helper functions
# ============================================================

def require_section(
    config: dict[str, Any],
    section: str,
) -> dict[str, Any]:
    """Return a required YAML section."""

    if section not in config:
        raise KeyError(
            f"Missing required YAML section: {section}"
        )

    value = config[section]

    if not isinstance(value, dict):
        raise TypeError(
            f"YAML section '{section}' must be a dictionary."
        )

    return value


def require_keys(
    dictionary: dict[str, Any],
    required_keys: set[str],
    section_name: str,
) -> None:
    """Fail if required keys are missing."""

    missing = required_keys.difference(dictionary)

    if missing:
        raise KeyError(
            f"Missing keys in {section_name}: {sorted(missing)}"
        )


def require_scalar(
    value: Any,
    parameter_name: str,
) -> Any:
    """Require one final scalar value rather than a diagnostic list."""

    if isinstance(value, list):
        if len(value) != 1:
            raise ValueError(
                f"{parameter_name} must contain one final value, "
                f"but received: {value}"
            )

        return value[0]

    return value


def assert_path_inside_project(
    path: Path,
    project_root: Path,
    label: str,
) -> None:
    """Fail when a path is outside the active project."""

    resolved_path = Path(path).resolve()
    resolved_root = Path(project_root).resolve()

    try:
        resolved_path.relative_to(resolved_root)
    except ValueError as exc:
        raise ValueError(
            f"{label} is outside the active project.\n"
            f"Path: {resolved_path}\n"
            f"Project root: {resolved_root}"
        ) from exc


def ensure_parent(path: Path) -> None:
    """Create an output parent directory."""

    Path(path).parent.mkdir(
        parents=True,
        exist_ok=True,
    )


def yes_response(text: str) -> bool:
    """Recognize common affirmative responses."""

    return text.strip().lower() in {
        "y",
        "yes",
        "ν",
        "ναι",
    }


def safe_resolution_label(value: float) -> str:
    """Convert 0.8 to 0_8 for filenames and obs keys."""

    return str(float(value)).replace(".", "_")


def parse_int_list(
    user_text: str,
    default: int,
) -> list[int]:
    """Parse one integer or comma-separated integers."""

    if user_text.strip() == "":
        return [int(default)]

    values = [
        int(item.strip())
        for item in user_text.split(",")
        if item.strip() != ""
    ]

    if not values:
        raise ValueError("No valid integer values were supplied.")

    return list(dict.fromkeys(values))


def parse_float_list(
    user_text: str,
    default: float,
) -> list[float]:
    """Parse one float or comma-separated floats."""

    if user_text.strip() == "":
        return [float(default)]

    values = [
        float(item.strip())
        for item in user_text.split(",")
        if item.strip() != ""
    ]

    if not values:
        raise ValueError("No valid float values were supplied.")

    return list(dict.fromkeys(values))


def save_current_figure(
    output_path: Path,
    dpi: int = 300,
    show: bool = True,
) -> None:
    """Save, optionally show and close the current figure."""

    ensure_parent(output_path)

    plt.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
    )

    if show:
        plt.show()

    plt.close("all")


# ============================================================
# 3. Paired validation
# ============================================================

def validate_paired_objects(
    rna: sc.AnnData,
    atac: sc.AnnData,
    multiome: sc.AnnData,
    stage: str,
) -> pd.DataFrame:
    """Validate identical paired cells and ordering."""

    same_cell_count = (
        rna.n_obs
        == atac.n_obs
        == multiome.n_obs
    )

    same_barcode_set = (
        set(rna.obs_names)
        == set(atac.obs_names)
        == set(multiome.obs_names)
    )

    same_barcode_order = (
        rna.obs_names.equals(atac.obs_names)
        and rna.obs_names.equals(multiome.obs_names)
    )

    rna_duplicates = int(
        rna.obs_names.duplicated().sum()
    )

    atac_duplicates = int(
        atac.obs_names.duplicated().sum()
    )

    multiome_duplicates = int(
        multiome.obs_names.duplicated().sum()
    )

    pairing_confirmed = (
        same_cell_count
        and same_barcode_set
        and same_barcode_order
        and rna_duplicates == 0
        and atac_duplicates == 0
        and multiome_duplicates == 0
    )

    table = pd.DataFrame(
        [
            {
                "stage": stage,
                "rna_n_cells": int(rna.n_obs),
                "atac_n_cells": int(atac.n_obs),
                "multiome_n_cells": int(multiome.n_obs),
                "same_cell_count": bool(same_cell_count),
                "same_barcode_set": bool(same_barcode_set),
                "same_barcode_order": bool(same_barcode_order),
                "rna_duplicate_barcodes": rna_duplicates,
                "atac_duplicate_barcodes": atac_duplicates,
                "multiome_duplicate_barcodes": multiome_duplicates,
                "pairing_confirmed": bool(pairing_confirmed),
            }
        ]
    )

    if not pairing_confirmed:
        raise RuntimeError(
            f"Paired validation failed at stage '{stage}'.\n"
            f"{table.to_string(index=False)}"
        )

    return table


# ============================================================
# 4. Step 02A handoff validation
# ============================================================

def load_selected_threshold_row(
    selected_thresholds_path: Path,
) -> pd.Series:
    """Load the single selected row created by paired Step 02A."""

    if not selected_thresholds_path.exists():
        raise FileNotFoundError(
            "Step 02A selected-threshold table not found:\n"
            f"{selected_thresholds_path}"
        )

    selected = pd.read_csv(selected_thresholds_path)

    if selected.shape[0] != 1:
        raise ValueError(
            "Expected exactly one selected-threshold row.\n"
            f"Rows found: {selected.shape[0]}\n"
            f"File: {selected_thresholds_path}"
        )

    required_columns = {
        "mito_prefix",
        "min_genes_per_cell",
        "min_cells_per_gene",
        "max_pct_mito",
        "max_total_counts",
        "doublet_score_threshold",
        "filter_predicted_doublets",
        "yaml_updated",
    }

    missing = required_columns.difference(selected.columns)

    if missing:
        raise KeyError(
            "Step 02A selected-threshold table is missing columns:\n"
            f"{sorted(missing)}"
        )

    row = selected.iloc[0]

    if not bool(row["yaml_updated"]):
        raise ValueError(
            "Step 02A reports yaml_updated=False.\n"
            "The final thresholds must be written to preprocessing_rna.yaml "
            "before Step 02B."
        )

    return row


def compare_yaml_and_step02a(
    rna_cfg: dict[str, Any],
    doublet_cfg: dict[str, Any],
    selected_row: pd.Series,
) -> pd.DataFrame:
    """Compare final YAML values with the Step 02A handoff."""

    yaml_values = {
        "mito_prefix": str(
            require_scalar(
                rna_cfg["mito_prefix"],
                "rna_preprocessing.mito_prefix",
            )
        ),
        "min_genes_per_cell": int(
            require_scalar(
                rna_cfg["min_genes_per_cell"],
                "rna_preprocessing.min_genes_per_cell",
            )
        ),
        "min_cells_per_gene": int(
            require_scalar(
                rna_cfg["min_cells_per_gene"],
                "rna_preprocessing.min_cells_per_gene",
            )
        ),
        "max_pct_mito": float(
            require_scalar(
                rna_cfg["max_pct_mito"],
                "rna_preprocessing.max_pct_mito",
            )
        ),
        "max_total_counts": float(
            require_scalar(
                rna_cfg["max_total_counts"],
                "rna_preprocessing.max_total_counts",
            )
        ),
        "doublet_score_threshold": float(
            require_scalar(
                doublet_cfg["score_threshold"],
                "doublet_detection.score_threshold",
            )
        ),
        "filter_predicted_doublets": bool(
            doublet_cfg["filter_predicted_doublets"]
        ),
    }

    step02a_values = {
        "mito_prefix": str(selected_row["mito_prefix"]),
        "min_genes_per_cell": int(
            selected_row["min_genes_per_cell"]
        ),
        "min_cells_per_gene": int(
            selected_row["min_cells_per_gene"]
        ),
        "max_pct_mito": float(
            selected_row["max_pct_mito"]
        ),
        "max_total_counts": float(
            selected_row["max_total_counts"]
        ),
        "doublet_score_threshold": float(
            selected_row["doublet_score_threshold"]
        ),
        "filter_predicted_doublets": bool(
            selected_row["filter_predicted_doublets"]
        ),
    }

    rows = []

    for parameter, yaml_value in yaml_values.items():
        selected_value = step02a_values[parameter]

        if isinstance(yaml_value, float):
            matches = bool(
                np.isclose(
                    yaml_value,
                    float(selected_value),
                    rtol=0,
                    atol=1e-10,
                )
            )
        else:
            matches = yaml_value == selected_value

        rows.append(
            {
                "parameter": parameter,
                "yaml_value": yaml_value,
                "step02a_value": selected_value,
                "matches": matches,
            }
        )

    comparison = pd.DataFrame(rows)

    if not bool(comparison["matches"].all()):
        raise ValueError(
            "Step 02A handoff does not match preprocessing_rna.yaml.\n\n"
            f"{comparison.to_string(index=False)}"
        )

    return comparison


def load_scrublet_scores(
    scores_path: Path,
    expected_obs_names: pd.Index,
) -> pd.DataFrame:
    """Load and align Step 02A Scrublet scores to raw paired cells."""

    if not scores_path.exists():
        raise FileNotFoundError(
            "Step 02A Scrublet score table not found:\n"
            f"{scores_path}"
        )

    scores = pd.read_csv(
        scores_path,
        dtype={"cell_id": str},
    )

    required_columns = {
        "cell_id",
        "doublet_score",
    }

    missing = required_columns.difference(scores.columns)

    if missing:
        raise KeyError(
            "Scrublet score table is missing columns:\n"
            f"{sorted(missing)}"
        )

    if scores["cell_id"].duplicated().any():
        duplicate_count = int(
            scores["cell_id"].duplicated().sum()
        )

        raise ValueError(
            "Scrublet score table has duplicated cell IDs.\n"
            f"Duplicates: {duplicate_count}"
        )

    expected = pd.Index(
        expected_obs_names.astype(str),
        name="cell_id",
    )

    observed = pd.Index(
        scores["cell_id"].astype(str),
        name="cell_id",
    )

    missing_cells = expected.difference(observed)
    extra_cells = observed.difference(expected)

    if len(missing_cells) > 0 or len(extra_cells) > 0:
        raise ValueError(
            "Scrublet scores do not match the raw paired cells.\n"
            f"Missing cells: {len(missing_cells)}\n"
            f"Extra cells: {len(extra_cells)}"
        )

    scores = (
        scores
        .set_index("cell_id")
        .loc[expected]
        .copy()
    )

    if not scores.index.equals(expected):
        raise RuntimeError(
            "Scrublet scores could not be aligned to raw cell order."
        )

    return scores


# ============================================================
# 5. QC and clustering summaries
# ============================================================

def summarize_qc(
    adata: sc.AnnData,
    stage: str,
) -> pd.DataFrame:
    """Summarize main RNA QC metrics."""

    rows = []

    for metric in [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
    ]:
        values = adata.obs[metric].astype(float)

        rows.append(
            {
                "stage": stage,
                "metric": metric,
                "n_cells": int(values.shape[0]),
                "mean": float(values.mean()),
                "median": float(values.median()),
                "std": float(values.std()),
                "min": float(values.min()),
                "q01": float(values.quantile(0.01)),
                "q05": float(values.quantile(0.05)),
                "q25": float(values.quantile(0.25)),
                "q75": float(values.quantile(0.75)),
                "q95": float(values.quantile(0.95)),
                "q99": float(values.quantile(0.99)),
                "max": float(values.max()),
            }
        )

    return pd.DataFrame(rows)


def summarize_grid_result(
    adata: sc.AnnData,
    cluster_key: str,
    n_pcs: int,
    n_neighbors: int,
    resolution: float,
) -> dict[str, Any]:
    """Summarize one clustering-grid combination."""

    counts = (
        adata.obs[cluster_key]
        .astype(str)
        .value_counts()
    )

    return {
        "n_pcs": int(n_pcs),
        "n_neighbors": int(n_neighbors),
        "leiden_resolution": float(resolution),
        "cluster_key": cluster_key,
        "n_clusters": int(counts.shape[0]),
        "min_cluster_size": int(counts.min()),
        "median_cluster_size": float(counts.median()),
        "max_cluster_size": int(counts.max()),
        "clusters_lt_20": int((counts < 20).sum()),
        "clusters_lt_50": int((counts < 50).sum()),
        "clusters_lt_100": int((counts < 100).sum()),
    }


def summarize_cluster_sizes(
    adata: sc.AnnData,
    cluster_key: str,
    n_pcs: int,
    n_neighbors: int,
    resolution: float,
) -> pd.DataFrame:
    """Create detailed cluster sizes for one parameter combination."""

    table = (
        adata.obs[cluster_key]
        .astype(str)
        .value_counts()
        .rename_axis("cluster")
        .reset_index(name="n_cells")
    )

    table["fraction_cells"] = (
        table["n_cells"] / adata.n_obs
    )

    table["n_pcs"] = int(n_pcs)
    table["n_neighbors"] = int(n_neighbors)
    table["leiden_resolution"] = float(resolution)
    table["cluster_key"] = cluster_key

    table["_sort"] = pd.to_numeric(
        table["cluster"],
        errors="coerce",
    )

    table = (
        table
        .sort_values(["_sort", "cluster"])
        .drop(columns="_sort")
        .reset_index(drop=True)
    )

    return table


def build_pca_variance_table(
    adata: sc.AnnData,
) -> pd.DataFrame:
    """Build PCA variance and cumulative-variance table."""

    if "pca" not in adata.uns:
        raise KeyError("PCA metadata missing from adata.uns['pca'].")

    variance = np.asarray(
        adata.uns["pca"]["variance"],
        dtype=float,
    )

    variance_ratio = np.asarray(
        adata.uns["pca"]["variance_ratio"],
        dtype=float,
    )

    return pd.DataFrame(
        {
            "pc": np.arange(1, len(variance) + 1),
            "variance": variance,
            "variance_ratio": variance_ratio,
            "cumulative_variance_ratio": np.cumsum(
                variance_ratio
            ),
        }
    )


def backup_and_update_clustering_yaml(
    config: dict[str, Any],
    config_path: Path,
    n_pcs: int,
    n_neighbors: int,
    resolution: float,
) -> tuple[Path, str]:
    """Backup YAML and write final clustering parameters."""

    timestamp = datetime.now().strftime(
        "%Y%m%d_%H%M%S"
    )

    backup_path = config_path.with_name(
        f"{config_path.stem}"
        f".before_step02B_clustering_{timestamp}"
        f"{config_path.suffix}"
    )

    shutil.copy2(
        config_path,
        backup_path,
    )

    leiden_key = (
        f"leiden_{safe_resolution_label(resolution)}"
    )

    dim_cfg = config["dimensionality_reduction"]

    dim_cfg["n_pcs"] = int(n_pcs)
    dim_cfg["n_neighbors"] = int(n_neighbors)
    dim_cfg["leiden_resolution"] = float(resolution)
    dim_cfg["leiden_key"] = leiden_key

    with open(
        config_path,
        "w",
        encoding="utf-8",
    ) as handle:
        yaml.safe_dump(
            config,
            handle,
            sort_keys=False,
            allow_unicode=True,
        )

    return backup_path, leiden_key


# ============================================================
# 6. Main workflow
# ============================================================

def main() -> None:
    # --------------------------------------------------------
    # 6.1 Active project
    # --------------------------------------------------------

    paths = ProjectPaths(
        refactor_root=refactor_dir
    )

    paths.ensure_output_dirs()

    expected_project_id = "pbmc_10k_multiome_v1"

    print("\n" + "=" * 80)
    print("RUN 02B — PAIRED RNA PREPROCESSING AND CLUSTERING")
    print("=" * 80)

    print("\nScript:")
    print(script_path)

    print("\nRefactored root:")
    print(refactor_dir)

    print("\nActive project:")
    print(paths.project_id)

    print("\nProject root:")
    print(paths.project_root)

    if paths.project_id != expected_project_id:
        raise RuntimeError(
            "Wrong active project.\n"
            f"Expected: {expected_project_id}\n"
            f"Observed: {paths.project_id}"
        )

    # --------------------------------------------------------
    # 6.2 Load YAML
    # --------------------------------------------------------

    config_path = (
        paths.configs
        / "preprocessing_rna.yaml"
    )

    if not config_path.exists():
        raise FileNotFoundError(
            "preprocessing_rna.yaml not found:\n"
            f"{config_path}"
        )

    config = load_yaml(config_path)

    if config is None:
        raise ValueError(
            "preprocessing_rna.yaml is empty or invalid."
        )

    rna_cfg = require_section(
        config,
        "rna_preprocessing",
    )

    doublet_cfg = require_section(
        config,
        "doublet_detection",
    )

    norm_cfg = require_section(
        config,
        "normalization",
    )

    hvg_cfg = require_section(
        config,
        "highly_variable_genes",
    )

    dim_cfg = require_section(
        config,
        "dimensionality_reduction",
    )

    qc_diag_cfg = require_section(
        config,
        "qc_diagnostics",
    )

    qc_outputs_cfg = require_section(
        qc_diag_cfg,
        "outputs",
    )

    outputs_cfg = require_section(
        config,
        "outputs",
    )

    require_keys(
        rna_cfg,
        {
            "mito_prefix",
            "min_genes_per_cell",
            "min_cells_per_gene",
            "max_pct_mito",
            "max_total_counts",
        },
        "rna_preprocessing",
    )

    require_keys(
        doublet_cfg,
        {
            "enabled",
            "score_threshold",
            "filter_predicted_doublets",
        },
        "doublet_detection",
    )

    require_keys(
        norm_cfg,
        {
            "target_sum",
            "log1p",
        },
        "normalization",
    )

    require_keys(
        hvg_cfg,
        {
            "n_top_genes",
            "flavor",
        },
        "highly_variable_genes",
    )

    require_keys(
        dim_cfg,
        {
            "scale_max_value",
            "pca_solver",
            "n_comps",
            "n_pcs",
            "n_neighbors",
            "leiden_resolution",
            "random_state",
        },
        "dimensionality_reduction",
    )

    if bool(norm_cfg["log1p"]) is not True:
        raise ValueError(
            "This Step 02B requires normalization.log1p: true."
        )

    if bool(doublet_cfg["enabled"]) and (
        doublet_cfg["score_threshold"] is None
    ):
        raise ValueError(
            "doublet_detection.score_threshold is null.\n"
            "Run or finalize Step 02A first."
        )

    # --------------------------------------------------------
    # 6.3 Resolve input paths
    # --------------------------------------------------------

    rna_raw_path = (
        paths.processed_rna
        / "pbmc_multiome_rna_raw.h5ad"
    )

    atac_raw_path = (
        paths.project_root
        / "data"
        / "processed"
        / "atac"
        / "pbmc_multiome_atac_raw.h5ad"
    )

    multiome_raw_path = (
        paths.project_root
        / "data"
        / "processed"
        / "multiome"
        / "pbmc_multiome_raw.h5ad"
    )

    selected_thresholds_path = (
        paths.prepare_output_path(
            qc_outputs_cfg["selected_thresholds"],
            label="step02A_selected_thresholds",
        )
    )

    scrublet_scores_path = (
        paths.prepare_output_path(
            qc_outputs_cfg["doublet_scores_table"],
            label="step02A_doublet_scores",
        )
    )

    input_paths = {
        "rna_raw": rna_raw_path,
        "atac_raw": atac_raw_path,
        "multiome_raw": multiome_raw_path,
        "selected_thresholds": selected_thresholds_path,
        "scrublet_scores": scrublet_scores_path,
        "preprocessing_config": config_path,
    }

    for label, path in input_paths.items():
        assert_path_inside_project(
            path,
            paths.project_root,
            label,
        )

        if not path.exists():
            raise FileNotFoundError(
                f"Required input not found ({label}):\n{path}"
            )

    # --------------------------------------------------------
    # 6.4 Resolve output paths
    # --------------------------------------------------------

    required_output_keys = {
        "paired_rna_cell_filtered_h5ad",
        "paired_atac_after_rna_qc_h5ad",
        "paired_multiome_after_rna_qc_h5ad",
        "rna_filtered_h5ad",
        "rna_hvg_h5ad",
        "rna_clustered_h5ad",
        "rna_qc_table",
        "rna_filtering_summary",
        "paired_filtering_summary",
        "pairing_validation",
        "retained_cell_ids",
        "removed_cell_ids",
        "doublets_removed",
        "rna_hvg_summary",
        "rna_cluster_summary",
        "rna_pca_elbow",
        "rna_pca_variance_table",
        "rna_pca_qc_plot",
        "rna_umap_qc_plot",
        "rna_umap_clusters",
        "rna_leiden_diagnostics",
        "rna_umap_diagnostics_dir",
        "manifest",
    }

    missing_output_keys = (
        required_output_keys
        .difference(outputs_cfg)
    )

    if missing_output_keys:
        raise KeyError(
            "Missing Step 02B output keys in preprocessing_rna.yaml:\n"
            f"{sorted(missing_output_keys)}"
        )

    output_paths = {
        key: paths.prepare_output_path(
            outputs_cfg[key],
            label=key,
        )
        for key in required_output_keys
    }

    clustering_grid_summary_out = (
        paths.prepare_output_path(
            "reports/tables/rna/step02B/"
            "rna_clustering_parameter_grid_summary.csv",
            label="clustering_grid_summary",
        )
    )

    clustering_grid_sizes_out = (
        paths.prepare_output_path(
            "reports/tables/rna/step02B/"
            "rna_clustering_parameter_grid_cluster_sizes.csv",
            label="clustering_grid_sizes",
        )
    )

    selected_clustering_out = (
        paths.prepare_output_path(
            "reports/tables/rna/step02B/"
            "rna_selected_clustering_parameters.csv",
            label="selected_clustering_parameters",
        )
    )

    threshold_validation_out = (
        paths.prepare_output_path(
            "reports/tables/rna/step02B/"
            "step02A_to_step02B_threshold_validation.csv",
            label="threshold_validation",
        )
    )

    all_output_paths = {
        **output_paths,
        "clustering_grid_summary": clustering_grid_summary_out,
        "clustering_grid_sizes": clustering_grid_sizes_out,
        "selected_clustering": selected_clustering_out,
        "threshold_validation": threshold_validation_out,
    }

    for label, path in all_output_paths.items():
        assert_path_inside_project(
            path,
            paths.project_root,
            label,
        )

        if path.suffix:
            ensure_parent(path)
        else:
            path.mkdir(
                parents=True,
                exist_ok=True,
            )

    output_paths[
        "rna_umap_diagnostics_dir"
    ].mkdir(
        parents=True,
        exist_ok=True,
    )

    print("\nResolved inputs:")
    for label, path in input_paths.items():
        print(f"{label}: {path}")

    print("\nResolved outputs:")
    for label, path in all_output_paths.items():
        print(f"{label}: {path}")

    # --------------------------------------------------------
    # 6.5 Load paired raw objects
    # --------------------------------------------------------

    print("\nLoading raw paired objects...")

    rna = sc.read_h5ad(rna_raw_path)
    atac = sc.read_h5ad(atac_raw_path)
    multiome = sc.read_h5ad(multiome_raw_path)

    print("\nRaw RNA:")
    print(rna)

    print("\nRaw ATAC:")
    print(atac)

    print("\nRaw Multiome:")
    print(multiome)

    pairing_before = validate_paired_objects(
        rna,
        atac,
        multiome,
        stage="before_step02B_filtering",
    )

    print("\nPairing validation before filtering:")
    print(pairing_before.to_string(index=False))

    # --------------------------------------------------------
    # 6.6 Validate Step 02A handoff
    # --------------------------------------------------------

    selected_row = load_selected_threshold_row(
        selected_thresholds_path
    )

    threshold_validation = (
        compare_yaml_and_step02a(
            rna_cfg=rna_cfg,
            doublet_cfg=doublet_cfg,
            selected_row=selected_row,
        )
    )

    write_csv(
        threshold_validation,
        threshold_validation_out,
        index=False,
    )

    print("\nStep 02A → Step 02B threshold validation:")
    print(threshold_validation.to_string(index=False))

    scrublet_scores = load_scrublet_scores(
        scrublet_scores_path,
        expected_obs_names=rna.obs_names,
    )

    # Attach Step 02A Scrublet evidence to all paired objects.
    doublet_scores = (
        scrublet_scores["doublet_score"]
        .astype(float)
        .to_numpy()
    )

    selected_doublet_threshold = float(
        selected_row["doublet_score_threshold"]
    )

    predicted_doublet = (
        doublet_scores
        >= selected_doublet_threshold
    )

    for adata in [rna, atac, multiome]:
        adata.obs["doublet_score"] = doublet_scores
        adata.obs["predicted_doublet"] = predicted_doublet

    # --------------------------------------------------------
    # 6.7 Calculate RNA QC metrics
    # --------------------------------------------------------

    selected_mito_prefix = str(
        selected_row["mito_prefix"]
    )

    rna_qc = calculate_rna_qc_metrics(
        adata=rna,
        mito_prefix=selected_mito_prefix,
    )

    # calculate_rna_qc_metrics may return a copy.
    rna = rna_qc

    # Restore Scrublet columns explicitly.
    rna.obs["doublet_score"] = doublet_scores
    rna.obs["predicted_doublet"] = predicted_doublet

    qc_before = summarize_qc(
        rna,
        stage="before_filtering",
    )

    # --------------------------------------------------------
    # 6.8 Build one final shared cell mask
    # --------------------------------------------------------

    min_genes = int(
        selected_row["min_genes_per_cell"]
    )

    min_cells_per_gene = int(
        selected_row["min_cells_per_gene"]
    )

    max_pct_mito = float(
        selected_row["max_pct_mito"]
    )

    max_total_counts = float(
        selected_row["max_total_counts"]
    )

    min_genes_pass = (
        rna.obs["n_genes_by_counts"]
        >= min_genes
    )

    mito_pass = (
        rna.obs["pct_counts_mt"]
        < max_pct_mito
    )

    total_counts_pass = (
        rna.obs["total_counts"]
        < max_total_counts
    )

    rna_qc_pass = (
        min_genes_pass
        & mito_pass
        & total_counts_pass
    )

    if bool(
        selected_row["filter_predicted_doublets"]
    ):
        singlet_pass = ~rna.obs[
            "predicted_doublet"
        ].astype(bool)
    else:
        singlet_pass = pd.Series(
            True,
            index=rna.obs_names,
        )

    final_keep_mask = (
        rna_qc_pass
        & singlet_pass
    )

    if int(final_keep_mask.sum()) == 0:
        raise RuntimeError(
            "Final Step 02B cell mask retained zero cells."
        )

    # Add audit columns to all objects before subsetting.
    for adata in [rna, atac, multiome]:
        adata.obs["pass_min_genes"] = (
            min_genes_pass.to_numpy()
        )
        adata.obs["pass_max_pct_mito"] = (
            mito_pass.to_numpy()
        )
        adata.obs["pass_max_total_counts"] = (
            total_counts_pass.to_numpy()
        )
        adata.obs["pass_rna_qc"] = (
            rna_qc_pass.to_numpy()
        )
        adata.obs["pass_singlet_filter"] = (
            singlet_pass.to_numpy()
        )
        adata.obs["retained_step02B"] = (
            final_keep_mask.to_numpy()
        )

    retained_cell_ids = pd.DataFrame(
        {
            "cell_id": rna.obs_names[
                final_keep_mask.to_numpy()
            ].astype(str)
        }
    )

    removed_obs = rna.obs.loc[
        ~final_keep_mask
    ].copy()

    removed_cell_ids = (
        removed_obs[
            [
                "total_counts",
                "n_genes_by_counts",
                "pct_counts_mt",
                "doublet_score",
                "predicted_doublet",
                "pass_min_genes",
                "pass_max_pct_mito",
                "pass_max_total_counts",
                "pass_rna_qc",
                "pass_singlet_filter",
            ]
        ]
        .reset_index(names="cell_id")
    )

    doublets_removed = (
        rna.obs.loc[
            rna.obs["predicted_doublet"].astype(bool)
            & ~final_keep_mask,
            [
                "total_counts",
                "n_genes_by_counts",
                "pct_counts_mt",
                "doublet_score",
                "predicted_doublet",
                "pass_rna_qc",
            ],
        ]
        .reset_index(names="cell_id")
    )

    # --------------------------------------------------------
    # 6.9 Apply exactly the same cell mask to all modalities
    # --------------------------------------------------------

    keep_array = final_keep_mask.to_numpy()

    rna_cell_filtered = (
        rna[keep_array, :]
        .copy()
    )

    atac_after_rna_qc = (
        atac[keep_array, :]
        .copy()
    )

    multiome_after_rna_qc = (
        multiome[keep_array, :]
        .copy()
    )

    pairing_after = validate_paired_objects(
        rna_cell_filtered,
        atac_after_rna_qc,
        multiome_after_rna_qc,
        stage="after_step02B_cell_filtering",
    )

    print("\nPairing validation after filtering:")
    print(pairing_after.to_string(index=False))

    # --------------------------------------------------------
    # 6.10 RNA gene filtering
    # --------------------------------------------------------

    genes_before = int(
        rna_cell_filtered.n_vars
    )

    detected_cells_per_gene = np.asarray(
        (rna_cell_filtered.X > 0).sum(axis=0)
    ).ravel()

    gene_keep_mask = (
        detected_cells_per_gene
        >= min_cells_per_gene
    )

    if int(gene_keep_mask.sum()) == 0:
        raise RuntimeError(
            "Gene filtering retained zero RNA genes."
        )

    rna_filtered = (
        rna_cell_filtered[
            :,
            gene_keep_mask,
        ]
        .copy()
    )

    genes_after = int(
        rna_filtered.n_vars
    )

    qc_after = summarize_qc(
        rna_filtered,
        stage="after_cell_and_gene_filtering",
    )

    # --------------------------------------------------------
    # 6.11 Save synchronized paired filtered objects
    # --------------------------------------------------------

    rna_cell_filtered.write_h5ad(
        output_paths[
            "paired_rna_cell_filtered_h5ad"
        ]
    )

    atac_after_rna_qc.write_h5ad(
        output_paths[
            "paired_atac_after_rna_qc_h5ad"
        ]
    )

    multiome_after_rna_qc.write_h5ad(
        output_paths[
            "paired_multiome_after_rna_qc_h5ad"
        ]
    )

    rna_filtered.write_h5ad(
        output_paths["rna_filtered_h5ad"]
    )

    write_csv(
        retained_cell_ids,
        output_paths["retained_cell_ids"],
        index=False,
    )

    write_csv(
        removed_cell_ids,
        output_paths["removed_cell_ids"],
        index=False,
    )

    write_csv(
        doublets_removed,
        output_paths["doublets_removed"],
        index=False,
    )

    pairing_validation = pd.concat(
        [
            pairing_before,
            pairing_after,
        ],
        ignore_index=True,
    )

    write_csv(
        pairing_validation,
        output_paths["pairing_validation"],
        index=False,
    )

    # --------------------------------------------------------
    # 6.12 Filtering summaries
    # --------------------------------------------------------

    filtering_summary = pd.DataFrame(
        [
            {
                "cells_before": int(rna.n_obs),
                "cells_passing_rna_qc": int(
                    rna_qc_pass.sum()
                ),
                "predicted_doublets": int(
                    predicted_doublet.sum()
                ),
                "cells_retained_final": int(
                    final_keep_mask.sum()
                ),
                "cells_removed_final": int(
                    (~final_keep_mask).sum()
                ),
                "cell_retention_fraction": float(
                    final_keep_mask.mean()
                ),
                "genes_before": genes_before,
                "genes_after": genes_after,
                "genes_removed": (
                    genes_before - genes_after
                ),
                "gene_retention_fraction": float(
                    genes_after / genes_before
                ),
                "min_genes_per_cell": min_genes,
                "min_cells_per_gene": min_cells_per_gene,
                "max_pct_mito": max_pct_mito,
                "max_total_counts": max_total_counts,
                "doublet_score_threshold": (
                    selected_doublet_threshold
                ),
                "filter_predicted_doublets": bool(
                    selected_row[
                        "filter_predicted_doublets"
                    ]
                ),
            }
        ]
    )

    paired_filtering_summary = pd.DataFrame(
        [
            {
                "modality": "RNA_cell_filtered",
                "n_cells": int(
                    rna_cell_filtered.n_obs
                ),
                "n_features": int(
                    rna_cell_filtered.n_vars
                ),
            },
            {
                "modality": "ATAC_after_RNA_QC",
                "n_cells": int(
                    atac_after_rna_qc.n_obs
                ),
                "n_features": int(
                    atac_after_rna_qc.n_vars
                ),
            },
            {
                "modality": "Multiome_after_RNA_QC",
                "n_cells": int(
                    multiome_after_rna_qc.n_obs
                ),
                "n_features": int(
                    multiome_after_rna_qc.n_vars
                ),
            },
            {
                "modality": "RNA_gene_filtered",
                "n_cells": int(
                    rna_filtered.n_obs
                ),
                "n_features": int(
                    rna_filtered.n_vars
                ),
            },
        ]
    )

    write_csv(
        filtering_summary,
        output_paths["rna_filtering_summary"],
        index=False,
    )

    write_csv(
        paired_filtering_summary,
        output_paths["paired_filtering_summary"],
        index=False,
    )

    write_csv(
        pd.concat(
            [qc_before, qc_after],
            ignore_index=True,
        ),
        output_paths["rna_qc_table"],
        index=False,
    )

    print("\nFiltering summary:")
    print(filtering_summary.to_string(index=False))

    # --------------------------------------------------------
    # 6.13 HVG selection from raw counts
    # --------------------------------------------------------

    print("\nCalculating highly variable genes from raw counts...")

    rna_processed_full = (
        rna_filtered.copy()
    )

    sc.pp.highly_variable_genes(
        rna_processed_full,
        n_top_genes=int(
            hvg_cfg["n_top_genes"]
        ),
        flavor=str(
            hvg_cfg["flavor"]
        ),
        subset=False,
        inplace=True,
    )

    if "highly_variable" not in rna_processed_full.var:
        raise RuntimeError(
            "HVG calculation did not create var['highly_variable']."
        )

    n_hvgs = int(
        rna_processed_full.var[
            "highly_variable"
        ].sum()
    )

    # --------------------------------------------------------
    # 6.14 Normalize and log-transform full filtered gene matrix
    # --------------------------------------------------------

    print("\nNormalizing RNA and applying log1p...")

    sc.pp.normalize_total(
        rna_processed_full,
        target_sum=float(
            norm_cfg["target_sum"]
        ),
    )

    sc.pp.log1p(
        rna_processed_full
    )

    # Store the complete log-normalized gene matrix for Step 03A.
    # Keep the complete log-normalized RNA object before HVG subsetting.
    rna_log_full = rna_processed_full.copy()

# PCA/clustering object contains only highly variable genes.
    rna_hvg = (
        rna_log_full[
            :,
            rna_log_full.var["highly_variable"].to_numpy(),
            ]
                .copy()
)

# Store the complete log-normalized gene matrix in .raw.
# The .raw setter requires an AnnData object, not an anndata.Raw object.
    rna_hvg.raw = rna_log_full
    print("\nRNA HVG/raw validation:")
    print(f"RNA HVG shape: {rna_hvg.shape}")
    print(f"RNA raw shape: {rna_hvg.raw.shape}")
    print(f"Same cells: {rna_hvg.n_obs == rna_hvg.raw.n_obs}")
    print(
    "Raw has more or equal genes than HVG matrix:",
    rna_hvg.raw.n_vars >= rna_hvg.n_vars,
)

    if rna_hvg.n_obs != rna_hvg.raw.n_obs:
        raise RuntimeError(
        "rna_hvg and rna_hvg.raw do not contain the same cells."
    )

    if rna_hvg.raw.n_vars < rna_hvg.n_vars:
        raise RuntimeError(
        "rna_hvg.raw contains fewer genes than the HVG matrix."
    )

    hvg_summary = pd.DataFrame(
        [
            {
                "n_cells": int(rna_hvg.n_obs),
                "genes_after_gene_filtering": int(
                    rna_processed_full.n_vars
                ),
                "n_hvgs": n_hvgs,
                "hvg_flavor": str(
                    hvg_cfg["flavor"]
                ),
                "requested_n_top_genes": int(
                    hvg_cfg["n_top_genes"]
                ),
                "target_sum": float(
                    norm_cfg["target_sum"]
                ),
                "log1p": bool(
                    norm_cfg["log1p"]
                ),
                "raw_contains_full_log_normalized_genes": True,
            }
        ]
    )

    rna_hvg.write_h5ad(
        output_paths["rna_hvg_h5ad"]
    )

    write_csv(
        hvg_summary,
        output_paths["rna_hvg_summary"],
        index=False,
    )

    print("\nRNA HVG object:")
    print(rna_hvg)

    print("\nrna_hvg.raw:")
    print(rna_hvg.raw)

    # --------------------------------------------------------
    # 6.15 Scale and PCA
    # --------------------------------------------------------

    print("\nScaling HVGs and running PCA...")

    rna_pca = rna_hvg.copy()

    sc.pp.scale(
        rna_pca,
        max_value=float(
            dim_cfg["scale_max_value"]
        ),
    )

    sc.tl.pca(
        rna_pca,
        n_comps=int(
            dim_cfg["n_comps"]
        ),
        svd_solver=str(
            dim_cfg["pca_solver"]
        ),
        random_state=int(
            dim_cfg["random_state"]
        ),
    )

    pca_variance_table = (
        build_pca_variance_table(
            rna_pca
        )
    )

    write_csv(
        pca_variance_table,
        output_paths["rna_pca_variance_table"],
        index=False,
    )

    sc.pl.pca_variance_ratio(
        rna_pca,
        n_pcs=int(
            dim_cfg["n_comps"]
        ),
        log=True,
        show=False,
    )

    save_current_figure(
        output_paths["rna_pca_elbow"],
        dpi=300,
        show=True,
    )

    # PCA colored by QC metrics.
    sc.pl.pca(
        rna_pca,
        color=[
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
            "doublet_score",
        ],
        show=False,
    )

    save_current_figure(
        output_paths["rna_pca_qc_plot"],
        dpi=300,
        show=True,
    )

    # --------------------------------------------------------
    # 6.16 Interactive clustering grid
    # --------------------------------------------------------

    print("\n" + "=" * 80)
    print("INTERACTIVE PCA / NEIGHBORS / LEIDEN GRID")
    print("=" * 80)

    print("\nInspect the PCA elbow plot.")
    print("You may provide one value or comma-separated values.")
    print("Example n_pcs: 20,30,40")
    print("Example n_neighbors: 10,15,20")
    print("Example Leiden resolution: 0.4,0.6,0.8,1.0")

    n_pcs_grid = parse_int_list(
        input(
            "\nn_pcs candidates "
            f"[default: {dim_cfg['n_pcs']}]: "
        ),
        default=int(dim_cfg["n_pcs"]),
    )

    n_neighbors_grid = parse_int_list(
        input(
            "n_neighbors candidates "
            f"[default: {dim_cfg['n_neighbors']}]: "
        ),
        default=int(dim_cfg["n_neighbors"]),
    )

    resolution_grid = parse_float_list(
        input(
            "Leiden resolution candidates "
            f"[default: {dim_cfg['leiden_resolution']}]: "
        ),
        default=float(
            dim_cfg["leiden_resolution"]
        ),
    )

    max_available_pcs = int(
        rna_pca.obsm["X_pca"].shape[1]
    )

    invalid_pcs = [
        value
        for value in n_pcs_grid
        if value < 2 or value > max_available_pcs
    ]

    if invalid_pcs:
        raise ValueError(
            f"Invalid n_pcs values: {invalid_pcs}\n"
            f"Available PCA dimensions: {max_available_pcs}"
        )

    print("\nCandidate grid:")
    print(f"n_pcs: {n_pcs_grid}")
    print(f"n_neighbors: {n_neighbors_grid}")
    print(f"Leiden resolutions: {resolution_grid}")

    grid_summary_rows = []
    grid_size_tables = []

    for (
        n_pcs,
        n_neighbors,
        resolution,
    ) in product(
        n_pcs_grid,
        n_neighbors_grid,
        resolution_grid,
    ):
        resolution_label = (
            safe_resolution_label(
                resolution
            )
        )

        diagnostic_key = (
            f"leiden_pcs{n_pcs}"
            f"_neighbors{n_neighbors}"
            f"_res{resolution_label}"
        )

        print(
            "\nRunning grid combination:"
            f" PCs={n_pcs},"
            f" neighbors={n_neighbors},"
            f" resolution={resolution}"
        )

        grid_adata = rna_pca.copy()

        sc.pp.neighbors(
            grid_adata,
            n_neighbors=int(n_neighbors),
            n_pcs=int(n_pcs),
            random_state=int(
                dim_cfg["random_state"]
            ),
        )

        sc.tl.umap(
            grid_adata,
            random_state=int(
                dim_cfg["random_state"]
            ),
        )

        sc.tl.leiden(
            grid_adata,
            resolution=float(resolution),
            key_added=diagnostic_key,
            random_state=int(
                dim_cfg["random_state"]
            ),
        )

        grid_summary_rows.append(
            summarize_grid_result(
                grid_adata,
                cluster_key=diagnostic_key,
                n_pcs=n_pcs,
                n_neighbors=n_neighbors,
                resolution=resolution,
            )
        )

        grid_size_tables.append(
            summarize_cluster_sizes(
                grid_adata,
                cluster_key=diagnostic_key,
                n_pcs=n_pcs,
                n_neighbors=n_neighbors,
                resolution=resolution,
            )
        )

        diagnostic_umap_out = (
            output_paths[
                "rna_umap_diagnostics_dir"
            ]
            / (
                f"umap_pcs{n_pcs}"
                f"_neighbors{n_neighbors}"
                f"_res{resolution_label}.png"
            )
        )

        sc.pl.umap(
            grid_adata,
            color=diagnostic_key,
            legend_loc="on data",
            legend_fontsize=6,
            legend_fontweight="normal",
            frameon=False,
            title=(
                f"PCs={n_pcs}, "
                f"neighbors={n_neighbors}, "
                f"Leiden={resolution}"
            ),
            show=False,
        )

        save_current_figure(
            diagnostic_umap_out,
            dpi=250,
            show=True,
        )

    grid_summary = pd.DataFrame(
        grid_summary_rows
    )

    grid_cluster_sizes = pd.concat(
        grid_size_tables,
        ignore_index=True,
    )

    write_csv(
        grid_summary,
        clustering_grid_summary_out,
        index=False,
    )

    write_csv(
        grid_cluster_sizes,
        clustering_grid_sizes_out,
        index=False,
    )

    print("\nClustering-grid summary:")
    print(
        grid_summary[
            [
                "n_pcs",
                "n_neighbors",
                "leiden_resolution",
                "n_clusters",
                "min_cluster_size",
                "clusters_lt_20",
                "clusters_lt_50",
                "clusters_lt_100",
            ]
        ].to_string(index=False)
    )

    # --------------------------------------------------------
    # 6.17 Final clustering selection
    # --------------------------------------------------------

    print("\n" + "=" * 80)
    print("FINAL CLUSTERING PARAMETER SELECTION")
    print("=" * 80)

    selected_n_pcs = int(
        input(
            "\nFinal selected n_pcs: "
        ).strip()
    )

    selected_n_neighbors = int(
        input(
            "Final selected n_neighbors: "
        ).strip()
    )

    selected_resolution = float(
        input(
            "Final selected Leiden resolution: "
        ).strip()
    )

    # --------------------------------------------------------
# Validate or generate the final selected combination
# --------------------------------------------------------

    combination_mask = (
        (grid_summary["n_pcs"].astype(int) == selected_n_pcs)
        & (
        grid_summary["n_neighbors"].astype(int)
        == selected_n_neighbors
    )
        & np.isclose(
        grid_summary["leiden_resolution"].astype(float),
        selected_resolution,
        rtol=0.0,
        atol=1e-8,
    )
)

    combination_exists = bool(
        combination_mask.any()
)

    if combination_exists:
        print(
        "\nThe final selected combination already exists "
        "in the diagnostic grid."
    )

    else:
        print(
        "\nWARNING: The final selected combination was not "
        "included in the original diagnostic grid."
        )

        print(
        "Running the selected combination now before "
        "locking it into preprocessing_rna.yaml:"
        )

        print(
        f"n_pcs={selected_n_pcs}, "
        f"n_neighbors={selected_n_neighbors}, "
        f"Leiden resolution={selected_resolution}"
        )

    selected_resolution_label = safe_resolution_label(
        selected_resolution
    )

    selected_diagnostic_key = (
        f"leiden_pcs{selected_n_pcs}"
        f"_neighbors{selected_n_neighbors}"
        f"_res{selected_resolution_label}"
    )

    selected_diagnostic_adata = rna_pca.copy()

    sc.pp.neighbors(
        selected_diagnostic_adata,
        n_neighbors=selected_n_neighbors,
        n_pcs=selected_n_pcs,
        random_state=int(
            dim_cfg["random_state"]
        ),
    )

    sc.tl.umap(
        selected_diagnostic_adata,
        random_state=int(
            dim_cfg["random_state"]
        ),
    )

    sc.tl.leiden(
        selected_diagnostic_adata,
        resolution=selected_resolution,
        key_added=selected_diagnostic_key,
        random_state=int(
            dim_cfg["random_state"]
        ),
    )

    # Add the missing combination to the summary table.
    selected_summary_row = summarize_grid_result(
        selected_diagnostic_adata,
        cluster_key=selected_diagnostic_key,
        n_pcs=selected_n_pcs,
        n_neighbors=selected_n_neighbors,
        resolution=selected_resolution,
    )

    grid_summary = pd.concat(
        [
            grid_summary,
            pd.DataFrame([selected_summary_row]),
        ],
        ignore_index=True,
    )

    # Add detailed cluster sizes.
    selected_cluster_sizes = summarize_cluster_sizes(
        selected_diagnostic_adata,
        cluster_key=selected_diagnostic_key,
        n_pcs=selected_n_pcs,
        n_neighbors=selected_n_neighbors,
        resolution=selected_resolution,
    )

    grid_cluster_sizes = pd.concat(
        [
            grid_cluster_sizes,
            selected_cluster_sizes,
        ],
        ignore_index=True,
    )

    # Save the newly evaluated UMAP.
    selected_diagnostic_umap_out = (
        output_paths["rna_umap_diagnostics_dir"]
        / (
            f"umap_pcs{selected_n_pcs}"
            f"_neighbors{selected_n_neighbors}"
            f"_res{selected_resolution_label}.png"
        )
    )

    sc.pl.umap(
        selected_diagnostic_adata,
        color=selected_diagnostic_key,
        legend_loc="on data",
        legend_fontsize=6,
        legend_fontweight="normal",
        frameon=False,
        title=(
            f"PCs={selected_n_pcs}, "
            f"neighbors={selected_n_neighbors}, "
            f"Leiden={selected_resolution}"
        ),
        show=False,
    )

    save_current_figure(
        selected_diagnostic_umap_out,
        dpi=250,
        show=True,
    )

    # Rewrite diagnostics so they include the selected combination.
    write_csv(
        grid_summary,
        clustering_grid_summary_out,
        index=False,
    )

    write_csv(
        grid_cluster_sizes,
        clustering_grid_sizes_out,
        index=False,
    )

    print(
        "\nThe selected combination was evaluated and added "
        "to the diagnostic outputs."
    )

    print("\nSelected combination summary:")
    print(
        pd.DataFrame(
            [selected_summary_row]
        ).to_string(index=False)
    )

    update_yaml = input(
        "\nUpdate preprocessing_rna.yaml with final clustering "
        "parameters? Type yes/no: "
    )

    if not yes_response(update_yaml):
        raise RuntimeError(
            "Step 02B stopped because final clustering parameters "
            "were not written to preprocessing_rna.yaml."
        )

    yaml_backup, final_leiden_key = (
        backup_and_update_clustering_yaml(
            config=config,
            config_path=config_path,
            n_pcs=selected_n_pcs,
            n_neighbors=selected_n_neighbors,
            resolution=selected_resolution,
        )
    )

    selected_clustering = pd.DataFrame(
        [
            {
                "n_pcs": selected_n_pcs,
                "n_neighbors": selected_n_neighbors,
                "leiden_resolution": (
                    selected_resolution
                ),
                "leiden_key": final_leiden_key,
                "selection_source": (
                    "interactive Step 02B grid review"
                ),
                "yaml_updated": True,
                "yaml_path": str(config_path),
                "yaml_backup": str(yaml_backup),
            }
        ]
    )

    write_csv(
        selected_clustering,
        selected_clustering_out,
        index=False,
    )

    # --------------------------------------------------------
    # 6.18 Final selected neighbors / UMAP / Leiden
    # --------------------------------------------------------

    print("\nRunning final selected clustering...")

    rna_clustered = rna_pca.copy()

    sc.pp.neighbors(
        rna_clustered,
        n_neighbors=selected_n_neighbors,
        n_pcs=selected_n_pcs,
        random_state=int(
            dim_cfg["random_state"]
        ),
    )

    sc.tl.umap(
        rna_clustered,
        random_state=int(
            dim_cfg["random_state"]
        ),
    )

    sc.tl.leiden(
        rna_clustered,
        resolution=selected_resolution,
        key_added=final_leiden_key,
        random_state=int(
            dim_cfg["random_state"]
        ),
    )

    cluster_summary = summarize_cluster_sizes(
        rna_clustered,
        cluster_key=final_leiden_key,
        n_pcs=selected_n_pcs,
        n_neighbors=selected_n_neighbors,
        resolution=selected_resolution,
    )

    write_csv(
        cluster_summary,
        output_paths["rna_cluster_summary"],
        index=False,
    )

    write_csv(
        grid_summary,
        output_paths["rna_leiden_diagnostics"],
        index=False,
    )

    rna_clustered.write_h5ad(
        output_paths["rna_clustered_h5ad"]
    )

    # Final cluster UMAP.
    sc.pl.umap(
        rna_clustered,
        color=final_leiden_key,
        legend_loc="on data",
        legend_fontsize=6,
        legend_fontweight="normal",
        frameon=False,
        title=(
            f"Final RNA clusters: {final_leiden_key}"
        ),
        show=False,
    )

    save_current_figure(
        output_paths["rna_umap_clusters"],
        dpi=300,
        show=True,
    )

    # Final UMAP QC.
    sc.pl.umap(
        rna_clustered,
        color=[
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
            "doublet_score",
        ],
        frameon=False,
        show=False,
    )

    save_current_figure(
        output_paths["rna_umap_qc_plot"],
        dpi=300,
        show=True,
    )

    # --------------------------------------------------------
    # 6.19 Manifest
    # --------------------------------------------------------

    manifest = {
        "step": "02B",
        "script": script_path.name,
        "description": (
            "Paired synchronized RNA QC and doublet filtering, "
            "RNA normalization/HVG/PCA, interactive clustering "
            "diagnostics and final RNA Leiden clustering."
        ),
        "timestamp": datetime.now().isoformat(
            timespec="seconds"
        ),
        "project": paths.as_dict(),
        "inputs": {
            label: str(path)
            for label, path in input_paths.items()
        },
        "selected_rna_qc": {
            "mito_prefix": selected_mito_prefix,
            "min_genes_per_cell": min_genes,
            "min_cells_per_gene": min_cells_per_gene,
            "max_pct_mito": max_pct_mito,
            "max_total_counts": max_total_counts,
        },
        "selected_doublet_filter": {
            "score_threshold": (
                selected_doublet_threshold
            ),
            "filter_predicted_doublets": bool(
                selected_row[
                    "filter_predicted_doublets"
                ]
            ),
            "predicted_doublets": int(
                predicted_doublet.sum()
            ),
        },
        "selected_clustering": {
            "n_pcs": selected_n_pcs,
            "n_neighbors": selected_n_neighbors,
            "leiden_resolution": (
                selected_resolution
            ),
            "leiden_key": final_leiden_key,
        },
        "pairing": {
            "before_filtering_confirmed": bool(
                pairing_before[
                    "pairing_confirmed"
                ].iloc[0]
            ),
            "after_filtering_confirmed": bool(
                pairing_after[
                    "pairing_confirmed"
                ].iloc[0]
            ),
            "same_shared_cell_mask_applied": True,
        },
        "summary": {
            "cells_before": int(rna.n_obs),
            "cells_after": int(
                rna_clustered.n_obs
            ),
            "cells_removed": int(
                rna.n_obs - rna_clustered.n_obs
            ),
            "genes_raw": int(rna.n_vars),
            "genes_after_filtering": int(
                rna_processed_full.n_vars
            ),
            "n_hvgs": int(rna_hvg.n_vars),
            "n_clusters": int(
                cluster_summary.shape[0]
            ),
        },
        "outputs": {
            label: str(path)
            for label, path in all_output_paths.items()
        },
        "validation": {
            "step02A_yaml_match": bool(
                threshold_validation[
                    "matches"
                ].all()
            ),
            "scrublet_rerun": False,
            "atac_peak_filtering_performed": False,
            "rna_raw_contains_full_log_normalized_genes": True,
            "final_clustered_rna_has_umap": bool(
                "X_umap"
                in rna_clustered.obsm
            ),
            "final_leiden_key_present": bool(
                final_leiden_key
                in rna_clustered.obs.columns
            ),
        },
    }

    write_json(
        manifest,
        output_paths["manifest"],
    )

    # --------------------------------------------------------
    # 6.20 Final report
    # --------------------------------------------------------

    print("\n" + "=" * 80)
    print("STEP 02B COMPLETED SUCCESSFULLY")
    print("=" * 80)

    print(f"\nCells before: {rna.n_obs}")
    print(
        f"Cells after paired filtering: "
        f"{rna_clustered.n_obs}"
    )
    print(
        f"Predicted doublets removed: "
        f"{doublets_removed.shape[0]}"
    )
    print(
        f"Genes after gene filtering: "
        f"{rna_processed_full.n_vars}"
    )
    print(f"HVGs: {rna_hvg.n_vars}")
    print(
        f"Final Leiden key: "
        f"{final_leiden_key}"
    )
    print(
        f"Final clusters: "
        f"{cluster_summary.shape[0]}"
    )

    print("\nFinal clustered RNA object:")
    print(
        output_paths["rna_clustered_h5ad"]
    )

    print("\nSynchronized ATAC object:")
    print(
        output_paths[
            "paired_atac_after_rna_qc_h5ad"
        ]
    )

    print("\nSynchronized Multiome object:")
    print(
        output_paths[
            "paired_multiome_after_rna_qc_h5ad"
        ]
    )

    print(
        "\nNext step: inspect all PCA/UMAP/grid outputs. "
        "Only after confirming the final clustering should "
        "Step 03A annotation be run."
    )


if __name__ == "__main__":
    main()