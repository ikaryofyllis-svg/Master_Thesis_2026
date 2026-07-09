"""
Run 02A — Interactive paired RNA QC and doublet diagnostics.

This diagnostic step:

1. Loads the raw paired RNA, ATAC and full Multiome objects.
2. Confirms that their cell barcodes are identical and ordered equally.
3. Calculates RNA QC metrics from raw RNA counts.
4. Runs Scrublet doublet detection on raw RNA counts.
5. Shows and saves RNA QC and doublet diagnostic plots.
6. Interactively explores candidate QC thresholds.
7. Allows selection of one final threshold combination.
8. Allows automatic or manual Scrublet threshold selection.
9. Creates a backup and updates preprocessing_rna.yaml.
10. Saves selected thresholds and a reproducibility manifest.

No cells or genes are removed in Step 02A.
Final synchronized filtering occurs in Step 02B.
"""

from __future__ import annotations

from datetime import datetime
from itertools import product
from pathlib import Path
from typing import Any
import importlib.util
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import yaml


# ============================================================
# 0. Locate Refactored root
# ============================================================

# Script location:
# Refactored/scripts/paired_scripts/
# run_02A_paired_rna_qc_diagnostics.py
#
# parents[0] = paired_scripts
# parents[1] = scripts
# parents[2] = Refactored

script_path = Path(__file__).resolve()
refactor_dir = script_path.parents[2]
src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))


# ============================================================
# 1. Import project helpers
# ============================================================

from pbmcgrn.config import load_yaml
from pbmcgrn.io import write_csv, write_json
from pbmcgrn.paths import ProjectPaths
from pbmcgrn.preprocessing.rna import calculate_rna_qc_metrics


# ============================================================
# 2. General helper functions
# ============================================================

def require_config_sections(
    config: dict[str, Any],
    required_sections: set[str],
) -> None:
    """Fail if required YAML sections are missing."""

    missing = required_sections.difference(config)

    if missing:
        raise KeyError(
            "Missing preprocessing_rna.yaml sections: "
            f"{sorted(missing)}"
        )


def require_keys(
    dictionary: dict[str, Any],
    required_keys: set[str],
    *,
    section_name: str,
) -> None:
    """Fail if required keys are missing from one YAML section."""

    missing = required_keys.difference(dictionary)

    if missing:
        raise KeyError(
            f"Missing keys in {section_name}: "
            f"{sorted(missing)}"
        )


def as_scalar_default(value: Any, parameter_name: str) -> Any:
    """
    Return a scalar default from YAML.

    A one-value list is accepted.
    Multi-value lists are ambiguous for the initial prompt.
    """

    if isinstance(value, list):
        if len(value) != 1:
            raise ValueError(
                f"{parameter_name} must initially be a scalar "
                f"or one-item list, but received: {value}"
            )

        return value[0]

    return value


def parse_number_list(
    user_text: str,
    default_value: Any,
    value_type: type,
) -> list[Any]:
    """Parse one number or comma-separated values."""

    if user_text.strip() == "":
        user_text = str(default_value)

    values = [
        item.strip()
        for item in user_text.split(",")
        if item.strip() != ""
    ]

    if not values:
        raise ValueError("No valid threshold values were supplied.")

    parsed = [value_type(item) for item in values]

    # Preserve order while removing duplicate values.
    return list(dict.fromkeys(parsed))


def ask_final_value(
    prompt_text: str,
    default_value: Any,
    value_type: type,
) -> Any:
    """Ask for one final scalar value."""

    user_text = input(
        f"{prompt_text} final value "
        f"[default: {default_value}]: "
    ).strip()

    if user_text == "":
        return value_type(default_value)

    return value_type(user_text)


def yes_response(user_text: str) -> bool:
    """Interpret common yes responses."""

    return user_text.strip().lower() in {
        "y",
        "yes",
        "ν",
        "ναι",
    }


def assert_raw_count_matrix(adata: sc.AnnData) -> None:
    """
    Basic check that RNA appears to contain non-negative raw counts.

    This is a sanity check, not mathematical proof that the matrix
    has never been normalized.
    """

    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("RNA AnnData is empty.")

    x = adata.X

    if sparse.issparse(x):
        values = x.data
    else:
        values = np.asarray(x).ravel()

    if values.size == 0:
        raise ValueError("RNA count matrix contains no non-zero values.")

    sample_size = min(values.size, 100_000)

    sample = np.asarray(
        values[:sample_size],
        dtype=float,
    )

    if np.any(~np.isfinite(sample)):
        raise ValueError(
            "RNA matrix contains non-finite values."
        )

    if np.any(sample < 0):
        raise ValueError(
            "RNA matrix contains negative values and does not "
            "look like raw counts."
        )

    integer_like_fraction = float(
        np.mean(np.isclose(sample, np.round(sample)))
    )

    print("\nRaw-count sanity check:")
    print(
        "Integer-like sampled values:",
        f"{integer_like_fraction:.4f}",
    )

    if integer_like_fraction < 0.99:
        raise ValueError(
            "RNA matrix does not look like raw integer counts. "
            "Scrublet and RNA QC should use unnormalized counts."
        )


def validate_paired_objects(
    rna: sc.AnnData,
    atac: sc.AnnData,
    multiome: sc.AnnData,
) -> pd.DataFrame:
    """Validate synchronized cell barcodes across three objects."""

    same_cell_count = (
        rna.n_obs
        == atac.n_obs
        == multiome.n_obs
    )

    same_barcode_order = (
        rna.obs_names.equals(atac.obs_names)
        and rna.obs_names.equals(multiome.obs_names)
    )

    same_barcode_set = (
        set(rna.obs_names)
        == set(atac.obs_names)
        == set(multiome.obs_names)
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

    no_duplicates = (
        rna_duplicates == 0
        and atac_duplicates == 0
        and multiome_duplicates == 0
    )

    pairing_confirmed = (
        same_cell_count
        and same_barcode_order
        and same_barcode_set
        and no_duplicates
    )

    summary = pd.DataFrame(
        [
            {
                "rna_n_cells": int(rna.n_obs),
                "atac_n_cells": int(atac.n_obs),
                "multiome_n_cells": int(
                    multiome.n_obs
                ),
                "same_cell_count": bool(
                    same_cell_count
                ),
                "same_barcode_set": bool(
                    same_barcode_set
                ),
                "same_barcode_order": bool(
                    same_barcode_order
                ),
                "rna_duplicate_barcodes": (
                    rna_duplicates
                ),
                "atac_duplicate_barcodes": (
                    atac_duplicates
                ),
                "multiome_duplicate_barcodes": (
                    multiome_duplicates
                ),
                "pairing_confirmed": bool(
                    pairing_confirmed
                ),
            }
        ]
    )

    if not pairing_confirmed:
        raise RuntimeError(
            "Paired input validation failed.\n"
            f"{summary.to_string(index=False)}"
        )

    return summary


# ============================================================
# 3. QC summary helpers
# ============================================================

def build_metric_summary(
    adata: sc.AnnData,
) -> pd.DataFrame:
    """Create descriptive statistics for main RNA QC metrics."""

    metrics = [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
    ]

    return (
        adata.obs[metrics]
        .describe()
        .T
        .reset_index(names="metric")
    )


def build_metric_quantiles(
    adata: sc.AnnData,
    quantiles: list[float],
) -> pd.DataFrame:
    """Create quantile table for main RNA QC metrics."""

    metrics = [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
    ]

    table = (
        adata.obs[metrics]
        .quantile(quantiles)
        .T
    )

    table.columns = [
        f"q{int(float(q) * 100):02d}"
        for q in quantiles
    ]

    return table.reset_index(names="metric")


def show_cell_metric_statistics(
    adata: sc.AnnData,
    metric_name: str,
) -> None:
    """Print descriptive statistics before threshold input."""

    series = adata.obs[metric_name]

    quantiles = series.quantile(
        [
            0.01,
            0.05,
            0.10,
            0.25,
            0.50,
            0.75,
            0.90,
            0.95,
            0.99,
        ]
    )

    rows = [
        {"statistic": "count", "value": series.count()},
        {"statistic": "mean", "value": series.mean()},
        {"statistic": "std", "value": series.std()},
        {"statistic": "min", "value": series.min()},
    ]

    for quantile, value in quantiles.items():
        rows.append(
            {
                "statistic": f"q{int(quantile * 100):02d}",
                "value": value,
            }
        )

    rows.append(
        {
            "statistic": "max",
            "value": series.max(),
        }
    )

    print("\n" + "=" * 80)
    print(f"STATISTICS FOR {metric_name}")
    print("=" * 80)
    print(pd.DataFrame(rows).to_string(index=False))


def detected_cells_per_gene(
    adata: sc.AnnData,
    cell_mask: np.ndarray | None = None,
) -> np.ndarray:
    """Count cells with non-zero expression for every gene."""

    matrix = adata.X

    if cell_mask is not None:
        matrix = matrix[cell_mask, :]

    if sparse.issparse(matrix):
        return np.asarray(
            matrix.getnnz(axis=0)
        ).ravel()

    return np.asarray(
        (matrix > 0).sum(axis=0)
    ).ravel()


def show_gene_detection_statistics(
    adata: sc.AnnData,
) -> None:
    """Print gene-detection statistics and common retention cutoffs."""

    detected = pd.Series(
        detected_cells_per_gene(adata),
        index=adata.var_names,
        name="detected_cells_per_gene",
    )

    quantiles = detected.quantile(
        [
            0.01,
            0.05,
            0.10,
            0.25,
            0.50,
            0.75,
            0.90,
            0.95,
            0.99,
        ]
    )

    summary_rows = [
        {
            "statistic": "count_genes",
            "value": detected.count(),
        },
        {
            "statistic": "mean_detected_cells",
            "value": detected.mean(),
        },
        {
            "statistic": "std",
            "value": detected.std(),
        },
        {
            "statistic": "min",
            "value": detected.min(),
        },
    ]

    for quantile, value in quantiles.items():
        summary_rows.append(
            {
                "statistic": (
                    f"q{int(quantile * 100):02d}"
                ),
                "value": value,
            }
        )

    summary_rows.append(
        {
            "statistic": "max",
            "value": detected.max(),
        }
    )

    common_thresholds = [
        1,
        2,
        3,
        5,
        10,
        20,
        50,
    ]

    retention = pd.DataFrame(
        {
            "min_cells_per_gene": common_thresholds,
            "genes_retained": [
                int((detected >= value).sum())
                for value in common_thresholds
            ],
        }
    )

    retention["genes_removed"] = (
        adata.n_vars
        - retention["genes_retained"]
    )

    retention["gene_retention_fraction"] = (
        retention["genes_retained"]
        / adata.n_vars
    )

    print("\n" + "=" * 80)
    print("STATISTICS FOR min_cells_per_gene")
    print("=" * 80)
    print(
        pd.DataFrame(summary_rows).to_string(
            index=False
        )
    )

    print("\nGene retention for common thresholds:")
    print(retention.to_string(index=False))


# ============================================================
# 4. Candidate threshold grid
# ============================================================

def build_threshold_candidate_summary(
    adata: sc.AnnData,
    min_genes_values: list[int],
    min_cells_gene_values: list[int],
    max_pct_mito_values: list[float],
    max_total_counts_values: list[float],
    selected_doublet_threshold: float,
) -> pd.DataFrame:
    """Evaluate all RNA QC threshold combinations."""

    rows: list[dict[str, Any]] = []

    obs = adata.obs

    predicted_doublet_mask = (
        obs["doublet_score"]
        >= selected_doublet_threshold
    )

    n_cells_before = int(adata.n_obs)
    n_genes_before = int(adata.n_vars)

    for (
        min_genes,
        min_cells_gene,
        max_pct_mito,
        max_total_counts,
    ) in product(
        min_genes_values,
        min_cells_gene_values,
        max_pct_mito_values,
        max_total_counts_values,
    ):
        min_genes_mask = (
            obs["n_genes_by_counts"]
            >= min_genes
        )

        mito_mask = (
            obs["pct_counts_mt"]
            < max_pct_mito
        )

        total_counts_mask = (
            obs["total_counts"]
            < max_total_counts
        )

        qc_cell_mask = (
            min_genes_mask
            & mito_mask
            & total_counts_mask
        )

        qc_and_singlet_mask = (
            qc_cell_mask
            & ~predicted_doublet_mask
        )

        gene_detection_counts = (
            detected_cells_per_gene(
                adata,
                qc_and_singlet_mask.to_numpy(),
            )
        )

        n_genes_after = int(
            (
                gene_detection_counts
                >= min_cells_gene
            ).sum()
        )

        n_cells_after_qc = int(
            qc_cell_mask.sum()
        )

        n_cells_after_qc_and_doublets = int(
            qc_and_singlet_mask.sum()
        )

        rows.append(
            {
                "min_genes_per_cell": int(
                    min_genes
                ),
                "min_cells_per_gene": int(
                    min_cells_gene
                ),
                "max_pct_mito": float(
                    max_pct_mito
                ),
                "max_total_counts": float(
                    max_total_counts
                ),
                "doublet_score_threshold": float(
                    selected_doublet_threshold
                ),
                "cells_before": n_cells_before,
                "cells_after_rna_qc": (
                    n_cells_after_qc
                ),
                "cells_after_qc_and_doublet_filter": (
                    n_cells_after_qc_and_doublets
                ),
                "cells_removed_by_rna_qc": (
                    n_cells_before
                    - n_cells_after_qc
                ),
                "predicted_doublets_total": int(
                    predicted_doublet_mask.sum()
                ),
                "predicted_doublets_among_qc_cells": int(
                    (
                        predicted_doublet_mask
                        & qc_cell_mask
                    ).sum()
                ),
                "final_cells_removed": (
                    n_cells_before
                    - n_cells_after_qc_and_doublets
                ),
                "final_cell_retention_fraction": (
                    n_cells_after_qc_and_doublets
                    / n_cells_before
                ),
                "genes_before": n_genes_before,
                "genes_after": n_genes_after,
                "genes_removed": (
                    n_genes_before
                    - n_genes_after
                ),
                "gene_retention_fraction": (
                    n_genes_after
                    / n_genes_before
                ),
                "removed_by_min_genes_only": int(
                    (~min_genes_mask).sum()
                ),
                "removed_by_mito_only": int(
                    (~mito_mask).sum()
                ),
                "removed_by_total_counts_only": int(
                    (~total_counts_mask).sum()
                ),
            }
        )

    result = pd.DataFrame(rows)

    result = result.sort_values(
        by=[
            "final_cell_retention_fraction",
            "gene_retention_fraction",
        ],
        ascending=False,
    ).reset_index(drop=True)

    result.insert(
        0,
        "candidate_id",
        range(1, result.shape[0] + 1),
    )

    return result


def choose_final_candidate(
    candidate_summary: pd.DataFrame,
) -> pd.Series | None:
    """Choose one candidate row or switch to manual mode."""

    review_columns = [
        "candidate_id",
        "min_genes_per_cell",
        "min_cells_per_gene",
        "max_pct_mito",
        "max_total_counts",
        "doublet_score_threshold",
        "cells_after_rna_qc",
        "cells_after_qc_and_doublet_filter",
        "final_cells_removed",
        "final_cell_retention_fraction",
        "genes_after",
        "genes_removed",
        "gene_retention_fraction",
    ]

    print("\n" + "=" * 80)
    print("FINAL THRESHOLD CANDIDATE SELECTION")
    print("=" * 80)

    print(
        candidate_summary[
            review_columns
        ].head(30).to_string(index=False)
    )

    selected_text = input(
        "\nFinal candidate_id or 'manual' "
        "[default: 1]: "
    ).strip()

    if selected_text == "":
        selected_text = "1"

    if selected_text.lower() == "manual":
        return None

    selected_id = int(selected_text)

    matches = candidate_summary.loc[
        candidate_summary["candidate_id"]
        == selected_id
    ]

    if matches.empty:
        raise ValueError(
            f"candidate_id not found: {selected_id}"
        )

    return matches.iloc[0]


# ============================================================
# 5. Plot functions
# ============================================================

def save_and_show_figure(
    output_path: Path,
) -> None:
    """Save current Matplotlib figure, display it, then close."""

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()
    plt.close()


def plot_rna_qc_diagnostics(
    adata: sc.AnnData,
    outputs: dict[str, Path],
    threshold_candidates: dict[str, list[Any]],
) -> None:
    """Save and display RNA QC plots."""

    # --------------------------------------------------------
    # Violin plots
    # --------------------------------------------------------

    sc.pl.violin(
        adata,
        keys=[
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
        ],
        jitter=0.25,
        multi_panel=True,
        show=False,
    )

    save_and_show_figure(
        outputs["qc_violin_plot"]
    )

    # --------------------------------------------------------
    # Histograms
    # --------------------------------------------------------

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(18, 5),
    )

    axes[0].hist(
        adata.obs["n_genes_by_counts"],
        bins=80,
        color="steelblue",
        alpha=0.85,
    )

    for value in threshold_candidates[
        "min_genes_per_cell"
    ]:
        axes[0].axvline(
            value,
            linestyle="--",
            linewidth=1.5,
            label=str(value),
        )

    axes[0].set_title("Detected genes per cell")
    axes[0].set_xlabel("n_genes_by_counts")
    axes[0].set_ylabel("Cells")
    axes[0].legend(title="min genes")

    axes[1].hist(
        adata.obs["pct_counts_mt"],
        bins=80,
        color="darkorange",
        alpha=0.85,
    )

    for value in threshold_candidates[
        "max_pct_mito"
    ]:
        axes[1].axvline(
            value,
            linestyle="--",
            linewidth=1.5,
            label=str(value),
        )

    axes[1].set_title("Mitochondrial percentage")
    axes[1].set_xlabel("pct_counts_mt")
    axes[1].set_ylabel("Cells")
    axes[1].legend(title="max mito")

    axes[2].hist(
        adata.obs["total_counts"],
        bins=80,
        color="seagreen",
        alpha=0.85,
    )

    for value in threshold_candidates[
        "max_total_counts"
    ]:
        axes[2].axvline(
            value,
            linestyle="--",
            linewidth=1.5,
            label=str(value),
        )

    axes[2].set_title("Total RNA counts per cell")
    axes[2].set_xlabel("total_counts")
    axes[2].set_ylabel("Cells")
    axes[2].legend(title="max counts")

    save_and_show_figure(
        outputs["qc_histogram_plot"]
    )

    # --------------------------------------------------------
    # Counts versus mitochondrial percentage
    # --------------------------------------------------------

    plt.figure(figsize=(8, 6))

    plt.scatter(
        adata.obs["total_counts"],
        adata.obs["pct_counts_mt"],
        s=6,
        alpha=0.35,
        color="slateblue",
    )

    for value in threshold_candidates[
        "max_pct_mito"
    ]:
        plt.axhline(
            value,
            linestyle="--",
            linewidth=1.2,
            label=f"mito={value}",
        )

    for value in threshold_candidates[
        "max_total_counts"
    ]:
        plt.axvline(
            value,
            linestyle=":",
            linewidth=1.2,
            label=f"counts={value}",
        )

    plt.xlabel("total_counts")
    plt.ylabel("pct_counts_mt")
    plt.title("RNA counts versus mitochondrial percentage")
    plt.legend(fontsize=8)

    save_and_show_figure(
        outputs["qc_counts_mito_scatter"]
    )

    # --------------------------------------------------------
    # Counts versus genes
    # --------------------------------------------------------

    plt.figure(figsize=(8, 6))

    plt.scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        s=6,
        alpha=0.35,
        color="teal",
    )

    for value in threshold_candidates[
        "min_genes_per_cell"
    ]:
        plt.axhline(
            value,
            linestyle="--",
            linewidth=1.2,
            label=f"genes={value}",
        )

    for value in threshold_candidates[
        "max_total_counts"
    ]:
        plt.axvline(
            value,
            linestyle=":",
            linewidth=1.2,
            label=f"counts={value}",
        )

    plt.xlabel("total_counts")
    plt.ylabel("n_genes_by_counts")
    plt.title("RNA counts versus detected genes")
    plt.legend(fontsize=8)

    save_and_show_figure(
        outputs["qc_counts_genes_scatter"]
    )


def plot_doublet_diagnostics(
    adata: sc.AnnData,
    outputs: dict[str, Path],
    threshold: float,
) -> None:
    """Save and display Scrublet diagnostic plots."""

    observed_scores = adata.obs[
        "doublet_score"
    ].astype(float)

    simulated_scores = np.asarray(
        adata.uns.get(
            "scrublet",
            {},
        ).get(
            "doublet_scores_sim",
            [],
        ),
        dtype=float,
    )

    # --------------------------------------------------------
    # Observed and simulated score distributions
    # --------------------------------------------------------

    plt.figure(figsize=(9, 6))

    plt.hist(
        observed_scores,
        bins=80,
        alpha=0.75,
        density=True,
        label="Observed cells",
        color="steelblue",
    )

    if simulated_scores.size > 0:
        plt.hist(
            simulated_scores,
            bins=80,
            alpha=0.45,
            density=True,
            label="Simulated doublets",
            color="darkorange",
        )

    plt.axvline(
        threshold,
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"threshold={threshold:.4f}",
    )

    plt.xlabel("Scrublet doublet score")
    plt.ylabel("Density")
    plt.title("Scrublet score distribution")
    plt.legend()

    save_and_show_figure(
        outputs["doublet_score_histogram"]
    )

    # --------------------------------------------------------
    # Doublet score versus RNA counts
    # --------------------------------------------------------

    plt.figure(figsize=(8, 6))

    plt.scatter(
        adata.obs["total_counts"],
        observed_scores,
        s=7,
        alpha=0.4,
        c=np.where(
            observed_scores >= threshold,
            "firebrick",
            "steelblue",
        ),
    )

    plt.axhline(
        threshold,
        color="black",
        linestyle="--",
        linewidth=1.5,
    )

    plt.xlabel("total_counts")
    plt.ylabel("doublet_score")
    plt.title("Scrublet score versus RNA counts")

    save_and_show_figure(
        outputs[
            "doublet_score_counts_scatter"
        ]
    )

    # --------------------------------------------------------
    # Doublet score versus detected genes
    # --------------------------------------------------------

    plt.figure(figsize=(8, 6))

    plt.scatter(
        adata.obs["n_genes_by_counts"],
        observed_scores,
        s=7,
        alpha=0.4,
        c=np.where(
            observed_scores >= threshold,
            "firebrick",
            "teal",
        ),
    )

    plt.axhline(
        threshold,
        color="black",
        linestyle="--",
        linewidth=1.5,
    )

    plt.xlabel("n_genes_by_counts")
    plt.ylabel("doublet_score")
    plt.title(
        "Scrublet score versus detected genes"
    )

    save_and_show_figure(
        outputs[
            "doublet_score_genes_scatter"
        ]
    )


# ============================================================
# 6. Scrublet functions
# ============================================================

def run_scrublet(
    adata: sc.AnnData,
    config: dict[str, Any],
) -> tuple[sc.AnnData, float, str]:
    """Run Scrublet on a copy of the raw RNA count matrix."""

    if importlib.util.find_spec("scrublet") is None:
        raise ImportError(
            "The 'scrublet' package is not installed.\n"
            "Install it once in the active environment with:\n"
            "python -m pip install scrublet"
        )

    result = adata.copy()

    configured_threshold = config.get(
        "score_threshold"
    )

    automatic_mode = (
        configured_threshold is None
    )

    print("\nRunning Scrublet on raw RNA counts...")

    sc.pp.scrublet(
        result,
        sim_doublet_ratio=float(
            config["sim_doublet_ratio"]
        ),
        expected_doublet_rate=float(
            config["expected_doublet_rate"]
        ),
        n_prin_comps=int(
            config["n_prin_comps"]
        ),
        threshold=(
            None
            if automatic_mode
            else float(configured_threshold)
        ),
        random_state=int(
            config["random_state"]
        ),
        copy=False,
        verbose=True,
    )

    if "doublet_score" not in result.obs:
        raise RuntimeError(
            "Scrublet did not create "
            "rna.obs['doublet_score']."
        )

    if "predicted_doublet" not in result.obs:
        raise RuntimeError(
            "Scrublet did not create "
            "rna.obs['predicted_doublet']."
        )

    scrublet_uns = result.uns.get(
        "scrublet",
        {}
    )

    detected_threshold = scrublet_uns.get(
        "threshold"
    )

    threshold_source = (
        "scanpy_scrublet_automatic"
        if automatic_mode
        else "preprocessing_rna_yaml"
    )

    if detected_threshold is None:
        predicted = result.obs[
            "predicted_doublet"
        ].astype(bool)

        if predicted.any():
            detected_threshold = float(
                result.obs.loc[
                    predicted,
                    "doublet_score",
                ].min()
            )

            threshold_source += (
                "_approximated_from_minimum_"
                "predicted_score"
            )

        elif configured_threshold is not None:
            detected_threshold = float(
                configured_threshold
            )

        else:
            detected_threshold = 0.25
            threshold_source = (
                "fallback_0.25_no_automatic_"
                "threshold_available"
            )

    detected_threshold = float(
        detected_threshold
    )

    # Recalculate explicitly so the stored labels correspond
    # exactly to the threshold recorded by this script.
    result.obs["predicted_doublet"] = (
        result.obs["doublet_score"]
        >= detected_threshold
    )

    return (
        result,
        detected_threshold,
        threshold_source,
    )


def choose_doublet_threshold(
    adata: sc.AnnData,
    automatic_threshold: float,
) -> tuple[float, str]:
    """Interactively accept or replace the Scrublet threshold."""

    scores = adata.obs[
        "doublet_score"
    ].astype(float)

    automatic_doublets = (
        scores >= automatic_threshold
    )

    print("\n" + "=" * 80)
    print("SCRUBLET THRESHOLD REVIEW")
    print("=" * 80)

    print(
        "Automatic/current threshold:",
        automatic_threshold,
    )

    print(
        "Predicted doublets:",
        int(automatic_doublets.sum()),
    )

    print(
        "Predicted doublet fraction:",
        float(automatic_doublets.mean()),
    )

    print("\nDoublet score quantiles:")

    print(
        scores.quantile(
            [
                0.01,
                0.05,
                0.10,
                0.25,
                0.50,
                0.75,
                0.90,
                0.95,
                0.99,
            ]
        )
    )

    user_text = input(
        "\nDoublet score threshold "
        f"[Enter keeps {automatic_threshold:.6f}]: "
    ).strip()

    if user_text == "":
        return (
            automatic_threshold,
            "accepted_automatic_scrublet_threshold",
        )

    return (
        float(user_text),
        "manual_interactive_doublet_threshold",
    )


# ============================================================
# 7. Main workflow
# ============================================================

def main() -> None:
    # --------------------------------------------------------
    # 7.1 Initialize project
    # --------------------------------------------------------

    paths = ProjectPaths(
        refactor_root=refactor_dir
    )

    paths.ensure_output_dirs()

    expected_project_id = (
        "pbmc_10k_multiome_v1"
    )

    print("\n" + "=" * 80)
    print(
        "RUN 02A — PAIRED INTERACTIVE RNA QC "
        "AND DOUBLET DIAGNOSTICS"
    )
    print("=" * 80)

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
    # 7.2 Load preprocessing YAML
    # --------------------------------------------------------

    config_path = (
        paths.configs
        / "preprocessing_rna.yaml"
    )

    if not config_path.exists():
        raise FileNotFoundError(
            "RNA preprocessing config not found:\n"
            f"{config_path}"
        )

    config = load_yaml(config_path)

    if config is None:
        raise ValueError(
            "preprocessing_rna.yaml was loaded as None.\n"
            "The YAML file is probably empty or was not saved.\n"
            f"Config path: {config_path}"
        )

    require_config_sections(
        config,
        {
            "rna_preprocessing",
            "doublet_detection",
            "qc_diagnostics",
        },
    )

    rna_cfg = config["rna_preprocessing"]

    doublet_cfg = config[
        "doublet_detection"
    ]

    qc_diag_cfg = config[
        "qc_diagnostics"
    ]

    diagnostic_outputs_cfg = (
        qc_diag_cfg["outputs"]
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
        section_name="rna_preprocessing",
    )

    require_keys(
        doublet_cfg,
        {
            "enabled",
            "expected_doublet_rate",
            "sim_doublet_ratio",
            "n_prin_comps",
            "random_state",
            "score_threshold",
            "filter_predicted_doublets",
        },
        section_name="doublet_detection",
    )

    print("\nLoaded config:")
    print(config_path)
    print("\nConfig path:")
    print(config_path)

    print("\nConfig exists:")
    print(config_path.exists())

    print("\nRaw YAML contents:")
    print(config_path.read_text(encoding="utf-8"))

    print("\nLoaded config object:")
    print(config)
    print(type(config))

    # --------------------------------------------------------
    # 7.3 Resolve outputs
    # --------------------------------------------------------

    outputs = {
        key: paths.prepare_output_path(
            value,
            label=key,
        )
        for key, value
        in diagnostic_outputs_cfg.items()
    }

    required_output_keys = {
        "raw_qc_summary",
        "raw_qc_quantiles",
        "threshold_candidate_summary",
        "selected_thresholds",
        "qc_violin_plot",
        "qc_histogram_plot",
        "qc_counts_mito_scatter",
        "qc_counts_genes_scatter",
        "doublet_scores_table",
        "doublet_summary",
        "doublet_score_quantiles",
        "doublet_score_histogram",
        "doublet_score_counts_scatter",
        "doublet_score_genes_scatter",
        "pairing_validation",
        "manifest",
    }

    missing_outputs = (
        required_output_keys
        - set(outputs)
    )

    if missing_outputs:
        raise KeyError(
            "Missing qc_diagnostics output keys: "
            f"{sorted(missing_outputs)}"
        )

    print("\nResolved outputs:")

    for key, path in outputs.items():
        print(f"{key}: {path}")

    # --------------------------------------------------------
    # 7.4 Load paired raw objects
    # --------------------------------------------------------

    rna_raw_path = (
        paths.processed_rna
        / "pbmc_multiome_rna_raw.h5ad"
    )

    atac_raw_path = (
        paths.processed_atac
        / "pbmc_multiome_atac_raw.h5ad"
    )

    multiome_raw_path = (
        paths.processed_multiome
        / "pbmc_multiome_raw.h5ad"
    )

    input_paths = {
        "rna_raw": rna_raw_path,
        "atac_raw": atac_raw_path,
        "multiome_raw": multiome_raw_path,
    }

    for label, path in input_paths.items():
        print(f"\n{label}:")
        print(path)

        if not path.exists():
            raise FileNotFoundError(
                f"{label} not found:\n{path}"
            )

    rna_raw = sc.read_h5ad(rna_raw_path)
    atac_raw = sc.read_h5ad(atac_raw_path)
    multiome_raw = sc.read_h5ad(
        multiome_raw_path
    )

    print("\nRaw RNA:")
    print(rna_raw)

    print("\nRaw ATAC:")
    print(atac_raw)

    print("\nRaw Multiome:")
    print(multiome_raw)

    # --------------------------------------------------------
    # 7.5 Validate pairing and raw counts
    # --------------------------------------------------------

    pairing_summary = (
        validate_paired_objects(
            rna=rna_raw,
            atac=atac_raw,
            multiome=multiome_raw,
        )
    )

    write_csv(
        pairing_summary,
        outputs["pairing_validation"],
        index=False,
    )

    print("\nPaired input validation:")
    print(
        pairing_summary.to_string(
            index=False
        )
    )

    assert_raw_count_matrix(rna_raw)

    # --------------------------------------------------------
    # 7.6 Calculate RNA QC metrics
    # --------------------------------------------------------

    default_mito_prefix = str(
        as_scalar_default(
            rna_cfg["mito_prefix"],
            "mito_prefix",
        )
    )

    mito_prefix_input = input(
        "\nmito_prefix "
        f"[default: {default_mito_prefix}]: "
    ).strip()

    mito_prefix = (
        mito_prefix_input
        if mito_prefix_input
        else default_mito_prefix
    )

    print("\nCalculating RNA QC metrics...")

    rna_qc = calculate_rna_qc_metrics(
        adata=rna_raw,
        mito_prefix=mito_prefix,
    )

    required_qc_columns = {
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
    }

    missing_qc_columns = (
        required_qc_columns
        - set(rna_qc.obs.columns)
    )

    if missing_qc_columns:
        raise RuntimeError(
            "RNA QC calculation did not create: "
            f"{sorted(missing_qc_columns)}"
        )

    if "mt" not in rna_qc.var:
        raise RuntimeError(
            "RNA QC calculation did not create "
            "rna.var['mt']."
        )

    n_mito_genes = int(
        rna_qc.var["mt"].sum()
    )

    print(
        "\nMitochondrial genes detected:",
        n_mito_genes,
    )

    if n_mito_genes == 0:
        raise RuntimeError(
            f"No mitochondrial genes matched "
            f"prefix '{mito_prefix}'."
        )

    qc_summary = build_metric_summary(
        rna_qc
    )

    qc_quantiles = build_metric_quantiles(
        rna_qc,
        [
            float(value)
            for value
            in qc_diag_cfg["quantiles"]
        ],
    )

    write_csv(
        qc_summary,
        outputs["raw_qc_summary"],
        index=False,
    )

    write_csv(
        qc_quantiles,
        outputs["raw_qc_quantiles"],
        index=False,
    )

    print("\nRaw RNA QC summary:")
    print(qc_summary.to_string(index=False))

    print("\nRaw RNA QC quantiles:")
    print(qc_quantiles.to_string(index=False))

    # --------------------------------------------------------
    # 7.7 Run Scrublet
    # --------------------------------------------------------

    if bool(doublet_cfg["enabled"]):
        (
            rna_qc,
            automatic_doublet_threshold,
            initial_threshold_source,
        ) = run_scrublet(
            rna_qc,
            doublet_cfg,
        )
    else:
        print(
            "\nScrublet is disabled in YAML."
        )

        rna_qc.obs["doublet_score"] = 0.0
        rna_qc.obs[
            "predicted_doublet"
        ] = False

        automatic_doublet_threshold = float(
            "inf"
        )

        initial_threshold_source = (
            "doublet_detection_disabled"
        )

    # Initial plot, so the user can inspect before choosing.
    plot_doublet_diagnostics(
        rna_qc,
        outputs,
        automatic_doublet_threshold,
    )

    if bool(doublet_cfg["enabled"]):
        (
            selected_doublet_threshold,
            selected_doublet_source,
        ) = choose_doublet_threshold(
            rna_qc,
            automatic_doublet_threshold,
        )
    else:
        selected_doublet_threshold = float(
            "inf"
        )

        selected_doublet_source = (
            "doublet_detection_disabled"
        )

    rna_qc.obs["predicted_doublet"] = (
        rna_qc.obs["doublet_score"]
        >= selected_doublet_threshold
    )

    # Recreate plots with the final selected threshold.
    plot_doublet_diagnostics(
        rna_qc,
        outputs,
        selected_doublet_threshold,
    )

    doublet_scores_table = (
        rna_qc.obs[
            [
                "total_counts",
                "n_genes_by_counts",
                "pct_counts_mt",
                "doublet_score",
                "predicted_doublet",
            ]
        ]
        .copy()
        .reset_index(names="cell_id")
    )

    write_csv(
        doublet_scores_table,
        outputs["doublet_scores_table"],
        index=False,
    )

    doublet_score_quantiles = (
        rna_qc.obs["doublet_score"]
        .quantile(
            [
                0.01,
                0.05,
                0.10,
                0.25,
                0.50,
                0.75,
                0.90,
                0.95,
                0.99,
            ]
        )
        .rename_axis("quantile")
        .reset_index(name="doublet_score")
    )

    write_csv(
        doublet_score_quantiles,
        outputs["doublet_score_quantiles"],
        index=False,
    )

    doublet_summary = pd.DataFrame(
        [
            {
                "n_cells": int(rna_qc.n_obs),
                "expected_doublet_rate": float(
                    doublet_cfg[
                        "expected_doublet_rate"
                    ]
                ),
                "automatic_threshold": float(
                    automatic_doublet_threshold
                ),
                "selected_threshold": float(
                    selected_doublet_threshold
                ),
                "initial_threshold_source": (
                    initial_threshold_source
                ),
                "selected_threshold_source": (
                    selected_doublet_source
                ),
                "predicted_doublets": int(
                    rna_qc.obs[
                        "predicted_doublet"
                    ].sum()
                ),
                "predicted_doublet_fraction": float(
                    rna_qc.obs[
                        "predicted_doublet"
                    ].mean()
                ),
                "filter_predicted_doublets": bool(
                    doublet_cfg[
                        "filter_predicted_doublets"
                    ]
                ),
            }
        ]
    )

    write_csv(
        doublet_summary,
        outputs["doublet_summary"],
        index=False,
    )

    print("\nFinal Scrublet summary:")
    print(
        doublet_summary.to_string(
            index=False
        )
    )

    # --------------------------------------------------------
    # 7.8 Ask for candidate QC thresholds
    # --------------------------------------------------------

    default_min_genes = int(
        as_scalar_default(
            rna_cfg["min_genes_per_cell"],
            "min_genes_per_cell",
        )
    )

    default_min_cells_gene = int(
        as_scalar_default(
            rna_cfg["min_cells_per_gene"],
            "min_cells_per_gene",
        )
    )

    default_max_pct_mito = float(
        as_scalar_default(
            rna_cfg["max_pct_mito"],
            "max_pct_mito",
        )
    )

    default_max_total_counts = float(
        as_scalar_default(
            rna_cfg["max_total_counts"],
            "max_total_counts",
        )
    )

    print("\n" + "=" * 80)
    print("INTERACTIVE RNA QC THRESHOLD EXPLORATION")
    print("=" * 80)

    print(
        "Give one value or comma-separated values, "
        "for example: 200,300,500"
    )

    show_cell_metric_statistics(
        rna_qc,
        "n_genes_by_counts",
    )

    min_genes_values = parse_number_list(
        input(
            "\nmin_genes_per_cell "
            f"[default: {default_min_genes}]: "
        ),
        default_min_genes,
        int,
    )

    show_gene_detection_statistics(
        rna_qc
    )

    min_cells_gene_values = parse_number_list(
        input(
            "\nmin_cells_per_gene "
            f"[default: {default_min_cells_gene}]: "
        ),
        default_min_cells_gene,
        int,
    )

    show_cell_metric_statistics(
        rna_qc,
        "pct_counts_mt",
    )

    max_pct_mito_values = parse_number_list(
        input(
            "\nmax_pct_mito "
            f"[default: {default_max_pct_mito}]: "
        ),
        default_max_pct_mito,
        float,
    )

    show_cell_metric_statistics(
        rna_qc,
        "total_counts",
    )

    max_total_counts_values = parse_number_list(
        input(
            "\nmax_total_counts "
            f"[default: {default_max_total_counts}]: "
        ),
        default_max_total_counts,
        float,
    )

    threshold_candidates = {
        "min_genes_per_cell": (
            min_genes_values
        ),
        "min_cells_per_gene": (
            min_cells_gene_values
        ),
        "max_pct_mito": (
            max_pct_mito_values
        ),
        "max_total_counts": (
            max_total_counts_values
        ),
    }

    print("\nCandidate values:")

    for key, value in threshold_candidates.items():
        print(f"{key}: {value}")

    # --------------------------------------------------------
    # 7.9 Show and save RNA QC plots
    # --------------------------------------------------------

    plot_rna_qc_diagnostics(
        rna_qc,
        outputs,
        threshold_candidates,
    )

    # --------------------------------------------------------
    # 7.10 Build candidate grid
    # --------------------------------------------------------

    candidate_summary = (
        build_threshold_candidate_summary(
            adata=rna_qc,
            min_genes_values=(
                min_genes_values
            ),
            min_cells_gene_values=(
                min_cells_gene_values
            ),
            max_pct_mito_values=(
                max_pct_mito_values
            ),
            max_total_counts_values=(
                max_total_counts_values
            ),
            selected_doublet_threshold=(
                selected_doublet_threshold
            ),
        )
    )

    write_csv(
        candidate_summary,
        outputs[
            "threshold_candidate_summary"
        ],
        index=False,
    )

    review_columns = [
        "candidate_id",
        "min_genes_per_cell",
        "min_cells_per_gene",
        "max_pct_mito",
        "max_total_counts",
        "cells_after_rna_qc",
        "cells_after_qc_and_doublet_filter",
        "final_cells_removed",
        "final_cell_retention_fraction",
        "genes_after",
        "genes_removed",
        "gene_retention_fraction",
    ]

    print("\nTop threshold combinations:")
    print(
        candidate_summary[
            review_columns
        ].head(30).to_string(index=False)
    )

    # --------------------------------------------------------
    # 7.11 Select final candidate
    # --------------------------------------------------------

    selected_candidate = (
        choose_final_candidate(
            candidate_summary
        )
    )

    final_mito_prefix_input = input(
        "\nmito_prefix final value "
        f"[default: {mito_prefix}]: "
    ).strip()

    final_mito_prefix = (
        final_mito_prefix_input
        if final_mito_prefix_input
        else mito_prefix
    )

    if selected_candidate is not None:
        final_candidate_id = int(
            selected_candidate["candidate_id"]
        )

        final_min_genes = int(
            selected_candidate[
                "min_genes_per_cell"
            ]
        )

        final_min_cells_gene = int(
            selected_candidate[
                "min_cells_per_gene"
            ]
        )

        final_max_pct_mito = float(
            selected_candidate[
                "max_pct_mito"
            ]
        )

        final_max_total_counts = float(
            selected_candidate[
                "max_total_counts"
            ]
        )

        print("\nSelected candidate:")
        print(
            selected_candidate[
                review_columns
            ].to_string()
        )

    else:
        final_candidate_id = None

        final_min_genes = ask_final_value(
            "min_genes_per_cell",
            default_min_genes,
            int,
        )

        final_min_cells_gene = (
            ask_final_value(
                "min_cells_per_gene",
                default_min_cells_gene,
                int,
            )
        )

        final_max_pct_mito = (
            ask_final_value(
                "max_pct_mito",
                default_max_pct_mito,
                float,
            )
        )

        final_max_total_counts = (
            ask_final_value(
                "max_total_counts",
                default_max_total_counts,
                float,
            )
        )

    # --------------------------------------------------------
    # 7.12 Final predicted doublet summary
    # --------------------------------------------------------

    final_doublet_mask = (
        rna_qc.obs["doublet_score"]
        >= selected_doublet_threshold
    )

    final_qc_mask = (
        (
            rna_qc.obs[
                "n_genes_by_counts"
            ]
            >= final_min_genes
        )
        & (
            rna_qc.obs[
                "pct_counts_mt"
            ]
            < final_max_pct_mito
        )
        & (
            rna_qc.obs[
                "total_counts"
            ]
            < final_max_total_counts
        )
    )

    if bool(
        doublet_cfg[
            "filter_predicted_doublets"
        ]
    ):
        final_keep_mask = (
            final_qc_mask
            & ~final_doublet_mask
        )
    else:
        final_keep_mask = final_qc_mask

    print("\n" + "=" * 80)
    print("FINAL SELECTION PREVIEW")
    print("=" * 80)

    print("Cells before:", rna_qc.n_obs)

    print(
        "Cells passing RNA QC:",
        int(final_qc_mask.sum()),
    )

    print(
        "Predicted doublets:",
        int(final_doublet_mask.sum()),
    )

    print(
        "Cells retained after configured filters:",
        int(final_keep_mask.sum()),
    )

    print(
        "Retention fraction:",
        float(final_keep_mask.mean()),
    )

    # --------------------------------------------------------
    # 7.13 Update YAML
    # --------------------------------------------------------

    update_text = input(
        "\nUpdate preprocessing_rna.yaml with "
        "these final values? [y/N]: "
    )

    yaml_updated = yes_response(update_text)

    backup_config_path: Path | None = None

    if yaml_updated:
        timestamp = datetime.now().strftime(
            "%Y%m%d_%H%M%S"
        )

        backup_config_path = (
            config_path.with_suffix(
                f".before_step02A_update_"
                f"{timestamp}.yaml"
            )
        )

        shutil.copy2(
            config_path,
            backup_config_path,
        )

        config[
            "rna_preprocessing"
        ]["mito_prefix"] = (
            final_mito_prefix
        )

        config[
            "rna_preprocessing"
        ]["min_genes_per_cell"] = (
            int(final_min_genes)
        )

        config[
            "rna_preprocessing"
        ]["min_cells_per_gene"] = (
            int(final_min_cells_gene)
        )

        config[
            "rna_preprocessing"
        ]["max_pct_mito"] = (
            float(final_max_pct_mito)
        )

        config[
            "rna_preprocessing"
        ]["max_total_counts"] = (
            float(final_max_total_counts)
        )

        config[
            "doublet_detection"
        ]["score_threshold"] = (
            float(selected_doublet_threshold)
        )

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

        print("\nUpdated YAML:")
        print(config_path)

        print("\nYAML backup:")
        print(backup_config_path)

    else:
        print(
            "\nYAML update skipped. "
            "Step 02B should not be run until "
            "the final values are stored in YAML."
        )

    # --------------------------------------------------------
    # 7.14 Save selected threshold handoff table
    # --------------------------------------------------------

    selected_thresholds = pd.DataFrame(
        [
            {
                "candidate_id": (
                    final_candidate_id
                ),
                "mito_prefix": (
                    final_mito_prefix
                ),
                "min_genes_per_cell": int(
                    final_min_genes
                ),
                "min_cells_per_gene": int(
                    final_min_cells_gene
                ),
                "max_pct_mito": float(
                    final_max_pct_mito
                ),
                "max_total_counts": float(
                    final_max_total_counts
                ),
                "doublet_detection_enabled": bool(
                    doublet_cfg["enabled"]
                ),
                "doublet_score_threshold": float(
                    selected_doublet_threshold
                ),
                "doublet_threshold_source": (
                    selected_doublet_source
                ),
                "filter_predicted_doublets": bool(
                    doublet_cfg[
                        "filter_predicted_doublets"
                    ]
                ),
                "predicted_doublets": int(
                    final_doublet_mask.sum()
                ),
                "cells_before": int(
                    rna_qc.n_obs
                ),
                "cells_passing_rna_qc": int(
                    final_qc_mask.sum()
                ),
                "cells_retained_final_preview": int(
                    final_keep_mask.sum()
                ),
                "final_retention_fraction": float(
                    final_keep_mask.mean()
                ),
                "yaml_updated": bool(
                    yaml_updated
                ),
                "yaml_path": str(
                    config_path
                ),
                "yaml_backup": (
                    str(backup_config_path)
                    if backup_config_path
                    else None
                ),
                "source_script": (
                    script_path.name
                ),
            }
        ]
    )

    write_csv(
        selected_thresholds,
        outputs["selected_thresholds"],
        index=False,
    )

    print("\nSelected thresholds:")
    print(
        selected_thresholds.to_string(
            index=False
        )
    )

    print("\nSaved selected thresholds:")
    print(outputs["selected_thresholds"])

    # --------------------------------------------------------
    # 7.15 Save manifest
    # --------------------------------------------------------

    manifest = {
        "step": "02A",
        "script": script_path.name,
        "description": (
            "Interactive paired RNA QC and Scrublet "
            "doublet diagnostics without filtering."
        ),
        "project": paths.as_dict(),
        "config": str(config_path),
        "inputs": {
            key: str(value)
            for key, value in input_paths.items()
        },
        "pairing": pairing_summary.iloc[
            0
        ].to_dict(),
        "rna_shape": {
            "n_cells": int(rna_qc.n_obs),
            "n_genes": int(rna_qc.n_vars),
        },
        "mitochondrial_genes_detected": (
            n_mito_genes
        ),
        "selected_qc_thresholds": {
            "mito_prefix": (
                final_mito_prefix
            ),
            "min_genes_per_cell": int(
                final_min_genes
            ),
            "min_cells_per_gene": int(
                final_min_cells_gene
            ),
            "max_pct_mito": float(
                final_max_pct_mito
            ),
            "max_total_counts": float(
                final_max_total_counts
            ),
        },
        "doublet_detection": {
            "enabled": bool(
                doublet_cfg["enabled"]
            ),
            "expected_doublet_rate": float(
                doublet_cfg[
                    "expected_doublet_rate"
                ]
            ),
            "automatic_threshold": float(
                automatic_doublet_threshold
            ),
            "selected_threshold": float(
                selected_doublet_threshold
            ),
            "selected_threshold_source": (
                selected_doublet_source
            ),
            "predicted_doublets": int(
                final_doublet_mask.sum()
            ),
            "predicted_doublet_fraction": float(
                final_doublet_mask.mean()
            ),
            "filter_predicted_doublets": bool(
                doublet_cfg[
                    "filter_predicted_doublets"
                ]
            ),
        },
        "final_preview": {
            "cells_before": int(
                rna_qc.n_obs
            ),
            "cells_passing_rna_qc": int(
                final_qc_mask.sum()
            ),
            "cells_retained_after_all_filters": int(
                final_keep_mask.sum()
            ),
            "retention_fraction": float(
                final_keep_mask.mean()
            ),
        },
        "yaml_updated": bool(yaml_updated),
        "yaml_backup": (
            str(backup_config_path)
            if backup_config_path
            else None
        ),
        "outputs": {
            key: str(value)
            for key, value in outputs.items()
        },
    }

    write_json(
        manifest,
        outputs["manifest"],
    )

    print("\nSaved manifest:")
    print(outputs["manifest"])

    print("\n" + "=" * 80)
    print("RUN 02A COMPLETED SUCCESSFULLY")
    print("=" * 80)

    print(
        "No RNA, ATAC or Multiome cells "
        "were removed in this diagnostic step."
    )

    print(
        "Final synchronized filtering will "
        "be performed in Step 02B."
    )


# ============================================================
# 8. Entry point
# ============================================================

if __name__ == "__main__":
    main()