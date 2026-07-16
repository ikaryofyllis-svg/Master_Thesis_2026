# ============================================================
# 05B_1_paired_prepare_atac_peak_detection_features.py
#
# Purpose
# -------
# Validate and freeze the official ATAC GRN input and calculate
# diagnostic peak-detection statistics:
#
# - globally across all approved cells,
# - within each cell_type_broad,
# - at thresholds 1%, 2%, 5%, 10%, 15%, and 20%.
#
# This script:
# 1. Loads the final post-ATAC-QC ATAC object.
# 2. Loads the official full-gene RNA GRN object.
# 3. Validates exact RNA-ATAC cell identity and order.
# 4. Validates required ATAC QC metadata.
# 5. Validates the ATAC matrix.
# 6. Parses genomic peak coordinates.
# 7. Freezes a byte-identical official ATAC GRN input.
# 8. Calculates global peak-detection metrics.
# 9. Calculates peak-detection metrics per cell_type_broad.
# 10. Produces diagnostic threshold summaries and plots.
#
# This script does NOT:
# - exclude peaks,
# - apply a final 10% peak threshold,
# - create cell-type-specific retained peak sets,
# - create BED files for motif scanning,
# - scan TF motifs,
# - create TF promoters,
# - construct TF-TF adjacency matrices,
# - calculate regulatory potential,
# - select final TF nodes.
#
# Final peak selection is deferred to 05B.2.
# ============================================================


# ============================================================
# 1. Imports
# ============================================================

from __future__ import annotations

import hashlib
import json
import re
import shutil
import sys
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
import yaml


# ============================================================
# 2. Resolve Refactored repository
# ============================================================

script_path = Path(__file__).resolve()

# Expected:
#
# Refactored/
# └── scripts/
#     └── paired_scripts/
#         └── 05B_1_paired_prepare_atac_peak_detection_features.py

refactor_dir = script_path.parents[2]

print("\nScript path:")
print(script_path)

print("\nRefactored repository:")
print(refactor_dir)

if not refactor_dir.is_dir():
    raise FileNotFoundError(
        "Refactored repository was not found:\n"
        f"{refactor_dir}"
    )


# ============================================================
# 3. Enable local imports
# ============================================================

src_dir = refactor_dir / "src"

if not src_dir.is_dir():
    raise FileNotFoundError(
        "The src directory was not found:\n"
        f"{src_dir}"
    )

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

from pbmcgrn.paths import ProjectPaths


# ============================================================
# 4. Validate active project
# ============================================================

paths = ProjectPaths(
    refactor_root=refactor_dir,
)

expected_project_id = "pbmc_10k_multiome_v1"

print("\nActive project ID:")
print(paths.project_id)

print("\nActive project root:")
print(paths.project_root)

if paths.project_id != expected_project_id:
    raise RuntimeError(
        "Wrong active project.\n"
        f"Expected: {expected_project_id}\n"
        f"Observed: {paths.project_id}\n"
        f"Project root: {paths.project_root}"
    )

if not paths.project_root.is_dir():
    raise FileNotFoundError(
        "Active project root does not exist:\n"
        f"{paths.project_root}"
    )


# ============================================================
# 5. Configuration
# ============================================================

expected_n_cells = 7_566

cell_type_column = "cell_type_broad"

detection_thresholds = [
    0.01,
    0.02,
    0.05,
    0.10,
    0.15,
    0.20,
]

integer_tolerance = 1e-6

expected_genome_build = "GRCh38"

required_obs_columns = [
    "cell_type_broad",
    "atac_fragments",
    "atac_frip",
    "atac_tss_fragment_fraction",
    "atac_mitochondrial_fraction",
    "atac_duplicate_fraction",
    "atac_lowmapq_fraction",
    "matrix_total_counts",
    "matrix_n_accessible_peaks",
    "pass_final_atac_qc",
]

required_var_columns = [
    "matrix_n_cells_by_counts",
    "matrix_total_counts",
]

standard_chromosomes = {
    *(f"chr{i}" for i in range(1, 23)),
    "chrX",
    "chrY",
    "chrM",
}

peak_patterns = [
    re.compile(
        r"^(?P<chrom>[^:]+):"
        r"(?P<start>\d+)-"
        r"(?P<end>\d+)$"
    ),
    re.compile(
        r"^(?P<chrom>chr[^-]+)-"
        r"(?P<start>\d+)-"
        r"(?P<end>\d+)$"
    ),
    re.compile(
        r"^(?P<chrom>chr[^_]+)_"
        r"(?P<start>\d+)_"
        r"(?P<end>\d+)$"
    ),
]


# ============================================================
# 6. Input paths
# ============================================================

atac_input_path = (
    paths.project_root
    / "data"
    / "processed"
    / "final_after_atac_qc"
    / "pbmc_multiome_atac_final_after_atac_qc.h5ad"
)

rna_reference_path = (
    paths.project_root
    / "data"
    / "processed"
    / "grn"
    / "rna"
    / "pbmc_multiome_rna_full_genes_final_cells.h5ad"
)


# ============================================================
# 7. Output directories
# ============================================================

official_atac_dir = (
    paths.project_root
    / "data"
    / "processed"
    / "grn"
    / "atac"
)

table_output_dir = (
    paths.project_root
    / "reports"
    / "tables"
    / "grn"
    / "step05B_1"
)

figure_output_dir = (
    paths.project_root
    / "reports"
    / "figures"
    / "grn"
    / "step05B_1"
)

summary_output_dir = (
    paths.project_root
    / "reports"
    / "summaries"
    / "grn"
    / "step05B_1"
)

for directory in (
    official_atac_dir,
    table_output_dir,
    figure_output_dir,
    summary_output_dir,
):
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )


# ============================================================
# 8. Output paths
# ============================================================

official_atac_path = (
    official_atac_dir
    / "pbmc_multiome_atac_final_cells_grn_input.h5ad"
)

input_audit_path = (
    table_output_dir
    / "05B_1_atac_input_audit.tsv"
)

cell_alignment_path = (
    table_output_dir
    / "05B_1_atac_rna_cell_alignment.tsv"
)

celltype_overview_path = (
    table_output_dir
    / "05B_1_cell_type_broad_overview.tsv"
)

peak_coordinate_path = (
    table_output_dir
    / "05B_1_atac_peak_coordinates.tsv.gz"
)

global_metrics_path = (
    table_output_dir
    / "05B_1_atac_peak_global_detection_metrics.tsv.gz"
)

celltype_metrics_path = (
    table_output_dir
    / "05B_1_atac_peak_detection_metrics_by_cell_type_broad.tsv.gz"
)

global_threshold_summary_path = (
    table_output_dir
    / "05B_1_atac_peak_threshold_summary_global.tsv"
)

celltype_threshold_summary_path = (
    table_output_dir
    / "05B_1_atac_peak_threshold_summary_by_cell_type_broad.tsv"
)

threshold_wide_path = (
    table_output_dir
    / "05B_1_atac_peak_retention_wide_by_cell_type_broad.tsv"
)

retention_curve_path = (
    figure_output_dir
    / "05B_1_atac_peak_retention_curve_by_cell_type_broad.png"
)

retention_heatmap_path = (
    figure_output_dir
    / "05B_1_atac_peak_retention_heatmap_by_cell_type_broad.png"
)

global_distribution_path = (
    figure_output_dir
    / "05B_1_atac_peak_global_detection_distribution.png"
)

summary_yaml_path = (
    summary_output_dir
    / "05B_1_atac_peak_detection_summary.yaml"
)

run_metadata_path = (
    summary_output_dir
    / "05B_1_run_metadata.json"
)


# ============================================================
# 9. Generic helper functions
# ============================================================

def utc_timestamp() -> str:
    """Return current UTC timestamp."""
    return datetime.now(timezone.utc).isoformat()


def validate_file(
    path: Path,
    label: str,
) -> None:
    """Require a non-empty regular file."""
    if not path.exists():
        raise FileNotFoundError(
            f"{label} was not found:\n{path}"
        )

    if not path.is_file():
        raise ValueError(
            f"{label} is not a regular file:\n{path}"
        )

    if path.stat().st_size == 0:
        raise ValueError(
            f"{label} is empty:\n{path}"
        )


def sha256_file(
    path: Path,
    chunk_size: int = 1024 * 1024,
) -> str:
    """Calculate SHA-256 incrementally."""
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)

            if not chunk:
                break

            digest.update(chunk)

    return digest.hexdigest()


def close_backed_anndata(
    adata: ad.AnnData,
) -> None:
    """Close an AnnData backing file safely."""
    if adata.isbacked and adata.file is not None:
        adata.file.close()


def json_safe(value: Any) -> Any:
    """Convert values recursively into JSON/YAML-safe objects."""
    if value is None:
        return None

    if isinstance(value, Path):
        return str(value)

    if isinstance(value, np.generic):
        return value.item()

    if isinstance(value, np.ndarray):
        return [
            json_safe(item)
            for item in value.tolist()
        ]

    if isinstance(value, pd.Index):
        return [
            json_safe(item)
            for item in value.tolist()
        ]

    if isinstance(value, pd.Series):
        return [
            json_safe(item)
            for item in value.tolist()
        ]

    if isinstance(value, dict):
        return {
            str(key): json_safe(item)
            for key, item in value.items()
        }

    if isinstance(value, (list, tuple, set)):
        return [
            json_safe(item)
            for item in value
        ]

    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        missing = False

    if isinstance(missing, (bool, np.bool_)) and missing:
        return None

    return value


def get_package_version(
    package_name: str,
) -> str | None:
    """Return package version if available."""
    try:
        return version(package_name)
    except Exception:
        return None


def normalize_chromosome(
    value: str,
) -> str:
    """Normalize chromosome labels to chr-prefixed names."""
    value = str(value).strip()

    if value.startswith("chr"):
        return value

    return f"chr{value}"


def parse_peak_name(
    peak_name: str,
) -> tuple[str, int, int] | None:
    """Parse common genomic peak identifier formats."""
    peak_name = str(peak_name).strip()

    for pattern in peak_patterns:
        match = pattern.match(peak_name)

        if match is None:
            continue

        return (
            normalize_chromosome(
                match.group("chrom")
            ),
            int(match.group("start")),
            int(match.group("end")),
        )

    return None


def calculate_peak_metrics(
    matrix: Any,
) -> dict[str, np.ndarray]:
    """
    Calculate peak metrics across cells.

    Parameters
    ----------
    matrix
        Cells × peaks matrix.

    Returns
    -------
    dict
        Arrays with one value per peak.
    """
    n_cells = int(matrix.shape[0])

    if n_cells == 0:
        raise ValueError(
            "Cannot calculate peak metrics for zero cells."
        )

    if sp.issparse(matrix):
        matrix = matrix.tocsr()

        n_cells_accessible = np.asarray(
            matrix.getnnz(axis=0)
        ).ravel().astype(np.int64)

        total_counts = np.asarray(
            matrix.sum(axis=0)
        ).ravel().astype(np.float64)

    else:
        matrix = np.asarray(matrix)

        n_cells_accessible = np.count_nonzero(
            matrix > 0,
            axis=0,
        ).astype(np.int64)

        total_counts = np.sum(
            matrix,
            axis=0,
            dtype=np.float64,
        )

    fraction_cells_accessible = (
        n_cells_accessible.astype(np.float64)
        / n_cells
    )

    mean_counts_all_cells = (
        total_counts
        / n_cells
    )

    mean_counts_accessible_cells = np.divide(
        total_counts,
        n_cells_accessible,
        out=np.zeros_like(
            total_counts,
            dtype=np.float64,
        ),
        where=n_cells_accessible > 0,
    )

    return {
        "n_cells": np.full(
            matrix.shape[1],
            n_cells,
            dtype=np.int64,
        ),
        "n_cells_accessible": n_cells_accessible,
        "fraction_cells_accessible": (
            fraction_cells_accessible
        ),
        "total_counts": total_counts,
        "mean_counts_all_cells": (
            mean_counts_all_cells
        ),
        "mean_counts_accessible_cells": (
            mean_counts_accessible_cells
        ),
    }


def build_threshold_summary(
    detection_fractions: pd.Series,
    n_accessible_cells: pd.Series,
    thresholds: list[float],
    group_label: str,
    group_value: str,
    n_cells_in_group: int,
) -> pd.DataFrame:
    """
    Summarize peak presence and retention within one group.

    Peak present:
        accessible in at least one cell.

    Peak passing:
        accessible in at least threshold fraction of cells.
    """
    fractions = pd.to_numeric(
        detection_fractions,
        errors="raise",
    )

    accessible_counts = pd.to_numeric(
        n_accessible_cells,
        errors="raise",
    )

    n_total_peaks = int(
        len(fractions)
    )

    present_mask = (
        accessible_counts > 0
    )

    n_peaks_present = int(
        present_mask.sum()
    )

    records: list[dict[str, Any]] = []

    for threshold in thresholds:
        passing_mask = (
            fractions >= threshold
        )

        n_peaks_passing = int(
            passing_mask.sum()
        )

        minimum_required_cells = int(
            np.ceil(
                threshold
                * n_cells_in_group
            )
        )

        records.append(
            {
                group_label: group_value,
                "n_cells_in_group": (
                    n_cells_in_group
                ),
                "threshold": threshold,
                "threshold_percent": (
                    100.0 * threshold
                ),
                "minimum_required_accessible_cells": (
                    minimum_required_cells
                ),
                "n_total_peaks": n_total_peaks,
                "n_peaks_present": (
                    n_peaks_present
                ),
                "n_peaks_present_but_below_threshold": int(
                    (
                        present_mask
                        & ~passing_mask
                    ).sum()
                ),
                "n_peaks_passing_threshold": (
                    n_peaks_passing
                ),
                "pct_peaks_passing_of_total": (
                    100.0
                    * n_peaks_passing
                    / n_total_peaks
                    if n_total_peaks > 0
                    else np.nan
                ),
                "pct_peaks_passing_of_present": (
                    100.0
                    * n_peaks_passing
                    / n_peaks_present
                    if n_peaks_present > 0
                    else np.nan
                ),
            }
        )

    return pd.DataFrame(records)


# ============================================================
# 10. Validate inputs
# ============================================================

required_inputs = {
    "Final post-QC ATAC object": (
        atac_input_path
    ),
    "Official full-gene RNA object": (
        rna_reference_path
    ),
}

print("\nValidating required inputs:")

for label, path in required_inputs.items():
    validate_file(
        path=path,
        label=label,
    )

    print(f"- {label}: OK")
    print(f"  {path}")


# ============================================================
# 11. Load ATAC and RNA
# ============================================================

run_started_at = utc_timestamp()

atac = ad.read_h5ad(
    atac_input_path,
    backed="r",
)

rna = ad.read_h5ad(
    rna_reference_path,
    backed="r",
)

print("\nLoaded objects:")
print(
    f"- ATAC: {atac.n_obs:,} cells × "
    f"{atac.n_vars:,} peaks"
)
print(
    f"- RNA: {rna.n_obs:,} cells × "
    f"{rna.n_vars:,} genes"
)


# ============================================================
# 12. Validate dimensions and identifiers
# ============================================================

if atac.n_obs != expected_n_cells:
    raise ValueError(
        "Unexpected number of ATAC cells.\n"
        f"Expected: {expected_n_cells:,}\n"
        f"Observed: {atac.n_obs:,}"
    )

if rna.n_obs != expected_n_cells:
    raise ValueError(
        "Unexpected number of RNA cells.\n"
        f"Expected: {expected_n_cells:,}\n"
        f"Observed: {rna.n_obs:,}"
    )

if atac.n_vars == 0:
    raise ValueError(
        "ATAC object contains zero peaks."
    )

if not atac.obs_names.is_unique:
    raise ValueError(
        "ATAC cell identifiers are not unique."
    )

if not atac.var_names.is_unique:
    raise ValueError(
        "ATAC peak identifiers are not unique."
    )

if not rna.obs_names.is_unique:
    raise ValueError(
        "RNA cell identifiers are not unique."
    )


# ============================================================
# 13. Validate exact RNA-ATAC cell alignment
# ============================================================

atac_cells = pd.Index(
    atac.obs_names.astype(str),
    name="cell_id",
)

rna_cells = pd.Index(
    rna.obs_names.astype(str),
    name="cell_id",
)

same_cell_set = (
    set(atac_cells)
    == set(rna_cells)
)

same_cell_order = (
    atac_cells.equals(rna_cells)
)

if not same_cell_set:
    raise ValueError(
        "RNA and ATAC do not contain the same cell set."
    )

if not same_cell_order:
    mismatch_positions = np.flatnonzero(
        atac_cells.to_numpy()
        != rna_cells.to_numpy()
    )

    raise ValueError(
        "RNA and ATAC cell order differs.\n"
        f"Mismatched positions: "
        f"{len(mismatch_positions):,}\n"
        f"Examples: "
        f"{mismatch_positions[:20].tolist()}"
    )

cell_alignment_df = pd.DataFrame(
    {
        "cell_position": np.arange(
            atac.n_obs,
            dtype=np.int64,
        ),
        "atac_cell_id": atac_cells,
        "rna_cell_id": rna_cells,
        "identical": True,
    }
)


# ============================================================
# 14. Validate required metadata
# ============================================================

missing_obs_columns = [
    column
    for column in required_obs_columns
    if column not in atac.obs.columns
]

if missing_obs_columns:
    raise KeyError(
        "Required ATAC obs columns are missing:\n"
        f"{missing_obs_columns}"
    )

missing_var_columns = [
    column
    for column in required_var_columns
    if column not in atac.var.columns
]

if missing_var_columns:
    raise KeyError(
        "Required ATAC var columns are missing:\n"
        f"{missing_var_columns}"
    )

if cell_type_column not in rna.obs.columns:
    raise KeyError(
        f"RNA is missing {cell_type_column!r}."
    )

atac_cell_types = (
    atac.obs[cell_type_column]
    .astype("string")
)

rna_cell_types = (
    rna.obs[cell_type_column]
    .astype("string")
)

if atac_cell_types.isna().any():
    raise ValueError(
        "ATAC cell_type_broad contains missing values."
    )

if rna_cell_types.isna().any():
    raise ValueError(
        "RNA cell_type_broad contains missing values."
    )

same_cell_type_annotations = np.array_equal(
    atac_cell_types.to_numpy(
        dtype=str
    ),
    rna_cell_types.to_numpy(
        dtype=str
    ),
)

if not same_cell_type_annotations:
    raise ValueError(
        "RNA and ATAC cell_type_broad annotations differ."
    )

pass_final_qc = (
    atac.obs["pass_final_atac_qc"]
    .astype(bool)
)

if not pass_final_qc.all():
    raise ValueError(
        "The final ATAC object contains cells with "
        "pass_final_atac_qc == False."
    )


# ============================================================
# 15. Create cell-type overview
# ============================================================

cell_type_values = sorted(
    atac_cell_types.unique().tolist()
)

celltype_overview_df = (
    atac_cell_types
    .value_counts(
        dropna=False
    )
    .rename_axis(
        cell_type_column
    )
    .reset_index(
        name="n_cells"
    )
)

celltype_overview_df[
    "fraction_of_all_cells"
] = (
    celltype_overview_df["n_cells"]
    / atac.n_obs
)

celltype_overview_df[
    "percent_of_all_cells"
] = (
    100.0
    * celltype_overview_df[
        "fraction_of_all_cells"
    ]
)

print("\nCell-type overview:")
print(
    celltype_overview_df.to_string(
        index=False
    )
)


# ============================================================
# 16. Determine ATAC counts source
# ============================================================

if "counts" in atac.layers:
    counts_source = "layers/counts"
else:
    counts_source = "X"

if counts_source == "layers/counts":
    matrix_object = atac.layers["counts"]
else:
    matrix_object = atac.X

print("\nATAC matrix:")
print(f"- Source: {counts_source}")
print(f"- Python type: {type(matrix_object)}")
print(f"- Shape: {matrix_object.shape}")
print(f"- Dtype: {matrix_object.dtype}")


# ============================================================
# 17. Load full ATAC matrix
# ============================================================

# The post-QC PBMC ATAC object should remain sparse.
# It is loaded into memory because peak detection must be
# calculated repeatedly across cell-type row subsets.

if counts_source == "layers/counts":
    atac_matrix = atac.layers["counts"][:]
else:
    atac_matrix = atac.X[:]

if sp.issparse(atac_matrix):
    atac_matrix = atac_matrix.tocsr()
else:
    atac_matrix = np.asarray(
        atac_matrix
    )

print("\nLoaded ATAC matrix:")
print(
    f"- Shape: {atac_matrix.shape[0]:,} cells × "
    f"{atac_matrix.shape[1]:,} peaks"
)


# ============================================================
# 18. Validate ATAC matrix values
# ============================================================

if sp.issparse(atac_matrix):
    stored_values = np.asarray(
        atac_matrix.data
    )
    total_nonzero = int(
        atac_matrix.nnz
    )
else:
    stored_values = np.asarray(
        atac_matrix
    ).ravel()
    total_nonzero = int(
        np.count_nonzero(
            atac_matrix
        )
    )

if stored_values.size > 0:
    if not np.issubdtype(
        stored_values.dtype,
        np.number,
    ):
        raise TypeError(
            "ATAC matrix is not numeric."
        )

    if np.isnan(stored_values).any():
        raise ValueError(
            "ATAC matrix contains NaN values."
        )

    if np.isinf(stored_values).any():
        raise ValueError(
            "ATAC matrix contains infinite values."
        )

    if (stored_values < 0).any():
        raise ValueError(
            "ATAC matrix contains negative values."
        )

    integer_like = bool(
        np.all(
            np.abs(
                stored_values
                - np.round(stored_values)
            )
            <= integer_tolerance
        )
    )

    minimum_stored_value = float(
        stored_values.min()
    )

    maximum_stored_value = float(
        stored_values.max()
    )

else:
    integer_like = True
    minimum_stored_value = 0.0
    maximum_stored_value = 0.0

total_entries = int(
    atac.n_obs
    * atac.n_vars
)

matrix_density = (
    total_nonzero
    / total_entries
)

matrix_sparsity = (
    1.0
    - matrix_density
)

print("\nATAC matrix validation:")
print(f"- Integer-like: {integer_like}")
print(
    f"- Minimum stored value: "
    f"{minimum_stored_value}"
)
print(
    f"- Maximum stored value: "
    f"{maximum_stored_value}"
)
print(f"- Density: {matrix_density:.8%}")
print(f"- Sparsity: {matrix_sparsity:.8%}")

if not integer_like:
    print(
        "\nWARNING:"
        "\n- ATAC values are non-negative and finite,"
        "\n  but not integer-like."
        "\n- Detection based on >0 remains valid."
    )


# ============================================================
# 19. Verify stored matrix-derived metrics
# ============================================================

if sp.issparse(atac_matrix):
    recomputed_cell_total_counts = np.asarray(
        atac_matrix.sum(axis=1)
    ).ravel()

    recomputed_cell_n_peaks = np.asarray(
        atac_matrix.getnnz(axis=1)
    ).ravel()

    recomputed_peak_n_cells = np.asarray(
        atac_matrix.getnnz(axis=0)
    ).ravel()

    recomputed_peak_total_counts = np.asarray(
        atac_matrix.sum(axis=0)
    ).ravel()

else:
    recomputed_cell_total_counts = np.sum(
        atac_matrix,
        axis=1,
    )

    recomputed_cell_n_peaks = np.count_nonzero(
        atac_matrix > 0,
        axis=1,
    )

    recomputed_peak_n_cells = np.count_nonzero(
        atac_matrix > 0,
        axis=0,
    )

    recomputed_peak_total_counts = np.sum(
        atac_matrix,
        axis=0,
    )

cell_total_counts_match = np.allclose(
    recomputed_cell_total_counts,
    pd.to_numeric(
        atac.obs["matrix_total_counts"],
        errors="raise",
    ).to_numpy(),
)

cell_n_peaks_match = np.array_equal(
    recomputed_cell_n_peaks.astype(
        np.int64
    ),
    pd.to_numeric(
        atac.obs[
            "matrix_n_accessible_peaks"
        ],
        errors="raise",
    )
    .to_numpy(
        dtype=np.int64
    ),
)

peak_n_cells_match = np.array_equal(
    recomputed_peak_n_cells.astype(
        np.int64
    ),
    pd.to_numeric(
        atac.var[
            "matrix_n_cells_by_counts"
        ],
        errors="raise",
    )
    .to_numpy(
        dtype=np.int64
    ),
)

peak_total_counts_match = np.allclose(
    recomputed_peak_total_counts,
    pd.to_numeric(
        atac.var["matrix_total_counts"],
        errors="raise",
    ).to_numpy(),
)

# Cell-level metrics must remain valid after cell subsetting,
# because they describe each retained cell independently.
if not cell_total_counts_match:
    raise ValueError(
        "Stored cell-level matrix_total_counts do not match "
        "the current final ATAC matrix."
    )

if not cell_n_peaks_match:
    raise ValueError(
        "Stored cell-level matrix_n_accessible_peaks do not "
        "match the current final ATAC matrix."
    )


# Peak-level metrics were calculated before synchronized cell
# filtering in Step 04A. AnnData subsetting preserved those var
# columns without recomputing them, so mismatch is expected.
stale_peak_metrics_detected = not (
    peak_n_cells_match
    and peak_total_counts_match
)

if stale_peak_metrics_detected:
    print(
        "\nWARNING:"
        "\n- Stored peak-level ATAC metrics do not match the "
        "current post-QC cell matrix."
        "\n- This is expected because matrix_n_cells_by_counts "
        "and matrix_total_counts were calculated before final "
        "cell filtering in Step 04A."
        "\n- Recomputed peak metrics from the current 7,566-cell "
        "matrix will be authoritative for the GRN pipeline."
    )

print("\nStored metric verification:")
print(
    f"- Cell total counts: "
    f"{cell_total_counts_match}"
)
print(
    f"- Cell accessible peaks: "
    f"{cell_n_peaks_match}"
)
print(
    f"- Peak accessible cells: "
    f"{peak_n_cells_match}"
)
print(
    f"- Peak total counts: "
    f"{peak_total_counts_match}"
)


# ============================================================
# 20. Parse genomic peak coordinates
# ============================================================

peak_records: list[dict[str, Any]] = []

for peak_index, peak_id in enumerate(
    atac.var_names.astype(str)
):
    parsed = parse_peak_name(
        peak_id
    )

    if parsed is None:
        chrom = pd.NA
        start = np.nan
        end = np.nan
        parsed_successfully = False
    else:
        chrom, start, end = parsed
        parsed_successfully = True

    peak_records.append(
        {
            "peak_index": peak_index,
            "peak_id": peak_id,
            "chrom": chrom,
            "start": start,
            "end": end,
            "parsed_successfully": (
                parsed_successfully
            ),
        }
    )

peak_coordinate_df = pd.DataFrame(
    peak_records
)

peak_coordinate_df[
    "coordinate_valid"
] = (
    peak_coordinate_df[
        "parsed_successfully"
    ]
    & peak_coordinate_df["chrom"].notna()
    & peak_coordinate_df["start"].notna()
    & peak_coordinate_df["end"].notna()
    & (
        peak_coordinate_df["start"]
        >= 0
    )
    & (
        peak_coordinate_df["end"]
        > peak_coordinate_df["start"]
    )
)

peak_coordinate_df[
    "peak_width"
] = (
    peak_coordinate_df["end"]
    - peak_coordinate_df["start"]
)

peak_coordinate_df[
    "standard_chromosome"
] = peak_coordinate_df[
    "chrom"
].isin(
    standard_chromosomes
)

peak_coordinate_df[
    "canonical_peak_id"
] = np.where(
    peak_coordinate_df[
        "coordinate_valid"
    ],
    (
        peak_coordinate_df[
            "chrom"
        ].astype(str)
        + ":"
        + peak_coordinate_df[
            "start"
        ].astype("Int64").astype(str)
        + "-"
        + peak_coordinate_df[
            "end"
        ].astype("Int64").astype(str)
    ),
    pd.NA,
)

peak_coordinate_df[
    "duplicate_coordinate"
] = (
    peak_coordinate_df[
        "canonical_peak_id"
    ].notna()
    & peak_coordinate_df[
        "canonical_peak_id"
    ].duplicated(
        keep=False
    )
)

n_invalid_coordinates = int(
    (
        ~peak_coordinate_df[
            "coordinate_valid"
        ]
    ).sum()
)

if n_invalid_coordinates > 0:
    examples = (
        peak_coordinate_df.loc[
            ~peak_coordinate_df[
                "coordinate_valid"
            ]
        ]
        .head(20)
    )

    raise ValueError(
        "Some ATAC peak identifiers could not be parsed.\n"
        f"Invalid peaks: {n_invalid_coordinates:,}\n"
        f"{examples.to_string(index=False)}"
    )

print("\nPeak coordinates:")
print(f"- Total peaks: {atac.n_vars:,}")
print(
    "- Valid coordinates: "
    f"{peak_coordinate_df['coordinate_valid'].sum():,}"
)
print(
    "- Duplicate-coordinate features: "
    f"{peak_coordinate_df['duplicate_coordinate'].sum():,}"
)
print(
    "- Non-standard chromosome peaks: "
    f"{(~peak_coordinate_df['standard_chromosome']).sum():,}"
)


# ============================================================
# 21. Calculate global peak metrics
# ============================================================

global_arrays = calculate_peak_metrics(
    atac_matrix
)

global_metrics_df = (
    peak_coordinate_df[
        [
            "peak_index",
            "peak_id",
            "canonical_peak_id",
            "chrom",
            "start",
            "end",
            "peak_width",
            "standard_chromosome",
            "duplicate_coordinate",
        ]
    ]
    .copy()
)

global_metrics_df[
    "n_total_cells"
] = global_arrays["n_cells"]

global_metrics_df[
    "n_cells_accessible_global"
] = global_arrays[
    "n_cells_accessible"
]

global_metrics_df[
    "fraction_cells_accessible_global"
] = global_arrays[
    "fraction_cells_accessible"
]

global_metrics_df[
    "percent_cells_accessible_global"
] = (
    100.0
    * global_metrics_df[
        "fraction_cells_accessible_global"
    ]
)

global_metrics_df[
    "total_counts_global"
] = global_arrays[
    "total_counts"
]

global_metrics_df[
    "mean_counts_all_cells_global"
] = global_arrays[
    "mean_counts_all_cells"
]

global_metrics_df[
    "mean_counts_accessible_cells_global"
] = global_arrays[
    "mean_counts_accessible_cells"
]

global_metrics_df[
    "detected_in_any_cell_global"
] = (
    global_metrics_df[
        "n_cells_accessible_global"
    ]
    > 0
)

for threshold in detection_thresholds:
    label = int(
        round(
            threshold * 100
        )
    )

    global_metrics_df[
        f"pass_{label}pct_global"
    ] = (
        global_metrics_df[
            "fraction_cells_accessible_global"
        ]
        >= threshold
    )


# ============================================================
# 22. Calculate global threshold summary
# ============================================================

global_threshold_summary_df = (
    build_threshold_summary(
        detection_fractions=(
            global_metrics_df[
                "fraction_cells_accessible_global"
            ]
        ),
        n_accessible_cells=(
            global_metrics_df[
                "n_cells_accessible_global"
            ]
        ),
        thresholds=detection_thresholds,
        group_label="scope",
        group_value="all_cells_global",
        n_cells_in_group=atac.n_obs,
    )
)


# ============================================================
# 23. Calculate peak metrics per cell_type_broad
# ============================================================

cell_type_array = atac_cell_types.to_numpy(
    dtype=str
)

celltype_metric_frames: list[pd.DataFrame] = []
celltype_threshold_frames: list[pd.DataFrame] = []

for cell_type_value in cell_type_values:
    cell_mask = (
        cell_type_array
        == str(cell_type_value)
    )

    n_cells_in_group = int(
        cell_mask.sum()
    )

    if n_cells_in_group == 0:
        continue

    group_matrix = atac_matrix[
        cell_mask,
        :,
    ]

    group_arrays = calculate_peak_metrics(
        group_matrix
    )

    group_df = (
        peak_coordinate_df[
            [
                "peak_index",
                "peak_id",
                "canonical_peak_id",
                "chrom",
                "start",
                "end",
                "peak_width",
                "standard_chromosome",
                "duplicate_coordinate",
            ]
        ]
        .copy()
    )

    group_df[
        cell_type_column
    ] = str(cell_type_value)

    group_df[
        "n_cells_in_group"
    ] = group_arrays["n_cells"]

    group_df[
        "n_cells_accessible"
    ] = group_arrays[
        "n_cells_accessible"
    ]

    group_df[
        "fraction_cells_accessible"
    ] = group_arrays[
        "fraction_cells_accessible"
    ]

    group_df[
        "percent_cells_accessible"
    ] = (
        100.0
        * group_df[
            "fraction_cells_accessible"
        ]
    )

    group_df[
        "total_counts"
    ] = group_arrays[
        "total_counts"
    ]

    group_df[
        "mean_counts_all_cells"
    ] = group_arrays[
        "mean_counts_all_cells"
    ]

    group_df[
        "mean_counts_accessible_cells"
    ] = group_arrays[
        "mean_counts_accessible_cells"
    ]

    group_df[
        "detected_in_any_cell"
    ] = (
        group_df[
            "n_cells_accessible"
        ]
        > 0
    )

    for threshold in detection_thresholds:
        label = int(
            round(
                threshold * 100
            )
        )

        group_df[
            f"pass_{label}pct"
        ] = (
            group_df[
                "fraction_cells_accessible"
            ]
            >= threshold
        )

    threshold_df = build_threshold_summary(
        detection_fractions=(
            group_df[
                "fraction_cells_accessible"
            ]
        ),
        n_accessible_cells=(
            group_df[
                "n_cells_accessible"
            ]
        ),
        thresholds=detection_thresholds,
        group_label=cell_type_column,
        group_value=str(cell_type_value),
        n_cells_in_group=(
            n_cells_in_group
        ),
    )

    celltype_metric_frames.append(
        group_df
    )

    celltype_threshold_frames.append(
        threshold_df
    )

    print(
        f"\nCalculated peak metrics: "
        f"{cell_type_value}"
    )
    print(
        f"- Cells: {n_cells_in_group:,}"
    )
    print(
        "- Peaks detected in at least one cell: "
        f"{group_df['detected_in_any_cell'].sum():,}"
    )


celltype_metrics_df = pd.concat(
    celltype_metric_frames,
    axis=0,
    ignore_index=True,
)

celltype_threshold_summary_df = pd.concat(
    celltype_threshold_frames,
    axis=0,
    ignore_index=True,
)


# ============================================================
# 24. Create wide retention table
# ============================================================

threshold_wide_df = (
    celltype_threshold_summary_df
    .pivot(
        index=[
            cell_type_column,
            "n_cells_in_group",
            "n_total_peaks",
            "n_peaks_present",
        ],
        columns="threshold_percent",
        values="n_peaks_passing_threshold",
    )
    .reset_index()
)

threshold_wide_df.columns.name = None

threshold_wide_df = threshold_wide_df.rename(
    columns={
        1.0: "n_peaks_pass_1pct",
        2.0: "n_peaks_pass_2pct",
        5.0: "n_peaks_pass_5pct",
        10.0: "n_peaks_pass_10pct",
        15.0: "n_peaks_pass_15pct",
        20.0: "n_peaks_pass_20pct",
    }
)

for threshold in (
    1,
    2,
    5,
    10,
    15,
    20,
):
    count_column = (
        f"n_peaks_pass_{threshold}pct"
    )

    fraction_column = (
        f"pct_total_peaks_pass_{threshold}pct"
    )

    threshold_wide_df[
        fraction_column
    ] = (
        100.0
        * threshold_wide_df[
            count_column
        ]
        / threshold_wide_df[
            "n_total_peaks"
        ]
    )


# ============================================================
# 25. Retention curve
# ============================================================

sns.set_theme(
    style="whitegrid",
    context="notebook",
)

plt.figure(
    figsize=(11, 7)
)

sns.lineplot(
    data=celltype_threshold_summary_df,
    x="threshold_percent",
    y="pct_peaks_passing_of_total",
    hue=cell_type_column,
    marker="o",
)

plt.xlabel(
    "Peak detection threshold (%)"
)

plt.ylabel(
    "Peaks retained (% of total peak universe)"
)

plt.title(
    "Diagnostic ATAC peak retention by cell_type_broad"
)

plt.legend(
    title=cell_type_column,
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
)

plt.tight_layout()

plt.savefig(
    retention_curve_path,
    dpi=300,
    bbox_inches="tight",
)

plt.close()


# ============================================================
# 26. Retention heatmap
# ============================================================

retention_heatmap_df = (
    celltype_threshold_summary_df
    .pivot(
        index=cell_type_column,
        columns="threshold_percent",
        values="pct_peaks_passing_of_total",
    )
)

retention_heatmap_df = retention_heatmap_df[
    sorted(
        retention_heatmap_df.columns
    )
]

plt.figure(
    figsize=(9, 6)
)

sns.heatmap(
    retention_heatmap_df,
    cmap="viridis",
    annot=True,
    fmt=".1f",
    linewidths=0.4,
    cbar_kws={
        "label": "% of total peaks retained",
    },
)

plt.xlabel(
    "Detection threshold (%)"
)

plt.ylabel(
    cell_type_column
)

plt.title(
    "ATAC peak retention heatmap"
)

plt.tight_layout()

plt.savefig(
    retention_heatmap_path,
    dpi=300,
    bbox_inches="tight",
)

plt.close()


# ============================================================
# 27. Global detection distribution
# ============================================================

plt.figure(
    figsize=(8, 5)
)

sns.histplot(
    data=global_metrics_df,
    x="percent_cells_accessible_global",
    bins=100,
)

for threshold in detection_thresholds:
    plt.axvline(
        100.0 * threshold,
        linestyle="--",
        linewidth=1,
        alpha=0.7,
    )

plt.xlabel(
    "Cells accessible globally (%)"
)

plt.ylabel(
    "Number of peaks"
)

plt.title(
    "Global ATAC peak detection distribution"
)

plt.tight_layout()

plt.savefig(
    global_distribution_path,
    dpi=300,
    bbox_inches="tight",
)

plt.close()


# ============================================================
# 28. Create input audit table
# ============================================================

input_audit_df = pd.DataFrame(
    [
        {
            "audit_item": "atac_input_path",
            "observed_value": str(
                atac_input_path
            ),
            "status": "PASS",
        },
        {
            "audit_item": "rna_reference_path",
            "observed_value": str(
                rna_reference_path
            ),
            "status": "PASS",
        },
        {
            "audit_item": "atac_n_cells",
            "observed_value": atac.n_obs,
            "status": "PASS",
        },
        {
            "audit_item": "atac_n_peaks",
            "observed_value": atac.n_vars,
            "status": "PASS",
        },
        {
            "audit_item": "same_cell_set",
            "observed_value": same_cell_set,
            "status": "PASS",
        },
        {
            "audit_item": "same_cell_order",
            "observed_value": same_cell_order,
            "status": "PASS",
        },
        {
            "audit_item": (
                "same_cell_type_annotations"
            ),
            "observed_value": (
                same_cell_type_annotations
            ),
            "status": "PASS",
        },
        {
            "audit_item": (
                "all_cells_pass_final_atac_qc"
            ),
            "observed_value": (
                bool(pass_final_qc.all())
            ),
            "status": "PASS",
        },
        {
            "audit_item": "counts_source",
            "observed_value": counts_source,
            "status": "PASS",
        },
        {
            "audit_item": "matrix_integer_like",
            "observed_value": integer_like,
            "status": (
                "PASS"
                if integer_like
                else "WARNING"
            ),
        },
        {
            "audit_item": "matrix_sparsity",
            "observed_value": matrix_sparsity,
            "status": "PASS",
        },
        {
            "audit_item": (
                "stored_cell_total_counts_match"
            ),
            "observed_value": (
                cell_total_counts_match
            ),
            "status": "PASS",
        },
        {
            "audit_item": (
                "stored_cell_n_peaks_match"
            ),
            "observed_value": (
                cell_n_peaks_match
            ),
            "status": "PASS",
        },
        {
            "audit_item": (
                "stored_peak_n_cells_match"
            ),
            "observed_value": (
                peak_n_cells_match
            ),
            "status": "PASS",
        },
        {
            "audit_item": (
                "stored_peak_total_counts_match"
            ),
            "observed_value": (
                peak_total_counts_match
            ),
            "status": "PASS",
        },
        {
            "audit_item": (
                "invalid_peak_coordinates"
            ),
            "observed_value": (
                n_invalid_coordinates
            ),
            "status": "PASS",
        },
        {
            "audit_item": (
                "duplicate_coordinate_features"
            ),
            "observed_value": int(
                peak_coordinate_df[
                    "duplicate_coordinate"
                ].sum()
            ),
            "status": (
                "WARNING"
                if peak_coordinate_df[
                    "duplicate_coordinate"
                ].any()
                else "PASS"
            ),
        },
        {
            "audit_item": (
                "expected_genome_build"
            ),
            "observed_value": (
                expected_genome_build
            ),
            "status": "REFERENCE",
        },
    ]
)


# ============================================================
# 29. Write tabular outputs
# ============================================================

input_audit_df.to_csv(
    input_audit_path,
    sep="\t",
    index=False,
)

cell_alignment_df.to_csv(
    cell_alignment_path,
    sep="\t",
    index=False,
)

celltype_overview_df.to_csv(
    celltype_overview_path,
    sep="\t",
    index=False,
)

peak_coordinate_df.to_csv(
    peak_coordinate_path,
    sep="\t",
    index=False,
    compression="gzip",
)

global_metrics_df.to_csv(
    global_metrics_path,
    sep="\t",
    index=False,
    compression="gzip",
)

celltype_metrics_df.to_csv(
    celltype_metrics_path,
    sep="\t",
    index=False,
    compression="gzip",
)

global_threshold_summary_df.to_csv(
    global_threshold_summary_path,
    sep="\t",
    index=False,
)

celltype_threshold_summary_df.to_csv(
    celltype_threshold_summary_path,
    sep="\t",
    index=False,
)

threshold_wide_df.to_csv(
    threshold_wide_path,
    sep="\t",
    index=False,
)


# ============================================================
# 30. Close backed objects before freezing
# ============================================================

close_backed_anndata(
    atac
)

close_backed_anndata(
    rna
)


# ============================================================
# 31. Freeze official ATAC GRN input
# ============================================================

source_sha256 = sha256_file(
    atac_input_path
)

if official_atac_path.exists():
    existing_sha256 = sha256_file(
        official_atac_path
    )

    if existing_sha256 == source_sha256:
        copy_action = (
            "existing_identical_output_retained"
        )

    else:
        backup_timestamp = (
            datetime.now(timezone.utc)
            .strftime("%Y%m%dT%H%M%SZ")
        )

        backup_path = (
            official_atac_path.with_name(
                official_atac_path.stem
                + f".backup_{backup_timestamp}"
                + official_atac_path.suffix
            )
        )

        official_atac_path.replace(
            backup_path
        )

        shutil.copy2(
            atac_input_path,
            official_atac_path,
        )

        copy_action = (
            "different_existing_output_backed_up_and_replaced"
        )

else:
    shutil.copy2(
        atac_input_path,
        official_atac_path,
    )

    copy_action = (
        "new_official_copy_created"
    )

official_sha256 = sha256_file(
    official_atac_path
)

if official_sha256 != source_sha256:
    raise RuntimeError(
        "Official ATAC copy is not byte-identical "
        "to the validated input."
    )


# ============================================================
# 32. Summary YAML
# ============================================================

run_completed_at = utc_timestamp()

summary_data = {
    "step": "05B_1",
    "description": (
        "Validate and freeze ATAC GRN input and calculate "
        "cell-type-specific peak detection diagnostics"
    ),
    "status": "completed_successfully",

    "input": {
        "atac_path": atac_input_path,
        "rna_reference_path": (
            rna_reference_path
        ),
        "source_sha256": source_sha256,
    },

    "official_output": {
        "atac_path": official_atac_path,
        "sha256": official_sha256,
        "copy_action": copy_action,
    },

    "dimensions": {
        "n_cells": expected_n_cells,
        "n_peaks": len(
            peak_coordinate_df
        ),
        "n_cell_types": len(
            cell_type_values
        ),
    },

    "matrix": {
        "counts_source": counts_source,
        "integer_like": integer_like,
        "density": matrix_density,
        "sparsity": matrix_sparsity,
        "minimum_stored_value": (
            minimum_stored_value
        ),
        "maximum_stored_value": (
            maximum_stored_value
        ),
    },

    "thresholds": detection_thresholds,

    "important": {
        "peak_filtering_applied": False,
        "final_10pct_selection_applied": False,
        "motif_scanning_applied": False,
        "tf_tf_adjacency_created": False,
    },

    "next_step": (
        "05B.2 will retain peaks detected in at least "
        "10% of cells within each cell_type_broad and "
        "create cell-type-specific BED files and masks."
    ),

    "timestamps": {
        "started_at_utc": run_started_at,
        "completed_at_utc": (
            run_completed_at
        ),
    },
}

with summary_yaml_path.open(
    "w",
    encoding="utf-8",
) as handle:
    yaml.safe_dump(
        json_safe(summary_data),
        handle,
        sort_keys=False,
        allow_unicode=True,
    )


# ============================================================
# 33. Run metadata
# ============================================================

run_metadata = {
    "script_path": script_path,
    "script_name": script_path.name,
    "project_id": paths.project_id,
    "project_root": paths.project_root,
    "started_at_utc": run_started_at,
    "completed_at_utc": run_completed_at,

    "package_versions": {
        "anndata": get_package_version(
            "anndata"
        ),
        "numpy": get_package_version(
            "numpy"
        ),
        "pandas": get_package_version(
            "pandas"
        ),
        "scipy": get_package_version(
            "scipy"
        ),
        "matplotlib": get_package_version(
            "matplotlib"
        ),
        "seaborn": get_package_version(
            "seaborn"
        ),
        "pyyaml": get_package_version(
            "PyYAML"
        ),
    },

    "inputs": {
        "atac": {
            "path": atac_input_path,
            "sha256": source_sha256,
        },
        "rna_reference": {
            "path": rna_reference_path,
            "sha256": sha256_file(
                rna_reference_path
            ),
        },
    },

    "outputs": {
        "official_atac": official_atac_path,
        "input_audit": input_audit_path,
        "cell_alignment": (
            cell_alignment_path
        ),
        "celltype_overview": (
            celltype_overview_path
        ),
        "peak_coordinates": (
            peak_coordinate_path
        ),
        "global_metrics": (
            global_metrics_path
        ),
        "celltype_metrics": (
            celltype_metrics_path
        ),
        "global_threshold_summary": (
            global_threshold_summary_path
        ),
        "celltype_threshold_summary": (
            celltype_threshold_summary_path
        ),
        "threshold_wide": (
            threshold_wide_path
        ),
        "retention_curve": (
            retention_curve_path
        ),
        "retention_heatmap": (
            retention_heatmap_path
        ),
        "global_distribution": (
            global_distribution_path
        ),
        "summary_yaml": summary_yaml_path,
    },
}

with run_metadata_path.open(
    "w",
    encoding="utf-8",
) as handle:
    json.dump(
        json_safe(run_metadata),
        handle,
        indent=2,
        ensure_ascii=False,
    )


# ============================================================
# 34. Validate outputs
# ============================================================

required_outputs = {
    "Official ATAC GRN input": (
        official_atac_path
    ),
    "Input audit": input_audit_path,
    "Cell alignment": (
        cell_alignment_path
    ),
    "Cell-type overview": (
        celltype_overview_path
    ),
    "Peak coordinates": (
        peak_coordinate_path
    ),
    "Global peak metrics": (
        global_metrics_path
    ),
    "Cell-type peak metrics": (
        celltype_metrics_path
    ),
    "Global threshold summary": (
        global_threshold_summary_path
    ),
    "Cell-type threshold summary": (
        celltype_threshold_summary_path
    ),
    "Wide threshold table": (
        threshold_wide_path
    ),
    "Retention curve": (
        retention_curve_path
    ),
    "Retention heatmap": (
        retention_heatmap_path
    ),
    "Global detection distribution": (
        global_distribution_path
    ),
    "Summary YAML": summary_yaml_path,
    "Run metadata": run_metadata_path,
}

print("\nValidating outputs:")

for label, path in required_outputs.items():
    validate_file(
        path=path,
        label=label,
    )

    print(f"- {label}: OK")
    print(f"  {path}")


# ============================================================
# 35. Final report
# ============================================================

print("\n" + "=" * 92)
print("05B.1 COMPLETED SUCCESSFULLY")
print("=" * 92)

print("\nOfficial ATAC GRN input:")
print(official_atac_path)

print("\nInput dimensions:")
print(
    f"{expected_n_cells:,} cells × "
    f"{len(peak_coordinate_df):,} peaks"
)

print("\nRNA-ATAC synchronization:")
print(f"- Same cell set: {same_cell_set}")
print(f"- Same cell order: {same_cell_order}")
print(
    "- Same cell_type_broad annotations: "
    f"{same_cell_type_annotations}"
)

print("\nATAC matrix:")
print(f"- Counts source: {counts_source}")
print(f"- Integer-like: {integer_like}")
print(f"- Sparsity: {matrix_sparsity:.8%}")

print("\nDiagnostic thresholds:")
print(
    ", ".join(
        f"{int(value * 100)}%"
        for value in detection_thresholds
    )
)

print("\nPeak retention at 10%:")
ten_percent_summary = (
    celltype_threshold_summary_df.loc[
        celltype_threshold_summary_df[
            "threshold"
        ].eq(0.10),
        [
            cell_type_column,
            "n_cells_in_group",
            "minimum_required_accessible_cells",
            "n_total_peaks",
            "n_peaks_present",
            "n_peaks_passing_threshold",
            "pct_peaks_passing_of_total",
        ],
    ]
    .sort_values(
        cell_type_column
    )
)

print(
    ten_percent_summary.to_string(
        index=False
    )
)

print("\nImportant:")
print("- No peak was excluded.")
print("- The 10% values are diagnostic only.")
print("- No BED files were created.")
print("- No motif scanning was performed.")
print("- No TF-TF adjacency was created.")

print("\nNext stage:")
print(
    "05B.2 will apply the >=10% cell-type-specific "
    "peak threshold and create retained peak masks "
    "and BED files."
)