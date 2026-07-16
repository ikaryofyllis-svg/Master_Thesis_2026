# ============================================================
# 05A_1_paired_prepare_tf_detection_features.py
#
# Purpose
# -------
# Prepare transcription-factor expression and detection
# statistics for downstream PBMC Multiome GRN analysis.
#
# This script:
# 1. Loads the official full-gene RNA object containing only
#    the final 7,566 approved cells.
# 2. Loads the frozen Lambert human TF reference.
# 3. Maps Lambert TF symbols to the full RNA gene universe.
# 4. Retains all uniquely mapped TFs as the master TF universe.
# 5. Calculates global TF expression/detection metrics.
# 6. Calculates TF metrics within each cell_type_broad.
# 7. Produces diagnostic summaries for thresholds:
#    1%, 2%, 5%, 10%, 15%, and 20%.
# 8. Produces TF-retention plots.
#
# This script does NOT:
# - exclude TFs using a detection threshold,
# - create final graph nodes,
# - use GENCODE, TSS, or promoters,
# - process JASPAR motifs,
# - map ATAC peaks,
# - construct GRN adjacency matrices,
# - train a GNN.
#
# The threshold summaries are diagnostic only.
# Final TF selection will be performed in a later script.
# ============================================================


# ============================================================
# 1. Imports
# ============================================================

from __future__ import annotations

import hashlib
import json
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
# 2. Resolve refactor repository
# ============================================================

script_path = Path(__file__).resolve()

# Expected location:
#
# Refactored/
# └── scripts/
#     └── paired_scripts/
#         └── 05A_1_paired_prepare_tf_detection_features.py
#
# parents[0] -> paired_scripts
# parents[1] -> scripts
# parents[2] -> Refactored

refactor_dir = script_path.parents[2]

print("\nScript path:")
print(script_path)

print("\nRefactor repository:")
print(refactor_dir)

if not refactor_dir.is_dir():
    raise FileNotFoundError(
        "Refactor repository was not found:\n"
        f"{refactor_dir}"
    )


# ============================================================
# 3. Enable local pbmcgrn imports
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
        f"Expected project ID: {expected_project_id}\n"
        f"Observed project ID: {paths.project_id}\n"
        f"Observed project root: {paths.project_root}"
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

# Diagnostic thresholds only.
# No TF is excluded by this script.
detection_thresholds = [
    0.01,
    0.02,
    0.05,
    0.10,
    0.15,
    0.20,
]

expected_lambert_version = "v1.01"

# Used only as a warning-level reference check.
expected_lambert_tf_count = 1_639

integer_tolerance = 1e-6
# ============================================================
# 5A. Marker-gene diagnostic configuration
# ============================================================

marker_gene_sets: dict[str, list[str]] = {
    "T_cells_general": [
        "CD3D",
        "CD3E",
        "CD3G",
        "TRAC",
    ],

    "CD4_T_cells": [
        "CD4",
        "IL7R",
        "CCR7",
        "LTB",
        "TCF7",
        "LEF1",
        "MAL",
        "NOSIP",
    ],

    "CD8_T_cells": [
        "CD8A",
        "CD8B",
        "CCL5",
        "GZMK",
        "LINC02446",
    ],

    "Cytotoxic_T_NK": [
        "NKG7",
        "GNLY",
        "PRF1",
        "GZMB",
        "CTSW",
        "CCL5",
    ],

    "NK_cells": [
        "NKG7",
        "GNLY",
        "PRF1",
        "KLRD1",
        "KLRF1",
        "TYROBP",
        "FCGR3A",
    ],

    "B_cells": [
        "MS4A1",
        "CD79A",
        "CD79B",
        "BANK1",
        "CD74",
        "HLA-DRA",
        "CD37",
    ],

    "Naive_B_cells": [
        "TCL1A",
        "IGHD",
        "IL4R",
        "FCER2",
        "IGHM",
    ],

    "Memory_B_cells": [
        "CD27",
        "AIM2",
        "TNFRSF13B",
        "GPR183",
        "CD82",
    ],

    "Plasma_cells": [
        "MZB1",
        "JCHAIN",
        "XBP1",
        "SDC1",
        "IGHG1",
        "IGHG3",
    ],

    "Monocytes_general": [
        "LYZ",
        "LST1",
        "AIF1",
        "TYROBP",
        "FCER1G",
        "CTSS",
    ],

    "Monocytes_CD14": [
        "S100A8",
        "S100A9",
        "FCN1",
        "CD14",
        "CTSD",
        "VCAN",
    ],

    "Monocytes_FCGR3A": [
        "FCGR3A",
        "MS4A7",
        "LST1",
        "AIF1",
        "IFITM3",
        "LILRB1",
    ],

    "Dendritic_cells_CD1C": [
        "FCER1A",
        "CST3",
        "CLEC10A",
        "CD1C",
        "CD74",
        "HLA-DPA1",
    ],

    "pDC": [
        "GZMB",
        "IRF7",
        "IL3RA",
        "TCF4",
        "CLEC4C",
    ],

    "Platelets": [
        "PPBP",
        "PF4",
        "NRGN",
        "SDPR",
        "RGS18",
        "GP9",
    ],

    "Cycling_cells": [
        "MKI67",
        "TOP2A",
        "STMN1",
        "TYMS",
        "UBE2C",
    ],

    "Interferon_response": [
        "ISG15",
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "IFI6",
        "MX1",
    ],

    "Stress_response": [
        "FOS",
        "JUN",
        "JUNB",
        "DUSP1",
        "HSPA1A",
        "HSPA1B",
    ],
}

# ============================================================
# 6. Resolve input paths
# ============================================================

rna_input_path = (
    paths.project_root
    / "data"
    / "processed"
    / "grn"
    / "rna"
    / "pbmc_multiome_rna_full_genes_final_cells.h5ad"
)

tf_reference_dir = (
    refactor_dir
    / "data"
    / "Human_TF_list"
)

tf_names_path = (
    tf_reference_dir
    / "TF_names_v_1.01.txt"
)

tf_database_path = (
    tf_reference_dir
    / "DatabaseExtract_v_1.01.csv"
)

tf_manifest_path = (
    tf_reference_dir
    / "lambert_human_tf_v1.01_manifest.yaml"
)


# ============================================================
# 7. Resolve output directories
# ============================================================

output_table_dir = (
    paths.project_root
    / "reports"
    / "tables"
    / "grn"
    / "step05A_1"
)

output_figure_dir = (
    paths.project_root
    / "reports"
    / "figures"
    / "grn"
    / "step05A_1"
)

output_summary_dir = (
    paths.project_root
    / "reports"
    / "summaries"
    / "grn"
    / "step05A_1"
)

for directory in (
    output_table_dir,
    output_figure_dir,
    output_summary_dir,
):
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )


# ============================================================
# 8. Resolve output paths
# ============================================================

input_audit_path = (
    output_table_dir
    / "05A_1_input_audit.tsv"
)

tf_catalogue_clean_path = (
    output_table_dir
    / "05A_1_human_tf_catalogue_clean.tsv"
)

tf_mapping_path = (
    output_table_dir
    / "05A_1_tf_to_full_rna_mapping.tsv"
)

tf_master_universe_path = (
    output_table_dir
    / "05A_1_master_mapped_tf_universe.tsv"
)

global_metrics_path = (
    output_table_dir
    / "05A_1_tf_global_detection_metrics.tsv"
)

celltype_metrics_path = (
    output_table_dir
    / "05A_1_tf_detection_metrics_by_cell_type_broad.tsv"
)

global_threshold_summary_path = (
    output_table_dir
    / "05A_1_tf_threshold_summary_global.tsv"
)

celltype_threshold_summary_path = (
    output_table_dir
    / "05A_1_tf_threshold_summary_by_cell_type_broad.tsv"
)

celltype_overview_path = (
    output_table_dir
    / "05A_1_cell_type_broad_overview.tsv"
)

mapping_exclusion_report_path = (
    output_table_dir
    / "05A_1_tf_mapping_exclusion_report.tsv"
)

retention_curve_path = (
    output_figure_dir
    / "05A_1_tf_retention_curve_by_cell_type_broad.png"
)

retention_heatmap_path = (
    output_figure_dir
    / "05A_1_tf_retention_heatmap_by_cell_type_broad.png"
)

global_distribution_path = (
    output_figure_dir
    / "05A_1_tf_global_detection_distribution.png"
)

summary_yaml_path = (
    output_summary_dir
    / "05A_1_tf_detection_summary.yaml"
)

run_metadata_path = (
    output_summary_dir
    / "05A_1_run_metadata.json"
)


# ============================================================
# 8A. Marker-gene diagnostic output paths
# ============================================================

marker_mapping_path = (
    output_table_dir
    / "05A_1_marker_gene_to_full_rna_mapping.tsv"
)

marker_global_metrics_path = (
    output_table_dir
    / "05A_1_marker_gene_global_detection_metrics.tsv"
)

marker_celltype_metrics_path = (
    output_table_dir
    / "05A_1_marker_gene_detection_by_cell_type_broad.tsv"
)

marker_presence_matrix_path = (
    output_table_dir
    / "05A_1_marker_gene_presence_matrix.tsv"
)

marker_detection_matrix_path = (
    output_table_dir
    / "05A_1_marker_gene_percent_detection_matrix.tsv"
)

marker_cluster_summary_path = (
    output_table_dir
    / "05A_1_marker_gene_cluster_summary.tsv"
)

marker_threshold_summary_path = (
    output_table_dir
    / "05A_1_marker_gene_threshold_summary.tsv"
)

marker_detection_heatmap_path = (
    output_figure_dir
    / "05A_1_marker_gene_detection_heatmap.png"
)

marker_presence_heatmap_path = (
    output_figure_dir
    / "05A_1_marker_gene_presence_heatmap.png"
)

# ============================================================
# 9. Generic helper functions
# ============================================================

def utc_timestamp() -> str:
    """Return the current UTC time in ISO-8601 format."""
    return datetime.now(timezone.utc).isoformat()


def validate_file(
    path: Path,
    label: str,
) -> None:
    """Validate that a required file exists and is non-empty."""
    if not path.exists():
        raise FileNotFoundError(
            f"{label} was not found:\n"
            f"{path}"
        )

    if not path.is_file():
        raise ValueError(
            f"{label} is not a regular file:\n"
            f"{path}"
        )

    if path.stat().st_size == 0:
        raise ValueError(
            f"{label} is empty:\n"
            f"{path}"
        )


def sha256_file(
    path: Path,
    chunk_size: int = 1024 * 1024,
) -> str:
    """Calculate SHA-256 without loading the entire file."""
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
    """Recursively convert values into JSON/YAML-safe values."""
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


def calculate_sparse_gene_metrics(
    matrix: Any,
) -> dict[str, np.ndarray]:
    """
    Calculate gene-level metrics across matrix rows.

    Parameters
    ----------
    matrix
        Cells × genes expression matrix.

    Returns
    -------
    dict
        Arrays with one value per gene.
    """
    n_cells = int(matrix.shape[0])

    if n_cells == 0:
        raise ValueError(
            "Cannot calculate metrics for a matrix with zero cells."
        )

    if sp.issparse(matrix):
        matrix = matrix.tocsr()

        n_expressing = np.asarray(
            matrix.getnnz(axis=0)
        ).ravel().astype(np.int64)

        total_counts = np.asarray(
            matrix.sum(axis=0)
        ).ravel().astype(np.float64)

    else:
        dense_matrix = np.asarray(matrix)

        n_expressing = np.count_nonzero(
            dense_matrix > 0,
            axis=0,
        ).astype(np.int64)

        total_counts = np.sum(
            dense_matrix,
            axis=0,
            dtype=np.float64,
        )

    fraction_expressing = (
        n_expressing.astype(np.float64)
        / n_cells
    )

    mean_counts_all_cells = (
        total_counts
        / n_cells
    )

    mean_counts_expressing_cells = np.divide(
        total_counts,
        n_expressing,
        out=np.zeros_like(
            total_counts,
            dtype=np.float64,
        ),
        where=n_expressing > 0,
    )

    return {
        "n_cells": np.full(
            matrix.shape[1],
            n_cells,
            dtype=np.int64,
        ),
        "n_expressing_cells": n_expressing,
        "fraction_expressing": fraction_expressing,
        "total_counts": total_counts,
        "mean_counts_all_cells": mean_counts_all_cells,
        "mean_counts_expressing_cells": (
            mean_counts_expressing_cells
        ),
    }


def build_threshold_summary(
    detection_fractions: pd.Series,
    n_expressing_cells: pd.Series,
    thresholds: list[float],
    group_label: str,
    group_value: str,
) -> pd.DataFrame:
    """
    Summarize TF presence and retention within one group.

    TF present in group:
        detected in at least one cell.

    TF retained at threshold:
        detected in at least the specified fraction of cells.
    """

    fractions = pd.to_numeric(
        detection_fractions,
        errors="raise",
    )

    expressing_counts = pd.to_numeric(
        n_expressing_cells,
        errors="raise",
    )

    n_master_mapped_tfs = int(
        len(fractions)
    )

    present_mask = (
        expressing_counts > 0
    )

    n_tfs_present_in_group = int(
        present_mask.sum()
    )

    pct_tfs_present_of_master = (
        100.0
        * n_tfs_present_in_group
        / n_master_mapped_tfs
        if n_master_mapped_tfs > 0
        else np.nan
    )

    records: list[dict[str, Any]] = []

    for threshold in thresholds:
        passing_mask = (
            fractions >= threshold
        )

        n_tfs_passing_threshold = int(
            passing_mask.sum()
        )

        n_tfs_present_but_below_threshold = int(
            (
                present_mask
                & ~passing_mask
            ).sum()
        )

        pct_tfs_passing_of_master = (
            100.0
            * n_tfs_passing_threshold
            / n_master_mapped_tfs
            if n_master_mapped_tfs > 0
            else np.nan
        )

        pct_tfs_passing_of_present = (
            100.0
            * n_tfs_passing_threshold
            / n_tfs_present_in_group
            if n_tfs_present_in_group > 0
            else np.nan
        )

        records.append(
            {
                group_label: group_value,
                "threshold": float(threshold),
                "threshold_percent": float(
                    threshold * 100
                ),

                "n_master_mapped_tfs": (
                    n_master_mapped_tfs
                ),

                "n_tfs_present_in_group": (
                    n_tfs_present_in_group
                ),

                "pct_tfs_present_of_master": (
                    pct_tfs_present_of_master
                ),

                "n_tfs_present_but_below_threshold": (
                    n_tfs_present_but_below_threshold
                ),

                "n_tfs_passing_threshold": (
                    n_tfs_passing_threshold
                ),

                "pct_tfs_passing_of_master": (
                    pct_tfs_passing_of_master
                ),

                "pct_tfs_passing_of_present": (
                    pct_tfs_passing_of_present
                ),
            }
        )

    return pd.DataFrame(
        records
    )

# ============================================================
# 10. Validate required input files
# ============================================================

required_inputs = {
    "Official full-gene final-cell RNA H5AD": (
        rna_input_path
    ),
    "Lambert TF-name list": tf_names_path,
    "Lambert complete TF database": tf_database_path,
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
# 11. Load and validate Lambert manifest if available
# ============================================================

manifest_data: dict[str, Any] | None = None

if tf_manifest_path.exists():
    validate_file(
        path=tf_manifest_path,
        label="Lambert reference manifest",
    )

    with tf_manifest_path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        manifest_data = yaml.safe_load(handle)

    observed_manifest_version = str(
        manifest_data.get(
            "reference_version",
            "",
        )
    )

    if observed_manifest_version != expected_lambert_version:
        raise ValueError(
            "Unexpected Lambert reference version.\n"
            f"Expected: {expected_lambert_version}\n"
            f"Observed: {observed_manifest_version}"
        )

    print("\nLambert manifest validated:")
    print(f"- Version: {observed_manifest_version}")

else:
    print(
        "\nWARNING:"
        "\n- Lambert manifest was not found."
        "\n- The TF files will still be validated directly."
    )


# ============================================================
# 12. Load official RNA object
# ============================================================

rna = ad.read_h5ad(
    rna_input_path,
    backed="r",
)

if rna.n_obs != expected_n_cells:
    raise ValueError(
        "Unexpected number of RNA cells.\n"
        f"Expected: {expected_n_cells:,}\n"
        f"Observed: {rna.n_obs:,}"
    )

if not rna.obs_names.is_unique:
    raise ValueError(
        "RNA cell identifiers are not unique."
    )

if not rna.var_names.is_unique:
    raise ValueError(
        "RNA feature identifiers are not unique."
    )

if cell_type_column not in rna.obs.columns:
    raise KeyError(
        "Required cell-type annotation was not found.\n"
        f"Required column: {cell_type_column}\n"
        f"Available columns: {list(rna.obs.columns)}"
    )

n_missing_cell_types = int(
    rna.obs[cell_type_column]
    .isna()
    .sum()
)

fraction_missing_cell_types = float(
    rna.obs[cell_type_column]
    .isna()
    .mean()
)

if n_missing_cell_types > 0:
    raise ValueError(
        f"{cell_type_column!r} contains missing values.\n"
        f"Missing cells: {n_missing_cell_types:,}\n"
        f"Missing fraction: "
        f"{fraction_missing_cell_types:.4%}"
    )


# ============================================================
# 13. Determine counts source
# ============================================================

if "counts" in rna.layers:
    counts_source = "layers/counts"
else:
    counts_source = "X"

print("\nRNA input:")
print(
    f"- Shape: {rna.n_obs:,} cells × "
    f"{rna.n_vars:,} genes"
)
print(f"- Counts source: {counts_source}")
print(
    f"- Cell-type column: {cell_type_column}"
)


# ============================================================
# 14. Validate count matrix values
# ============================================================

# In backed mode, rna.X may be an AnnData sparse-dataset
# wrapper rather than a standard scipy sparse matrix.
# We therefore load only the mapped TF matrix later for full
# numeric validation. At this stage, we inspect the declared
# matrix dtype and shape without densifying the full RNA matrix.

if counts_source == "layers/counts":
    full_count_matrix = rna.layers["counts"]
else:
    full_count_matrix = rna.X

print("\nFull RNA count-matrix representation:")
print(f"- Python type: {type(full_count_matrix)}")
print(f"- Shape: {full_count_matrix.shape}")
print(f"- Dtype: {full_count_matrix.dtype}")

# The full 7,566 × 36,601 matrix is not loaded into memory here.
# Numeric validation will be performed after subsetting to the
# mapped TF columns.

integer_like = None

# ============================================================
# 15. Create cell-type overview
# ============================================================

cell_type_series = (
    rna.obs[cell_type_column]
    .astype("string")
)

cell_type_values = sorted(
    cell_type_series.unique().tolist()
)

celltype_overview_df = (
    cell_type_series
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
    / rna.n_obs
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
# 16. Determine RNA gene-symbol universe
# ============================================================

rna_var_names = pd.Index(
    rna.var_names.astype(str),
    name="rna_var_name",
)

candidate_symbol_columns = [
    "gene_symbol",
    "gene_symbols",
    "feature_name",
    "feature_names",
    "symbol",
    "gene_name",
    "gene_names",
]

observed_symbol_columns = [
    column
    for column in candidate_symbol_columns
    if column in rna.var.columns
]

if observed_symbol_columns:
    selected_symbol_column = (
        observed_symbol_columns[0]
    )

    rna_gene_symbols_original = (
        rna.var[selected_symbol_column]
        .astype("string")
        .str.strip()
    )

    gene_symbol_source = (
        f"rna.var[{selected_symbol_column!r}]"
    )

else:
    selected_symbol_column = None

    rna_gene_symbols_original = pd.Series(
        rna_var_names,
        index=rna_var_names,
        dtype="string",
    ).str.strip()

    gene_symbol_source = "rna.var_names"


rna_gene_table = pd.DataFrame(
    {
        "rna_var_name": rna_var_names.astype(str),
        "rna_gene_symbol_original": (
            rna_gene_symbols_original.to_numpy()
        ),
        "rna_position": np.arange(
            rna.n_vars,
            dtype=np.int64,
        ),
    }
)

rna_gene_table[
    "rna_gene_symbol_normalized"
] = (
    rna_gene_table[
        "rna_gene_symbol_original"
    ]
    .astype("string")
    .str.strip()
    .str.upper()
)

rna_gene_table[
    "rna_gene_symbol_missing"
] = (
    rna_gene_table[
        "rna_gene_symbol_normalized"
    ]
    .isna()
    | rna_gene_table[
        "rna_gene_symbol_normalized"
    ]
    .eq("")
)

rna_gene_table[
    "rna_symbol_mapping_multiplicity"
] = (
    rna_gene_table
    .groupby(
        "rna_gene_symbol_normalized",
        dropna=False,
    )[
        "rna_var_name"
    ]
    .transform("size")
)

if "gene_ids" in rna.var.columns:
    rna_gene_table["gene_id"] = (
        rna.var["gene_ids"]
        .astype("string")
        .to_numpy()
    )
else:
    rna_gene_table["gene_id"] = pd.NA


print("\nRNA gene-symbol universe:")
print(f"- Source: {gene_symbol_source}")
print(f"- RNA genes: {len(rna_gene_table):,}")
print(
    "- Missing normalized symbols: "
    f"{rna_gene_table['rna_gene_symbol_missing'].sum():,}"
)
print(
    "- Features participating in duplicated-symbol mappings: "
    f"{(rna_gene_table['rna_symbol_mapping_multiplicity'] > 1).sum():,}"
)


# ============================================================
# 17. Load and clean Lambert TF catalogue
# ============================================================

tf_names_raw = pd.read_csv(
    tf_names_path,
    header=None,
    names=["tf_symbol_original"],
    dtype="string",
)

tf_catalogue = tf_names_raw.copy()

tf_catalogue[
    "tf_symbol_original"
] = (
    tf_catalogue[
        "tf_symbol_original"
    ]
    .astype("string")
    .str.strip()
)

tf_catalogue[
    "tf_symbol_normalized"
] = (
    tf_catalogue[
        "tf_symbol_original"
    ]
    .str.upper()
)

missing_tf_symbols = (
    tf_catalogue[
        "tf_symbol_normalized"
    ]
    .isna()
    | tf_catalogue[
        "tf_symbol_normalized"
    ]
    .eq("")
)

if missing_tf_symbols.any():
    raise ValueError(
        "Lambert TF catalogue contains missing or empty "
        "TF symbols."
    )

duplicated_tf_symbols = (
    tf_catalogue[
        "tf_symbol_normalized"
    ]
    .duplicated(
        keep=False
    )
)

if duplicated_tf_symbols.any():
    duplicated_examples = (
        tf_catalogue.loc[
            duplicated_tf_symbols,
            "tf_symbol_original",
        ]
        .tolist()
    )

    raise ValueError(
        "Lambert TF catalogue contains duplicated symbols.\n"
        f"Examples: {duplicated_examples[:20]}"
    )

tf_catalogue = (
    tf_catalogue
    .sort_values(
        "tf_symbol_normalized"
    )
    .reset_index(
        drop=True
    )
)

tf_catalogue.insert(
    0,
    "catalogue_index",
    np.arange(
        len(tf_catalogue),
        dtype=np.int64,
    ),
)

tf_catalogue[
    "reference_name"
] = "Lambert_human_TFs"

tf_catalogue[
    "reference_version"
] = expected_lambert_version


if len(tf_catalogue) != expected_lambert_tf_count:
    print(
        "\nWARNING:"
        "\n- Lambert TF count differs from the expected "
        f"{expected_lambert_tf_count:,} TFs."
        f"\n- Observed: {len(tf_catalogue):,}"
    )

print("\nLambert TF catalogue:")
print(f"- Curated TFs: {len(tf_catalogue):,}")


# ============================================================
# 18. Build TF-to-RNA mapping
# ============================================================

tf_mapping = tf_catalogue.merge(
    rna_gene_table[
        [
            "rna_var_name",
            "rna_gene_symbol_original",
            "rna_gene_symbol_normalized",
            "rna_position",
            "gene_id",
            "rna_symbol_mapping_multiplicity",
        ]
    ],
    how="left",
    left_on="tf_symbol_normalized",
    right_on="rna_gene_symbol_normalized",
)

tf_mapping[
    "present_in_rna"
] = (
    tf_mapping[
        "rna_var_name"
    ]
    .notna()
)

tf_mapping["mapping_status"] = (
    "not_present_in_rna"
)

tf_mapping.loc[
    tf_mapping["present_in_rna"]
    & (
        tf_mapping[
            "rna_symbol_mapping_multiplicity"
        ]
        == 1
    ),
    "mapping_status",
] = "unique_match"

tf_mapping.loc[
    tf_mapping["present_in_rna"]
    & (
        tf_mapping[
            "rna_symbol_mapping_multiplicity"
        ]
        > 1
    ),
    "mapping_status",
] = "ambiguous_multiple_rna_features"


# ============================================================
# 19. Create master uniquely mapped TF universe
# ============================================================

master_tf_universe = (
    tf_mapping.loc[
        tf_mapping[
            "mapping_status"
        ]
        .eq("unique_match"),
        [
            "catalogue_index",
            "tf_symbol_original",
            "tf_symbol_normalized",
            "rna_var_name",
            "rna_gene_symbol_original",
            "rna_position",
            "gene_id",
            "reference_name",
            "reference_version",
        ],
    ]
    .sort_values(
        "tf_symbol_normalized"
    )
    .reset_index(
        drop=True
    )
)

master_tf_universe.insert(
    0,
    "master_tf_index",
    np.arange(
        len(master_tf_universe),
        dtype=np.int64,
    ),
)

if master_tf_universe.empty:
    raise ValueError(
        "No Lambert TFs were uniquely mapped to the RNA "
        "gene universe."
    )

tf_positions = (
    master_tf_universe[
        "rna_position"
    ]
    .to_numpy(
        dtype=np.int64
    )
)

print("\nTF-to-RNA mapping:")
print(f"- Lambert TFs: {len(tf_catalogue):,}")
print(
    "- Unique RNA matches: "
    f"{len(master_tf_universe):,}"
)
print(
    "- Absent from RNA: "
    f"{tf_mapping['mapping_status'].eq('not_present_in_rna').sum():,}"
)
print(
    "- Ambiguous RNA mappings: "
    f"{tf_mapping['mapping_status'].eq('ambiguous_multiple_rna_features').sum():,}"
)


# ============================================================
# 20. Load mapped TF count matrix into memory
# ============================================================

print("\nLoading mapped TF count matrix:")

if counts_source == "layers/counts":
    tf_count_matrix = (
        rna[
            :,
            tf_positions,
        ]
        .layers["counts"]
    )
else:
    tf_count_matrix = (
        rna[
            :,
            tf_positions,
        ]
        .X
    )

if sp.issparse(tf_count_matrix):
    tf_count_matrix = tf_count_matrix.tocsr()
else:
    tf_count_matrix = np.asarray(
        tf_count_matrix
    )

print(
    f"- Shape: {tf_count_matrix.shape[0]:,} cells × "
    f"{tf_count_matrix.shape[1]:,} TFs"
)


# ============================================================
# 20A. Validate mapped TF count matrix values
# ============================================================

if sp.issparse(tf_count_matrix):
    stored_values = np.asarray(
        tf_count_matrix.data
    )
else:
    stored_values = np.asarray(
        tf_count_matrix
    ).ravel()

if stored_values.size > 0:
    if not np.issubdtype(
        stored_values.dtype,
        np.number,
    ):
        raise TypeError(
            "Mapped TF count matrix does not contain a "
            "numeric dtype.\n"
            f"Observed dtype: {stored_values.dtype}"
        )

    if np.isnan(stored_values).any():
        raise ValueError(
            "Mapped TF count matrix contains NaN values."
        )

    if np.isinf(stored_values).any():
        raise ValueError(
            "Mapped TF count matrix contains infinite values."
        )

    if (stored_values < 0).any():
        raise ValueError(
            "Mapped TF count matrix contains negative values."
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

    min_stored_value = float(
        stored_values.min()
    )

    max_stored_value = float(
        stored_values.max()
    )

else:
    integer_like = True
    min_stored_value = 0.0
    max_stored_value = 0.0

print("\nMapped TF count-matrix validation:")
print(f"- Numeric dtype: {stored_values.dtype}")
print(f"- Integer-like: {integer_like}")
print(f"- Minimum stored value: {min_stored_value}")
print(f"- Maximum stored value: {max_stored_value}")
print("- Contains NaN: False")
print("- Contains Inf: False")
print("- Contains negative values: False")

if not integer_like:
    print(
        "\nWARNING:"
        "\n- The mapped TF matrix is non-negative and finite,"
        "\n  but its values are not integer-like."
        "\n- Detection based on >0 remains valid, but total"
        "\n  counts may represent normalized expression."
    )

# ============================================================
# 21. Calculate global TF metrics
# ============================================================

global_arrays = calculate_sparse_gene_metrics(
    tf_count_matrix
)

global_metrics_df = (
    master_tf_universe[
        [
            "master_tf_index",
            "tf_symbol_original",
            "tf_symbol_normalized",
            "rna_var_name",
            "gene_id",
        ]
    ]
    .copy()
)

global_metrics_df[
    "n_total_cells"
] = global_arrays["n_cells"]

global_metrics_df[
    "n_expressing_cells_global"
] = global_arrays[
    "n_expressing_cells"
]

global_metrics_df[
    "fraction_expressing_global"
] = global_arrays[
    "fraction_expressing"
]

global_metrics_df[
    "percent_expressing_global"
] = (
    100.0
    * global_metrics_df[
        "fraction_expressing_global"
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
    "mean_counts_expressing_cells_global"
] = global_arrays[
    "mean_counts_expressing_cells"
]


# ============================================================
# 22. Calculate TF metrics per cell_type_broad
# ============================================================

celltype_metric_frames: list[pd.DataFrame] = []

for cell_type_value in cell_type_values:
    cell_mask = (
        cell_type_series
        .eq(cell_type_value)
        .to_numpy()
    )

    n_group_cells = int(
        cell_mask.sum()
    )

    if n_group_cells == 0:
        raise RuntimeError(
            "Encountered an empty cell_type_broad group:\n"
            f"{cell_type_value}"
        )

    group_matrix = tf_count_matrix[
        cell_mask,
        :,
    ]

    group_arrays = calculate_sparse_gene_metrics(
        group_matrix
    )

    group_df = (
        master_tf_universe[
            [
                "master_tf_index",
                "tf_symbol_original",
                "tf_symbol_normalized",
                "rna_var_name",
                "gene_id",
            ]
        ]
        .copy()
    )

    group_df[
        cell_type_column
    ] = cell_type_value

    group_df[
        "n_cells_in_group"
    ] = group_arrays["n_cells"]

    group_df[
        "n_expressing_cells"
    ] = group_arrays[
        "n_expressing_cells"
    ]

    group_df[
        "fraction_expressing"
    ] = group_arrays[
        "fraction_expressing"
    ]

    group_df[
        "percent_expressing"
    ] = (
        100.0
        * group_df[
            "fraction_expressing"
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
        "mean_counts_expressing_cells"
    ] = group_arrays[
        "mean_counts_expressing_cells"
    ]

    celltype_metric_frames.append(
        group_df
    )

celltype_metrics_df = pd.concat(
    celltype_metric_frames,
    ignore_index=True,
)

celltype_metrics_df = celltype_metrics_df[
    [
        "master_tf_index",
        "tf_symbol_original",
        "tf_symbol_normalized",
        "rna_var_name",
        "gene_id",
        cell_type_column,
        "n_cells_in_group",
        "n_expressing_cells",
        "fraction_expressing",
        "percent_expressing",
        "total_counts",
        "mean_counts_all_cells",
        "mean_counts_expressing_cells",
    ]
]


# ============================================================
# 23. Validate global vs grouped totals
# ============================================================

grouped_count_check = (
    celltype_metrics_df
    .groupby(
        "master_tf_index",
        sort=True,
    )[
        [
            "n_expressing_cells",
            "total_counts",
        ]
    ]
    .sum()
    .reset_index()
)

global_count_check = (
    global_metrics_df[
        [
            "master_tf_index",
            "n_expressing_cells_global",
            "total_counts_global",
        ]
    ]
    .sort_values(
        "master_tf_index"
    )
    .reset_index(
        drop=True
    )
)

validation_df = global_count_check.merge(
    grouped_count_check,
    how="left",
    on="master_tf_index",
    validate="one_to_one",
)

if not np.array_equal(
    validation_df[
        "n_expressing_cells_global"
    ].to_numpy(),
    validation_df[
        "n_expressing_cells"
    ].to_numpy(),
):
    raise RuntimeError(
        "Global expressing-cell counts do not equal the sum "
        "across cell_type_broad groups."
    )

if not np.allclose(
    validation_df[
        "total_counts_global"
    ].to_numpy(),
    validation_df[
        "total_counts"
    ].to_numpy(),
    rtol=1e-8,
    atol=1e-8,
):
    raise RuntimeError(
        "Global TF total counts do not equal the sum across "
        "cell_type_broad groups."
    )


# ============================================================
# 24. Global threshold diagnostics
# ============================================================

global_threshold_summary_df = build_threshold_summary(
    detection_fractions=global_metrics_df[
        "fraction_expressing_global"
    ],
    n_expressing_cells=global_metrics_df[
        "n_expressing_cells_global"
    ],
    thresholds=detection_thresholds,
    group_label="scope",
    group_value="global_all_cells",
)


# ============================================================
# 25. Cell-type threshold diagnostics
# ============================================================

celltype_threshold_frames: list[pd.DataFrame] = []

for cell_type_value in cell_type_values:
    group_metrics = celltype_metrics_df.loc[
        celltype_metrics_df[
            cell_type_column
        ]
        .eq(cell_type_value)
    ]

    threshold_df = build_threshold_summary(
    detection_fractions=group_metrics[
        "fraction_expressing"
    ],
    n_expressing_cells=group_metrics[
        "n_expressing_cells"
    ],
    thresholds=detection_thresholds,
    group_label=cell_type_column,
    group_value=cell_type_value,
   )

    threshold_df.insert(
        1,
        "n_cells_in_group",
        int(
            group_metrics[
                "n_cells_in_group"
            ]
            .iloc[0]
        ),
    )

    celltype_threshold_frames.append(
        threshold_df
    )

celltype_threshold_summary_df = pd.concat(
    celltype_threshold_frames,
    ignore_index=True,
)


# ============================================================
# 26. Create mapping exclusion report
# ============================================================

mapping_exclusion_report_df = (
    tf_mapping.loc[
        ~tf_mapping[
            "mapping_status"
        ]
        .eq("unique_match")
    ]
    .copy()
)

mapping_exclusion_report_df[
    "exclusion_scope"
] = "mapping_only"

mapping_exclusion_report_df[
    "important_note"
] = (
    "No TF was excluded because of expression detection. "
    "Only absent or ambiguous RNA symbol mappings appear "
    "in this report."
)


# ============================================================
# 27. Create input audit
# ============================================================

input_audit_df = pd.DataFrame(
    [
        {
            "input": "full_gene_final_cell_rna",
            "path": str(rna_input_path),
            "n_obs": int(rna.n_obs),
            "n_vars": int(rna.n_vars),
            "file_size_bytes": int(
                rna_input_path.stat().st_size
            ),
            "sha256": sha256_file(
                rna_input_path
            ),
        },
        {
            "input": "lambert_tf_names",
            "path": str(tf_names_path),
            "n_obs": int(
                len(tf_catalogue)
            ),
            "n_vars": pd.NA,
            "file_size_bytes": int(
                tf_names_path.stat().st_size
            ),
            "sha256": sha256_file(
                tf_names_path
            ),
        },
        {
            "input": "lambert_complete_database",
            "path": str(tf_database_path),
            "n_obs": pd.NA,
            "n_vars": pd.NA,
            "file_size_bytes": int(
                tf_database_path.stat().st_size
            ),
            "sha256": sha256_file(
                tf_database_path
            ),
        },
    ]
)


# ============================================================
# 28. Save tables
# ============================================================

tf_catalogue.to_csv(
    tf_catalogue_clean_path,
    sep="\t",
    index=False,
)

tf_mapping.to_csv(
    tf_mapping_path,
    sep="\t",
    index=False,
)

master_tf_universe.to_csv(
    tf_master_universe_path,
    sep="\t",
    index=False,
)

global_metrics_df.to_csv(
    global_metrics_path,
    sep="\t",
    index=False,
)

celltype_metrics_df.to_csv(
    celltype_metrics_path,
    sep="\t",
    index=False,
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

celltype_overview_df.to_csv(
    celltype_overview_path,
    sep="\t",
    index=False,
)

mapping_exclusion_report_df.to_csv(
    mapping_exclusion_report_path,
    sep="\t",
    index=False,
)

input_audit_df.to_csv(
    input_audit_path,
    sep="\t",
    index=False,
)


# ============================================================
# 29. Plot styling
# ============================================================

sns.set_theme(
    style="whitegrid",
    context="notebook",
)

threshold_percent_values = [
    threshold * 100
    for threshold in detection_thresholds
]


# ============================================================
# 30. Plot TF retention curves by cell_type_broad
# ============================================================

plt.figure(
    figsize=(12, 8)
)

sns.lineplot(
    data=celltype_threshold_summary_df,
    x="threshold_percent",
    y="pct_tfs_passing_of_present",
    hue=cell_type_column,
    marker="o",
)
plt.xlabel(
    "TF detection threshold within cell type (%)"
)

plt.ylabel(
    "TFs retained among TFs present in cell type (%)"
)

plt.title(
    "TF retention across detection thresholds "
    "by cell_type_broad"
)

plt.xticks(
    threshold_percent_values
)

plt.ylim(
    0,
    100,
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
# 31. Plot TF-retention heatmap
# ============================================================

retention_heatmap_df = (
    celltype_threshold_summary_df
    .pivot(
        index=cell_type_column,
        columns="threshold_percent",
        values="pct_tfs_passing_of_present",
    )
    .reindex(
        index=cell_type_values
    )
)

plt.figure(
    figsize=(10, 7)
)

sns.heatmap(
    retention_heatmap_df,
    annot=True,
    fmt=".1f",
    cmap="viridis",
    vmin=0,
    vmax=100,
    cbar_kws={
        "label": "TFs retained among TFs present (%)",
    },
)

plt.xlabel(
    "Detection threshold within cell type (%)"
)

plt.ylabel(
    cell_type_column
)

plt.title(
    "Mapped TF retention by cell type and threshold"
)

plt.tight_layout()

plt.savefig(
    retention_heatmap_path,
    dpi=300,
    bbox_inches="tight",
)

plt.close()


# ============================================================
# 32. Plot global TF detection distribution
# ============================================================

plt.figure(
    figsize=(11, 7)
)

sns.histplot(
    data=global_metrics_df,
    x="percent_expressing_global",
    bins=50,
    color="#3B82F6",
)

for threshold_percent in threshold_percent_values:
    plt.axvline(
        threshold_percent,
        linestyle="--",
        linewidth=1.2,
        alpha=0.8,
        label=f"{threshold_percent:g}%",
    )

plt.xlabel(
    "Cells expressing TF globally (%)"
)

plt.ylabel(
    "Number of mapped TFs"
)

plt.title(
    "Global TF detection distribution"
)

plt.legend(
    title="Diagnostic thresholds"
)

plt.tight_layout()

plt.savefig(
    global_distribution_path,
    dpi=300,
    bbox_inches="tight",
)

plt.close()


# ============================================================
# 33. Build summary statistics
# ============================================================

n_lambert_tfs = int(
    len(tf_catalogue)
)

n_unique_mapped_tfs = int(
    len(master_tf_universe)
)

n_absent_tfs = int(
    tf_mapping[
        "mapping_status"
    ]
    .eq("not_present_in_rna")
    .sum()
)

n_ambiguous_tfs = int(
    tf_mapping[
        "mapping_status"
    ]
    .eq(
        "ambiguous_multiple_rna_features"
    )
    .sum()
)

n_zero_detected_global = int(
    global_metrics_df[
        "n_expressing_cells_global"
    ]
    .eq(0)
    .sum()
)

global_detection_quantiles = (
    global_metrics_df[
        "fraction_expressing_global"
    ]
    .quantile(
        [
            0.00,
            0.25,
            0.50,
            0.75,
            0.90,
            1.00,
        ]
    )
    .to_dict()
)


# ============================================================
# 34. Save YAML summary
# ============================================================

summary = {
    "step": "05A.1",
    "script": script_path.name,
    "created_at_utc": utc_timestamp(),
    "project_id": paths.project_id,
    "objective": (
        "Map the Lambert human TF catalogue to the official "
        "full-gene PBMC RNA object and calculate global and "
        "cell_type_broad expression-detection statistics."
    ),
    "important_policy": {
        "all_uniquely_mapped_tfs_retained": True,
        "expression_threshold_filtering_applied": False,
        "thresholds_are_diagnostic_only": True,
        "final_tf_selection_deferred": True,
    },
    "inputs": {
        "rna_h5ad": str(
            rna_input_path
        ),
        "lambert_tf_names": str(
            tf_names_path
        ),
        "lambert_database": str(
            tf_database_path
        ),
        "lambert_manifest": (
            str(tf_manifest_path)
            if tf_manifest_path.exists()
            else None
        ),
    },
    "rna": {
        "n_cells": int(
            rna.n_obs
        ),
        "n_genes": int(
            rna.n_vars
        ),
        "counts_source": counts_source,
        "counts_integer_like": bool(
            integer_like
        ),
        "gene_symbol_source": (
            gene_symbol_source
        ),
        "cell_type_column": (
            cell_type_column
        ),
        "n_cell_types": int(
            len(cell_type_values)
        ),
        "cell_types": cell_type_values,
    },
    "tf_mapping": {
        "n_lambert_tfs": n_lambert_tfs,
        "n_unique_mapped_tfs": (
            n_unique_mapped_tfs
        ),
        "n_absent_from_rna": (
            n_absent_tfs
        ),
        "n_ambiguous_rna_mappings": (
            n_ambiguous_tfs
        ),
        "mapped_fraction_of_catalogue": (
            n_unique_mapped_tfs
            / n_lambert_tfs
        ),
    },
    "detection": {
        "n_zero_detected_global": (
            n_zero_detected_global
        ),
        "diagnostic_thresholds": (
            detection_thresholds
        ),
        "diagnostic_thresholds_percent": (
            threshold_percent_values
        ),
        "global_detection_fraction_quantiles": {
            str(key): float(value)
            for key, value
            in global_detection_quantiles.items()
        },
    },
    "outputs": {
        "tf_mapping": str(
            tf_mapping_path
        ),
        "master_tf_universe": str(
            tf_master_universe_path
        ),
        "global_metrics": str(
            global_metrics_path
        ),
        "celltype_metrics": str(
            celltype_metrics_path
        ),
        "global_threshold_summary": str(
            global_threshold_summary_path
        ),
        "celltype_threshold_summary": str(
            celltype_threshold_summary_path
        ),
        "retention_curve": str(
            retention_curve_path
        ),
        "retention_heatmap": str(
            retention_heatmap_path
        ),
    },
    "scope_exclusions": [
        "No TF expression filtering",
        "No final TF node selection",
        "No GENCODE processing",
        "No promoter or TSS construction",
        "No JASPAR motif analysis",
        "No ATAC peak processing",
        "No adjacency construction",
        "No GNN training",
    ],
}

with summary_yaml_path.open(
    "w",
    encoding="utf-8",
) as handle:
    yaml.safe_dump(
        json_safe(summary),
        handle,
        sort_keys=False,
        allow_unicode=True,
    )


# ============================================================
# 35. Save run metadata
# ============================================================

run_metadata = {
    "script": script_path.name,
    "script_path": str(script_path),
    "completed_at_utc": utc_timestamp(),
    "project_id": paths.project_id,
    "project_root": str(paths.project_root),
    "python_version": sys.version,
    "numpy_version": np.__version__,
    "pandas_version": pd.__version__,
    "scipy_version": version("scipy"),
    "anndata_version": version("anndata"),
    "matplotlib_version": version("matplotlib"),
    "seaborn_version": version("seaborn"),
    "n_lambert_tfs": n_lambert_tfs,
    "n_unique_mapped_tfs": (
        n_unique_mapped_tfs
    ),
    "n_cell_types": int(
        len(cell_type_values)
    ),
    "status": "completed",
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
# 36. Validate outputs
# ============================================================

required_outputs = {
    "Input audit": input_audit_path,
    "Clean TF catalogue": tf_catalogue_clean_path,
    "TF-to-RNA mapping": tf_mapping_path,
    "Master mapped TF universe": (
        tf_master_universe_path
    ),
    "Global TF metrics": global_metrics_path,
    "Cell-type TF metrics": celltype_metrics_path,
    "Global threshold summary": (
        global_threshold_summary_path
    ),
    "Cell-type threshold summary": (
        celltype_threshold_summary_path
    ),
    "Cell-type overview": celltype_overview_path,
    "Mapping exclusion report": (
        mapping_exclusion_report_path
    ),
    "TF retention curve": retention_curve_path,
    "TF retention heatmap": retention_heatmap_path,
    "Global TF detection distribution": (
        global_distribution_path
    ),
    "Summary YAML": summary_yaml_path,
    "Run metadata JSON": run_metadata_path,
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
# 37. Close RNA object
# ============================================================

close_backed_anndata(
    rna
)


# ============================================================
# 38. Final report
# ============================================================

print("\n" + "=" * 72)
print("05A.1 COMPLETED SUCCESSFULLY")
print("=" * 72)

print("\nRNA input:")
print(
    f"{expected_n_cells:,} cells × "
    f"{len(rna_gene_table):,} genes"
)

print("\nLambert catalogue:")
print(f"- Curated TFs: {n_lambert_tfs:,}")
print(
    f"- Uniquely mapped TFs: "
    f"{n_unique_mapped_tfs:,}"
)
print(
    f"- Absent from RNA: "
    f"{n_absent_tfs:,}"
)
print(
    f"- Ambiguous mappings: "
    f"{n_ambiguous_tfs:,}"
)

print("\nCell-type analysis:")
print(
    f"- Metadata column: "
    f"{cell_type_column}"
)
print(
    f"- Broad cell types: "
    f"{len(cell_type_values):,}"
)

print("\nDiagnostic thresholds:")
print(
    ", ".join(
        f"{threshold * 100:g}%"
        for threshold
        in detection_thresholds
    )
)

print("\nImportant:")
print(
    "- No mapped TF was excluded using expression detection."
)
print(
    "- Threshold summaries are diagnostic only."
)
print(
    "- Final TF selection is deferred to 05B.1."
)

print("\nPrimary outputs:")
print(tf_master_universe_path)
print(global_metrics_path)
print(celltype_metrics_path)
print(celltype_threshold_summary_path)

print("\nNext stage:")
print(
    "Inspect the retention summaries and plots before "
    "configuring global or per-cell_type_broad selection "
    "in 05B.1."
)

# ============================================================
# 39. Display TF presence per cell_type_broad
# ============================================================

tf_presence_overview = (
    celltype_threshold_summary_df.loc[
        celltype_threshold_summary_df[
            "threshold"
        ].eq(
            detection_thresholds[0]
        ),
        [
            "cell_type_broad",
            "n_cells_in_group",
            "n_master_mapped_tfs",
            "n_tfs_present_in_group",
            "pct_tfs_present_of_master",
        ],
    ]
    .drop_duplicates(
        subset=[
            "cell_type_broad",
        ]
    )
    .sort_values(
        "n_tfs_present_in_group",
        ascending=False,
    )
    .reset_index(
        drop=True
    )
)

print("\n" + "=" * 90)
print("TFs PRESENT IN EACH CELL_TYPE_BROAD")
print("=" * 90)

print(
    tf_presence_overview.to_string(
        index=False
    )
)


# ============================================================
# 40. Display 5% and 10% TF-retention comparison
# ============================================================

comparison_5_10 = (
    celltype_threshold_summary_df.loc[
        celltype_threshold_summary_df[
            "threshold"
        ].isin(
            [0.05, 0.10]
        ),
        [
            "cell_type_broad",
            "n_cells_in_group",
            "threshold_percent",
            "n_master_mapped_tfs",
            "n_tfs_present_in_group",
            "pct_tfs_present_of_master",
            "n_tfs_present_but_below_threshold",
            "n_tfs_passing_threshold",
            "pct_tfs_passing_of_master",
            "pct_tfs_passing_of_present",
        ],
    ]
    .sort_values(
        [
            "cell_type_broad",
            "threshold_percent",
        ]
    )
    .reset_index(
        drop=True
    )
)

print("\n" + "=" * 110)
print("TF RETENTION COMPARISON: 5% VS 10%")
print("=" * 110)

print(
    comparison_5_10.to_string(
        index=False
    )
)


# ============================================================
# 41. Display wide 5% vs 10% retention table
# ============================================================

retention_wide = (
    celltype_threshold_summary_df.loc[
        celltype_threshold_summary_df[
            "threshold"
        ].isin(
            [0.05, 0.10]
        )
    ]
    .pivot(
        index=[
            "cell_type_broad",
            "n_cells_in_group",
            "n_tfs_present_in_group",
        ],
        columns="threshold_percent",
        values="pct_tfs_passing_of_present",
    )
    .rename(
        columns={
            5.0: "pct_present_tfs_retained_at_5pct",
            10.0: "pct_present_tfs_retained_at_10pct",
        }
    )
    .reset_index()
)

retention_wide[
    "retention_loss_5_to_10_percentage_points"
] = (
    retention_wide[
        "pct_present_tfs_retained_at_5pct"
    ]
    - retention_wide[
        "pct_present_tfs_retained_at_10pct"
    ]
)

retention_wide = (
    retention_wide
    .sort_values(
        "retention_loss_5_to_10_percentage_points",
        ascending=False,
    )
    .reset_index(
        drop=True
    )
)

print("\n" + "=" * 110)
print("WIDE TF RETENTION TABLE")
print("=" * 110)

print(
    retention_wide.to_string(
        index=False
    )
)


# ============================================================
# Marker diagnostic 1. Build marker-gene catalogue
# ============================================================

marker_records: list[dict[str, Any]] = []

for marker_group, genes in marker_gene_sets.items():
    for marker_order, gene_symbol in enumerate(genes):
        marker_records.append(
            {
                "marker_group": marker_group,
                "marker_order": marker_order,
                "marker_gene_original": gene_symbol,
                "marker_gene_normalized": (
                    str(gene_symbol)
                    .strip()
                    .upper()
                ),
            }
        )

marker_catalogue = pd.DataFrame(marker_records)

if marker_catalogue.empty:
    raise ValueError(
        "The marker-gene catalogue is empty."
    )

# The same gene can legitimately belong to multiple marker groups,
# for example CCL5, GZMB, LST1, and FCGR3A.
marker_catalogue[
    "marker_group_gene_id"
] = (
    marker_catalogue["marker_group"]
    + "::"
    + marker_catalogue["marker_gene_normalized"]
)

if marker_catalogue[
    "marker_group_gene_id"
].duplicated().any():
    duplicated_rows = marker_catalogue.loc[
        marker_catalogue[
            "marker_group_gene_id"
        ].duplicated(
            keep=False
        )
    ]

    raise ValueError(
        "A marker gene is duplicated inside the same marker "
        "group.\n"
        f"{duplicated_rows.to_string(index=False)}"
    )


# ============================================================
# Marker diagnostic 2. Map markers to full RNA gene universe
# ============================================================

marker_mapping = marker_catalogue.merge(
    rna_gene_table[
        [
            "rna_var_name",
            "rna_gene_symbol_original",
            "rna_gene_symbol_normalized",
            "rna_position",
            "gene_id",
            "rna_symbol_mapping_multiplicity",
        ]
    ],
    how="left",
    left_on="marker_gene_normalized",
    right_on="rna_gene_symbol_normalized",
)

marker_mapping[
    "present_in_rna_gene_universe"
] = marker_mapping[
    "rna_var_name"
].notna()

marker_mapping["mapping_status"] = "not_present_in_rna"

marker_mapping.loc[
    marker_mapping[
        "present_in_rna_gene_universe"
    ]
    & (
        marker_mapping[
            "rna_symbol_mapping_multiplicity"
        ]
        == 1
    ),
    "mapping_status",
] = "unique_match"

marker_mapping.loc[
    marker_mapping[
        "present_in_rna_gene_universe"
    ]
    & (
        marker_mapping[
            "rna_symbol_mapping_multiplicity"
        ]
        > 1
    ),
    "mapping_status",
] = "ambiguous_multiple_rna_features"


# ============================================================
# Marker diagnostic 3. Create unique mapped marker universe
# ============================================================

# A gene may occur in several marker groups.
# Expression is loaded only once per unique gene.

unique_marker_universe = (
    marker_mapping.loc[
        marker_mapping["mapping_status"].eq(
            "unique_match"
        ),
        [
            "marker_gene_normalized",
            "rna_var_name",
            "rna_gene_symbol_original",
            "rna_position",
            "gene_id",
        ],
    ]
    .drop_duplicates(
        subset=["marker_gene_normalized"]
    )
    .sort_values(
        "marker_gene_normalized"
    )
    .reset_index(
        drop=True
    )
)

unique_marker_universe.insert(
    0,
    "marker_gene_index",
    np.arange(
        len(unique_marker_universe),
        dtype=np.int64,
    ),
)

if unique_marker_universe.empty:
    raise ValueError(
        "None of the marker genes mapped uniquely to the "
        "full RNA gene universe."
    )

marker_positions = (
    unique_marker_universe["rna_position"]
    .to_numpy(
        dtype=np.int64
    )
)

print("\nMarker-gene mapping:")
print(
    f"- Marker-group entries: "
    f"{len(marker_catalogue):,}"
)
print(
    f"- Unique requested marker genes: "
    f"{marker_catalogue['marker_gene_normalized'].nunique():,}"
)
print(
    f"- Unique mapped marker genes: "
    f"{len(unique_marker_universe):,}"
)
print(
    "- Marker entries absent from RNA: "
    f"{marker_mapping['mapping_status'].eq('not_present_in_rna').sum():,}"
)
print(
    "- Ambiguous marker entries: "
    f"{marker_mapping['mapping_status'].eq('ambiguous_multiple_rna_features').sum():,}"
)


# ============================================================
# Marker diagnostic 4. Load marker count matrix
# ============================================================

if counts_source == "layers/counts":
    marker_count_matrix = (
        rna[
            :,
            marker_positions,
        ]
        .layers["counts"]
    )
else:
    marker_count_matrix = (
        rna[
            :,
            marker_positions,
        ]
        .X
    )

if sp.issparse(marker_count_matrix):
    marker_count_matrix = marker_count_matrix.tocsr()
else:
    marker_count_matrix = np.asarray(
        marker_count_matrix
    )

print(
    "\nMarker count matrix:"
    f"\n- Shape: {marker_count_matrix.shape[0]:,} cells × "
    f"{marker_count_matrix.shape[1]:,} marker genes"
)


# ============================================================
# Marker diagnostic 5. Calculate global marker metrics
# ============================================================

marker_global_arrays = calculate_sparse_gene_metrics(
    marker_count_matrix
)

marker_global_metrics_df = (
    unique_marker_universe[
        [
            "marker_gene_index",
            "marker_gene_normalized",
            "rna_var_name",
            "gene_id",
        ]
    ]
    .copy()
)

marker_global_metrics_df[
    "n_total_cells"
] = marker_global_arrays["n_cells"]

marker_global_metrics_df[
    "n_expressing_cells_global"
] = marker_global_arrays[
    "n_expressing_cells"
]

marker_global_metrics_df[
    "fraction_expressing_global"
] = marker_global_arrays[
    "fraction_expressing"
]

marker_global_metrics_df[
    "percent_expressing_global"
] = (
    100.0
    * marker_global_metrics_df[
        "fraction_expressing_global"
    ]
)

marker_global_metrics_df[
    "total_counts_global"
] = marker_global_arrays[
    "total_counts"
]

marker_global_metrics_df[
    "mean_counts_all_cells_global"
] = marker_global_arrays[
    "mean_counts_all_cells"
]

marker_global_metrics_df[
    "mean_counts_expressing_cells_global"
] = marker_global_arrays[
    "mean_counts_expressing_cells"
]


# ============================================================
# Marker diagnostic 6. Calculate metrics per cell_type_broad
# ============================================================

marker_celltype_frames: list[pd.DataFrame] = []

cell_type_array = cell_type_series.to_numpy()

for cell_type_value in cell_type_values:
    cell_mask = (
        cell_type_array == cell_type_value
    )

    n_cells_in_group = int(
        np.sum(cell_mask)
    )

    if n_cells_in_group == 0:
        continue

    group_matrix = marker_count_matrix[
        cell_mask,
        :,
    ]

    group_arrays = calculate_sparse_gene_metrics(
        group_matrix
    )

    group_df = (
        unique_marker_universe[
            [
                "marker_gene_index",
                "marker_gene_normalized",
                "rna_var_name",
                "gene_id",
            ]
        ]
        .copy()
    )

    group_df[cell_type_column] = cell_type_value

    group_df[
        "n_cells_in_group"
    ] = group_arrays["n_cells"]

    group_df[
        "n_expressing_cells"
    ] = group_arrays[
        "n_expressing_cells"
    ]

    group_df[
        "fraction_expressing"
    ] = group_arrays[
        "fraction_expressing"
    ]

    group_df[
        "percent_expressing"
    ] = (
        100.0
        * group_df["fraction_expressing"]
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
        "mean_counts_expressing_cells"
    ] = group_arrays[
        "mean_counts_expressing_cells"
    ]

    group_df[
        "detected_in_cluster"
    ] = (
        group_df["n_expressing_cells"] > 0
    )

    for threshold in detection_thresholds:
        threshold_label = (
            f"detected_at_least_"
            f"{int(round(threshold * 100))}pct"
        )

        group_df[threshold_label] = (
            group_df["fraction_expressing"]
            >= threshold
        )

    marker_celltype_frames.append(
        group_df
    )

marker_celltype_metrics_df = pd.concat(
    marker_celltype_frames,
    axis=0,
    ignore_index=True,
)

if marker_celltype_metrics_df.empty:
    raise ValueError(
        "No marker metrics were calculated by cell type."
    )
    
# ============================================================
# Marker diagnostic 7. Attach marker-group annotations
# ============================================================

marker_group_annotations = (
    marker_catalogue[
        [
            "marker_group",
            "marker_order",
            "marker_gene_original",
            "marker_gene_normalized",
        ]
    ]
    .copy()
)

marker_celltype_annotated_df = (
    marker_group_annotations.merge(
        marker_celltype_metrics_df,
        how="left",
        on="marker_gene_normalized",
        validate="many_to_many",
    )
)

marker_global_annotated_df = (
    marker_group_annotations.merge(
        marker_global_metrics_df,
        how="left",
        on="marker_gene_normalized",
        validate="many_to_one",
    )
)

marker_celltype_annotated_df = (
    marker_celltype_annotated_df.sort_values(
        [
            "marker_group",
            "marker_order",
            cell_type_column,
        ]
    )
    .reset_index(
        drop=True
    )
)

# ============================================================
# Marker diagnostic 8. Create detection matrices
# ============================================================

marker_presence_matrix_df = (
    marker_celltype_metrics_df
    .pivot(
        index="marker_gene_normalized",
        columns=cell_type_column,
        values="detected_in_cluster",
    )
    .fillna(False)
    .astype(bool)
)

marker_percent_matrix_df = (
    marker_celltype_metrics_df
    .pivot(
        index="marker_gene_normalized",
        columns=cell_type_column,
        values="percent_expressing",
    )
    .fillna(0.0)
)

# Keep marker genes in biological-group order rather than
# alphabetical order.
ordered_marker_genes = (
    marker_catalogue[
        "marker_gene_normalized"
    ]
    .drop_duplicates()
    .tolist()
)

ordered_marker_genes_present = [
    gene
    for gene in ordered_marker_genes
    if gene in marker_percent_matrix_df.index
]

marker_presence_matrix_df = (
    marker_presence_matrix_df.loc[
        ordered_marker_genes_present
    ]
)

marker_percent_matrix_df = (
    marker_percent_matrix_df.loc[
        ordered_marker_genes_present
    ]
)

# ============================================================
# Marker diagnostic 9. Summarize clusters per marker gene
# ============================================================

marker_cluster_summary_records: list[dict[str, Any]] = []

for marker_group, marker_order, marker_gene in (
    marker_catalogue[
        [
            "marker_group",
            "marker_order",
            "marker_gene_normalized",
        ]
    ]
    .itertuples(
        index=False,
        name=None,
    )
):
    mapping_rows = marker_mapping.loc[
        (
            marker_mapping["marker_group"]
            == marker_group
        )
        & (
            marker_mapping[
                "marker_gene_normalized"
            ]
            == marker_gene
        )
    ]

    mapping_statuses = (
        mapping_rows["mapping_status"]
        .dropna()
        .unique()
        .tolist()
    )

    mapping_status = (
        mapping_statuses[0]
        if len(mapping_statuses) == 1
        else ";".join(mapping_statuses)
    )

    gene_metrics = marker_celltype_metrics_df.loc[
        marker_celltype_metrics_df[
            "marker_gene_normalized"
        ].eq(marker_gene)
    ].copy()

    detected_metrics = (
        gene_metrics.loc[
            gene_metrics[
                "detected_in_cluster"
            ]
        ]
        .sort_values(
            "percent_expressing",
            ascending=False,
        )
    )

    clusters_detected = (
        detected_metrics[
            cell_type_column
        ]
        .astype(str)
        .tolist()
    )

    clusters_at_1pct = (
        gene_metrics.loc[
            gene_metrics[
                "fraction_expressing"
            ]
            >= 0.01,
            cell_type_column,
        ]
        .astype(str)
        .tolist()
    )

    clusters_at_5pct = (
        gene_metrics.loc[
            gene_metrics[
                "fraction_expressing"
            ]
            >= 0.05,
            cell_type_column,
        ]
        .astype(str)
        .tolist()
    )

    clusters_at_10pct = (
        gene_metrics.loc[
            gene_metrics[
                "fraction_expressing"
            ]
            >= 0.10,
            cell_type_column,
        ]
        .astype(str)
        .tolist()
    )

    if detected_metrics.empty:
        highest_cluster = None
        highest_percent = 0.0
        highest_n_expressing = 0
    else:
        highest_row = detected_metrics.iloc[0]

        highest_cluster = str(
            highest_row[cell_type_column]
        )

        highest_percent = float(
            highest_row["percent_expressing"]
        )

        highest_n_expressing = int(
            highest_row["n_expressing_cells"]
        )

    marker_cluster_summary_records.append(
        {
            "marker_group": marker_group,
            "marker_order": marker_order,
            "marker_gene": marker_gene,
            "mapping_status": mapping_status,

            "present_in_rna_gene_universe": (
                mapping_status == "unique_match"
            ),

            "detected_in_any_cell": (
                len(clusters_detected) > 0
            ),

            "n_clusters_detected": len(
                clusters_detected
            ),

            "clusters_detected_any_cell": ";".join(
                clusters_detected
            ),

            "n_clusters_at_least_1pct": len(
                clusters_at_1pct
            ),

            "clusters_at_least_1pct": ";".join(
                clusters_at_1pct
            ),

            "n_clusters_at_least_5pct": len(
                clusters_at_5pct
            ),

            "clusters_at_least_5pct": ";".join(
                clusters_at_5pct
            ),

            "n_clusters_at_least_10pct": len(
                clusters_at_10pct
            ),

            "clusters_at_least_10pct": ";".join(
                clusters_at_10pct
            ),

            "highest_detection_cluster": (
                highest_cluster
            ),

            "highest_percent_expressing": (
                highest_percent
            ),

            "highest_n_expressing_cells": (
                highest_n_expressing
            ),
        }
    )

marker_cluster_summary_df = pd.DataFrame(
    marker_cluster_summary_records
)


# ============================================================
# Marker diagnostic 10. Threshold summary
# ============================================================

marker_threshold_records: list[dict[str, Any]] = []

for (
    marker_group,
    cell_type_value,
), group_df in marker_celltype_annotated_df.groupby(
    [
        "marker_group",
        cell_type_column,
    ],
    observed=True,
    sort=True,
):
    n_markers_requested = int(
        group_df[
            "marker_gene_normalized"
        ].nunique()
    )

    n_markers_mapped = int(
        group_df.loc[
            group_df[
                "marker_gene_index"
            ].notna(),
            "marker_gene_normalized",
        ]
        .nunique()
    )

    n_markers_detected = int(
        group_df.loc[
            group_df[
                "n_expressing_cells"
            ].fillna(0)
            > 0,
            "marker_gene_normalized",
        ]
        .nunique()
    )

    for threshold in detection_thresholds:
        passing_genes = (
            group_df.loc[
                group_df[
                    "fraction_expressing"
                ].fillna(0.0)
                >= threshold,
                "marker_gene_normalized",
            ]
            .drop_duplicates()
            .tolist()
        )

        marker_threshold_records.append(
            {
                "marker_group": marker_group,
                cell_type_column: cell_type_value,
                "threshold": threshold,
                "threshold_percent": (
                    100.0 * threshold
                ),
                "n_markers_requested": (
                    n_markers_requested
                ),
                "n_markers_mapped": (
                    n_markers_mapped
                ),
                "n_markers_detected_any_cell": (
                    n_markers_detected
                ),
                "n_markers_passing_threshold": (
                    len(passing_genes)
                ),
                "fraction_requested_passing": (
                    len(passing_genes)
                    / n_markers_requested
                    if n_markers_requested > 0
                    else np.nan
                ),
                "passing_marker_genes": ";".join(
                    passing_genes
                ),
            }
        )

marker_threshold_summary_df = pd.DataFrame(
    marker_threshold_records
)


# ============================================================
# Marker diagnostic 11. Detection heatmap
# ============================================================

heatmap_height = max(
    10.0,
    0.25 * len(marker_percent_matrix_df),
)

heatmap_width = max(
    10.0,
    0.9 * len(marker_percent_matrix_df.columns),
)

plt.figure(
    figsize=(
        heatmap_width,
        heatmap_height,
    )
)

sns.heatmap(
    marker_percent_matrix_df,
    cmap="viridis",
    vmin=0.0,
    vmax=min(
        100.0,
        max(
            20.0,
            float(
                marker_percent_matrix_df
                .to_numpy()
                .max()
            ),
        ),
    ),
    linewidths=0.15,
    linecolor="white",
    cbar_kws={
        "label": "Percent of cells expressing marker",
    },
)

plt.title(
    "Marker-gene RNA detection by cell_type_broad"
)

plt.xlabel(
    "cell_type_broad"
)

plt.ylabel(
    "Marker gene"
)

plt.xticks(
    rotation=45,
    ha="right",
)

plt.tight_layout()

plt.savefig(
    marker_detection_heatmap_path,
    dpi=300,
    bbox_inches="tight",
)

plt.close()

# ============================================================
# Marker diagnostic 12. Binary presence heatmap
# ============================================================

plt.figure(
    figsize=(
        heatmap_width,
        heatmap_height,
    )
)

sns.heatmap(
    marker_presence_matrix_df.astype(int),
    cmap=[
        "#f2f2f2",
        "#2166ac",
    ],
    vmin=0,
    vmax=1,
    linewidths=0.15,
    linecolor="white",
    cbar_kws={
        "ticks": [0, 1],
        "label": "Detected in at least one cell",
    },
)

plt.title(
    "Marker-gene presence by cell_type_broad"
)

plt.xlabel(
    "cell_type_broad"
)

plt.ylabel(
    "Marker gene"
)

plt.xticks(
    rotation=45,
    ha="right",
)

plt.tight_layout()

plt.savefig(
    marker_presence_heatmap_path,
    dpi=300,
    bbox_inches="tight",
)

plt.close()

# ============================================================
# Marker diagnostic 13. Write outputs
# ============================================================

marker_mapping.to_csv(
    marker_mapping_path,
    sep="\t",
    index=False,
)

marker_global_annotated_df.to_csv(
    marker_global_metrics_path,
    sep="\t",
    index=False,
)

marker_celltype_annotated_df.to_csv(
    marker_celltype_metrics_path,
    sep="\t",
    index=False,
)

marker_presence_matrix_df.astype(int).to_csv(
    marker_presence_matrix_path,
    sep="\t",
    index=True,
)

marker_percent_matrix_df.to_csv(
    marker_detection_matrix_path,
    sep="\t",
    index=True,
)

marker_cluster_summary_df.to_csv(
    marker_cluster_summary_path,
    sep="\t",
    index=False,
)

marker_threshold_summary_df.to_csv(
    marker_threshold_summary_path,
    sep="\t",
    index=False,
)
# ============================================================
# Marker diagnostic 14. Print important audit findings
# ============================================================

absent_marker_rows = (
    marker_cluster_summary_df.loc[
        ~marker_cluster_summary_df[
            "present_in_rna_gene_universe"
        ]
    ]
)

undetected_marker_rows = (
    marker_cluster_summary_df.loc[
        marker_cluster_summary_df[
            "present_in_rna_gene_universe"
        ]
        & ~marker_cluster_summary_df[
            "detected_in_any_cell"
        ]
    ]
)

print("\nMarker-gene audit summary:")
print(
    "- Requested marker-group entries: "
    f"{len(marker_catalogue):,}"
)
print(
    "- Unique requested genes: "
    f"{marker_catalogue['marker_gene_normalized'].nunique():,}"
)
print(
    "- Unique mapped genes: "
    f"{len(unique_marker_universe):,}"
)
print(
    "- Marker-group entries absent from RNA universe: "
    f"{len(absent_marker_rows):,}"
)
print(
    "- Mapped marker-group entries with zero detection: "
    f"{len(undetected_marker_rows):,}"
)

if not absent_marker_rows.empty:
    print("\nMarkers absent from RNA gene universe:")
    print(
        absent_marker_rows[
            [
                "marker_group",
                "marker_gene",
                "mapping_status",
            ]
        ].to_string(
            index=False
        )
    )

if not undetected_marker_rows.empty:
    print("\nMapped markers not detected in any cell:")
    print(
        undetected_marker_rows[
            [
                "marker_group",
                "marker_gene",
            ]
        ].to_string(
            index=False
        )
    )

print("\nMarker detection outputs:")
print(f"- Mapping: {marker_mapping_path}")
print(
    "- Cell-type metrics: "
    f"{marker_celltype_metrics_path}"
)
print(
    "- Cluster summary: "
    f"{marker_cluster_summary_path}"
)
print(
    "- Percent matrix: "
    f"{marker_detection_matrix_path}"
)
print(
    "- Detection heatmap: "
    f"{marker_detection_heatmap_path}"
)



