# ============================================================
# 05A_paired_build_full_gene_rna_final_cells.py
#
# Purpose
# -------
# Create the official full-gene RNA object for downstream
# PBMC Multiome GRN analysis.
#
# The output contains:
# - all RNA genes from the original full-gene RNA object,
# - only the 7,566 cells approved after RNA + ATAC QC,
# - cells in exactly the same order as the final RNA, ATAC
#   and Multiome objects,
# - authoritative final RNA metadata,
# - all raw-RNA metadata with the prefix "raw_rna__".
#
# This script does NOT:
# - download or filter transcription factors,
# - apply gene detection thresholds,
# - define promoters or TSS coordinates,
# - use GENCODE or JASPAR,
# - perform motif analysis,
# - construct an adjacency matrix,
# - train a graph neural network.
# ============================================================


# ============================================================
# 1. Imports
# ============================================================

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml


# ============================================================
# 2. Resolve refactor repository
# ============================================================

script_path = Path(__file__).resolve()

# Expected script location:
#
# Refactored/
# └── scripts/
#     └── paired_scripts/
#         └── 05A_paired_build_full_gene_rna_final_cells.py
#
# parents[0] -> paired_scripts
# parents[1] -> scripts
# parents[2] -> Refactored

refactor_dir = script_path.parents[2]

print("\nScript path:")
print(script_path)

print("\nRefactor repository:")
print(refactor_dir)


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
# 5. Resolve official input paths
# ============================================================

processed_dir = (
    paths.project_root
    / "data"
    / "processed"
)

final_filtered_dir = (
    processed_dir
    / "final_after_atac_qc"
)

rna_raw_dir = (
    processed_dir
    / "rna"
)

rna_final_path = (
    final_filtered_dir
    / "pbmc_multiome_rna_final_after_atac_qc.h5ad"
)

atac_final_path = (
    final_filtered_dir
    / "pbmc_multiome_atac_final_after_atac_qc.h5ad"
)

multiome_final_path = (
    final_filtered_dir
    / "pbmc_multiome_final_after_atac_qc.h5ad"
)

rna_raw_path = (
    rna_raw_dir
    / "pbmc_multiome_rna_raw.h5ad"
)


# ============================================================
# 6. Expected final-object dimensions
# ============================================================

expected_n_cells = 7_566
expected_rna_hvg_n_vars = 3_000
expected_atac_n_vars = 143_887
expected_multiome_n_vars = 180_488


# ============================================================
# 7. Resolve output directories and paths
# ============================================================

grn_rna_output_dir = (
    processed_dir
    / "grn"
    / "rna"
)

output_table_dir = (
    paths.project_root
    / "reports"
    / "tables"
    / "grn"
    / "step05A"
)

output_summary_dir = (
    paths.project_root
    / "reports"
    / "summaries"
    / "grn"
    / "step05A"
)

for directory in (
    grn_rna_output_dir,
    output_table_dir,
    output_summary_dir,
):
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )


rna_grn_output_path = (
    grn_rna_output_dir
    / "pbmc_multiome_rna_full_genes_final_cells.h5ad"
)

input_audit_path = (
    output_table_dir
    / "05A_full_gene_rna_input_audit.tsv"
)

cell_sync_audit_path = (
    output_table_dir
    / "05A_cell_synchronization_audit.tsv"
)

metadata_audit_path = (
    output_table_dir
    / "05A_obs_metadata_audit.tsv"
)

feature_audit_path = (
    output_table_dir
    / "05A_full_gene_rna_feature_audit.tsv"
)

summary_yaml_path = (
    output_summary_dir
    / "05A_full_gene_rna_summary.yaml"
)

run_metadata_path = (
    output_summary_dir
    / "05A_run_metadata.json"
)


# ============================================================
# 8. Runtime configuration
# ============================================================

# False protects an existing official output from being
# overwritten accidentally.
overwrite_output = False

# Prefix applied to every raw-RNA obs column.
raw_obs_prefix = "raw_rna__"

# We expect raw counts in the full-gene RNA object.
# Integer-like validation is reported but is not a hard failure,
# because some pipelines store counts using floating-point dtype.
integer_tolerance = 1e-6


# ============================================================
# 9. Generic helper functions
# ============================================================

def utc_timestamp() -> str:
    """Return the current time as ISO-8601 UTC."""
    return datetime.now(timezone.utc).isoformat()


def validate_file(
    path: Path,
    label: str,
    allow_empty: bool = False,
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

    if not allow_empty and path.stat().st_size == 0:
        raise ValueError(
            f"{label} is empty:\n"
            f"{path}"
        )


def sha256_file(
    path: Path,
    chunk_size: int = 1024 * 1024,
) -> str:
    """Calculate a SHA-256 checksum without loading the file."""
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)

            if not chunk:
                break

            digest.update(chunk)

    return digest.hexdigest()


def close_backed_anndata(adata: ad.AnnData) -> None:
    """Close an AnnData backing file safely."""
    if adata.isbacked and adata.file is not None:
        adata.file.close()


def json_safe(value: Any) -> Any:
    """
    Recursively convert values into JSON/YAML-safe Python
    objects.

    Handles:
    - pathlib.Path,
    - NumPy scalar values,
    - NumPy arrays,
    - pandas Index and Series,
    - dictionaries,
    - lists and tuples,
    - scalar missing values.
    """

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

    # Only scalar values should reach pd.isna().
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        missing = False

    if isinstance(missing, (bool, np.bool_)) and missing:
        return None

    return value


def audit_anndata(
    adata: ad.AnnData,
    object_name: str,
    input_path: Path,
    expected_obs: int | None = None,
    expected_vars: int | None = None,
) -> dict[str, Any]:
    """Audit an AnnData object without modifying it."""
    observed_obs = int(adata.n_obs)
    observed_vars = int(adata.n_vars)

    if expected_obs is not None and observed_obs != expected_obs:
        raise ValueError(
            f"{object_name}: unexpected number of cells.\n"
            f"Expected: {expected_obs:,}\n"
            f"Observed: {observed_obs:,}"
        )

    if expected_vars is not None and observed_vars != expected_vars:
        raise ValueError(
            f"{object_name}: unexpected number of features.\n"
            f"Expected: {expected_vars:,}\n"
            f"Observed: {observed_vars:,}"
        )

    if not adata.obs_names.is_unique:
        raise ValueError(
            f"{object_name}: cell identifiers are not unique."
        )

    if not adata.var_names.is_unique:
        raise ValueError(
            f"{object_name}: feature identifiers are not unique."
        )

    return {
        "object": object_name,
        "path": str(input_path),
        "n_obs": observed_obs,
        "n_vars": observed_vars,
        "obs_names_unique": bool(adata.obs_names.is_unique),
        "var_names_unique": bool(adata.var_names.is_unique),
        "n_obs_columns": int(len(adata.obs.columns)),
        "n_var_columns": int(len(adata.var.columns)),
        "n_layers": int(len(adata.layers.keys())),
        "is_backed": bool(adata.isbacked),
        "file_size_bytes": int(input_path.stat().st_size),
        "sha256": sha256_file(input_path),
    }


def matrix_values(
    matrix: Any,
) -> np.ndarray:
    """
    Return stored matrix values.

    For sparse matrices, only non-zero stored values are
    returned. For dense matrices, all values are flattened.
    """
    if sp.issparse(matrix):
        return np.asarray(matrix.data)

    return np.asarray(matrix).ravel()


def audit_count_matrix(
    matrix: Any,
    matrix_name: str,
    integer_tolerance_value: float = 1e-6,
) -> dict[str, Any]:
    """
    Audit a candidate count matrix.

    Zero entries in a sparse matrix are implicit and are not
    treated as missing values.
    """
    values = matrix_values(matrix)

    is_sparse = bool(sp.issparse(matrix))
    matrix_format = (
        matrix.getformat()
        if is_sparse
        else "dense"
    )

    if values.size == 0:
        min_value = 0.0
        max_value = 0.0
        n_negative = 0
        n_nan = 0
        n_inf = 0
        integer_like = True
    else:
        finite_mask = np.isfinite(values)

        n_nan = int(np.isnan(values).sum())
        n_inf = int(np.isinf(values).sum())
        n_negative = int((values < 0).sum())

        finite_values = values[finite_mask]

        if finite_values.size == 0:
            min_value = None
            max_value = None
            integer_like = False
        else:
            min_value = float(finite_values.min())
            max_value = float(finite_values.max())

            integer_like = bool(
                np.all(
                    np.abs(
                        finite_values
                        - np.round(finite_values)
                    )
                    <= integer_tolerance_value
                )
            )

    if n_nan > 0:
        raise ValueError(
            f"{matrix_name} contains {n_nan:,} NaN values."
        )

    if n_inf > 0:
        raise ValueError(
            f"{matrix_name} contains {n_inf:,} infinite values."
        )

    if n_negative > 0:
        raise ValueError(
            f"{matrix_name} contains {n_negative:,} negative values."
        )

    return {
        "matrix_name": matrix_name,
        "shape": list(matrix.shape),
        "dtype": str(matrix.dtype),
        "is_sparse": is_sparse,
        "matrix_format": matrix_format,
        "n_stored_values": int(values.size),
        "min_stored_value": min_value,
        "max_stored_value": max_value,
        "contains_nan": bool(n_nan > 0),
        "contains_inf": bool(n_inf > 0),
        "contains_negative": bool(n_negative > 0),
        "integer_like": integer_like,
    }


def series_values_equal(
    left: pd.Series,
    right: pd.Series,
) -> pd.Series:
    """
    Compare two aligned Series while treating missing values
    in both Series as equal.
    """
    left_object = left.astype("object")
    right_object = right.astype("object")

    both_missing = (
        pd.isna(left_object)
        & pd.isna(right_object)
    )

    equal_non_missing = (
        left_object.eq(right_object)
        .fillna(False)
    )

    return (
        both_missing
        | equal_non_missing
    )


# ============================================================
# 10. Validate required input files
# ============================================================

required_inputs = {
    "Final RNA H5AD": rna_final_path,
    "Final ATAC H5AD": atac_final_path,
    "Final Multiome H5AD": multiome_final_path,
    "Full-gene raw RNA H5AD": rna_raw_path,
}

print("\nValidating required project inputs:")

for label, path in required_inputs.items():
    validate_file(
        path=path,
        label=label,
    )

    print(
        f"- {label}: OK "
        f"({path.stat().st_size:,} bytes)"
    )


# ============================================================
# 11. Protect existing official output
# ============================================================

if rna_grn_output_path.exists():
    if not overwrite_output:
        raise FileExistsError(
            "The official GRN RNA output already exists and "
            "overwrite_output=False.\n"
            f"Output path:\n{rna_grn_output_path}\n\n"
            "Inspect or remove the existing file before rerunning, "
            "or set overwrite_output=True intentionally."
        )

    if not rna_grn_output_path.is_file():
        raise ValueError(
            "The output path exists but is not a regular file:\n"
            f"{rna_grn_output_path}"
        )

    print("\nExisting output will be overwritten:")
    print(rna_grn_output_path)


# ============================================================
# 12. Load AnnData objects in backed mode
# ============================================================

rna_final = ad.read_h5ad(
    rna_final_path,
    backed="r",
)

atac_final = ad.read_h5ad(
    atac_final_path,
    backed="r",
)

multiome_final = ad.read_h5ad(
    multiome_final_path,
    backed="r",
)

rna_raw = ad.read_h5ad(
    rna_raw_path,
    backed="r",
)


# ============================================================
# 13. Audit input objects
# ============================================================

input_audit_records = [
    audit_anndata(
        adata=rna_final,
        object_name="RNA_final_HVG",
        input_path=rna_final_path,
        expected_obs=expected_n_cells,
        expected_vars=expected_rna_hvg_n_vars,
    ),
    audit_anndata(
        adata=atac_final,
        object_name="ATAC_final",
        input_path=atac_final_path,
        expected_obs=expected_n_cells,
        expected_vars=expected_atac_n_vars,
    ),
    audit_anndata(
        adata=multiome_final,
        object_name="Multiome_final",
        input_path=multiome_final_path,
        expected_obs=expected_n_cells,
        expected_vars=expected_multiome_n_vars,
    ),
    audit_anndata(
        adata=rna_raw,
        object_name="RNA_raw_full_gene",
        input_path=rna_raw_path,
        expected_obs=None,
        expected_vars=None,
    ),
]

input_audit_df = pd.DataFrame(
    input_audit_records
)

print("\nInput object audit:")
print(
    input_audit_df.to_string(
        index=False
    )
)


# ============================================================
# 14. Define canonical final-cell order
# ============================================================

final_cells = pd.Index(
    rna_final.obs_names.astype(str),
    name="cell_barcode",
)

atac_cells = pd.Index(
    atac_final.obs_names.astype(str),
    name="cell_barcode",
)

multiome_cells = pd.Index(
    multiome_final.obs_names.astype(str),
    name="cell_barcode",
)

raw_rna_cells = pd.Index(
    rna_raw.obs_names.astype(str),
    name="cell_barcode",
)


# ============================================================
# 15. Validate synchronization of final objects
# ============================================================

rna_atac_same_membership = (
    len(final_cells.difference(atac_cells)) == 0
    and len(atac_cells.difference(final_cells)) == 0
)

rna_atac_same_order = final_cells.equals(
    atac_cells
)

rna_multiome_same_membership = (
    len(final_cells.difference(multiome_cells)) == 0
    and len(multiome_cells.difference(final_cells)) == 0
)

rna_multiome_same_order = final_cells.equals(
    multiome_cells
)

if not rna_atac_same_membership:
    raise ValueError(
        "Final RNA and final ATAC do not contain the same "
        "cell barcode set."
    )

if not rna_atac_same_order:
    raise ValueError(
        "Final RNA and final ATAC contain the same cells but "
        "not in exactly the same order."
    )

if not rna_multiome_same_membership:
    raise ValueError(
        "Final RNA and final Multiome do not contain the same "
        "cell barcode set."
    )

if not rna_multiome_same_order:
    raise ValueError(
        "Final RNA and final Multiome contain the same cells "
        "but not in exactly the same order."
    )

if not raw_rna_cells.is_unique:
    raise ValueError(
        "Full-gene raw RNA cell identifiers are not unique."
    )


# ============================================================
# 16. Validate final cells inside full-gene raw RNA
# ============================================================

missing_final_cells_in_raw = (
    final_cells.difference(raw_rna_cells)
)

extra_raw_cells = (
    raw_rna_cells.difference(final_cells)
)

if len(missing_final_cells_in_raw) > 0:
    raise ValueError(
        "Some approved final cells are absent from the "
        "full-gene raw RNA object.\n"
        f"Number missing: {len(missing_final_cells_in_raw):,}\n"
        f"Examples: {missing_final_cells_in_raw[:20].tolist()}"
    )

cell_sync_records = [
    {
        "comparison": "RNA_final_vs_ATAC_final",
        "n_reference_cells": int(len(final_cells)),
        "n_compared_cells": int(len(atac_cells)),
        "same_membership": rna_atac_same_membership,
        "same_order": rna_atac_same_order,
        "n_missing_from_compared": int(
            len(final_cells.difference(atac_cells))
        ),
        "n_extra_in_compared": int(
            len(atac_cells.difference(final_cells))
        ),
    },
    {
        "comparison": "RNA_final_vs_Multiome_final",
        "n_reference_cells": int(len(final_cells)),
        "n_compared_cells": int(len(multiome_cells)),
        "same_membership": rna_multiome_same_membership,
        "same_order": rna_multiome_same_order,
        "n_missing_from_compared": int(
            len(final_cells.difference(multiome_cells))
        ),
        "n_extra_in_compared": int(
            len(multiome_cells.difference(final_cells))
        ),
    },
    {
        "comparison": "RNA_final_vs_RNA_raw_full_gene",
        "n_reference_cells": int(len(final_cells)),
        "n_compared_cells": int(len(raw_rna_cells)),
        "same_membership": bool(
            len(missing_final_cells_in_raw) == 0
            and len(extra_raw_cells) == 0
        ),
        "same_order": final_cells.equals(raw_rna_cells),
        "n_missing_from_compared": int(
            len(missing_final_cells_in_raw)
        ),
        "n_extra_in_compared": int(
            len(extra_raw_cells)
        ),
    },
]

cell_sync_df = pd.DataFrame(
    cell_sync_records
)

print("\nCell synchronization audit:")
print(
    cell_sync_df.to_string(
        index=False
    )
)

print("\nFinal approved cells:")
print(f"- Number of final cells: {len(final_cells):,}")
print("- All final cells are present in raw RNA: True")
print(
    "- Additional pre-QC cells in raw RNA: "
    f"{len(extra_raw_cells):,}"
)


# ============================================================
# 17. Resolve final-cell positions inside raw RNA
# ============================================================

raw_cell_positions = pd.Series(
    data=np.arange(
        rna_raw.n_obs,
        dtype=np.int64,
    ),
    index=raw_rna_cells,
    dtype="int64",
)

final_positions_in_raw = (
    raw_cell_positions
    .loc[final_cells]
    .to_numpy(dtype=np.int64)
)

if len(final_positions_in_raw) != expected_n_cells:
    raise RuntimeError(
        "Unexpected number of final-cell positions resolved "
        "inside the raw RNA object."
    )

if len(np.unique(final_positions_in_raw)) != expected_n_cells:
    raise RuntimeError(
        "Resolved raw-RNA cell positions are not unique."
    )


# ============================================================
# 18. Load final-cell full-gene RNA subset into memory
# ============================================================

print("\nCreating full-gene RNA subset:")
print(
    f"- Raw RNA input shape: "
    f"{rna_raw.n_obs:,} × {rna_raw.n_vars:,}"
)
print(
    f"- Selecting approved cells: "
    f"{len(final_positions_in_raw):,}"
)
print("- Retaining all RNA genes.")
print("- Applying canonical final RNA cell order.")

# Integer positional indexing is used because the raw object is
# backed. The position vector is already ordered according to
# final_cells.
rna_grn = (
    rna_raw[
        final_positions_in_raw,
        :,
    ]
    .to_memory()
)

if rna_grn.is_view:
    rna_grn = rna_grn.copy()


# ============================================================
# 19. Validate subset before metadata replacement
# ============================================================

if rna_grn.n_obs != expected_n_cells:
    raise ValueError(
        "The subset RNA object contains an unexpected number "
        "of cells.\n"
        f"Expected: {expected_n_cells:,}\n"
        f"Observed: {rna_grn.n_obs:,}"
    )

if rna_grn.n_vars != rna_raw.n_vars:
    raise ValueError(
        "The subset RNA object did not retain all RNA genes.\n"
        f"Raw RNA genes: {rna_raw.n_vars:,}\n"
        f"Subset RNA genes: {rna_grn.n_vars:,}"
    )

if not pd.Index(
    rna_grn.obs_names.astype(str)
).equals(final_cells):
    raise ValueError(
        "The subset RNA cell order does not match the "
        "canonical final-cell order."
    )

if not pd.Index(
    rna_grn.var_names.astype(str)
).equals(
    pd.Index(rna_raw.var_names.astype(str))
):
    raise ValueError(
        "The subset RNA gene order differs from the full-gene "
        "raw RNA object."
    )


# ============================================================
# 20. Prepare final and raw RNA metadata
# ============================================================

final_obs = (
    rna_final.obs
    .copy()
)

final_obs.index = final_obs.index.astype(str)
final_obs = final_obs.loc[final_cells].copy()

raw_obs_selected = (
    rna_raw.obs
    .iloc[final_positions_in_raw]
    .copy()
)

raw_obs_selected.index = (
    raw_obs_selected.index.astype(str)
)

if not pd.Index(
    raw_obs_selected.index
).equals(final_cells):
    raise ValueError(
        "Selected raw-RNA metadata is not aligned to the "
        "canonical final-cell order."
    )


# ============================================================
# 21. Audit shared metadata columns
# ============================================================

shared_obs_columns = sorted(
    set(final_obs.columns)
    .intersection(raw_obs_selected.columns)
)

metadata_audit_records: list[dict[str, Any]] = []

for column in shared_obs_columns:
    equality_mask = series_values_equal(
        final_obs[column],
        raw_obs_selected[column],
    )

    n_equal = int(equality_mask.sum())
    n_different = int((~equality_mask).sum())

    metadata_audit_records.append(
        {
            "column": column,
            "present_in_final_obs": True,
            "present_in_raw_obs": True,
            "n_cells": int(expected_n_cells),
            "n_equal_values": n_equal,
            "n_different_values": n_different,
            "fraction_equal": float(
                n_equal / expected_n_cells
            ),
            "authoritative_source": "final_rna_obs",
            "raw_copy_output_column": (
                f"{raw_obs_prefix}{column}"
            ),
        }
    )

final_only_columns = sorted(
    set(final_obs.columns)
    .difference(raw_obs_selected.columns)
)

for column in final_only_columns:
    metadata_audit_records.append(
        {
            "column": column,
            "present_in_final_obs": True,
            "present_in_raw_obs": False,
            "n_cells": int(expected_n_cells),
            "n_equal_values": pd.NA,
            "n_different_values": pd.NA,
            "fraction_equal": pd.NA,
            "authoritative_source": "final_rna_obs",
            "raw_copy_output_column": pd.NA,
        }
    )

raw_only_columns = sorted(
    set(raw_obs_selected.columns)
    .difference(final_obs.columns)
)

for column in raw_only_columns:
    metadata_audit_records.append(
        {
            "column": column,
            "present_in_final_obs": False,
            "present_in_raw_obs": True,
            "n_cells": int(expected_n_cells),
            "n_equal_values": pd.NA,
            "n_different_values": pd.NA,
            "fraction_equal": pd.NA,
            "authoritative_source": "raw_rna_obs_only",
            "raw_copy_output_column": (
                f"{raw_obs_prefix}{column}"
            ),
        }
    )

metadata_audit_df = pd.DataFrame(
    metadata_audit_records
)

if metadata_audit_df.empty:
    metadata_audit_df = pd.DataFrame(
        columns=[
            "column",
            "present_in_final_obs",
            "present_in_raw_obs",
            "n_cells",
            "n_equal_values",
            "n_different_values",
            "fraction_equal",
            "authoritative_source",
            "raw_copy_output_column",
        ]
    )

print("\nMetadata audit summary:")
print(
    f"- Final RNA obs columns: {len(final_obs.columns):,}"
)
print(
    f"- Raw RNA obs columns: "
    f"{len(raw_obs_selected.columns):,}"
)
print(
    f"- Shared columns: {len(shared_obs_columns):,}"
)
print(
    f"- Final-only columns: {len(final_only_columns):,}"
)
print(
    f"- Raw-only columns: {len(raw_only_columns):,}"
)

shared_columns_with_differences = (
    metadata_audit_df.loc[
        metadata_audit_df[
            "n_different_values"
        ]
        .fillna(0)
        .astype(int)
        > 0,
        "column",
    ]
    .astype(str)
    .tolist()
)

print(
    "- Shared columns with at least one differing value: "
    f"{len(shared_columns_with_differences):,}"
)

if shared_columns_with_differences:
    print("- Examples:")
    for column in shared_columns_with_differences[:20]:
        print(f"  - {column}")


# ============================================================
# 22. Build combined output metadata
# ============================================================

# Every raw-RNA metadata column is retained with a prefix.
# This avoids collisions and makes provenance explicit.
raw_obs_prefixed = raw_obs_selected.copy()

raw_obs_prefixed.columns = [
    f"{raw_obs_prefix}{column}"
    for column in raw_obs_prefixed.columns
]

duplicated_output_columns = (
    set(final_obs.columns)
    .intersection(raw_obs_prefixed.columns)
)

if duplicated_output_columns:
    raise ValueError(
        "Metadata column collision after applying raw-RNA "
        "prefix.\n"
        f"Colliding columns: "
        f"{sorted(duplicated_output_columns)[:20]}"
    )

combined_obs = pd.concat(
    [
        final_obs,
        raw_obs_prefixed,
    ],
    axis=1,
)

combined_obs.index = final_cells

if combined_obs.columns.duplicated().any():
    duplicated_columns = (
        combined_obs.columns[
            combined_obs.columns.duplicated(
                keep=False
            )
        ]
        .tolist()
    )

    raise ValueError(
        "The combined output metadata contains duplicated "
        "column names.\n"
        f"Examples: {duplicated_columns[:20]}"
    )

if not pd.Index(
    combined_obs.index.astype(str)
).equals(final_cells):
    raise ValueError(
        "Combined output metadata is not aligned to the "
        "canonical final-cell order."
    )

rna_grn.obs = combined_obs


# ============================================================
# 23. Check required cell-type annotation
# ============================================================

cell_type_column = "cell_type_broad"

cell_type_column_present = (
    cell_type_column
    in rna_grn.obs.columns
)

if not cell_type_column_present:
    print(
        "\nWARNING:"
        f"\n- {cell_type_column!r} is not present in final RNA obs."
        "\n- 05A can still create the RNA object."
        "\n- 05A.1 cell-type-aware TF filtering will require this "
        "column or an explicitly configured replacement."
    )

    n_missing_cell_type = None
    fraction_missing_cell_type = None
    n_cell_types = None
else:
    n_missing_cell_type = int(
        rna_grn.obs[
            cell_type_column
        ]
        .isna()
        .sum()
    )

    fraction_missing_cell_type = float(
        rna_grn.obs[
            cell_type_column
        ]
        .isna()
        .mean()
    )

    n_cell_types = int(
        rna_grn.obs[
            cell_type_column
        ]
        .dropna()
        .nunique()
    )

    print("\nCell-type annotation audit:")
    print(
        f"- Column: {cell_type_column}"
    )
    print(
        f"- Missing values: "
        f"{n_missing_cell_type:,}"
    )
    print(
        f"- Missing fraction: "
        f"{fraction_missing_cell_type:.4%}"
    )
    print(
        f"- Non-missing broad cell types: "
        f"{n_cell_types:,}"
    )


# ============================================================
# 24. Determine and audit counts source
# ============================================================

if "counts" in rna_grn.layers:
    counts_source = "layers/counts"
    counts_matrix = rna_grn.layers["counts"]
else:
    counts_source = "X"
    counts_matrix = rna_grn.X

counts_audit = audit_count_matrix(
    matrix=counts_matrix,
    matrix_name=counts_source,
    integer_tolerance_value=integer_tolerance,
)

print("\nCount-matrix audit:")
for key, value in counts_audit.items():
    print(f"- {key}: {value}")

if not counts_audit["integer_like"]:
    print(
        "\nWARNING:"
        "\n- The selected counts source is non-negative and finite,"
        "\n  but its stored values are not integer-like."
        "\n- Confirm whether the full-gene RNA object contains raw"
        "\n  UMI counts or normalized expression."
        "\n- Detection based on >0 may still be valid, but total"
        "\n  counts would not represent raw UMI totals."
    )


# ============================================================
# 25. Audit full RNA feature universe
# ============================================================

candidate_gene_symbol_columns = [
    "gene_symbol",
    "gene_symbols",
    "feature_name",
    "feature_names",
    "symbol",
    "gene_name",
    "gene_names",
]

observed_gene_symbol_columns = [
    column
    for column in candidate_gene_symbol_columns
    if column in rna_grn.var.columns
]

var_names_as_string = pd.Index(
    rna_grn.var_names.astype(str)
)

ensembl_like_fraction = float(
    pd.Series(
        var_names_as_string,
        dtype="string",
    )
    .str.match(
        r"^ENS[A-Z]*G\d+(?:\.\d+)?$",
        na=False,
    )
    .mean()
)

if ensembl_like_fraction >= 0.80:
    inferred_var_identifier_type = "ensembl_gene_id"
elif ensembl_like_fraction <= 0.05:
    inferred_var_identifier_type = "likely_gene_symbol_or_other"
else:
    inferred_var_identifier_type = "mixed_or_uncertain"

selected_symbol_column = (
    observed_gene_symbol_columns[0]
    if observed_gene_symbol_columns
    else None
)

if selected_symbol_column is not None:
    symbol_series = (
        rna_grn.var[
            selected_symbol_column
        ]
        .astype("string")
        .str.strip()
    )

    n_missing_gene_symbols = int(
        symbol_series.isna().sum()
        + symbol_series.eq("").sum()
    )

    non_missing_symbols = symbol_series[
        symbol_series.notna()
        & symbol_series.ne("")
    ]

    n_duplicated_gene_symbols = int(
        non_missing_symbols
        .duplicated()
        .sum()
    )
else:
    n_missing_gene_symbols = None
    n_duplicated_gene_symbols = None

feature_audit_record = {
    "n_genes": int(rna_grn.n_vars),
    "var_names_unique": bool(
        rna_grn.var_names.is_unique
    ),
    "n_var_columns": int(
        len(rna_grn.var.columns)
    ),
    "var_columns": "|".join(
        map(str, rna_grn.var.columns)
    ),
    "inferred_var_identifier_type": (
        inferred_var_identifier_type
    ),
    "ensembl_like_fraction": (
        ensembl_like_fraction
    ),
    "candidate_symbol_columns_found": "|".join(
        observed_gene_symbol_columns
    ),
    "selected_symbol_column_for_audit": (
        selected_symbol_column
    ),
    "n_missing_gene_symbols": (
        n_missing_gene_symbols
    ),
    "n_additional_duplicated_gene_symbols": (
        n_duplicated_gene_symbols
    ),
    "counts_source": counts_source,
    "counts_dtype": counts_audit["dtype"],
    "counts_is_sparse": counts_audit["is_sparse"],
    "counts_matrix_format": counts_audit[
        "matrix_format"
    ],
    "counts_integer_like": counts_audit[
        "integer_like"
    ],
    "counts_min_stored_value": counts_audit[
        "min_stored_value"
    ],
    "counts_max_stored_value": counts_audit[
        "max_stored_value"
    ],
}

feature_audit_df = pd.DataFrame(
    [feature_audit_record]
)

print("\nFull-gene RNA feature audit:")
print(
    feature_audit_df.to_string(
        index=False
    )
)


# ============================================================
# 26. Add provenance to output AnnData
# ============================================================

rna_grn.uns["grn_rna_provenance"] = {
    "created_at_utc": utc_timestamp(),
    "script_name": script_path.name,
    "project_id": paths.project_id,
    "purpose": (
        "Official full-gene RNA object restricted to the "
        "final cells approved after RNA and ATAC QC."
    ),
    "canonical_cell_order_source": str(
        rna_final_path
    ),
    "expression_and_gene_universe_source": str(
        rna_raw_path
    ),
    "final_atac_validation_source": str(
        atac_final_path
    ),
    "final_multiome_validation_source": str(
        multiome_final_path
    ),
    "n_final_cells": int(
        rna_grn.n_obs
    ),
    "n_full_rna_genes": int(
        rna_grn.n_vars
    ),
    "n_pre_qc_raw_rna_cells_excluded": int(
        len(extra_raw_cells)
    ),
    "counts_source": counts_source,
    "counts_integer_like": bool(
        counts_audit["integer_like"]
    ),
    "final_obs_is_authoritative": True,
    "raw_obs_prefix": raw_obs_prefix,
    "n_final_obs_columns": int(
        len(final_obs.columns)
    ),
    "n_raw_obs_columns_retained": int(
        len(raw_obs_selected.columns)
    ),
    "n_combined_obs_columns": int(
        len(combined_obs.columns)
    ),
    "cell_type_column": (
        cell_type_column
        if cell_type_column_present
        else None
    ),
    "cell_type_missing_count": (
        n_missing_cell_type
    ),
    "important_scope_note": (
        "No transcription-factor filtering, gene-detection "
        "filtering, genomic annotation, promoter definition, "
        "motif analysis, adjacency construction or GNN "
        "training was performed in this script."
    ),
}


# ============================================================
# 27. Final in-memory validations
# ============================================================

if rna_grn.n_obs != expected_n_cells:
    raise RuntimeError(
        "Final in-memory output has an incorrect cell count."
    )

if rna_grn.n_vars != rna_raw.n_vars:
    raise RuntimeError(
        "Final in-memory output has an incorrect gene count."
    )

if not pd.Index(
    rna_grn.obs_names.astype(str)
).equals(final_cells):
    raise RuntimeError(
        "Final in-memory output has an incorrect cell order."
    )

if not pd.Index(
    rna_grn.var_names.astype(str)
).equals(
    pd.Index(rna_raw.var_names.astype(str))
):
    raise RuntimeError(
        "Final in-memory output has an incorrect gene order."
    )

if not rna_grn.obs_names.is_unique:
    raise RuntimeError(
        "Final in-memory output has duplicated cell names."
    )

if not rna_grn.var_names.is_unique:
    raise RuntimeError(
        "Final in-memory output has duplicated feature names."
    )


# ============================================================
# 28. Write official H5AD output
# ============================================================

print("\nWriting official GRN RNA object:")
print(rna_grn_output_path)

temporary_output_path = (
    rna_grn_output_path.with_suffix(
        rna_grn_output_path.suffix + ".tmp"
    )
)

if temporary_output_path.exists():
    temporary_output_path.unlink()

rna_grn.write_h5ad(
    temporary_output_path,
    compression="gzip",
)

validate_file(
    path=temporary_output_path,
    label="Temporary full-gene final-cell RNA H5AD",
)

temporary_output_sha256 = sha256_file(
    temporary_output_path
)

print("\nTemporary output SHA-256:")
print(temporary_output_sha256)

# ============================================================
# 29. Reopen and validate serialized H5AD
# ============================================================

print("\nReopening written H5AD for serialization validation.")

rna_grn_check = ad.read_h5ad(
    temporary_output_path,
    backed="r",
)

try:
    if rna_grn_check.n_obs != expected_n_cells:
        raise RuntimeError(
            "Serialized output contains an incorrect number "
            "of cells."
        )

    if rna_grn_check.n_vars != rna_raw.n_vars:
        raise RuntimeError(
            "Serialized output contains an incorrect number "
            "of genes."
        )

    if not pd.Index(
        rna_grn_check.obs_names.astype(str)
    ).equals(final_cells):
        raise RuntimeError(
            "Serialized output cell order differs from the "
            "canonical final-cell order."
        )

    if not pd.Index(
        rna_grn_check.var_names.astype(str)
    ).equals(
        pd.Index(rna_raw.var_names.astype(str))
    ):
        raise RuntimeError(
            "Serialized output gene order differs from the "
            "full-gene raw RNA object."
        )

    if list(rna_grn_check.obs.columns) != list(
        combined_obs.columns
    ):
        raise RuntimeError(
            "Serialized output metadata columns differ from "
            "the constructed metadata schema."
        )

    if (
        "grn_rna_provenance"
        not in rna_grn_check.uns
    ):
        raise RuntimeError(
            "Serialized output is missing "
            "uns['grn_rna_provenance']."
        )

finally:
    close_backed_anndata(
        rna_grn_check
    )


# ============================================================
# 30. Save audit tables
# ============================================================

input_audit_df.to_csv(
    input_audit_path,
    sep="\t",
    index=False,
)

cell_sync_df.to_csv(
    cell_sync_audit_path,
    sep="\t",
    index=False,
)

metadata_audit_df.to_csv(
    metadata_audit_path,
    sep="\t",
    index=False,
)

feature_audit_df.to_csv(
    feature_audit_path,
    sep="\t",
    index=False,
)


# ============================================================
# 31. Save YAML summary
# ============================================================

summary = {
    "step": "05A",
    "script": script_path.name,
    "created_at_utc": utc_timestamp(),
    "project_id": paths.project_id,
    "objective": (
        "Create the official full-gene RNA object containing "
        "only final cells approved after RNA and ATAC QC."
    ),
    "inputs": {
        "rna_final_hvg": str(
            rna_final_path
        ),
        "atac_final": str(
            atac_final_path
        ),
        "multiome_final": str(
            multiome_final_path
        ),
        "rna_raw_full_gene": str(
            rna_raw_path
        ),
    },
    "input_checksums": {
        "rna_final_hvg_sha256": sha256_file(
            rna_final_path
        ),
        "atac_final_sha256": sha256_file(
            atac_final_path
        ),
        "multiome_final_sha256": sha256_file(
            multiome_final_path
        ),
        "rna_raw_full_gene_sha256": sha256_file(
            rna_raw_path
        ),
    },
    "cell_universe": {
        "canonical_source": str(
            rna_final_path
        ),
        "n_final_cells": int(
            len(final_cells)
        ),
        "n_raw_rna_cells": int(
            len(raw_rna_cells)
        ),
        "n_extra_pre_qc_raw_rna_cells": int(
            len(extra_raw_cells)
        ),
        "n_final_cells_missing_from_raw": int(
            len(missing_final_cells_in_raw)
        ),
        "rna_atac_same_membership": bool(
            rna_atac_same_membership
        ),
        "rna_atac_same_order": bool(
            rna_atac_same_order
        ),
        "rna_multiome_same_membership": bool(
            rna_multiome_same_membership
        ),
        "rna_multiome_same_order": bool(
            rna_multiome_same_order
        ),
    },
    "gene_universe": {
        "source": str(
            rna_raw_path
        ),
        "n_genes": int(
            rna_grn.n_vars
        ),
        "all_raw_rna_genes_retained": True,
        "gene_order_retained": True,
        "inferred_var_identifier_type": (
            inferred_var_identifier_type
        ),
        "candidate_symbol_columns": (
            observed_gene_symbol_columns
        ),
    },
    "metadata": {
        "final_rna_obs_authoritative": True,
        "raw_rna_obs_retained": True,
        "raw_rna_obs_prefix": raw_obs_prefix,
        "n_final_obs_columns": int(
            len(final_obs.columns)
        ),
        "n_raw_obs_columns": int(
            len(raw_obs_selected.columns)
        ),
        "n_combined_obs_columns": int(
            len(combined_obs.columns)
        ),
        "n_shared_columns": int(
            len(shared_obs_columns)
        ),
        "shared_columns_with_differences": (
            shared_columns_with_differences
        ),
    },
    "counts": {
        key: json_safe(value)
        for key, value in counts_audit.items()
    },
    "cell_type_annotation": {
        "column": (
            cell_type_column
            if cell_type_column_present
            else None
        ),
        "column_present": bool(
            cell_type_column_present
        ),
        "n_missing": (
            n_missing_cell_type
        ),
        "fraction_missing": (
            fraction_missing_cell_type
        ),
        "n_non_missing_categories": (
            n_cell_types
        ),
    },
    "output": {
        "rna_full_genes_final_cells_h5ad": str(
            rna_grn_output_path
        ),
        "shape": [
            int(rna_grn.n_obs),
            int(rna_grn.n_vars),
        ],
        "sha256": temporary_output_sha256,
    },
    "reports": {
        "input_audit": str(
            input_audit_path
        ),
        "cell_synchronization_audit": str(
            cell_sync_audit_path
        ),
        "metadata_audit": str(
            metadata_audit_path
        ),
        "feature_audit": str(
            feature_audit_path
        ),
    },
    "scope_exclusions": [
        "No TF catalogue processing",
        "No TF filtering",
        "No RNA detection filtering",
        "No GENCODE annotation",
        "No TSS or promoter construction",
        "No JASPAR motif processing",
        "No peak-to-gene mapping",
        "No adjacency matrix construction",
        "No GNN training",
    ],
}

with summary_yaml_path.open(
    "w",
    encoding="utf-8",
) as handle:
    yaml.safe_dump(
        summary,
        handle,
        sort_keys=False,
        allow_unicode=True,
    )


# ============================================================
# 32. Save run metadata
# ============================================================

run_metadata = {
    "created_at_utc": utc_timestamp(),
    "script": script_path.name,
    "script_path": str(script_path),
    "project_id": paths.project_id,
    "project_root": str(paths.project_root),
    "python_version": sys.version,
    "numpy_version": np.__version__,
    "pandas_version": pd.__version__,
    "scipy_version": __import__("scipy").__version__,
    "anndata_version": ad.__version__,
    "output_path": str(rna_grn_output_path),
    "output_sha256": temporary_output_sha256,
    "script_status": "completed",
}

with run_metadata_path.open(
    "w",
    encoding="utf-8",
) as handle:
    json.dump(
        run_metadata,
        handle,
        indent=2,
        ensure_ascii=False,
    )
# ============================================================
# Promote validated temporary H5AD to official output
# ============================================================

if rna_grn_output_path.exists():
    if overwrite_output:
        rna_grn_output_path.unlink()
    else:
        raise FileExistsError(
            "The official output appeared during execution:\n"
            f"{rna_grn_output_path}"
        )

temporary_output_path.replace(
    rna_grn_output_path
)

validate_file(
    path=rna_grn_output_path,
    label="Official full-gene final-cell RNA H5AD",
)

print("\nTemporary H5AD promoted to official output:")
print(rna_grn_output_path)

# ============================================================
# 33. Validate all written outputs
# ============================================================

required_outputs = {
    "Official GRN RNA H5AD": rna_grn_output_path,
    "Input audit": input_audit_path,
    "Cell synchronization audit": cell_sync_audit_path,
    "Metadata audit": metadata_audit_path,
    "Feature audit": feature_audit_path,
    "Summary YAML": summary_yaml_path,
    "Run metadata JSON": run_metadata_path,
}

print("\nValidating written outputs:")

for label, path in required_outputs.items():
    validate_file(
        path=path,
        label=label,
        allow_empty=False,
    )

    print(f"- {label}: OK")
    print(f"  {path}")


# ============================================================
# 34. Close backed input objects
# ============================================================

for backed_object in (
    rna_final,
    atac_final,
    multiome_final,
    rna_raw,
):
    close_backed_anndata(
        backed_object
    )


# ============================================================
# 35. Final report
# ============================================================

print("\n" + "=" * 70)
print("05A COMPLETED SUCCESSFULLY")
print("=" * 70)

print("\nOfficial output:")
print(rna_grn_output_path)

print("\nOutput shape:")
print(
    f"{rna_grn.n_obs:,} cells × "
    f"{rna_grn.n_vars:,} genes"
)

print("\nCell synchronization:")
print("- Final RNA vs ATAC: identical membership and order")
print("- Final RNA vs Multiome: identical membership and order")
print("- All approved cells found in full-gene raw RNA")

print("\nGene universe:")
print("- All genes from the full-gene raw RNA object retained")
print("- No HVG filtering applied")
print("- No TF filtering applied")

print("\nMetadata:")
print(
    f"- Final authoritative columns: "
    f"{len(final_obs.columns):,}"
)
print(
    f"- Raw RNA columns retained with prefix "
    f"{raw_obs_prefix!r}: "
    f"{len(raw_obs_selected.columns):,}"
)
print(
    f"- Total output obs columns: "
    f"{len(combined_obs.columns):,}"
)

print("\nCounts source:")
print(counts_source)

print("\nNext official input:")
print(
    "Use this H5AD as the only RNA expression input "
    "for 05A.1 TF feature preparation."
)