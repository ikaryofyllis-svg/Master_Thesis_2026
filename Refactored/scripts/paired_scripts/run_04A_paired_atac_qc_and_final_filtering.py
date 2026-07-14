# -*- coding: utf-8 -*-
"""
Step 04A — paired scATAC QC and synchronized final filtering.

Current stage:
- Locate the Refactored root.
- Validate the active project.
- Resolve selected project inputs.
- Resolve external 10x raw-data directory.
- Check that required paths exist.

No AnnData loading.
No filtering.
No TF-IDF, LSI, UMAP, or Leiden.
"""

from pathlib import Path
import sys


import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib.pyplot as plt
import seaborn as sns

# ============================================================
# 1. Find Refactored root
# ============================================================

def find_refactored_root(start_path: Path) -> Path:
    """
    Find the parent directory containing both src/ and projects/.

    This is safer than using a fixed parents[1] or parents[2],
    because it also works if the script is moved to another
    subdirectory inside Refactored.
    """

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
        "Expected a parent directory containing both:\n"
        "- src/\n"
        "- projects/"
    )


try:
    script_path = Path(__file__).resolve()
    refactor_dir = find_refactored_root(script_path.parent)
except NameError:
    # Interactive fallback, for example when running in a notebook
    script_path = Path("run_04A_atac_qc_and_final_filtering.py")
    refactor_dir = find_refactored_root(Path.cwd())


print("\nScript path:")
print(script_path)

print("\nRefactored root:")
print(refactor_dir)


# ============================================================
# 2. Enable local pbmcgrn imports
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
# 3. Validate active project
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
# 4. Resolve selected project inputs
# ============================================================

selected_dir = (
    paths.project_root
    / "data"
    / "processed"
    / "selected_for_scatac_analysis"
)

rna_selected_path = (
    selected_dir
    / "pbmc_multiome_rna_selected_for_scatac_analysis.h5ad"
)

atac_selected_path = (
    selected_dir
    / "pbmc_multiome_atac_selected_for_scatac_analysis.h5ad"
)

multiome_selected_path = (
    selected_dir
    / "pbmc_multiome_selected_for_scatac_analysis.h5ad"
)

selection_mask_path = (
    paths.project_root
    / "reports"
    / "tables"
    / "rna"
    / "step03R"
    / "selected_cells_for_scatac_analysis_mask.csv"
)


# ============================================================
# 5. Resolve external 10x raw-data directory
# ============================================================

# The original 10x files are outside Refactored and are treated
# as read-only source data.
thesis_root = refactor_dir.parent

external_raw_dir = (
    thesis_root
    / "data"
    / "raw"
    / "PBMC from a Healthy Donor - Granulocytes Removed Through Cell Sorting (10k)"
)


# ============================================================
# 6. Resolve future output directories
# ============================================================

final_output_dir = (
    paths.project_root
    / "data"
    / "processed"
    / "final_after_atac_qc"
)

tables_output_dir = (
    paths.project_root
    / "reports"
    / "tables"
    / "atac"
    / "step04A"
)

figures_output_dir = (
    paths.project_root
    / "reports"
    / "figures"
    / "atac"
    / "step04A"
)

manifest_output_dir = (
    paths.project_root
    / "reports"
    / "manifests"
)


# ============================================================
# 7. Validate required existing paths
# ============================================================

required_existing_paths = {
    "selected input directory": selected_dir,
    "selected RNA": rna_selected_path,
    "selected ATAC": atac_selected_path,
    "selected Multiome": multiome_selected_path,
    "selection mask": selection_mask_path,
    "external 10x raw directory": external_raw_dir,
}

missing_paths = []

print("\nRequired input paths:")

for label, path in required_existing_paths.items():
    exists = path.exists()

    print(f"\n{label}:")
    print(path)
    print(f"Exists: {exists}")

    if not exists:
        missing_paths.append((label, path))

if missing_paths:
    formatted_missing = "\n\n".join(
        f"{label}:\n{path}"
        for label, path in missing_paths
    )

    raise FileNotFoundError(
        "One or more required paths were not found:\n\n"
        f"{formatted_missing}"
    )


# ============================================================
# 8. Validate that outputs remain inside active project
# ============================================================

def assert_path_inside_project(
    path: Path,
    project_root: Path,
    label: str,
) -> None:
    """Fail if a planned output path is outside the active project."""

    resolved_path = path.resolve()
    resolved_project_root = project_root.resolve()

    try:
        resolved_path.relative_to(resolved_project_root)
    except ValueError as exc:
        raise ValueError(
            f"{label} is outside the active project.\n"
            f"Path: {resolved_path}\n"
            f"Project root: {resolved_project_root}"
        ) from exc


planned_output_dirs = {
    "final output directory": final_output_dir,
    "tables output directory": tables_output_dir,
    "figures output directory": figures_output_dir,
    "manifest output directory": manifest_output_dir,
}

print("\nPlanned output directories:")

for label, path in planned_output_dirs.items():
    assert_path_inside_project(
        path=path,
        project_root=paths.project_root,
        label=label,
    )

    print(f"\n{label}:")
    print(path)


# ============================================================
# 9. Final architecture check
# ============================================================

print("\nArchitecture validation completed successfully.")

print("\nRules confirmed:")
print("- External 10x raw data are read-only inputs.")
print("- Selected AnnData objects are project-local inputs.")
print("- All future outputs remain inside the active project.")
print("- No AnnData objects were loaded.")
print("- No filtering or preprocessing was performed.")
# ============================================================
# 10. Audit external 10x raw files
# ============================================================

def format_size_bytes(size_bytes: int) -> str:
    """Return a human-readable file size."""

    size = float(size_bytes)

    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024 or unit == "TB":
            return f"{size:.2f} {unit}"

        size /= 1024

    return f"{size_bytes} B"


print("\nExternal 10x raw files:")

external_files = sorted(
    path
    for path in external_raw_dir.rglob("*")
    if path.is_file()
)

if not external_files:
    raise FileNotFoundError(
        "The external 10x raw directory exists, "
        "but it contains no files:\n"
        f"{external_raw_dir}"
    )

for path in external_files:
    relative_path = path.relative_to(external_raw_dir)

    print(
        f"- {relative_path} "
        f"({format_size_bytes(path.stat().st_size)})"
    )


# ============================================================
# 11. Identify required 10x files by filename patterns
# ============================================================

def find_files_by_patterns(
    files: list[Path],
    patterns: list[str],
) -> list[Path]:
    """
    Return files whose lowercase filename contains
    every pattern in at least one pattern group.
    """

    matches = []

    for path in files:
        filename = path.name.lower()

        if all(pattern.lower() in filename for pattern in patterns):
            matches.append(path)

    return matches


fragment_candidates = [
    path
    for path in external_files
    if (
        path.name.lower().endswith(".tsv.gz")
        and "fragment" in path.name.lower()
    )
]

fragment_index_candidates = [
    path
    for path in external_files
    if (
        path.name.lower().endswith(".tbi")
        and "fragment" in path.name.lower()
    )
]

peak_bed_candidates = [
    path
    for path in external_files
    if (
        path.suffix.lower() == ".bed"
        and "peak" in path.name.lower()
    )
]

per_barcode_metrics_candidates = [
    path
    for path in external_files
    if (
        path.suffix.lower() == ".csv"
        and "barcode" in path.name.lower()
        and "metric" in path.name.lower()
    )
]

summary_metrics_candidates = [
    path
    for path in external_files
    if (
        path.suffix.lower() == ".csv"
        and "summary" in path.name.lower()
        and "metric" in path.name.lower()
        and "barcode" not in path.name.lower()
    )
]


candidate_groups = {
    "ATAC fragments TSV.GZ": fragment_candidates,
    "ATAC fragments index": fragment_index_candidates,
    "ATAC peak BED": peak_bed_candidates,
    "per-barcode metrics CSV": per_barcode_metrics_candidates,
    "summary metrics CSV": summary_metrics_candidates,
}

print("\nDetected required-file candidates:")

for label, candidates in candidate_groups.items():
    print(f"\n{label}:")

    if not candidates:
        print("  No candidate found.")
        continue

    for candidate in candidates:
        print(f"  {candidate}")


# ============================================================
# 12. Require exactly one essential candidate
# ============================================================

essential_candidate_groups = {
    "ATAC fragments TSV.GZ": fragment_candidates,
    "ATAC fragments index": fragment_index_candidates,
    "ATAC peak BED": peak_bed_candidates,
    "per-barcode metrics CSV": per_barcode_metrics_candidates,
}

resolved_external_inputs = {}

for label, candidates in essential_candidate_groups.items():
    if len(candidates) == 0:
        raise FileNotFoundError(
            f"No candidate was found for: {label}\n"
            f"Searched inside:\n{external_raw_dir}"
        )

    if len(candidates) > 1:
        formatted_candidates = "\n".join(
            str(path)
            for path in candidates
        )

        raise RuntimeError(
            f"Multiple candidates were found for: {label}\n\n"
            f"{formatted_candidates}\n\n"
            "Use explicit filenames instead of automatic detection."
        )

    resolved_external_inputs[label] = candidates[0]


print("\nResolved essential external inputs:")

for label, path in resolved_external_inputs.items():
    print(f"\n{label}:")
    print(path)
    print(f"Size: {format_size_bytes(path.stat().st_size)}")


# ============================================================
# 13. Validate fragments/index pairing
# ============================================================

fragments_path = resolved_external_inputs[
    "ATAC fragments TSV.GZ"
]

fragments_index_path = resolved_external_inputs[
    "ATAC fragments index"
]

expected_index_names = {
    fragments_path.name + ".tbi",
    fragments_path.with_suffix("").name + ".tbi",
}

if fragments_index_path.name not in expected_index_names:
    print(
        "\nWarning: fragments index filename does not exactly "
        "match the detected fragments filename."
    )
    print(f"Fragments: {fragments_path.name}")
    print(f"Index:     {fragments_index_path.name}")

print("\nExternal raw-file audit completed successfully.")
print("No large raw file was opened or parsed.")


# ============================================================
# 14. Resolve explicit external 10x input files
# ============================================================

fragments_path = (
    external_raw_dir
    / "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
)

fragments_index_path = (
    external_raw_dir
    / "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi"
)

peaks_bed_path = (
    external_raw_dir
    / "pbmc_granulocyte_sorted_10k_atac_peaks.bed"
)

peak_annotation_path = (
    external_raw_dir
    / "pbmc_granulocyte_sorted_10k_atac_peak_annotation.tsv"
)

per_barcode_metrics_path = (
    external_raw_dir
    / "pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv"
)

summary_metrics_path = (
    external_raw_dir
    / "pbmc_granulocyte_sorted_10k_summary.csv"
)

filtered_matrix_h5_path = (
    external_raw_dir
    / "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
)


explicit_external_inputs = {
    "fragments": fragments_path,
    "fragments index": fragments_index_path,
    "peaks BED": peaks_bed_path,
    "peak annotation": peak_annotation_path,
    "per-barcode metrics": per_barcode_metrics_path,
    "summary metrics": summary_metrics_path,
    "filtered feature matrix": filtered_matrix_h5_path,
}

print("\nExplicit external input validation:")

for label, path in explicit_external_inputs.items():
    print(f"\n{label}:")
    print(path)
    print(f"Exists: {path.exists()}")

    if not path.exists():
        raise FileNotFoundError(
            f"Required external input was not found:\n"
            f"{label}\n{path}"
        )
    

    # ============================================================
# 15. Inspect per-barcode metrics structure
# ============================================================

import pandas as pd


print("\nReading only the first 5 rows of per-barcode metrics...")

per_barcode_preview = pd.read_csv(
    per_barcode_metrics_path,
    nrows=5,
)

print("\nPer-barcode metrics shape preview:")
print(per_barcode_preview.shape)

print("\nPer-barcode metrics columns:")
for column in per_barcode_preview.columns:
    print(f"- {column}")

print("\nFirst 5 rows:")
print(
    per_barcode_preview.to_string(
        index=False
    )
)

# ============================================================
# 16. Inspect run summary metrics
# ============================================================

summary_metrics = pd.read_csv(
    summary_metrics_path,
)

print("\nRun summary metrics shape:")
print(summary_metrics.shape)

print("\nRun summary metrics columns:")
for column in summary_metrics.columns:
    print(f"- {column}")

print("\nRun summary metrics:")
print(
    summary_metrics.to_string(
        index=False
    )
)




# ============================================================
# 17. Load selected barcode and ATAC metric columns
# ============================================================

import scanpy as sc


per_barcode_columns = [
    "barcode",
    "gex_barcode",
    "atac_barcode",
    "is_cell",
    "excluded_reason",
    "atac_raw_reads",
    "atac_unmapped_reads",
    "atac_lowmapq",
    "atac_dup_reads",
    "atac_chimeric_reads",
    "atac_mitochondrial_reads",
    "atac_fragments",
    "atac_TSS_fragments",
    "atac_peak_region_fragments",
    "atac_peak_region_cutsites",
]

print(
    "\nLoading only barcode and ATAC-specific columns "
    "from per-barcode metrics..."
)

per_barcode_metrics = pd.read_csv(
    per_barcode_metrics_path,
    usecols=per_barcode_columns,
)

print("\nLoaded per-barcode metrics:")
print(per_barcode_metrics.shape)

print("\nMemory usage:")
print(
    f"{per_barcode_metrics.memory_usage(deep=True).sum() / 1024**2:.2f} MB"
)


# ============================================================
# 18. Basic barcode-table validation
# ============================================================

barcode_columns = [
    "barcode",
    "gex_barcode",
    "atac_barcode",
]

for column in barcode_columns:
    per_barcode_metrics[column] = (
        per_barcode_metrics[column]
        .astype(str)
        .str.strip()
    )

    n_missing = int(
        per_barcode_metrics[column]
        .isin(["", "nan", "None"])
        .sum()
    )

    n_duplicated = int(
        per_barcode_metrics[column]
        .duplicated()
        .sum()
    )

    print(f"\nBarcode column: {column}")
    print(
        f"Unique values: "
        f"{per_barcode_metrics[column].nunique()}"
    )
    print(f"Missing-like values: {n_missing}")
    print(f"Duplicated values: {n_duplicated}")



    # ============================================================
# 19. Load selected ATAC for barcode comparison
# ============================================================

print("\nLoading selected ATAC AnnData...")

atac = sc.read_h5ad(
    atac_selected_path,
)

print("\nSelected ATAC:")
print(atac)

selected_atac_barcodes = pd.Index(
    atac.obs_names.astype(str),
    name="selected_atac_barcode",
)

if not selected_atac_barcodes.is_unique:
    duplicated = selected_atac_barcodes[
        selected_atac_barcodes.duplicated()
    ].unique()

    raise ValueError(
        "Selected ATAC obs_names contain duplicated barcodes.\n"
        f"Examples:\n{duplicated[:10].tolist()}"
    )

print("\nSelected ATAC barcode examples:")
print(selected_atac_barcodes[:5].tolist())


# ============================================================
# 20. Identify matching 10x barcode column
# ============================================================

barcode_overlap_rows = []

selected_barcode_set = set(
    selected_atac_barcodes
)

for column in barcode_columns:
    metric_barcode_set = set(
        per_barcode_metrics[column]
    )

    n_matching = len(
        selected_barcode_set.intersection(
            metric_barcode_set
        )
    )

    n_missing_from_metrics = (
        len(selected_barcode_set)
        - n_matching
    )

    overlap_fraction = (
        n_matching
        / len(selected_barcode_set)
    )

    barcode_overlap_rows.append(
        {
            "candidate_column": column,
            "selected_atac_cells": len(
                selected_atac_barcodes
            ),
            "matching_barcodes": n_matching,
            "missing_from_metrics": (
                n_missing_from_metrics
            ),
            "overlap_fraction": overlap_fraction,
        }
    )


barcode_overlap = pd.DataFrame(
    barcode_overlap_rows
).sort_values(
    "overlap_fraction",
    ascending=False,
)

print("\nBarcode overlap with selected ATAC:")
print(
    barcode_overlap.to_string(
        index=False
    )
)


best_barcode_column = str(
    barcode_overlap.iloc[0][
        "candidate_column"
    ]
)

best_overlap_fraction = float(
    barcode_overlap.iloc[0][
        "overlap_fraction"
    ]
)

print("\nBest matching barcode column:")
print(best_barcode_column)

print("\nBest overlap fraction:")
print(f"{best_overlap_fraction:.4f}")


if best_overlap_fraction < 1.0:
    raise ValueError(
        "No per-barcode metrics column matched all "
        "selected ATAC obs_names.\n\n"
        f"Best column: {best_barcode_column}\n"
        f"Overlap fraction: {best_overlap_fraction:.4f}"
    )


print(
    "\nAll selected ATAC barcodes were found "
    f"in '{best_barcode_column}'."
)














# ============================================================
# 21. Prepare per-barcode ATAC QC table
# ============================================================

qc_columns = [
    "barcode",
    "is_cell",
    "excluded_reason",
    "atac_raw_reads",
    "atac_unmapped_reads",
    "atac_lowmapq",
    "atac_dup_reads",
    "atac_chimeric_reads",
    "atac_mitochondrial_reads",
    "atac_fragments",
    "atac_TSS_fragments",
    "atac_peak_region_fragments",
    "atac_peak_region_cutsites",
]

per_barcode_qc = (
    per_barcode_metrics[qc_columns]
    .copy()
)

per_barcode_qc["barcode"] = (
    per_barcode_qc["barcode"]
    .astype(str)
    .str.strip()
)

if per_barcode_qc["barcode"].duplicated().any():
    duplicated_barcodes = (
        per_barcode_qc.loc[
            per_barcode_qc["barcode"].duplicated(keep=False),
            "barcode",
        ]
        .unique()
        .tolist()
    )

    raise ValueError(
        "Duplicate barcodes were found in the per-barcode QC table.\n"
        f"Examples:\n{duplicated_barcodes[:10]}"
    )

per_barcode_qc = (
    per_barcode_qc
    .set_index("barcode")
)
# ============================================================
# 22. Compute derived ATAC QC fractions
# ============================================================

def safe_fraction(
    numerator: pd.Series,
    denominator: pd.Series,
) -> pd.Series:
    """Compute a fraction and return NaN where denominator is zero."""

    denominator = denominator.astype(float)
    numerator = numerator.astype(float)

    result = numerator / denominator.replace(0, np.nan)

    return result


per_barcode_qc["atac_frip"] = safe_fraction(
    per_barcode_qc["atac_peak_region_fragments"],
    per_barcode_qc["atac_fragments"],
)

per_barcode_qc["atac_tss_fragment_fraction"] = safe_fraction(
    per_barcode_qc["atac_TSS_fragments"],
    per_barcode_qc["atac_fragments"],
)

per_barcode_qc["atac_mitochondrial_fraction"] = safe_fraction(
    per_barcode_qc["atac_mitochondrial_reads"],
    per_barcode_qc["atac_raw_reads"],
)

per_barcode_qc["atac_duplicate_fraction"] = safe_fraction(
    per_barcode_qc["atac_dup_reads"],
    per_barcode_qc["atac_raw_reads"],
)

per_barcode_qc["atac_lowmapq_fraction"] = safe_fraction(
    per_barcode_qc["atac_lowmapq"],
    per_barcode_qc["atac_raw_reads"],
)


# ============================================================
# 23. Align selected QC rows to ATAC barcode order
# ============================================================

atac_barcodes = pd.Index(
    atac.obs_names.astype(str),
    name="barcode",
)

per_barcode_qc.index = pd.Index(
    per_barcode_qc.index.astype(str),
    name="barcode",
)

selected_qc = per_barcode_qc.reindex(
    atac_barcodes
).copy()

# Check that reindexing did not introduce missing barcode rows
missing_barcode_rows = selected_qc.isna().all(axis=1)

if missing_barcode_rows.any():
    missing_examples = (
        selected_qc.index[missing_barcode_rows]
        .tolist()[:10]
    )

    raise ValueError(
        "Some selected ATAC barcodes were not found "
        "in the QC metrics table.\n"
        f"Examples: {missing_examples}"
    )

# Force exactly the same index object/order as atac.obs_names
selected_qc.index = atac.obs_names.copy()

if not np.array_equal(
    selected_qc.index.astype(str),
    atac.obs_names.astype(str),
):
    raise ValueError(
        "Selected QC rows are not aligned with "
        "ATAC obs_names after reindexing."
    )

print(
    "\nSelected QC table successfully aligned "
    "to ATAC cell order."
)

overlapping_columns = [
    column
    for column in selected_qc.columns
    if column in atac.obs.columns
]

if overlapping_columns:
    raise ValueError(
        "The following QC columns already exist in atac.obs:\n"
        f"{overlapping_columns}"
    )

atac.obs = atac.obs.join(
    selected_qc,
    how="left",
)

print("\nJoined ATAC QC metadata into atac.obs.")

print("\nNew ATAC obs columns:")
for column in selected_qc.columns:
    print(f"- {column}")
    # ============================================================
# 24. Validate joined QC metadata
# ============================================================

qc_missingness = pd.DataFrame(
    {
        "column": selected_qc.columns,
        "missing_count": [
            int(atac.obs[column].isna().sum())
            for column in selected_qc.columns
        ],
        "missing_fraction": [
            float(atac.obs[column].isna().mean())
            for column in selected_qc.columns
        ],
    }
)

print("\nATAC QC metadata missingness:")
print(
    qc_missingness.to_string(
        index=False
    )
)

if (
    qc_missingness["missing_count"]
    .gt(0)
    .any()
):
    print(
        "\nWarning: some joined QC columns contain missing values."
    )
else:
    print(
        "\nAll selected cells have complete joined ATAC QC metadata."
    )


    # ============================================================
# 25. Validate 10x cell-calling status
# ============================================================

print("\n10x is_cell counts:")
print(
    atac.obs["is_cell"]
    .value_counts(dropna=False)
    .sort_index()
)

n_not_cell = int(
    (atac.obs["is_cell"] != 1).sum()
)

print(f"\nSelected cells with is_cell != 1: {n_not_cell}")

if n_not_cell > 0:
    examples = (
        atac.obs.loc[
            atac.obs["is_cell"] != 1,
            [
                "is_cell",
                "excluded_reason",
                "cell_type_broad",
            ],
        ]
        .head(10)
    )

    print("\nExamples:")
    print(examples)












    # ============================================================
# 26. Overall ATAC QC summary
# ============================================================

qc_metric_columns = [
    "atac_raw_reads",
    "atac_fragments",
    "atac_frip",
    "atac_tss_fragment_fraction",
    "atac_mitochondrial_fraction",
    "atac_duplicate_fraction",
    "atac_lowmapq_fraction",
]

overall_qc_summary = (
    atac.obs[qc_metric_columns]
    .describe(
        percentiles=[
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
    .T
)


# ============================================================
# 27. ATAC QC summary by broad cell type
# ============================================================

primary_cell_type_key = "cell_type_broad"

if primary_cell_type_key not in atac.obs.columns:
    raise KeyError(
        f"Missing required obs column: {primary_cell_type_key}"
    )

cell_type_counts = (
    atac.obs[primary_cell_type_key]
    .value_counts(dropna=False)
    .rename_axis(primary_cell_type_key)
    .reset_index(name="n_cells")
)

cell_type_counts["fraction"] = (
    cell_type_counts["n_cells"]
    / atac.n_obs
)

print("\nSelected cells by broad cell type:")
print(
    cell_type_counts.to_string(
        index=False
    )
)

qc_by_cell_type = (
    atac.obs
    .groupby(
        primary_cell_type_key,
        observed=True,
    )[qc_metric_columns]
    .agg(
        ["count", "median", "mean", "min", "max"]
    )
)

print("\nATAC QC by broad cell type:")
print(
    qc_by_cell_type.to_string()
)






# ============================================================
# 28. Compute matrix-derived ATAC QC metrics
# ============================================================

from scipy import sparse


X = atac.X

if sparse.issparse(X):
    X = X.tocsr()
else:
    X = sparse.csr_matrix(X)

atac.obs["matrix_total_counts"] = (
    np.asarray(X.sum(axis=1))
    .ravel()
)

atac.obs["matrix_n_accessible_peaks"] = (
    np.asarray((X > 0).sum(axis=1))
    .ravel()
)

atac.var["matrix_n_cells_by_counts"] = (
    np.asarray((X > 0).sum(axis=0))
    .ravel()
)

atac.var["matrix_total_counts"] = (
    np.asarray(X.sum(axis=0))
    .ravel()
)

print("\nAdded matrix-derived ATAC QC metrics:")
print("- matrix_total_counts")
print("- matrix_n_accessible_peaks")
print("- matrix_n_cells_by_counts")
print("- matrix_total_counts")




# ============================================================
# 29. Build overall ATAC QC quantile table
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
]

cell_qc_metrics = [
    "atac_fragments",
    "matrix_total_counts",
    "matrix_n_accessible_peaks",
    "atac_frip",
    "atac_tss_fragment_fraction",
    "atac_mitochondrial_fraction",
    "atac_duplicate_fraction",
    "atac_lowmapq_fraction",
]

overall_quantiles = (
    atac.obs[cell_qc_metrics]
    .quantile(quantile_levels)
)

overall_quantiles.index.name = "quantile"

print("\nOverall ATAC QC quantiles:")
print(
    overall_quantiles.to_string()
)





# ============================================================
# 30. Save overall quantile table
# ============================================================

tables_output_dir.mkdir(
    parents=True,
    exist_ok=True,
)

overall_quantiles_out = (
    tables_output_dir
    / "atac_qc_overall_quantiles_before_filtering.csv"
)

overall_quantiles.reset_index().to_csv(
    overall_quantiles_out,
    index=False,
    encoding="utf-8-sig",
)

print("\nSaved overall quantiles:")
print(overall_quantiles_out)

# ============================================================
# 31. Build ATAC QC quantiles by broad cell type
# ============================================================

cell_type_quantile_tables = []

for cell_type, group in atac.obs.groupby(
    "cell_type_broad",
    observed=True,
):
    group_quantiles = (
        group[cell_qc_metrics]
        .quantile(quantile_levels)
        .reset_index()
        .rename(columns={"index": "quantile"})
    )

    group_quantiles.insert(
        0,
        "cell_type_broad",
        cell_type,
    )

    group_quantiles.insert(
        1,
        "n_cells",
        int(group.shape[0]),
    )

    cell_type_quantile_tables.append(
        group_quantiles
    )

quantiles_by_cell_type = pd.concat(
    cell_type_quantile_tables,
    axis=0,
    ignore_index=True,
)

print("\nATAC QC quantiles by broad cell type:")
print(
    quantiles_by_cell_type.head(20).to_string(
        index=False
    )
)
# ============================================================
# 32. Save quantiles by cell type
# ============================================================

quantiles_by_cell_type_out = (
    tables_output_dir
    / "atac_qc_quantiles_by_cell_type_before_filtering.csv"
)

quantiles_by_cell_type.to_csv(
    quantiles_by_cell_type_out,
    index=False,
    encoding="utf-8-sig",
)

print("\nSaved cell-type quantiles:")
print(quantiles_by_cell_type_out)








# ============================================================
# 33. Configure diagnostic plotting
# ============================================================

figures_output_dir.mkdir(
    parents=True,
    exist_ok=True,
)

sns.set_theme(
    style="whitegrid",
    context="notebook",
)


def save_show_close(
    output_path: Path,
    dpi: int = 300,
    show: bool = True,
) -> None:
    """Save, optionally display, and then close the figure."""

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
    )

    if show:
        plt.show()

    plt.close()


# ============================================================
# 34. Overall ATAC QC histograms
# ============================================================

histogram_metrics = {
    "atac_fragments": "ATAC fragments per cell",
    "matrix_n_accessible_peaks": "Accessible peaks per cell",
    "atac_frip": "FRiP-like fraction",
    "atac_mitochondrial_fraction": "Mitochondrial read fraction",
}

for metric, label in histogram_metrics.items():

    plt.figure(
        figsize=(7, 4)
    )

    sns.histplot(
        data=atac.obs,
        x=metric,
        bins=100,
    )

    plt.xlabel(label)
    plt.ylabel("Number of cells")

    plt.title(
        f"{label} before ATAC QC filtering"
    )

    output_path = (
        figures_output_dir
        / f"{metric}_histogram_before_filtering.png"
    )

    save_show_close(
        output_path,
        show=True,
    )

    print(f"Saved: {output_path}")
    
    
    
    
    
# ============================================================
# 35. Log10 histograms for count-like metrics
# ============================================================

log_histogram_metrics = {
    "atac_fragments": "ATAC fragments per cell",
    "matrix_total_counts": "ATAC matrix counts per cell",
    "matrix_n_accessible_peaks": "Accessible peaks per cell",
}

for metric, label in log_histogram_metrics.items():

    positive_values = atac.obs.loc[
        atac.obs[metric] > 0,
        metric,
    ]

    plt.figure(
        figsize=(7, 4)
    )

    sns.histplot(
        x=np.log10(positive_values),
        bins=100,
    )

    plt.xlabel(
        f"log10({label})"
    )

    plt.ylabel("Number of cells")

    plt.title(
        f"{label}, log10 scale"
    )

    output_path = (
        figures_output_dir
        / f"{metric}_log10_histogram_before_filtering.png"
    )

    save_show_close(
        output_path,
        show=True,
    )

    print(f"Saved: {output_path}")


# ============================================================
# Step 36: Violin plots of ATAC QC metrics by cell_type_broad
# ============================================================

import matplotlib.pyplot as plt
import seaborn as sns

step36_output_dir = (
    figures_output_dir
    / "step_36_violin_plots_by_cell_type"
)

step36_output_dir.mkdir(
    parents=True,
    exist_ok=True,
)


cell_type_col = "cell_type_broad"

violin_metrics = [
    "matrix_total_counts",
    "matrix_n_accessible_peaks",
    "atac_frip",
    "atac_tss_fragment_fraction",
    "atac_mitochondrial_fraction",
    "atac_duplicate_fraction",
    "atac_lowmapq_fraction",
]

missing_violin_columns = [
    column
    for column in [cell_type_col, *violin_metrics]
    if column not in atac.obs.columns
]

if missing_violin_columns:
    raise KeyError(
        "Step 36 cannot run because these atac.obs columns are missing: "
        f"{missing_violin_columns}"
    )

cell_type_order = (
    atac.obs[cell_type_col]
    .astype(str)
    .value_counts()
    .index
    .tolist()
)

print("\nStep 36")
print("Creating violin plots by cell_type_broad.")
print(f"Cell-type order, from largest to smallest group: {cell_type_order}")


for metric in violin_metrics:

    plot_df = (
        atac.obs[
            [
                cell_type_col,
                metric,
            ]
        ]
        .dropna()
        .copy()
    )

    plot_df[cell_type_col] = (
        plot_df[cell_type_col]
        .astype(str)
    )

    figure_width = max(
        10,
        1.25 * len(cell_type_order),
    )

    plt.figure(
        figsize=(figure_width, 6)
    )

    sns.violinplot(
        data=plot_df,
        x=cell_type_col,
        y=metric,
        order=cell_type_order,
        inner="quartile",
        cut=0,
        density_norm="width",
        linewidth=0.8,
        color="steelblue",
    )

    plt.xlabel("Broad cell type")
    plt.ylabel(metric)

    plt.title(
        f"{metric} by cell_type_broad "
        "before ATAC QC filtering"
    )

    plt.xticks(
        rotation=45,
        ha="right",
    )

    output_path = (
        step36_output_dir
        / f"{metric}_violin_by_cell_type_before_filtering.png"
    )

    save_show_close(
        output_path,
        show=True,
    )

    print(f"Saved: {output_path}")
    
    
    
    
# ============================================================
# Step 37: Scatter plots for fragments, accessible peaks, FRiP
# ============================================================

step37_output_dir = (
    figures_output_dir
    / "step_37_qc_scatter_plots"
)

step37_output_dir.mkdir(
    parents=True,
    exist_ok=True,
)

required_scatter_columns = [
    "cell_type_broad",
    "atac_fragments",
    "matrix_n_accessible_peaks",
    "atac_frip",
]

missing_scatter_columns = [
    column
    for column in required_scatter_columns
    if column not in atac.obs.columns
]

if missing_scatter_columns:
    raise KeyError(
        "Step 37 cannot run because these atac.obs "
        f"columns are missing: {missing_scatter_columns}"
    )

scatter_df = (
    atac.obs[required_scatter_columns]
    .dropna()
    .copy()
)

scatter_df["cell_type_broad"] = (
    scatter_df["cell_type_broad"]
    .astype(str)
)

scatter_df = scatter_df.loc[
    (scatter_df["atac_fragments"] > 0)
    & (scatter_df["matrix_n_accessible_peaks"] > 0)
].copy()

scatter_df["log10_atac_fragments"] = np.log10(
    scatter_df["atac_fragments"].astype(float)
)

scatter_df["log10_matrix_n_accessible_peaks"] = np.log10(
    scatter_df["matrix_n_accessible_peaks"].astype(float)
)

print("\nStep 37")
print("Prepared scatter-plot dataframe.")
print(f"Cells included: {scatter_df.shape[0]:,}")
print(f"Cells excluded because of missing or zero values: "
      f"{atac.n_obs - scatter_df.shape[0]:,}")


# ============================================================
# Step 37A: Fragments versus accessible peaks
# ============================================================

plt.figure(
    figsize=(8, 7)
)

sns.scatterplot(
    data=scatter_df,
    x="log10_atac_fragments",
    y="log10_matrix_n_accessible_peaks",
    hue="cell_type_broad",
    hue_order=cell_type_order,
    s=14,
    alpha=0.45,
    linewidth=0,
)

plt.xlabel(
    "log10(ATAC fragments per cell)"
)

plt.ylabel(
    "log10(accessible peaks per cell)"
)

plt.title(
    "ATAC fragments versus accessible peaks\n"
    "before ATAC QC filtering"
)

plt.legend(
    title="cell_type_broad",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0,
    markerscale=1.4,
)

output_path = (
    step37_output_dir
    / "fragments_vs_accessible_peaks_before_filtering.png"
)

save_show_close(
    output_path,
    show=True,
)

print(f"Saved: {output_path}")



# ============================================================
# Step 37B: Fragments versus FRiP
# ============================================================

plt.figure(
    figsize=(8, 7)
)

sns.scatterplot(
    data=scatter_df,
    x="log10_atac_fragments",
    y="atac_frip",
    hue="cell_type_broad",
    hue_order=cell_type_order,
    s=14,
    alpha=0.45,
    linewidth=0,
)

plt.xlabel(
    "log10(ATAC fragments per cell)"
)

plt.ylabel(
    "ATAC FRiP"
)

plt.title(
    "ATAC fragments versus FRiP\n"
    "before ATAC QC filtering"
)

plt.legend(
    title="cell_type_broad",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0,
    markerscale=1.4,
)

output_path = (
    step37_output_dir
    / "fragments_vs_frip_before_filtering.png"
)

save_show_close(
    output_path,
    show=True,
)

print(f"Saved: {output_path}")


# ============================================================
# Step 37C: Accessible peaks versus FRiP
# ============================================================

plt.figure(
    figsize=(8, 7)
)

sns.scatterplot(
    data=scatter_df,
    x="log10_matrix_n_accessible_peaks",
    y="atac_frip",
    hue="cell_type_broad",
    hue_order=cell_type_order,
    s=14,
    alpha=0.45,
    linewidth=0,
)

plt.xlabel(
    "log10(accessible peaks per cell)"
)

plt.ylabel(
    "ATAC FRiP"
)

plt.title(
    "Accessible peaks versus FRiP\n"
    "before ATAC QC filtering"
)

plt.legend(
    title="cell_type_broad",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0,
    markerscale=1.4,
)

output_path = (
    step37_output_dir
    / "accessible_peaks_vs_frip_before_filtering.png"
)

save_show_close(
    output_path,
    show=True,
)

print(f"Saved: {output_path}")



# ============================================================
# Step 38: Inspect candidate global QC thresholds
# ============================================================

threshold_metrics = [
    "atac_fragments",
    "matrix_n_accessible_peaks",
    "atac_frip",
]

threshold_quantiles = [
    0.00,
    0.01,
    0.02,
    0.05,
    0.10,
    0.25,
    0.50,
]

candidate_threshold_table = (
    atac.obs[threshold_metrics]
    .quantile(threshold_quantiles)
)

candidate_threshold_table.index.name = "quantile"

print("\nStep 38")
print("Candidate global threshold support table:")
print(
    candidate_threshold_table.to_string()
)

print(
    "\nNo filtering mask or subset was created."
)


# ============================================================
# Step 39: Interactively define candidate threshold sets
# ============================================================

def prompt_threshold(
    label: str,
    default_value: float,
    minimum: float = 0.0,
    maximum: float | None = None,
    integer: bool = False,
) -> int | float:
    """
    Prompt for a validated numeric threshold.

    Press Enter to keep the data-derived default.
    """

    while True:
        raw_value = input(
            f"{label} [{default_value}]: "
        ).strip()

        if raw_value == "":
            value = float(default_value)
        else:
            try:
                value = float(raw_value)
            except ValueError:
                print(
                    "Invalid input. Enter a number or press Enter."
                )
                continue

        if value < minimum:
            print(f"Value must be >= {minimum}.")
            continue

        if maximum is not None and value > maximum:
            print(f"Value must be <= {maximum}.")
            continue

        if integer:
            return int(round(value))

        return float(value)
    
conservative_defaults = {
    "min_atac_fragments": int(
        round(
            candidate_threshold_table.loc[
                0.05,
                "atac_fragments",
            ]
            / 100
        )
        * 100
    ),
    "min_accessible_peaks": int(
        round(
            candidate_threshold_table.loc[
                0.05,
                "matrix_n_accessible_peaks",
            ]
            / 100
        )
        * 100
    ),
    "min_atac_frip": round(
        float(
            candidate_threshold_table.loc[
                0.05,
                "atac_frip",
            ]
        ),
        2,
    ),
}

lenient_defaults = {
    "min_atac_fragments": int(
        round(
            candidate_threshold_table.loc[
                0.02,
                "atac_fragments",
            ]
            / 100
        )
        * 100
    ),
    "min_accessible_peaks": int(
        round(
            candidate_threshold_table.loc[
                0.02,
                "matrix_n_accessible_peaks",
            ]
            / 100
        )
        * 100
    ),
    "min_atac_frip": round(
        float(
            candidate_threshold_table.loc[
                0.02,
                "atac_frip",
            ]
        ),
        2,
    ),
}

print("\nStep 39")
print("Define the conservative candidate set.")
print("Press Enter to keep each data-derived suggestion.")

conservative_thresholds = {
    "min_atac_fragments": prompt_threshold(
        label="Conservative minimum ATAC fragments",
        default_value=conservative_defaults[
            "min_atac_fragments"
        ],
        integer=True,
    ),
    "min_accessible_peaks": prompt_threshold(
        label="Conservative minimum accessible peaks",
        default_value=conservative_defaults[
            "min_accessible_peaks"
        ],
        integer=True,
    ),
    "min_atac_frip": prompt_threshold(
        label="Conservative minimum ATAC FRiP",
        default_value=conservative_defaults[
            "min_atac_frip"
        ],
        minimum=0,
        maximum=1,
    ),
}

print("\nDefine the lenient candidate set.")
print("Press Enter to keep each data-derived suggestion.")

lenient_thresholds = {
    "min_atac_fragments": prompt_threshold(
        label="Lenient minimum ATAC fragments",
        default_value=lenient_defaults[
            "min_atac_fragments"
        ],
        integer=True,
    ),
    "min_accessible_peaks": prompt_threshold(
        label="Lenient minimum accessible peaks",
        default_value=lenient_defaults[
            "min_accessible_peaks"
        ],
        integer=True,
    ),
    "min_atac_frip": prompt_threshold(
        label="Lenient minimum ATAC FRiP",
        default_value=lenient_defaults[
            "min_atac_frip"
        ],
        minimum=0,
        maximum=1,
    ),
}

candidate_sets = {
    "conservative": conservative_thresholds,
    "lenient": lenient_thresholds,
}

print("\nCandidate threshold sets entered by the user:")

for candidate_name, thresholds in candidate_sets.items():
    print(f"\n{candidate_name}:")
    for threshold_name, threshold_value in thresholds.items():
        print(f"- {threshold_name}: {threshold_value}")

print("\nNo filtering mask or subset was created.")


# ============================================================
# Step 40: Compare interactive candidate threshold sets
# ============================================================

candidate_comparison_rows = []
candidate_masks = {}

for candidate_name, thresholds in candidate_sets.items():

    candidate_mask = (
        (
            atac.obs["atac_fragments"]
            >= thresholds["min_atac_fragments"]
        )
        &
        (
            atac.obs["matrix_n_accessible_peaks"]
            >= thresholds["min_accessible_peaks"]
        )
        &
        (
            atac.obs["atac_frip"]
            >= thresholds["min_atac_frip"]
        )
    )

    candidate_masks[candidate_name] = candidate_mask

    n_pass = int(candidate_mask.sum())
    n_fail = int((~candidate_mask).sum())

    candidate_comparison_rows.append(
        {
            "candidate_set": candidate_name,
            "min_atac_fragments": (
                thresholds["min_atac_fragments"]
            ),
            "min_accessible_peaks": (
                thresholds["min_accessible_peaks"]
            ),
            "min_atac_frip": (
                thresholds["min_atac_frip"]
            ),
            "n_pass": n_pass,
            "n_fail": n_fail,
            "retention_fraction": n_pass / atac.n_obs,
            "loss_fraction": n_fail / atac.n_obs,
        }
    )

candidate_comparison = pd.DataFrame(
    candidate_comparison_rows
)

print("\nStep 40")
print("Candidate threshold comparison:")

print(
    candidate_comparison.to_string(
        index=False,
        formatters={
            "retention_fraction": "{:.2%}".format,
            "loss_fraction": "{:.2%}".format,
        },
    )
)

print("\nNo AnnData subset was performed.")

# ============================================================
# Step 41: Select candidate threshold set
# ============================================================

valid_candidate_names = set(candidate_sets)

while True:
    selected_candidate_name = input(
        "\nSelect candidate set "
        "[conservative/lenient]: "
    ).strip().lower()

    if selected_candidate_name in valid_candidate_names:
        break

    print(
        "Invalid selection. Enter 'conservative' or 'lenient'."
    )


selected_candidate_thresholds = candidate_sets[
    selected_candidate_name
].copy()

candidate_pass_mask = candidate_masks[
    selected_candidate_name
].copy()

print("\nStep 41")
print(f"Selected candidate set: {selected_candidate_name}")

print("\nSelected candidate thresholds:")

for threshold_name, threshold_value in (
    selected_candidate_thresholds.items()
):
    print(f"- {threshold_name}: {threshold_value}")

print("\nSelection is not yet final.")
print("No AnnData subset was performed.")

# ============================================================
# Step 42: Preview retention by cell_type_broad
# ============================================================

retention_by_cell_type = (
    pd.DataFrame(
        {
            "cell_type_broad": (
                atac.obs["cell_type_broad"]
                .astype(str)
            ),
            "candidate_pass": candidate_pass_mask,
        },
        index=atac.obs_names,
    )
    .groupby(
        "cell_type_broad",
        observed=True,
    )
    .agg(
        n_before=("candidate_pass", "size"),
        n_after=("candidate_pass", "sum"),
    )
    .reset_index()
)

retention_by_cell_type["n_lost"] = (
    retention_by_cell_type["n_before"]
    - retention_by_cell_type["n_after"]
)

retention_by_cell_type["retention_fraction"] = (
    retention_by_cell_type["n_after"]
    / retention_by_cell_type["n_before"]
)

retention_by_cell_type["loss_fraction"] = (
    retention_by_cell_type["n_lost"]
    / retention_by_cell_type["n_before"]
)

retention_by_cell_type = (
    retention_by_cell_type
    .sort_values(
        "loss_fraction",
        ascending=False,
    )
)

print("\nStep 42")
print(
    f"Retention preview for selected set: "
    f"{selected_candidate_name}"
)

print(
    retention_by_cell_type.to_string(
        index=False,
        formatters={
            "retention_fraction": "{:.2%}".format,
            "loss_fraction": "{:.2%}".format,
        },
    )
)

loss_warning_threshold = 0.15

high_loss_groups = retention_by_cell_type.loc[
    retention_by_cell_type["loss_fraction"]
    > loss_warning_threshold
].copy()

if high_loss_groups.empty:
    print(
        "\nNo cell_type_broad group loses more than 15%."
    )
else:
    print(
        "\nWARNING: These cell_type_broad groups "
        "lose more than 15%:"
    )

    print(
        high_loss_groups.to_string(
            index=False,
            formatters={
                "retention_fraction": "{:.2%}".format,
                "loss_fraction": "{:.2%}".format,
            },
        )
    )

fail_fragments = (
    atac.obs["atac_fragments"]
    < selected_candidate_thresholds[
        "min_atac_fragments"
    ]
)

fail_accessible_peaks = (
    atac.obs["matrix_n_accessible_peaks"]
    < selected_candidate_thresholds[
        "min_accessible_peaks"
    ]
)

fail_frip = (
    atac.obs["atac_frip"]
    < selected_candidate_thresholds[
        "min_atac_frip"
    ]
)

fail_any = (
    fail_fragments
    | fail_accessible_peaks
    | fail_frip
)

failure_reason_df = pd.DataFrame(
    {
        "cell_type_broad": (
            atac.obs["cell_type_broad"]
            .astype(str)
        ),
        "fail_fragments": fail_fragments,
        "fail_accessible_peaks": fail_accessible_peaks,
        "fail_frip": fail_frip,
        "fail_any": fail_any,
    },
    index=atac.obs_names,
)

failure_by_cell_type = (
    failure_reason_df
    .groupby(
        "cell_type_broad",
        observed=True,
    )
    .agg(
        n_cells=("fail_any", "size"),
        n_fail_fragments=("fail_fragments", "sum"),
        n_fail_accessible_peaks=(
            "fail_accessible_peaks",
            "sum",
        ),
        n_fail_frip=("fail_frip", "sum"),
        n_fail_any=("fail_any", "sum"),
    )
    .reset_index()
)

failure_by_cell_type["fraction_fail_fragments"] = (
    failure_by_cell_type["n_fail_fragments"]
    / failure_by_cell_type["n_cells"]
)

failure_by_cell_type[
    "fraction_fail_accessible_peaks"
] = (
    failure_by_cell_type[
        "n_fail_accessible_peaks"
    ]
    / failure_by_cell_type["n_cells"]
)

failure_by_cell_type["fraction_fail_frip"] = (
    failure_by_cell_type["n_fail_frip"]
    / failure_by_cell_type["n_cells"]
)

failure_by_cell_type["fraction_fail_any"] = (
    failure_by_cell_type["n_fail_any"]
    / failure_by_cell_type["n_cells"]
)

failure_by_cell_type = (
    failure_by_cell_type
    .sort_values(
        "fraction_fail_any",
        ascending=False,
    )
)

print("\nFailure reasons by cell_type_broad:")

print(
    failure_by_cell_type.to_string(
        index=False,
        formatters={
            "fraction_fail_fragments": "{:.2%}".format,
            "fraction_fail_accessible_peaks": (
                "{:.2%}".format
            ),
            "fraction_fail_frip": "{:.2%}".format,
            "fraction_fail_any": "{:.2%}".format,
        },
    )
)

print("\nThis is report-only.")
print("No AnnData subset was performed.")

## ============================================================
# Step 43: Confirm selected thresholds and save YAML
# ============================================================

import yaml


print("\nStep 43")
print("Final threshold approval.")

print(f"\nSelected set: {selected_candidate_name}")

for threshold_name, threshold_value in (
    selected_candidate_thresholds.items()
):
    print(f"- {threshold_name}: {threshold_value}")

overall_n_pass = int(candidate_pass_mask.sum())
overall_n_fail = int((~candidate_pass_mask).sum())

print("\nFinal preview:")
print(f"- total cells: {atac.n_obs:,}")
print(f"- passing cells: {overall_n_pass:,}")
print(f"- failing cells: {overall_n_fail:,}")
print(
    f"- retention: "
    f"{overall_n_pass / atac.n_obs:.2%}"
)

while True:
    confirmation = input(
        "\nApprove and save these thresholds? [y/n]: "
    ).strip().lower()

    if confirmation in {"y", "yes"}:
        thresholds_approved = True
        break

    if confirmation in {"n", "no"}:
        thresholds_approved = False
        break

    print("Enter 'y' or 'n'.")


if not thresholds_approved:
    raise RuntimeError(
        "ATAC QC threshold approval was cancelled. "
        "No YAML was written and no filtering was performed."
    )
    

atac_qc_thresholds_path = (
    manifest_output_dir
    / "step04A_approved_atac_qc_thresholds.yaml"
)

manifest_output_dir.mkdir(
    parents=True,
    exist_ok=True,
)

approved_threshold_config = {
    "project_id": paths.project_id,
    "script": (
        "run_04A_paired_atac_qc_and_final_filtering.py"
    ),
    "selected_candidate_set": selected_candidate_name,
    "thresholds": selected_candidate_thresholds,
    "preview": {
        "n_cells_before": int(atac.n_obs),
        "n_cells_pass": overall_n_pass,
        "n_cells_fail": overall_n_fail,
        "retention_fraction": float(
            overall_n_pass / atac.n_obs
        ),
    },
    "warning_rules": {
        "max_allowed_cell_type_loss_fraction": 0.15,
        "warning_only": True,
    },
}

with atac_qc_thresholds_path.open(
    "w",
    encoding="utf-8",
) as file:
    yaml.safe_dump(
        approved_threshold_config,
        file,
        sort_keys=False,
        allow_unicode=True,
    )

print("\nApproved ATAC QC thresholds saved to:")
print(atac_qc_thresholds_path)

print("\nNo AnnData subset was performed.")    

# ============================================================
# Step 44: Load and validate approved ATAC QC thresholds
# ============================================================

import yaml


if not atac_qc_thresholds_path.is_file():
    raise FileNotFoundError(
        "Approved ATAC QC threshold YAML was not found:\n"
        f"{atac_qc_thresholds_path}"
    )


with atac_qc_thresholds_path.open(
    "r",
    encoding="utf-8",
) as file:
    loaded_atac_qc_config = yaml.safe_load(file)


if not isinstance(loaded_atac_qc_config, dict):
    raise TypeError(
        "The approved ATAC QC YAML must contain a dictionary."
    )


loaded_project_id = loaded_atac_qc_config.get(
    "project_id"
)

if loaded_project_id != paths.project_id:
    raise RuntimeError(
        "The approved ATAC QC YAML belongs to a different project.\n"
        f"Expected project ID: {paths.project_id}\n"
        f"Observed project ID: {loaded_project_id}"
    )


loaded_thresholds = loaded_atac_qc_config.get(
    "thresholds"
)

if not isinstance(loaded_thresholds, dict):
    raise TypeError(
        "Missing or invalid 'thresholds' section in YAML."
    )


required_threshold_keys = {
    "min_atac_fragments",
    "min_accessible_peaks",
    "min_atac_frip",
}

missing_threshold_keys = (
    required_threshold_keys
    - set(loaded_thresholds)
)

if missing_threshold_keys:
    raise KeyError(
        "Missing required threshold keys in YAML: "
        f"{sorted(missing_threshold_keys)}"
    )
    
    
# ============================================================
# Step 44 continued: Validate threshold values
# ============================================================

approved_atac_qc_thresholds = {
    "min_atac_fragments": int(
        loaded_thresholds["min_atac_fragments"]
    ),
    "min_accessible_peaks": int(
        loaded_thresholds["min_accessible_peaks"]
    ),
    "min_atac_frip": float(
        loaded_thresholds["min_atac_frip"]
    ),
}


if approved_atac_qc_thresholds[
    "min_atac_fragments"
] < 0:
    raise ValueError(
        "min_atac_fragments cannot be negative."
    )

if approved_atac_qc_thresholds[
    "min_accessible_peaks"
] < 0:
    raise ValueError(
        "min_accessible_peaks cannot be negative."
    )

if not (
    0.0
    <= approved_atac_qc_thresholds["min_atac_frip"]
    <= 1.0
):
    raise ValueError(
        "min_atac_frip must be between 0 and 1."
    )


print("\nStep 44")
print("Approved ATAC QC thresholds loaded from YAML:")

for threshold_name, threshold_value in (
    approved_atac_qc_thresholds.items()
):
    print(f"- {threshold_name}: {threshold_value}")

print("\nYAML project validation passed.")
print("No filtering mask was created.")
print("No AnnData subset was performed.")


# ============================================================
# Step 45: Build final ATAC QC mask and filter ATAC
# ============================================================

pass_atac_fragments_qc = (
    atac.obs["atac_fragments"]
    >= approved_atac_qc_thresholds["min_atac_fragments"]
)

pass_atac_accessible_peaks_qc = (
    atac.obs["matrix_n_accessible_peaks"]
    >= approved_atac_qc_thresholds["min_accessible_peaks"]
)

pass_atac_frip_qc = (
    atac.obs["atac_frip"]
    >= approved_atac_qc_thresholds["min_atac_frip"]
)

final_atac_qc_mask = (
    pass_atac_fragments_qc
    & pass_atac_accessible_peaks_qc
    & pass_atac_frip_qc
)

final_atac_qc_mask = pd.Series(
    final_atac_qc_mask.to_numpy(dtype=bool),
    index=atac.obs_names.copy(),
    name="pass_final_atac_qc",
)


if final_atac_qc_mask.isna().any():
    raise ValueError(
        "The final ATAC QC mask contains missing values."
    )

if not final_atac_qc_mask.index.equals(atac.obs_names):
    raise ValueError(
        "The final ATAC QC mask is not aligned with atac.obs_names."
    )

n_before = int(atac.n_obs)
n_pass = int(final_atac_qc_mask.sum())
n_fail = int((~final_atac_qc_mask).sum())

print("\nStep 45")
print("Final ATAC QC mask created.")
print(f"Cells before filtering: {n_before:,}")
print(f"Cells passing QC: {n_pass:,}")
print(f"Cells failing QC: {n_fail:,}")
print(f"Retention: {n_pass / n_before:.2%}")

atac.obs["pass_atac_fragments_qc"] = (
    pass_atac_fragments_qc
    .reindex(atac.obs_names)
    .to_numpy(dtype=bool)
)

atac.obs["pass_atac_accessible_peaks_qc"] = (
    pass_atac_accessible_peaks_qc
    .reindex(atac.obs_names)
    .to_numpy(dtype=bool)
)

atac.obs["pass_atac_frip_qc"] = (
    pass_atac_frip_qc
    .reindex(atac.obs_names)
    .to_numpy(dtype=bool)
)

atac.obs["pass_final_atac_qc"] = (
    final_atac_qc_mask.to_numpy(dtype=bool)
)


atac_filtered = atac[
    final_atac_qc_mask.to_numpy()
].copy()

print("\nFiltered ATAC object:")
print(atac_filtered)

print("\nATAC filtering audit:")
print(f"Original cells: {atac.n_obs:,}")
print(f"Filtered cells: {atac_filtered.n_obs:,}")
print(f"Removed cells: {atac.n_obs - atac_filtered.n_obs:,}")
print(f"Peaks retained: {atac_filtered.n_vars:,}")


# ============================================================
# Step 46: Load selected RNA and Multiome objects
# ============================================================

print("\nStep 46")
print("Loading selected RNA and Multiome objects...")

rna = sc.read_h5ad(
    rna_selected_path,
)

multiome = sc.read_h5ad(
    multiome_selected_path,
)

print("\nSelected RNA:")
print(rna)

print("\nSelected Multiome:")
print(multiome)


# ============================================================
# Step 46 continued: Validate synchronized cell barcodes
# ============================================================

if not rna.obs_names.is_unique:
    raise ValueError(
        "RNA obs_names contain duplicated barcodes."
    )

if not multiome.obs_names.is_unique:
    raise ValueError(
        "Multiome obs_names contain duplicated barcodes."
    )

if not atac.obs_names.equals(rna.obs_names):
    raise ValueError(
        "RNA obs_names do not exactly match ATAC obs_names "
        "in identity and order."
    )

if not atac.obs_names.equals(multiome.obs_names):
    raise ValueError(
        "Multiome obs_names do not exactly match ATAC obs_names "
        "in identity and order."
    )

print("\nBarcode synchronization validation passed.")
print(f"ATAC cells: {atac.n_obs:,}")
print(f"RNA cells: {rna.n_obs:,}")
print(f"Multiome cells: {multiome.n_obs:,}")

print("\nNo RNA or Multiome filtering was performed yet.")



# ============================================================
# Step 47: Apply synchronized final filtering
# ============================================================

print("\nStep 47")
print("Applying synchronized filtering to RNA, ATAC, and Multiome.")

synchronized_mask = (
    final_atac_qc_mask
    .reindex(atac.obs_names)
)

if synchronized_mask.isna().any():
    raise ValueError(
        "Synchronized QC mask contains missing values."
    )

if not synchronized_mask.index.equals(atac.obs_names):
    raise ValueError(
        "Synchronized QC mask is not aligned with ATAC obs_names."
    )
    
    
rna_filtered = rna[
    synchronized_mask.to_numpy(dtype=bool)
].copy()

atac_filtered = atac[
    synchronized_mask.to_numpy(dtype=bool)
].copy()

multiome_filtered = multiome[
    synchronized_mask.to_numpy(dtype=bool)
].copy()




expected_n_cells = int(
    synchronized_mask.sum()
)

if rna_filtered.n_obs != expected_n_cells:
    raise RuntimeError(
        "Filtered RNA cell count does not match the QC mask."
    )

if atac_filtered.n_obs != expected_n_cells:
    raise RuntimeError(
        "Filtered ATAC cell count does not match the QC mask."
    )

if multiome_filtered.n_obs != expected_n_cells:
    raise RuntimeError(
        "Filtered Multiome cell count does not match the QC mask."
    )

if not rna_filtered.obs_names.equals(
    atac_filtered.obs_names
):
    raise RuntimeError(
        "Filtered RNA and ATAC barcodes are not identical."
    )

if not rna_filtered.obs_names.equals(
    multiome_filtered.obs_names
):
    raise RuntimeError(
        "Filtered RNA and Multiome barcodes are not identical."
    )
    
    
    
    
print("\nSynchronized filtering audit:")
print(f"Cells before filtering: {atac.n_obs:,}")
print(f"Cells after filtering: {expected_n_cells:,}")
print(f"Cells removed: {atac.n_obs - expected_n_cells:,}")

print("\nFiltered objects:")
print(f"RNA: {rna_filtered.n_obs:,} × {rna_filtered.n_vars:,}")
print(f"ATAC: {atac_filtered.n_obs:,} × {atac_filtered.n_vars:,}")
print(
    f"Multiome: "
    f"{multiome_filtered.n_obs:,} × {multiome_filtered.n_vars:,}"
)

print("\nBarcode synchronization after filtering passed.")
print("No h5ad files have been saved yet.")


# ============================================================
# Step 48: Save final filtered RNA, ATAC, and Multiome
# ============================================================

import anndata as ad

ad.settings.allow_write_nullable_strings = True

final_output_dir.mkdir(
    parents=True,
    exist_ok=True,
)

rna_final_path = (
    final_output_dir
    / "pbmc_multiome_rna_final_after_atac_qc.h5ad"
)

atac_final_path = (
    final_output_dir
    / "pbmc_multiome_atac_final_after_atac_qc.h5ad"
)

multiome_final_path = (
    final_output_dir
    / "pbmc_multiome_final_after_atac_qc.h5ad"
)

print("\nStep 48")
print("Saving filtered RNA...")
rna_filtered.write_h5ad(
    str(rna_final_path),
    compression="gzip",
)

print("Saving filtered ATAC...")
atac_filtered.write_h5ad(
    str(atac_final_path),
    compression="gzip",
)

print("Saving filtered Multiome...")
multiome_filtered.write_h5ad(
    str(multiome_final_path),
    compression="gzip",
)

print("\nAll three files were saved.")


# ============================================================
# Step 49: Reload and confirm identical filtered cells
# ============================================================

rna_check = sc.read_h5ad(
    rna_final_path,
    backed="r",
)

atac_check = sc.read_h5ad(
    atac_final_path,
    backed="r",
)

multiome_check = sc.read_h5ad(
    multiome_final_path,
    backed="r",
)

try:
    print("\nStep 49")
    print(f"RNA cells: {rna_check.n_obs:,}")
    print(f"ATAC cells: {atac_check.n_obs:,}")
    print(f"Multiome cells: {multiome_check.n_obs:,}")

    if not rna_check.obs_names.equals(atac_check.obs_names):
        raise RuntimeError(
            "RNA and ATAC cells are not identical or not in the same order."
        )

    if not rna_check.obs_names.equals(multiome_check.obs_names):
        raise RuntimeError(
            "RNA and Multiome cells are not identical or not in the same order."
        )

    print("\nAll three saved objects contain identical cells.")
    print("Cell order is also identical.")

finally:
    rna_check.file.close()
    atac_check.file.close()
    multiome_check.file.close()