from pathlib import Path
import sys

import pandas as pd
import scanpy as sc
import yaml


# ============================================================
# 1. Configuration
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
    script_path = Path("run_03R_apply_selected_cluster_mask.py")
    refactor_dir = find_refactored_root(Path.cwd())

src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

from pbmcgrn.paths import ProjectPaths


paths = ProjectPaths(
    refactor_root=refactor_dir,
)

expected_project_id = "pbmc_10k_multiome_v1"

if paths.project_id != expected_project_id:
    raise RuntimeError(
        "Wrong active project.\n"
        f"Expected: {expected_project_id}\n"
        f"Observed: {paths.project_id}"
    )

project_dir = paths.project_root

preprocessing_config_path = (
    paths.configs
    / "preprocessing_rna.yaml"
)

if not preprocessing_config_path.exists():
    raise FileNotFoundError(
        "preprocessing_rna.yaml not found:\n"
        f"{preprocessing_config_path}"
    )

with open(
    preprocessing_config_path,
    "r",
    encoding="utf-8",
) as handle:
    preprocessing_config = yaml.safe_load(handle) or {}

cluster_key = (
    preprocessing_config
    .get("dimensionality_reduction", {})
    .get("leiden_key")
)

if not cluster_key:
    raise KeyError(
        "Missing dimensionality_reduction.leiden_key in:\n"
        f"{preprocessing_config_path}"
    )

annotation_path = (
    project_dir
    / "reports"
    / "tables"
    / "rna"
    / "step03A"
    / "rna_informed_annotation_template.csv"
)

rna_path = (
    project_dir
    / "data"
    / "processed"
    / "rna"
    / "pbmc_multiome_rna_clustered.h5ad"
)

atac_path = (
    project_dir
    / "data"
    / "processed"
    / "atac"
    / "pbmc_multiome_atac_after_rna_qc.h5ad"
)

multiome_path = (
    project_dir
    / "data"
    / "processed"
    / "multiome"
    / "pbmc_multiome_after_rna_qc.h5ad"
)

output_dir = (
    project_dir
    / "data"
    / "processed"
    / "selected_for_scatac_analysis"
)

output_dir.mkdir(
    parents=True,
    exist_ok=True,
)

rna_output = (
    output_dir
    / "pbmc_multiome_rna_selected_for_scatac_analysis.h5ad"
)

atac_output = (
    output_dir
    / "pbmc_multiome_atac_selected_for_scatac_analysis.h5ad"
)

multiome_output = (
    output_dir
    / "pbmc_multiome_selected_for_scatac_analysis.h5ad"
)

mask_output = (
    project_dir
    / "reports"
    / "tables"
    / "rna"
    / "step03R"
    / "selected_cells_for_scatac_analysis_mask.csv"
)

mask_output.parent.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 2. Validate input files
# ============================================================

input_files = {
    "annotation": annotation_path,
    "RNA": rna_path,
    "ATAC": atac_path,
    "Multiome": multiome_path,
}

for label, path in input_files.items():
    if not path.exists():
        raise FileNotFoundError(
            f"{label} input file was not found:\n{path}"
        )


# ============================================================
# 3. Read annotation table
# ============================================================

def read_annotation_csv(path: Path) -> pd.DataFrame:
    """Read Excel-saved CSV robustly across UTF-8, Greek ANSI, and Windows encodings."""

    encodings_to_try = [
        "utf-8-sig",
        "utf-8",
        "cp1253",      # Greek Windows / Excel ANSI
        "windows-1253",
        "cp1252",      # Western Windows fallback
        "latin1",      # last-resort fallback; never fails but may misread symbols
    ]

    last_error = None

    for encoding in encodings_to_try:
        try:
            df = pd.read_csv(
                path,
                sep=None,
                engine="python",
                encoding=encoding,
                dtype=str,
            )

            print(f"\nAnnotation CSV read successfully with encoding: {encoding}")

            return df

        except UnicodeDecodeError as error:
            last_error = error
            print(f"Failed reading annotation CSV with encoding: {encoding}")

    raise UnicodeDecodeError(
        last_error.encoding,
        last_error.object,
        last_error.start,
        last_error.end,
        (
            "Could not read annotation CSV with common encodings. "
            "Open the file in Excel and save it as 'CSV UTF-8'."
        ),
    )


annotation = read_annotation_csv(annotation_path)

annotation.columns = (
    annotation.columns
    .astype(str)
    .str.replace("\ufeff", "", regex=False)
    .str.strip()
)

for column in annotation.columns:
    annotation[column] = (
        annotation[column]
        .fillna("")
        .astype(str)
        .str.strip()
    )

required_columns = {
    "cluster",
    "candidate_cell_type_broad",
    "include_in_downstream_provisional",
}

missing_columns = required_columns.difference(
    annotation.columns
)

if missing_columns:
    raise ValueError(
        "Annotation table is missing required columns:\n"
        f"{sorted(missing_columns)}\n\n"
        f"Available columns:\n{annotation.columns.tolist()}"
    )


# ============================================================
# 4. Validate manual inclusion decisions
# ============================================================

annotation["cluster"] = (
    annotation["cluster"]
    .astype(str)
    .str.strip()
)

annotation["include_in_downstream_provisional"] = (
    annotation["include_in_downstream_provisional"]
    .fillna("")
    .astype(str)
    .str.strip()
    .str.lower()
)

allowed_values = {
    "yes",
    "no",
    "review",
}

observed_values = set(
    annotation["include_in_downstream_provisional"]
)

invalid_values = (
    observed_values
    - allowed_values
    - {""}
)

if invalid_values:
    raise ValueError(
        "Invalid values in "
        "'include_in_downstream_provisional':\n"
        f"{sorted(invalid_values)}\n\n"
        "Use only: yes, no, or review."
    )

empty_decisions = annotation.loc[
    annotation[
        "include_in_downstream_provisional"
    ].eq(""),
    "cluster",
].tolist()

if empty_decisions:
    raise ValueError(
        "No downstream decision was provided for clusters:\n"
        f"{empty_decisions}\n\n"
        "Fill each cluster with yes, no, or review."
    )

duplicated_clusters = annotation.loc[
    annotation["cluster"].duplicated(keep=False),
    "cluster",
].unique().tolist()

if duplicated_clusters:
    raise ValueError(
        "Duplicate annotation rows were found for clusters:\n"
        f"{duplicated_clusters}"
    )


# ============================================================
# 5. Select clusters marked yes
# ============================================================

selected_annotation = annotation.loc[
    annotation[
        "include_in_downstream_provisional"
    ].eq("yes")
].copy()

selected_clusters = set(
    selected_annotation["cluster"].astype(str)
)

if not selected_clusters:
    raise ValueError(
        "No clusters were marked as 'yes'."
    )

print("\nSelected clusters:")
print(
    selected_annotation[
        [
            "cluster",
            "candidate_cell_type_broad",
            "candidate_cell_type_detailed",
            "annotation_confidence",
            "include_in_downstream_provisional",
        ]
    ].to_string(index=False)
)


# ============================================================
# 6. Load paired objects
# ============================================================

print("\nLoading RNA...")
rna = sc.read_h5ad(rna_path)

print("Loading ATAC...")
atac = sc.read_h5ad(atac_path)

print("Loading Multiome...")
multiome = sc.read_h5ad(multiome_path)


# ============================================================
# 7. Validate paired barcodes before filtering
# ============================================================

print("\nCell counts before selection:")
print(f"RNA:      {rna.n_obs}")
print(f"ATAC:     {atac.n_obs}")
print(f"Multiome: {multiome.n_obs}")

if rna.n_obs != atac.n_obs:
    raise ValueError(
        "RNA and ATAC have different cell counts before masking."
    )

if rna.n_obs != multiome.n_obs:
    raise ValueError(
        "RNA and Multiome have different cell counts "
        "before masking."
    )

if not rna.obs_names.equals(atac.obs_names):
    raise ValueError(
        "RNA and ATAC obs_names are not identical "
        "and in the same order."
    )

if not rna.obs_names.equals(multiome.obs_names):
    raise ValueError(
        "RNA and Multiome obs_names are not identical "
        "and in the same order."
    )

if cluster_key not in rna.obs.columns:
    raise ValueError(
        f"Cluster key '{cluster_key}' is missing from RNA obs."
    )


# ============================================================
# 8. Validate selected clusters against RNA
# ============================================================

rna_clusters = set(
    rna.obs[cluster_key]
    .astype(str)
    .unique()
)

unknown_selected_clusters = (
    selected_clusters
    - rna_clusters
)

if unknown_selected_clusters:
    raise ValueError(
        "The annotation table selected clusters "
        "that do not exist in the RNA object:\n"
        f"{sorted(unknown_selected_clusters)}"
    )


# ============================================================
# 9. Build the shared cell mask
# ============================================================

rna_cluster_values = (
    rna.obs[cluster_key]
    .astype(str)
)

selected_mask = (
    rna_cluster_values
    .isin(selected_clusters)
)

n_selected = int(selected_mask.sum())
n_removed = int((~selected_mask).sum())

if n_selected == 0:
    raise ValueError(
        "The selected cluster mask contains zero cells."
    )

print("\nSelection summary:")
print(f"Selected cells: {n_selected}")
print(f"Removed cells:  {n_removed}")
print(
    f"Selected fraction: "
    f"{n_selected / rna.n_obs:.3f}"
)


# ============================================================
# 10. Save reusable barcode mask
# ============================================================

broad_mapping = (
    annotation
    .set_index("cluster")[
        "candidate_cell_type_broad"
    ]
    .to_dict()
)

detailed_mapping = (
    annotation
    .set_index("cluster")[
        "candidate_cell_type_detailed"
    ]
    .to_dict()
)

mask_table = pd.DataFrame(
    {
        "barcode": rna.obs_names.astype(str),
        "cluster": rna_cluster_values.to_numpy(),
        "candidate_cell_type_broad": (
            rna_cluster_values
            .map(broad_mapping)
            .to_numpy()
        ),
        "candidate_cell_type_detailed": (
            rna_cluster_values
            .map(detailed_mapping)
            .to_numpy()
        ),
        "selected_for_scatac_analysis": (
            selected_mask.to_numpy()
        ),
    }
)

mask_table.to_csv(
    mask_output,
    index=False,
    encoding="utf-8-sig",
)

print("\nSaved reusable barcode mask:")
print(mask_output)


# ============================================================
# 11. Apply exactly the same mask to all modalities
# ============================================================

rna_selected = rna[selected_mask.to_numpy()].copy()
atac_selected = atac[selected_mask.to_numpy()].copy()
multiome_selected = multiome[
    selected_mask.to_numpy()
].copy()


# ============================================================
# 12. Transfer RNA annotations to all selected objects
# ============================================================

selected_cluster_series = (
    rna_selected.obs[cluster_key]
    .astype(str)
)

rna_selected.obs[
    "cell_type_broad"
] = selected_cluster_series.map(
    broad_mapping
)

rna_selected.obs[
    "cell_type_detailed"
] = selected_cluster_series.map(
    detailed_mapping
)

for obj in [
    atac_selected,
    multiome_selected,
]:
    obj.obs[cluster_key] = (
        rna_selected.obs[cluster_key]
        .astype(str)
        .to_numpy()
    )

    obj.obs["cell_type_broad"] = (
        rna_selected.obs["cell_type_broad"]
        .astype(str)
        .to_numpy()
    )

    obj.obs["cell_type_detailed"] = (
        rna_selected.obs["cell_type_detailed"]
        .astype(str)
        .to_numpy()
    )

    obj.obs["selected_for_scatac_analysis"] = True

rna_selected.obs["selected_for_scatac_analysis"] = True


# ============================================================
# 13. Validate after masking
# ============================================================

if not rna_selected.obs_names.equals(
    atac_selected.obs_names
):
    raise ValueError(
        "Selected RNA and ATAC barcodes differ."
    )

if not rna_selected.obs_names.equals(
    multiome_selected.obs_names
):
    raise ValueError(
        "Selected RNA and Multiome barcodes differ."
    )

if not (
    rna_selected.n_obs
    == atac_selected.n_obs
    == multiome_selected.n_obs
):
    raise ValueError(
        "Selected objects have different cell counts."
    )

print("\nCell counts after selection:")
print(f"RNA:      {rna_selected.n_obs}")
print(f"ATAC:     {atac_selected.n_obs}")
print(f"Multiome: {multiome_selected.n_obs}")

print("\nSelected cells per cluster:")
print(
    rna_selected.obs[cluster_key]
    .astype(str)
    .value_counts()
    .sort_index()
)

print("\nSelected cells per broad cell type:")
print(
    rna_selected.obs["cell_type_broad"]
    .value_counts()
)


# ============================================================
# 14. Save selected copies
# ============================================================

rna_selected.write_h5ad(
    rna_output,
    compression="gzip",
)

atac_selected.write_h5ad(
    atac_output,
    compression="gzip",
)

multiome_selected.write_h5ad(
    multiome_output,
    compression="gzip",
)

print("\nSaved selected RNA:")
print(rna_output)

print("\nSaved selected ATAC:")
print(atac_output)

print("\nSaved selected Multiome:")
print(multiome_output)

print("\nStep 03R completed successfully.")
