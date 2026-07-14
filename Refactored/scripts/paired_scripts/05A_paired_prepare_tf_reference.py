# ============================================================
# 05A_paired_prepare_human_tf_reference.py
#
# Purpose
# -------
# 1. Validate the active paired PBMC Multiome project.
# 2. Validate the official selected RNA, ATAC and Multiome inputs.
# 3. Confirm that all three objects contain the same cells
#    in exactly the same order.
# 4. Download and validate the curated human TF catalogue.
# 5. Match curated human TF symbols against the RNA gene universe.
# 6. Create a stable, ordered TF node index for downstream GRN work.
#
# This script does NOT:
# - define promoters or TSS coordinates,
# - download JASPAR motifs,
# - scan ATAC peaks,
# - create TF-TF adjacency matrices,
# - train a neural network.
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
from urllib.request import Request, urlopen

import anndata as ad
import pandas as pd
import yaml




# ============================================================
# 2. Resolve refactor repository
# ============================================================

script_path = Path(__file__).resolve()

# Script location:
# Refactored/scripts/paired_scripts/05A_....py
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
# 5. Resolve official final downstream inputs
# ============================================================

final_filtered_dir = (
    paths.project_root
    / "data"
    / "processed"
    / "final_after_atac_qc"
)

rna_filtered_path = (
    final_filtered_dir
    / "pbmc_multiome_rna_final_after_atac_qc.h5ad"
)

atac_filtered_path = (
    final_filtered_dir
    / "pbmc_multiome_atac_final_after_atac_qc.h5ad"
)

multiome_filtered_path = (
    final_filtered_dir
    / "pbmc_multiome_final_after_atac_qc.h5ad"
)

expected_n_cells = 7_566
expected_rna_n_vars = 3_000
expected_atac_n_vars = 143_887
expected_multiome_n_vars = 180_488

# ============================================================
# 6. Resolve reference and output directories
# ============================================================

tf_reference_dir = (
    paths.project_root
    / "data"
    / "external"
    / "references"
    / "transcription_factors"
    / "human_tfs_lambert_v1.01"
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
    tf_reference_dir,
    output_table_dir,
    output_summary_dir,
):
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )


# ============================================================
# 7. Define official TF reference URLs
# ============================================================

tf_names_url = (
    "https://humantfs.ccbr.utoronto.ca/"
    "download/v_1.01/TF_names_v_1.01.txt"
)

tf_database_url = (
    "https://humantfs.ccbr.utoronto.ca/"
    "download/v_1.01/DatabaseExtract_v_1.01.csv"
)

tf_names_path = (
    tf_reference_dir
    / "TF_names_v_1.01.txt"
)

tf_database_path = (
    tf_reference_dir
    / "DatabaseExtract_v_1.01.csv"
)


# ============================================================
# 8. Define output paths
# ============================================================

input_audit_path = (
    output_table_dir
    / "05A_paired_input_audit.tsv"
)

tf_catalogue_clean_path = (
    output_table_dir
    / "05A_human_tf_catalogue_clean.tsv"
)

tf_rna_mapping_path = (
    output_table_dir
    / "05A_tf_to_rna_mapping.tsv"
)

tf_node_index_path = (
    output_table_dir
    / "05A_tf_node_index.tsv"
)

tf_exclusion_report_path = (
    output_table_dir
    / "05A_tf_exclusion_report.tsv"
)

summary_yaml_path = (
    output_summary_dir
    / "05A_tf_reference_summary.yaml"
)

run_metadata_path = (
    output_summary_dir
    / "05A_run_metadata.json"
)


# ============================================================
# 9. Generic validation helpers
# ============================================================

def validate_file(
    path: Path,
    label: str,
    allow_empty: bool = False,
) -> None:
    """
    Validate that a required file exists and is non-empty.
    """
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
    """
    Calculate a SHA-256 checksum without loading
    the entire file into memory.
    """
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)

            if not chunk:
                break

            digest.update(chunk)

    return digest.hexdigest()


def utc_timestamp() -> str:
    """
    Return an ISO-8601 UTC timestamp.
    """
    return datetime.now(timezone.utc).isoformat()




# ============================================================
# 10. Validate official project inputs
# ============================================================

required_project_files = {
    "Final filtered RNA H5AD": rna_filtered_path,
    "Final filtered ATAC H5AD": atac_filtered_path,
    "Final filtered Multiome H5AD": multiome_filtered_path,
}

print("\nValidating project inputs:")

for label, path in required_project_files.items():
    validate_file(
        path=path,
        label=label,
    )

    print(
        f"- {label}: OK "
        f"({path.stat().st_size:,} bytes)"
    )

# ============================================================
# 11. Lightweight AnnData loading
# ============================================================

rna = ad.read_h5ad(
    rna_filtered_path,
    backed="r",
)

atac = ad.read_h5ad(
    atac_filtered_path,
    backed="r",
)

multiome = ad.read_h5ad(
    multiome_filtered_path,
    backed="r",
)
# ============================================================
# 12. Audit object dimensions and identifiers
# ============================================================

def audit_anndata(
    adata: ad.AnnData,
    object_name: str,
    expected_obs: int,
    expected_vars: int,
) -> dict[str, object]:
    """
    Audit one AnnData object without modifying it.
    """
    observed_obs = int(adata.n_obs)
    observed_vars = int(adata.n_vars)

    obs_names_unique = bool(adata.obs_names.is_unique)
    var_names_unique = bool(adata.var_names.is_unique)

    if observed_obs != expected_obs:
        raise ValueError(
            f"{object_name}: unexpected number of cells.\n"
            f"Expected: {expected_obs:,}\n"
            f"Observed: {observed_obs:,}"
        )

    if observed_vars != expected_vars:
        raise ValueError(
            f"{object_name}: unexpected number of features.\n"
            f"Expected: {expected_vars:,}\n"
            f"Observed: {observed_vars:,}"
        )

    if not obs_names_unique:
        raise ValueError(
            f"{object_name}: cell identifiers are not unique."
        )

    if not var_names_unique:
        raise ValueError(
            f"{object_name}: feature identifiers are not unique."
        )

    return {
        "object": object_name,
        "n_obs": observed_obs,
        "n_vars": observed_vars,
        "obs_names_unique": obs_names_unique,
        "var_names_unique": var_names_unique,
        "backed_mode": bool(adata.isbacked),
    }


input_audit_records = [
    audit_anndata(
        adata=rna,
        object_name="RNA",
        expected_obs=expected_n_cells,
        expected_vars=expected_rna_n_vars,
    ),
    audit_anndata(
        adata=atac,
        object_name="ATAC",
        expected_obs=expected_n_cells,
        expected_vars=expected_atac_n_vars,
    ),
    audit_anndata(
        adata=multiome,
        object_name="Multiome",
        expected_obs=expected_n_cells,
        expected_vars=expected_multiome_n_vars,
    ),
]

input_audit_df = pd.DataFrame(
    input_audit_records
)

print("\nInput object audit:")
print(
    input_audit_df.to_string(
        index=False,
    )
)


# ============================================================
# 13. Validate paired-cell synchronization
# ============================================================

rna_cells = rna.obs_names.astype(str)
atac_cells = atac.obs_names.astype(str)
multiome_cells = multiome.obs_names.astype(str)

rna_atac_same_order = rna_cells.equals(
    atac_cells
)

rna_multiome_same_order = rna_cells.equals(
    multiome_cells
)

atac_multiome_same_order = atac_cells.equals(
    multiome_cells
)

if not rna_atac_same_order:
    raise ValueError(
        "RNA and ATAC cell identifiers are not identical "
        "in the same order."
    )

if not rna_multiome_same_order:
    raise ValueError(
        "RNA and Multiome cell identifiers are not identical "
        "in the same order."
    )

if not atac_multiome_same_order:
    raise ValueError(
        "ATAC and Multiome cell identifiers are not identical "
        "in the same order."
    )

print(
    "\nPaired-cell synchronization validated:"
)
print(
    f"- Identical cells: {len(rna_cells):,}"
)
print(
    "- RNA, ATAC and Multiome cell order: identical"
)


# ============================================================
# 14. Inspect RNA gene identifiers
# ============================================================

rna_var_names = pd.Index(
    rna.var_names.astype(str),
    name="rna_var_name",
)

print("\nRNA var columns:")

if len(rna.var.columns) == 0:
    print("- No additional RNA var columns")
else:
    for column in rna.var.columns:
        print(f"- {column}")

print("\nFirst 20 RNA var_names:")
print(rna_var_names[:20].tolist())


# ============================================================
# 15. Identify possible RNA gene-symbol columns
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

observed_symbol_columns = [
    column
    for column in candidate_gene_symbol_columns
    if column in rna.var.columns
]

print("\nCandidate RNA gene-symbol columns found:")
print(
    observed_symbol_columns
    if observed_symbol_columns
    else "None"
)


# ============================================================
# 16. Determine the RNA gene-symbol universe
# ============================================================

def looks_like_ensembl_gene_id(
    values: pd.Index,
    sample_size: int = 500,
) -> bool:
    """
    Heuristic test for Ensembl-style RNA var_names.
    """
    sample = (
        pd.Series(values[:sample_size], dtype="string")
        .dropna()
        .str.upper()
    )

    if sample.empty:
        return False

    fraction_ensembl = (
        sample
        .str.match(r"^ENSG\d+(\.\d+)?$")
        .mean()
    )

    return bool(fraction_ensembl >= 0.80)


rna_var_names_are_ensembl = looks_like_ensembl_gene_id(
    rna_var_names
)

if observed_symbol_columns:
    selected_gene_symbol_column = (
        observed_symbol_columns[0]
    )

    rna_gene_symbols = (
        rna.var[selected_gene_symbol_column]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    rna_gene_symbol_source = (
        f"rna.var[{selected_gene_symbol_column!r}]"
    )

elif not rna_var_names_are_ensembl:
    selected_gene_symbol_column = None

    rna_gene_symbols = pd.Series(
        rna_var_names,
        index=rna_var_names,
        dtype="string",
    ).str.strip().str.upper()

    rna_gene_symbol_source = "rna.var_names"

else:
    raise ValueError(
        "RNA var_names appear to contain Ensembl gene IDs, "
        "but no recognized gene-symbol column was found in rna.var.\n"
        "A valid Ensembl-to-HGNC mapping is required before TF matching.\n"
        f"Available RNA var columns: {list(rna.var.columns)}"
    )


rna_gene_table = pd.DataFrame(
    {
        "rna_var_name": rna_var_names.astype(str),
        "rna_gene_symbol": rna_gene_symbols.to_numpy(),
    }
)

rna_gene_table["rna_gene_symbol"] = (
    rna_gene_table["rna_gene_symbol"]
    .astype("string")
    .str.strip()
    .str.upper()
)

rna_gene_table["rna_gene_symbol_missing"] = (
    rna_gene_table["rna_gene_symbol"].isna()
    | rna_gene_table["rna_gene_symbol"].eq("")
)

rna_gene_table["rna_gene_symbol_duplicated"] = (
    rna_gene_table["rna_gene_symbol"]
    .duplicated(
        keep=False,
    )
)

n_missing_rna_symbols = int(
    rna_gene_table["rna_gene_symbol_missing"].sum()
)

n_duplicated_rna_symbols = int(
    rna_gene_table.loc[
        ~rna_gene_table["rna_gene_symbol_missing"],
        "rna_gene_symbol",
    ]
    .duplicated()
    .sum()
)

print("\nRNA gene-symbol audit:")
print(f"- Source: {rna_gene_symbol_source}")
print(f"- Total RNA features: {len(rna_gene_table):,}")
print(f"- Missing gene symbols: {n_missing_rna_symbols:,}")
print(
    "- Additional duplicated symbol entries: "
    f"{n_duplicated_rna_symbols:,}"
)

if n_missing_rna_symbols > 0:
    raise ValueError(
        "Missing RNA gene symbols were detected. "
        "Resolve them before creating the TF node universe."
    )

if n_duplicated_rna_symbols > 0:
    duplicated_symbols = (
        rna_gene_table.loc[
            rna_gene_table[
                "rna_gene_symbol_duplicated"
            ],
            "rna_gene_symbol",
        ]
        .dropna()
        .unique()
        .tolist()
    )

    raise ValueError(
        "Duplicated RNA gene symbols were detected.\n"
        "The TF-to-RNA mapping would not be one-to-one.\n"
        f"Examples: {duplicated_symbols[:20]}"
    )


# ============================================================
# 17. Download helper for external references
# ============================================================

def download_reference_file(
    url: str,
    output_path: Path,
    overwrite: bool = False,
    timeout_seconds: int = 120,
) -> None:
    """
    Download an external reference file safely.

    Existing non-empty files are retained unless overwrite=True.
    """
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if output_path.exists() and not overwrite:
        validate_file(
            path=output_path,
            label="Existing reference file",
        )

        print("\nReference file already exists; download skipped:")
        print(output_path)
        print(
            f"File size: {output_path.stat().st_size:,} bytes"
        )
        return

    temporary_path = output_path.with_suffix(
        output_path.suffix + ".part"
    )

    if temporary_path.exists():
        temporary_path.unlink()

    request = Request(
        url,
        headers={
            "User-Agent": (
                "pbmcgrn-reference-downloader/1.0"
            )
        },
    )

    print("\nDownloading reference:")
    print(url)

    print("Destination:")
    print(output_path)

    try:
        with urlopen(
            request,
            timeout=timeout_seconds,
        ) as response:
            content = response.read()

        if len(content) == 0:
            raise ValueError(
                "The remote server returned an empty file."
            )

        temporary_path.write_bytes(content)

        validate_file(
            path=temporary_path,
            label="Temporary downloaded reference",
        )

        temporary_path.replace(
            output_path
        )

    except Exception as exc:
        if temporary_path.exists():
            temporary_path.unlink()

        raise RuntimeError(
            "Failed to download the external TF reference.\n"
            f"URL: {url}\n"
            f"Destination: {output_path}"
        ) from exc

    validate_file(
        path=output_path,
        label="Downloaded TF reference",
    )

    print("Download completed successfully.")
    print(
        f"File size: {output_path.stat().st_size:,} bytes"
    )


# ============================================================
# 18. Download official human TF references
# ============================================================

download_reference_file(
    url=tf_names_url,
    output_path=tf_names_path,
    overwrite=False,
)

# The complete database is downloaded for provenance and
# future annotation. The simple official TF-name list is used
# as the authoritative node-catalogue input in this script.
download_reference_file(
    url=tf_database_url,
    output_path=tf_database_path,
    overwrite=False,
)


# ============================================================
# 19. Validate downloaded reference files
# ============================================================

validate_file(
    path=tf_names_path,
    label="Official human TF-name list",
)

validate_file(
    path=tf_database_path,
    label="Official complete human TF database",
)

tf_names_sha256 = sha256_file(
    tf_names_path
)

tf_database_sha256 = sha256_file(
    tf_database_path
)

print("\nReference checksums:")
print(f"- TF names SHA-256: {tf_names_sha256}")
print(f"- TF database SHA-256: {tf_database_sha256}")


# ============================================================
# 20. Load and clean official TF names
# ============================================================

tf_names_raw = pd.read_csv(
    tf_names_path,
    header=None,
    names=["tf_symbol"],
    dtype="string",
)

tf_catalogue = tf_names_raw.copy()

tf_catalogue["tf_symbol_original"] = (
    tf_catalogue["tf_symbol"]
)

tf_catalogue["tf_symbol"] = (
    tf_catalogue["tf_symbol"]
    .astype("string")
    .str.strip()
    .str.upper()
)

tf_catalogue["missing_tf_symbol"] = (
    tf_catalogue["tf_symbol"].isna()
    | tf_catalogue["tf_symbol"].eq("")
)

tf_catalogue["duplicated_tf_symbol"] = (
    tf_catalogue["tf_symbol"]
    .duplicated(
        keep=False,
    )
)

if tf_catalogue["missing_tf_symbol"].any():
    invalid_rows = tf_catalogue.loc[
        tf_catalogue["missing_tf_symbol"]
    ]

    raise ValueError(
        "The official TF catalogue contains missing symbols.\n"
        f"Invalid rows:\n{invalid_rows.head(20)}"
    )

if tf_catalogue["duplicated_tf_symbol"].any():
    duplicated_tfs = (
        tf_catalogue.loc[
            tf_catalogue["duplicated_tf_symbol"],
            "tf_symbol",
        ]
        .unique()
        .tolist()
    )

    raise ValueError(
        "The official TF catalogue contains duplicated symbols.\n"
        f"Examples: {duplicated_tfs[:30]}"
    )

tf_catalogue = (
    tf_catalogue
    .drop(
        columns=[
            "missing_tf_symbol",
            "duplicated_tf_symbol",
        ]
    )
    .sort_values(
        "tf_symbol"
    )
    .reset_index(
        drop=True
    )
)

tf_catalogue["catalogue_index"] = (
    range(len(tf_catalogue))
)

tf_catalogue["reference_name"] = (
    "Lambert_human_TFs"
)

tf_catalogue["reference_version"] = (
    "v1.01"
)

print("\nOfficial human TF catalogue:")
print(f"- Number of curated TFs: {len(tf_catalogue):,}")
print("\nFirst 20 TF symbols:")
print(
    tf_catalogue[
        ["catalogue_index", "tf_symbol"]
    ]
    .head(20)
    .to_string(index=False)
)


# ============================================================
# 21. Parse the complete database for audit only
# ============================================================

tf_database = pd.read_csv(
    tf_database_path,
    low_memory=False,
)

if tf_database.empty:
    raise ValueError(
        "The complete human TF database is empty."
    )

print("\nComplete TF database:")
print(f"- Shape: {tf_database.shape}")
print("- Columns:")

for column in tf_database.columns:
    print(f"  - {column}")


# ============================================================
# 22. Match official TFs to the RNA gene universe
# ============================================================

tf_rna_mapping = tf_catalogue.merge(
    rna_gene_table[
        [
            "rna_var_name",
            "rna_gene_symbol",
        ]
    ],
    how="left",
    left_on="tf_symbol",
    right_on="rna_gene_symbol",
    validate="one_to_one",
)

tf_rna_mapping["present_in_rna"] = (
    tf_rna_mapping["rna_var_name"].notna()
)

tf_rna_mapping["included_as_tf_node"] = (
    tf_rna_mapping["present_in_rna"]
)

tf_rna_mapping["exclusion_reason"] = pd.NA

tf_rna_mapping.loc[
    ~tf_rna_mapping["present_in_rna"],
    "exclusion_reason",
] = "not_present_in_current_rna_gene_universe"


n_catalogue_tfs = len(
    tf_rna_mapping
)

n_tfs_present_in_rna = int(
    tf_rna_mapping["present_in_rna"].sum()
)

n_tfs_absent_from_rna = int(
    (~tf_rna_mapping["present_in_rna"]).sum()
)

tf_retention_fraction = (
    n_tfs_present_in_rna
    / n_catalogue_tfs
)

print("\nTF-to-RNA mapping:")
print(f"- Catalogue TFs: {n_catalogue_tfs:,}")
print(
    f"- TFs present in RNA: "
    f"{n_tfs_present_in_rna:,}"
)
print(
    f"- TFs absent from RNA: "
    f"{n_tfs_absent_from_rna:,}"
)
print(
    f"- TF retention: "
    f"{tf_retention_fraction:.2%}"
)


# ============================================================
# 23. Create stable ordered TF node index
# ============================================================

# Alphabetical order is deterministic and independent of the
# original ordering of the downloaded catalogue or RNA object.
tf_node_index = (
    tf_rna_mapping.loc[
        tf_rna_mapping["included_as_tf_node"],
        [
            "tf_symbol",
            "rna_var_name",
            "catalogue_index",
            "reference_name",
            "reference_version",
        ],
    ]
    .sort_values(
        "tf_symbol"
    )
    .reset_index(
        drop=True
    )
)

tf_node_index.insert(
    0,
    "tf_index",
    range(len(tf_node_index)),
)

if tf_node_index.empty:
    raise ValueError(
        "No curated TFs matched the RNA gene universe."
    )

if tf_node_index["tf_symbol"].duplicated().any():
    raise ValueError(
        "Duplicated TF symbols were found in the final node index."
    )

if tf_node_index["rna_var_name"].duplicated().any():
    raise ValueError(
        "Duplicated RNA feature identifiers were found "
        "in the final TF node index."
    )

expected_tf_indices = list(
    range(len(tf_node_index))
)

observed_tf_indices = (
    tf_node_index["tf_index"]
    .astype(int)
    .tolist()
)

if observed_tf_indices != expected_tf_indices:
    raise RuntimeError(
        "The final tf_index is not contiguous from 0 to N-1."
    )

print("\nFinal preliminary TF node index:")
print(f"- Number of TF nodes: {len(tf_node_index):,}")
print("\nFirst 20 TF nodes:")
print(
    tf_node_index.head(20).to_string(
        index=False
    )
)


# ============================================================
# 24. Build exclusion report
# ============================================================

tf_exclusion_report = (
    tf_rna_mapping.loc[
        ~tf_rna_mapping["included_as_tf_node"],
        [
            "catalogue_index",
            "tf_symbol",
            "present_in_rna",
            "exclusion_reason",
            "reference_name",
            "reference_version",
        ],
    ]
    .sort_values(
        [
            "exclusion_reason",
            "tf_symbol",
        ]
    )
    .reset_index(
        drop=True
    )
)

print("\nTF exclusion report:")
print(
    tf_exclusion_report[
        "exclusion_reason"
    ]
    .value_counts(
        dropna=False
    )
    .to_string()
)


# ============================================================
# 25. Save tabular outputs
# ============================================================

input_audit_df.to_csv(
    input_audit_path,
    sep="\t",
    index=False,
)

tf_catalogue.to_csv(
    tf_catalogue_clean_path,
    sep="\t",
    index=False,
)

tf_rna_mapping.to_csv(
    tf_rna_mapping_path,
    sep="\t",
    index=False,
)

tf_node_index.to_csv(
    tf_node_index_path,
    sep="\t",
    index=False,
)

tf_exclusion_report.to_csv(
    tf_exclusion_report_path,
    sep="\t",
    index=False,
)


# ============================================================
# 26. Create YAML summary
# ============================================================

summary = {
    "step": "05A",
    "script_name": (
        "05A_paired_prepare_human_tf_reference.py"
    ),
    "created_at_utc": utc_timestamp(),
    "project": {
        "project_id": paths.project_id,
        "project_root": str(paths.project_root),
    },
    "inputs": {
        "rna_filtered_path": str(
            rna_filtered_path
        ),
        "atac_filtered_path": str(
            atac_filtered_path
        ),
        "multiome_filtered_path": str(
            multiome_filtered_path
        ),
},
    "expected_dimensions": {
        "n_cells": expected_n_cells,
        "rna_n_vars": expected_rna_n_vars,
        "atac_n_vars": expected_atac_n_vars,
        "multiome_n_vars": expected_multiome_n_vars,
    },
    "observed_dimensions": {
        "rna": {
            "n_obs": int(rna.n_obs),
            "n_vars": int(rna.n_vars),
        },
        "atac": {
            "n_obs": int(atac.n_obs),
            "n_vars": int(atac.n_vars),
        },
        "multiome": {
            "n_obs": int(multiome.n_obs),
            "n_vars": int(multiome.n_vars),
        },
    },
    "paired_cell_validation": {
        "rna_atac_same_order": (
            rna_atac_same_order
        ),
        "rna_multiome_same_order": (
            rna_multiome_same_order
        ),
        "atac_multiome_same_order": (
            atac_multiome_same_order
        ),
    },
    "rna_gene_universe": {
        "source": rna_gene_symbol_source,
        "var_names_look_like_ensembl": (
            rna_var_names_are_ensembl
        ),
        "n_rna_features": int(
            len(rna_gene_table)
        ),
        "n_missing_gene_symbols": (
            n_missing_rna_symbols
        ),
        "n_additional_duplicated_symbols": (
            n_duplicated_rna_symbols
        ),
    },
    "tf_reference": {
        "name": "Lambert human TF catalogue",
        "version": "v1.01",
        "tf_names_url": tf_names_url,
        "tf_database_url": tf_database_url,
        "tf_names_path": str(tf_names_path),
        "tf_database_path": str(
            tf_database_path
        ),
        "tf_names_sha256": tf_names_sha256,
        "tf_database_sha256": (
            tf_database_sha256
        ),
    },
    "tf_mapping": {
        "n_catalogue_tfs": n_catalogue_tfs,
        "n_tfs_present_in_rna": (
            n_tfs_present_in_rna
        ),
        "n_tfs_absent_from_rna": (
            n_tfs_absent_from_rna
        ),
        "tf_retention_fraction": float(
            tf_retention_fraction
        ),
        "n_preliminary_tf_nodes": int(
            len(tf_node_index)
        ),
    },
    "outputs": {
        "input_audit": str(
            input_audit_path
        ),
        "tf_catalogue_clean": str(
            tf_catalogue_clean_path
        ),
        "tf_rna_mapping": str(
            tf_rna_mapping_path
        ),
        "tf_node_index": str(
            tf_node_index_path
        ),
        "tf_exclusion_report": str(
            tf_exclusion_report_path
        ),
        "summary_yaml": str(
            summary_yaml_path
        ),
    },
    "important_note": (
        "The current TF node list is preliminary. "
        "It has been filtered by presence in the current "
        "RNA gene universe but has not yet been filtered "
        "using GENCODE annotation, TSS availability or "
        "JASPAR motif availability."
    ),
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
# 27. Save run metadata
# ============================================================

run_metadata = {
    "created_at_utc": utc_timestamp(),
    "python_version": sys.version,
    "pandas_version": pd.__version__,
    "anndata_version": ad.__version__,
    "project_id": paths.project_id,
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
# 28. Validate written outputs
# ============================================================

required_outputs = {
    "Input audit": input_audit_path,
    "Clean TF catalogue": tf_catalogue_clean_path,
    "TF-to-RNA mapping": tf_rna_mapping_path,
    "TF node index": tf_node_index_path,
    "TF exclusion report": tf_exclusion_report_path,
    "Summary YAML": summary_yaml_path,
    "Run metadata": run_metadata_path,
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
# 29. Close backed AnnData objects
# ============================================================

for adata in (
    rna,
    atac,
    multiome,
):
    if adata.isbacked:
        adata.file.close()


# ============================================================
# 30. Final report
# ============================================================

print("\n" + "=" * 60)
print("05A COMPLETED SUCCESSFULLY")
print("=" * 60)

print("\nCurated TF catalogue:")
print(f"{n_catalogue_tfs:,} TFs")

print("\nTFs present in current RNA gene universe:")
print(f"{n_tfs_present_in_rna:,}")

print("\nTFs excluded because they are absent from RNA:")
print(f"{n_tfs_absent_from_rna:,}")

print("\nPreliminary TF nodes:")
print(f"{len(tf_node_index):,}")

print("\nPrimary output:")
print(tf_node_index_path)

print("\nImportant:")
print(
    "This is a preliminary TF node universe. "
    "The next script will add genomic annotation, "
    "TSS coordinates and promoter intervals."
)