"""
Run 01 — Load and validate paired 10x Multiome PBMC data.

Workflow
--------
1. Read the active project from project_registry.yaml.
2. Load one 10x Multiome H5 containing Gene Expression and Peaks.
3. Validate the raw input and feature types.
4. Split the shared matrix into RNA and ATAC AnnData objects.
5. Confirm that RNA and ATAC contain exactly the same cell barcodes.
6. Save project-local raw H5AD objects, validation tables and manifest.

Raw inputs may be outside the active project.
All outputs must remain inside the active project.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc


# ============================================================
# 0. Locate Refactored root
# ============================================================

# Script location:
# Refactored/scripts/paired_scripts/run_01_paired_load_and_validate_multiome.py
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


# ============================================================
# 2. Helper functions
# ============================================================

def hash_barcodes(barcodes: pd.Index) -> str:
    """
    Return an order-sensitive SHA256 hash for a barcode index.

    Identical hashes mean that the barcode strings and their order
    are identical.
    """
    text = "\n".join(barcodes.astype(str))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def require_columns(
    dataframe: pd.DataFrame,
    required_columns: set[str],
    *,
    dataframe_name: str,
) -> None:
    """Raise an error when required columns are missing."""
    missing = required_columns.difference(dataframe.columns)

    if missing:
        raise KeyError(
            f"{dataframe_name} is missing required columns: "
            f"{sorted(missing)}"
        )


def confirm_output_inside_project(
    output_path: Path,
    project_root: Path,
    *,
    label: str,
) -> None:
    """Additional explicit safety check for project-local outputs."""
    output_path = output_path.resolve()
    project_root = project_root.resolve()

    try:
        output_path.relative_to(project_root)
    except ValueError as exc:
        raise ValueError(
            f"{label} is outside the active project.\n"
            f"Output: {output_path}\n"
            f"Project root: {project_root}"
        ) from exc

print("Script path:", script_path)
print("Refactored root:", refactor_dir)

registry_path = refactor_dir / "project_registry.yaml"

print("Registry path:", registry_path)
print("Registry exists:", registry_path.exists())

if registry_path.exists():
    print("\nRegistry contents:")
    print(registry_path.read_text(encoding="utf-8"))
# ============================================================
# 3. Main workflow
# ============================================================

def main() -> None:
    # --------------------------------------------------------
    # 3.1 Initialize active project
    # --------------------------------------------------------

    paths = ProjectPaths(refactor_root=refactor_dir)
    paths.ensure_output_dirs()

    expected_project_id = "pbmc_10k_multiome_v1"

    print("\n" + "=" * 70)
    print("RUN 01 — PAIRED MULTIOME LOAD AND VALIDATION")
    print("=" * 70)

    print("\nActive project:")
    print(paths.project_id)

    print("\nProject root:")
    print(paths.project_root)

    if paths.project_id != expected_project_id:
        raise RuntimeError(
            "Wrong active project.\n"
            f"Expected: {expected_project_id}\n"
            f"Observed: {paths.project_id}\n"
            "Update Refactored/project_registry.yaml before running."
        )

    # --------------------------------------------------------
    # 3.2 Load project-local YAML
    # --------------------------------------------------------

    config_path = paths.configs / "raw_inputs.yaml"

    if not config_path.exists():
        raise FileNotFoundError(
            f"Configuration file not found:\n{config_path}"
        )

    config = load_yaml(config_path)

    if "raw_inputs" not in config:
        raise KeyError("raw_inputs.yaml has no 'raw_inputs' section")

    if "outputs" not in config:
        raise KeyError("raw_inputs.yaml has no 'outputs' section")

    raw_inputs = config["raw_inputs"]
    outputs = config["outputs"]
    validation = config.get("validation", {})

    print("\nLoaded config:")
    print(config_path)

    # --------------------------------------------------------
    # 3.3 Resolve external/shared raw input
    # --------------------------------------------------------

    if "multiome_10x_h5" not in raw_inputs:
        raise KeyError(
            "raw_inputs.yaml is missing raw_inputs.multiome_10x_h5"
        )

    multiome_file = paths.resolve_project_path(
    raw_inputs["multiome_10x_h5"]
 )
    
    input_exists = multiome_file.exists()
    input_is_file = multiome_file.is_file()

    input_size_mb = (
        multiome_file.stat().st_size / 1_000_000
        if input_exists and input_is_file
        else None
    )

    file_check = pd.DataFrame(
        [
            {
                "input_name": "multiome_10x_h5",
                "configured_path": raw_inputs["multiome_10x_h5"],
                "resolved_path": str(multiome_file),
                "exists": input_exists,
                "is_file": input_is_file,
                "size_mb": input_size_mb,
            }
        ]
    )

    print("\nRaw input check:")
    print(file_check.to_string(index=False))

    if not input_exists:
        raise FileNotFoundError(
            f"Multiome input file not found:\n{multiome_file}"
        )

    if not input_is_file:
        raise RuntimeError(
            "The configured Multiome input is not a file:\n"
            f"{multiome_file}"
        )

    # --------------------------------------------------------
    # 3.4 Prepare safe project-local outputs
    # --------------------------------------------------------

    required_output_keys = {
        "raw_multiome_h5ad",
        "raw_rna_h5ad",
        "raw_atac_h5ad",
        "file_check_table",
        "feature_type_table",
        "pairing_summary_table",
        "barcode_preview_table",
        "manifest",
    }

    missing_output_keys = required_output_keys.difference(outputs)

    if missing_output_keys:
        raise KeyError(
            "Missing output keys in raw_inputs.yaml: "
            f"{sorted(missing_output_keys)}"
        )

    multiome_out = paths.prepare_output_path(
        outputs["raw_multiome_h5ad"],
        label="raw_multiome_h5ad",
    )

    rna_out = paths.prepare_output_path(
        outputs["raw_rna_h5ad"],
        label="raw_rna_h5ad",
    )

    atac_out = paths.prepare_output_path(
        outputs["raw_atac_h5ad"],
        label="raw_atac_h5ad",
    )

    file_check_out = paths.prepare_output_path(
        outputs["file_check_table"],
        label="file_check_table",
    )

    feature_type_out = paths.prepare_output_path(
        outputs["feature_type_table"],
        label="feature_type_table",
    )

    pairing_summary_out = paths.prepare_output_path(
        outputs["pairing_summary_table"],
        label="pairing_summary_table",
    )

    barcode_preview_out = paths.prepare_output_path(
        outputs["barcode_preview_table"],
        label="barcode_preview_table",
    )

    manifest_out = paths.prepare_output_path(
        outputs["manifest"],
        label="step01_manifest",
    )

    prepared_outputs = {
        "raw_multiome_h5ad": multiome_out,
        "raw_rna_h5ad": rna_out,
        "raw_atac_h5ad": atac_out,
        "file_check_table": file_check_out,
        "feature_type_table": feature_type_out,
        "pairing_summary_table": pairing_summary_out,
        "barcode_preview_table": barcode_preview_out,
        "manifest": manifest_out,
    }

    for label, output_path in prepared_outputs.items():
        confirm_output_inside_project(
            output_path,
            paths.project_root,
            label=label,
        )

    print("\nPrepared project-local outputs:")

    for label, output_path in prepared_outputs.items():
        print(f"{label}: {output_path}")

    write_csv(
        file_check,
        file_check_out,
        index=False,
    )

    # --------------------------------------------------------
    # 3.5 Load the full 10x Multiome H5
    # --------------------------------------------------------

    print("\nLoading paired 10x Multiome H5...")

    multiome = sc.read_10x_h5(
        multiome_file,
        gex_only=False,
    )

    print("\nRaw Multiome object:")
    print(multiome)

    require_columns(
        multiome.var,
        {"feature_types"},
        dataframe_name="multiome.var",
    )

    # Preserve the original feature names before making them unique.
    multiome.var["feature_name_original"] = (
        multiome.var_names.astype(str)
    )

    multiome.var_names_make_unique()

    # Add an explicit shared cell key before splitting.
    multiome.obs["cell_id"] = multiome.obs_names.astype(str)

    # --------------------------------------------------------
    # 3.6 Validate feature types
    # --------------------------------------------------------

    feature_type_series = (
        multiome.var["feature_types"].astype(str)
    )

    feature_type_counts = (
        feature_type_series
        .value_counts()
        .rename_axis("feature_type")
        .reset_index(name="n_features")
    )

    print("\nFeature-type counts:")
    print(feature_type_counts.to_string(index=False))

    required_feature_types = set(
        validation.get(
            "required_feature_types",
            ["Gene Expression", "Peaks"],
        )
    )

    observed_feature_types = set(
        feature_type_series.unique()
    )

    missing_feature_types = (
        required_feature_types - observed_feature_types
    )

    if missing_feature_types:
        write_csv(
            feature_type_counts,
            feature_type_out,
            index=False,
        )

        raise RuntimeError(
            "Input H5 is not a valid RNA + ATAC Multiome matrix.\n"
            f"Missing feature types: {sorted(missing_feature_types)}\n"
            f"Observed feature types: {sorted(observed_feature_types)}"
        )

    # --------------------------------------------------------
    # 3.7 Split into paired RNA and ATAC objects
    # --------------------------------------------------------

    rna_mask = feature_type_series == "Gene Expression"
    atac_mask = feature_type_series == "Peaks"

    print("\nSplitting modalities...")
    print("Gene Expression features:", int(rna_mask.sum()))
    print("Peak features:", int(atac_mask.sum()))

    rna = multiome[:, rna_mask].copy()
    atac = multiome[:, atac_mask].copy()

    # Explicit common join key.
    rna.obs["cell_id"] = rna.obs_names.astype(str)
    atac.obs["cell_id"] = atac.obs_names.astype(str)

    # Helpful modality metadata.
    rna.uns["modality"] = "RNA"
    atac.uns["modality"] = "ATAC"
    multiome.uns["dataset_type"] = "paired_10x_multiome"

    print("\nRNA object:")
    print(rna)

    print("\nATAC object:")
    print(atac)

    # --------------------------------------------------------
    # 3.8 Strict barcode validation
    # --------------------------------------------------------

    rna_barcodes = pd.Index(
        rna.obs_names.astype(str),
        name="cell_id",
    )

    atac_barcodes = pd.Index(
        atac.obs_names.astype(str),
        name="cell_id",
    )

    multiome_barcodes = pd.Index(
        multiome.obs_names.astype(str),
        name="cell_id",
    )

    rna_duplicate_count = int(
        rna_barcodes.duplicated().sum()
    )

    atac_duplicate_count = int(
        atac_barcodes.duplicated().sum()
    )

    multiome_duplicate_count = int(
        multiome_barcodes.duplicated().sum()
    )

    shared_barcodes = rna_barcodes.intersection(
        atac_barcodes
    )

    same_cell_count = (
        rna.n_obs == atac.n_obs == multiome.n_obs
    )

    same_barcode_set = (
        set(rna_barcodes)
        == set(atac_barcodes)
        == set(multiome_barcodes)
    )

    same_barcode_order = (
        rna_barcodes.equals(atac_barcodes)
        and rna_barcodes.equals(multiome_barcodes)
    )

    multiome_barcode_hash = hash_barcodes(
        multiome_barcodes
    )

    rna_barcode_hash = hash_barcodes(
        rna_barcodes
    )

    atac_barcode_hash = hash_barcodes(
        atac_barcodes
    )

    same_barcode_hash = (
        multiome_barcode_hash
        == rna_barcode_hash
        == atac_barcode_hash
    )

    require_identical_order = bool(
        validation.get(
            "require_identical_barcode_order",
            True,
        )
    )

    require_no_duplicates = bool(
        validation.get(
            "require_no_duplicate_barcodes",
            True,
        )
    )

    no_duplicate_barcodes = (
        multiome_duplicate_count == 0
        and rna_duplicate_count == 0
        and atac_duplicate_count == 0
    )

    paired_confirmed = (
        same_cell_count
        and same_barcode_set
        and same_barcode_hash
        and (
            same_barcode_order
            if require_identical_order
            else True
        )
        and (
            no_duplicate_barcodes
            if require_no_duplicates
            else True
        )
    )

    pairing_summary = pd.DataFrame(
        [
            {
                "multiome_n_cells": int(multiome.n_obs),
                "multiome_n_features": int(multiome.n_vars),
                "rna_n_cells": int(rna.n_obs),
                "rna_n_genes": int(rna.n_vars),
                "atac_n_cells": int(atac.n_obs),
                "atac_n_peaks": int(atac.n_vars),
                "shared_barcodes": int(
                    len(shared_barcodes)
                ),
                "rna_shared_fraction": (
                    len(shared_barcodes) / rna.n_obs
                ),
                "atac_shared_fraction": (
                    len(shared_barcodes) / atac.n_obs
                ),
                "multiome_duplicate_barcodes": (
                    multiome_duplicate_count
                ),
                "rna_duplicate_barcodes": (
                    rna_duplicate_count
                ),
                "atac_duplicate_barcodes": (
                    atac_duplicate_count
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
                "same_barcode_hash": bool(
                    same_barcode_hash
                ),
                "paired_confirmed": bool(
                    paired_confirmed
                ),
            }
        ]
    )

    print("\nPairing validation:")
    print(pairing_summary.to_string(index=False))

    # --------------------------------------------------------
    # 3.9 Barcode preview
    # --------------------------------------------------------

    preview_n = min(25, multiome.n_obs)

    barcode_preview = pd.DataFrame(
        {
            "position": range(preview_n),
            "multiome_barcode": (
                multiome_barcodes[:preview_n].tolist()
            ),
            "rna_barcode": (
                rna_barcodes[:preview_n].tolist()
            ),
            "atac_barcode": (
                atac_barcodes[:preview_n].tolist()
            ),
        }
    )

    barcode_preview["all_same_position"] = (
        (
            barcode_preview["multiome_barcode"]
            == barcode_preview["rna_barcode"]
        )
        & (
            barcode_preview["multiome_barcode"]
            == barcode_preview["atac_barcode"]
        )
    )

    # --------------------------------------------------------
    # 3.10 Save validation tables
    # --------------------------------------------------------

    write_csv(
        feature_type_counts,
        feature_type_out,
        index=False,
    )

    write_csv(
        pairing_summary,
        pairing_summary_out,
        index=False,
    )

    write_csv(
        barcode_preview,
        barcode_preview_out,
        index=False,
    )

    # Hard fail after preserving diagnostic tables.
    if not paired_confirmed:
        raise RuntimeError(
            "Paired Multiome validation failed.\n"
            f"Review the pairing table:\n{pairing_summary_out}"
        )

    # --------------------------------------------------------
    # 3.11 Save validated H5AD objects
    # --------------------------------------------------------

    print("\nPairing confirmed.")
    print("Saving raw H5AD objects...")

    multiome.write_h5ad(
        multiome_out,
        compression="gzip",
    )

    rna.write_h5ad(
        rna_out,
        compression="gzip",
    )

    atac.write_h5ad(
        atac_out,
        compression="gzip",
    )

    print("\nSaved full Multiome object:")
    print(multiome_out)

    print("\nSaved RNA object:")
    print(rna_out)

    print("\nSaved ATAC object:")
    print(atac_out)

    # --------------------------------------------------------
    # 3.12 Manifest
    # --------------------------------------------------------

    manifest = {
        "step": "01",
        "script": script_path.name,
        "description": (
            "Load one paired 10x Multiome H5, validate feature "
            "types and barcode alignment, split RNA and ATAC, "
            "and save project-aware raw AnnData objects."
        ),
        "project": paths.as_dict(),
        "config": str(config_path),
        "inputs": {
            "multiome_10x_h5": str(multiome_file),
            "multiome_10x_h5_size_mb": input_size_mb,
        },
        "validation_parameters": {
            "required_feature_types": sorted(
                required_feature_types
            ),
            "require_identical_barcode_order": (
                require_identical_order
            ),
            "require_no_duplicate_barcodes": (
                require_no_duplicates
            ),
        },
        "outputs": {
            "raw_multiome_h5ad": str(multiome_out),
            "raw_rna_h5ad": str(rna_out),
            "raw_atac_h5ad": str(atac_out),
            "file_check_table": str(file_check_out),
            "feature_type_table": str(feature_type_out),
            "pairing_summary_table": str(
                pairing_summary_out
            ),
            "barcode_preview_table": str(
                barcode_preview_out
            ),
            "manifest": str(manifest_out),
        },
        "summary": {
            "multiome_n_cells": int(multiome.n_obs),
            "multiome_n_features": int(multiome.n_vars),
            "rna_n_cells": int(rna.n_obs),
            "rna_n_genes": int(rna.n_vars),
            "atac_n_cells": int(atac.n_obs),
            "atac_n_peaks": int(atac.n_vars),
            "shared_barcodes": int(
                len(shared_barcodes)
            ),
            "paired_confirmed": bool(
                paired_confirmed
            ),
            "multiome_barcode_hash": (
                multiome_barcode_hash
            ),
            "rna_barcode_hash": rna_barcode_hash,
            "atac_barcode_hash": atac_barcode_hash,
        },
    }

    write_json(
        manifest,
        manifest_out,
    )

    print("\nSaved manifest:")
    print(manifest_out)

    print("\n" + "=" * 70)
    print("RUN 01 COMPLETED SUCCESSFULLY")
    print("=" * 70)
    print("Paired Multiome confirmed:", paired_confirmed)
    print("Shared cells:", len(shared_barcodes))
    print("RNA genes:", rna.n_vars)
    print("ATAC peaks:", atac.n_vars)
    print("Project root:", paths.project_root)


# ============================================================
# 4. Entry point
# ============================================================

if __name__ == "__main__":
    main()