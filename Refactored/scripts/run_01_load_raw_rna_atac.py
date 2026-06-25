"""Βήμα 01 — Load raw RNA/ATAC files και save project-aware raw h5ad objects"""

from pathlib import Path  # χειρισμός paths
import sys  # για import path setup

import pandas as pd  # για file check table
import scanpy as sc  # για read_10x_h5


# ============================================================
# 0. Locate Refactored root in a portable way
# ============================================================

refactor_dir = Path(__file__).resolve().parents[1]  # Refactored/

src_dir = refactor_dir / "src"  # Refactored/src

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))  # κάνει διαθέσιμο το pbmcgrn package


# ============================================================
# 1. Import helpers
# ============================================================

from pbmcgrn.config import load_yaml  # YAML loader
from pbmcgrn.io import write_csv, write_json  # save helpers
from pbmcgrn.paths import ProjectPaths  # project-aware paths


# ============================================================
# 2. Initialize project paths
# ============================================================

paths = ProjectPaths(refactor_root=refactor_dir)  # διαβάζει active project από project_registry.yaml
paths.ensure_output_dirs()  # δημιουργεί project folders αν λείπουν

print("\nActive project:")
print(paths.project_id)  # πρέπει να είναι pbmc_10k_v1

print("\nProject root:")
print(paths.project_root)  # Refactored/projects/pbmc_10k_v1

# ============================================================
# 3. Load raw input config
# ============================================================

config_path = paths.configs / "raw_inputs.yaml"  # project-local config
config = load_yaml(config_path)  # φορτώνει raw input paths

raw_inputs = config["raw_inputs"]  # raw RNA/ATAC input paths
outputs = config["outputs"]  # project output paths

print("\nLoaded raw input config:")
print(config_path)


# ============================================================
# 4. Resolve raw input paths
# ============================================================

rna_file = paths.resolve_project_path(raw_inputs["rna_10x_h5"])  # λύνει path από project_root
atac_file = paths.resolve_project_path(raw_inputs["atac_10x_h5"])  # λύνει path από project_root

print("\nRNA input:")
print(rna_file)
print("Exists:", rna_file.exists())  # πρέπει να είναι True

print("\nATAC input:")
print(atac_file)
print("Exists:", atac_file.exists())  # πρέπει να είναι True

if not rna_file.exists():
    raise FileNotFoundError(f"RNA input file not found:\n{rna_file}")  # σταματά αν λείπει RNA

if not atac_file.exists():
    raise FileNotFoundError(f"ATAC input file not found:\n{atac_file}")  # σταματά αν λείπει ATAC


# ============================================================
# 5. Check all raw input files from YAML
# ============================================================

file_check_rows = []  # λίστα για input file check table

for input_name, input_path_string in raw_inputs.items():
    input_path = paths.resolve_project_path(input_path_string)  # resolve από project_root
    exists = input_path.exists()  # check αν υπάρχει
    size_mb = input_path.stat().st_size / 1e6 if exists else None  # μέγεθος σε MB

    file_check_rows.append(
        {
            "input_name": input_name,
            "path": str(input_path),
            "exists": exists,
            "size_mb": size_mb,
        }
    )  # μία γραμμή ανά raw input

file_check = pd.DataFrame(file_check_rows)  # γίνεται table

print("\nRaw input file check:")
print(file_check)

missing_inputs = file_check.loc[~file_check["exists"], "input_name"].tolist()  # ποια λείπουν

if missing_inputs:
    raise FileNotFoundError(
        "Λείπουν raw input files:\n" + "\n".join(missing_inputs)
    )  # σταματά αν λείπει κάτι από YAML


# ============================================================
# 6. Prepare project-aware outputs
# ============================================================

rna_out = paths.prepare_output_path(
    outputs["raw_rna_h5ad"],
    label="raw_rna_h5ad",
)  # project-local RNA raw h5ad

atac_out = paths.prepare_output_path(
    outputs["raw_atac_h5ad"],
    label="raw_atac_h5ad",
)  # project-local ATAC raw h5ad

file_check_out = paths.prepare_output_path(
    outputs["file_check_table"],
    label="file_check_table",
)  # project-local file check CSV

manifest_out = paths.prepare_output_path(
    "reports/manifests/step01_load_raw_rna_atac_manifest.json",
    label="step01_manifest",
)  # project-local manifest


# ============================================================
# 7. Save raw input file check table
# ============================================================

write_csv(file_check, file_check_out, index=False)  # save CSV με raw input checks

print("\nSaved raw input file check table to:")
print(file_check_out)


# ============================================================
# 8. Load RNA 10x h5
# ============================================================

print("\nLoading RNA 10x h5...")

rna = sc.read_10x_h5(rna_file)  # φορτώνει RNA counts
rna.var_names_make_unique()  # κάνει μοναδικά gene names

print("\nRNA object:")
print(rna)

if "feature_types" in rna.var.columns:
    print("\nRNA feature types:")
    print(rna.var["feature_types"].value_counts())  # sanity check feature types


# ============================================================
# 9. Load ATAC 10x h5
# ============================================================

print("\nLoading ATAC 10x h5...")

atac = sc.read_10x_h5(
    atac_file,
    gex_only=False,
)  # για ATAC θέλουμε όλα τα features, όχι μόνο Gene Expression

atac.var_names_make_unique()  # κάνει μοναδικά peak names

print("\nATAC object before filtering:")
print(atac)

if "feature_types" in atac.var.columns:
    print("\nATAC feature types before filtering:")
    print(atac.var["feature_types"].value_counts())  # δείχνει Peaks κλπ

    peak_mask = atac.var["feature_types"] == "Peaks"  # κρατά μόνο Peaks

    print("\nFiltering ATAC to Peaks only...")
    print("Features before:", atac.n_vars)
    print("Peaks:", int(peak_mask.sum()))

    atac = atac[:, peak_mask].copy()  # κρατά μόνο peak matrix

    print("\nATAC object after peak filtering:")
    print(atac)
else:
    print("\nNo feature_types column found. Keeping all ATAC features.")  # fallback


# ============================================================
# 10. Raw object summary
# ============================================================

summary = {
    "rna_n_cells": int(rna.n_obs),
    "rna_n_genes": int(rna.n_vars),
    "atac_n_cells": int(atac.n_obs),
    "atac_n_peaks": int(atac.n_vars),
}  # βασικό summary για manifest

print("\nRaw object summary:")
for key, value in summary.items():
    print(f"{key}: {value}")  # εμφανίζει βασικά μεγέθη


# ============================================================
# 11. Save raw h5ad objects
# ============================================================

rna.write_h5ad(rna_out)  # σώζει RNA raw object μέσα στο active project
atac.write_h5ad(atac_out)  # σώζει ATAC raw object μέσα στο active project

print("\nSaved raw RNA h5ad to:")
print(rna_out)

print("\nSaved raw ATAC h5ad to:")
print(atac_out)


# ============================================================
# 12. Save manifest
# ============================================================

manifest = {
    "step": "01",
    "description": "Load raw RNA and ATAC 10x h5 files and save project-aware raw h5ad files",
    "project": paths.as_dict(),
    "config": str(config_path),
    "inputs": {
        key: str(paths.resolve_project_path(value))
        for key, value in raw_inputs.items()
    },
    "outputs": {
        "raw_rna_h5ad": str(rna_out),
        "raw_atac_h5ad": str(atac_out),
        "file_check_table": str(file_check_out),
    },
    "summary": summary,
}  # reproducibility manifest

write_json(manifest, manifest_out)  # save manifest JSON

print("\nSaved manifest to:")
print(manifest_out)


# ============================================================
# 13. Final summary
# ============================================================

print("\nDone.")
print("Το Step 01 ολοκληρώθηκε.")
print("Τα raw h5ad αρχεία γράφτηκαν μέσα στο active project workspace.")
print(f"Project root: {paths.project_root}")