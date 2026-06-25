"""Βήμα 00S — Safe copy υπαρχόντων Refactored outputs στο active project workspace"""

from pathlib import Path
import sys
import shutil
import json
import yaml


# ============================================================
# 0. Refactored root και import path
# ============================================================

refactor_dir = Path(
    r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project/Refactored"
)

src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))


# ============================================================
# 1. Load project-aware paths
# ============================================================

from pbmcgrn.paths import ProjectPaths


paths = ProjectPaths(refactor_root=refactor_dir)
paths.ensure_output_dirs()

print("Active project id:")
print(paths.project_id)

print("\nProject root:")
print(paths.project_root)


# ============================================================
# 2. Migration settings
# ============================================================
# overwrite_existing = False:
#   Δεν αντικαθιστά υπάρχοντα αρχεία στο project workspace.
#
# Αυτό είναι ασφαλές γιατί δεν σβήνει και δεν μετακινεί τίποτα.

overwrite_existing = False

old_locations = {
    "configs": refactor_dir / "configs",
    "data": refactor_dir / "data",
    "reports": refactor_dir / "reports",
    "models": refactor_dir / "models",
}

new_locations = {
    "configs": paths.configs,
    "data": paths.data,
    "reports": paths.reports,
    "models": paths.models,
}


# ============================================================
# 3. Helper copy function
# ============================================================

def copy_tree_safe(src_root, dst_root, overwrite=False):
    """
    Copy files recursively from src_root to dst_root.

    Does not delete source files.
    Skips existing files unless overwrite=True.
    """

    src_root = Path(src_root)
    dst_root = Path(dst_root)

    result = {
        "source": str(src_root),
        "destination": str(dst_root),
        "copied_files": [],
        "skipped_existing_files": [],
        "missing_source": False,
    }

    if not src_root.exists():
        result["missing_source"] = True
        return result

    for src_file in src_root.rglob("*"):
        if src_file.is_dir():
            continue

        relative_path = src_file.relative_to(src_root)
        dst_file = dst_root / relative_path

        dst_file.parent.mkdir(parents=True, exist_ok=True)

        if dst_file.exists() and not overwrite:
            result["skipped_existing_files"].append(str(dst_file))
            continue

        shutil.copy2(src_file, dst_file)
        result["copied_files"].append(str(dst_file))

    return result


# ============================================================
# 4. Copy existing global outputs/configs into active project
# ============================================================

migration_results = {}

for key in old_locations:
    src = old_locations[key]
    dst = new_locations[key]

    print(f"\nCopying {key}:")
    print("FROM:", src)
    print("TO:  ", dst)

    result = copy_tree_safe(
        src_root=src,
        dst_root=dst,
        overwrite=overwrite_existing,
    )

    migration_results[key] = result

    print("Copied files:", len(result["copied_files"]))
    print("Skipped existing:", len(result["skipped_existing_files"]))
    print("Missing source:", result["missing_source"])


# ============================================================
# 5. Fix project raw_inputs.yaml relative raw data paths
# ============================================================
# Πριν:
#   ../data/raw/...
#
# Αυτό ήταν σωστό όταν config ήταν στο Refactored/configs.
# Τώρα config είναι στο Refactored/projects/pbmc_10k_v1/configs.
#
# Από project root προς original thesis data/raw:
#   ../../../data/raw/...
#
# project_root:
#   Refactored/projects/pbmc_10k_v1
#
# ../../../data:
#   Msc_Thesis_Project/data

raw_inputs_yaml = paths.configs / "raw_inputs.yaml"

raw_path_replacements = {
    "../data/raw/": "../../../data/raw/",
    "..\\data\\raw\\": "../../../data/raw/",
}

raw_inputs_fix_status = {
    "file": str(raw_inputs_yaml),
    "exists": raw_inputs_yaml.exists(),
    "changed": False,
}

if raw_inputs_yaml.exists():
    with open(raw_inputs_yaml, "r", encoding="utf-8") as handle:
        raw_text = handle.read()

    fixed_text = raw_text

    for old, new in raw_path_replacements.items():
        fixed_text = fixed_text.replace(old, new)

    if fixed_text != raw_text:
        with open(raw_inputs_yaml, "w", encoding="utf-8") as handle:
            handle.write(fixed_text)

        raw_inputs_fix_status["changed"] = True

    print("\nChecked/fixed raw_inputs.yaml:")
    print(raw_inputs_yaml)
    print("Changed:", raw_inputs_fix_status["changed"])

else:
    print("\nraw_inputs.yaml not found in project configs.")


# ============================================================
# 6. Create project README note
# ============================================================

readme_out = paths.project_root / "README_project_workspace.md"

readme_text = f"""# {paths.project_id}

This folder contains dataset/project-specific inputs, processed outputs, reports, figures, manifests and models.

Shared code lives outside this folder:

- `Refactored/src/`
- `Refactored/scripts/`

Project-specific outputs live here:

- `data/processed/`
- `reports/tables/`
- `reports/figures/`
- `reports/manifests/`
- `models/`

This prevents outputs from different datasets or experiments from being mixed.
"""

if not readme_out.exists():
    readme_out.write_text(readme_text, encoding="utf-8")

print("\nProject README:")
print(readme_out)


# ============================================================
# 7. Save migration manifest
# ============================================================

manifest = {
    "step": "00S",
    "description": (
        "Safe copy of existing global Refactored outputs into active project workspace. "
        "No source files were deleted or moved."
    ),
    "active_project_id": paths.project_id,
    "project_root": str(paths.project_root),
    "overwrite_existing": overwrite_existing,
    "migration_results": migration_results,
    "raw_inputs_fix_status": raw_inputs_fix_status,
}

manifest_out = paths.manifests / "step00S_migrate_existing_outputs_manifest.json"

with open(manifest_out, "w", encoding="utf-8") as handle:
    json.dump(
        manifest,
        handle,
        indent=2,
        ensure_ascii=False,
    )

print("\nSaved migration manifest to:")
print(manifest_out)


# ============================================================
# 8. Final summary
# ============================================================

print("\nDone.")
print("Τα υπάρχοντα Refactored outputs αντιγράφηκαν στο active project workspace.")
print("Δεν διαγράφηκε τίποτα από τα παλιά folders.")
print("\nNew project workspace:")
print(paths.project_root)

"""Βήμα 00T — Έλεγχος ότι τα βασικά migrated outputs υπάρχουν στο active project workspace"""

from pathlib import Path
import sys
import pandas as pd


# ============================================================
# 0. Refactored root και import path
# ============================================================

refactor_dir = Path(
    r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project/Refactored"
)

src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))


# ============================================================
# 1. Load project-aware paths
# ============================================================

from pbmcgrn.paths import ProjectPaths


paths = ProjectPaths(refactor_root=refactor_dir)
paths.ensure_output_dirs()

print("Active project:")
print(paths.project_id)

print("\nProject root:")
print(paths.project_root)


# ============================================================
# 2. Expected migrated files
# ============================================================

expected_files = {
    "rna_raw": paths.processed_rna / "pbmc_rna_raw.h5ad",
    "rna_filtered": paths.processed_rna / "pbmc_rna_filtered.h5ad",
    "rna_hvg": paths.processed_rna / "pbmc_rna_hvg.h5ad",
    "rna_clustered": paths.processed_rna / "pbmc_rna_clustered.h5ad",
    "rna_celltypist_annotated": paths.processed_rna / "pbmc_rna_celltypist_annotated.h5ad",
    "raw_inputs_config": paths.configs / "raw_inputs.yaml",
    "annotation_backend_config": paths.configs / "annotation_backend.yaml",
}


# ============================================================
# 3. Check existence
# ============================================================

rows = []

for name, path in expected_files.items():
    rows.append(
        {
            "name": name,
            "path": str(path),
            "exists": path.exists(),
            "size_mb": path.stat().st_size / 1e6 if path.exists() else None,
        }
    )

check_table = pd.DataFrame(rows)

print("\nMigration check table:")
print(check_table.to_string(index=False))


# ============================================================
# 4. Save check table
# ============================================================

out = paths.tables / "migration_check_project_outputs.csv"
out.parent.mkdir(parents=True, exist_ok=True)

check_table.to_csv(out, index=False)

print("\nSaved check table to:")
print(out)

print("\nDone.")

