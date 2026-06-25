"""Βήμα 00R — Έλεγχος project-aware paths και active project workspace"""

from pathlib import Path
import sys
import json


# ============================================================
# 0. Refactored project root και import path
# ============================================================

refactor_dir = Path(
    r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project/Refactored"
)

src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))


# ============================================================
# 1. Load ProjectPaths
# ============================================================

from pbmcgrn.paths import ProjectPaths


paths = ProjectPaths(refactor_root=refactor_dir)
paths.ensure_output_dirs()


# ============================================================
# 2. Print active project paths
# ============================================================

print("Active project id:")
print(paths.project_id)

print("\nProject root:")
print(paths.project_root)

print("\nRNA processed directory:")
print(paths.processed_rna)

print("\nAnnotation tables directory:")
print(paths.tables_annotation)

print("\nRNA figures directory:")
print(paths.figures_rna)

print("\nManifests directory:")
print(paths.manifests)


# ============================================================
# 3. Save test manifest
# ============================================================

test_manifest = {
    "status": "success",
    "message": "Project-aware paths are working.",
    "paths": paths.as_dict(),
}

test_manifest_out = (
    paths.manifests / "test_00R_project_paths_manifest.json"
)

with open(test_manifest_out, "w", encoding="utf-8") as handle:
    json.dump(
        test_manifest,
        handle,
        indent=2,
        ensure_ascii=False,
    )

print("\nSaved test manifest to:")
print(test_manifest_out)

print("\nDone.")