# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 14:10:30 2026

@author: jkary
"""

"""Βήμα 00E — Έλεγχος ότι το refactored Python package φορτώνεται σωστά"""

from pathlib import Path
import sys


# ============================================================
# 0. Ορισμός Refactored project root
# ============================================================

refactor_dir = Path(
    r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project/Refactored"
)

src_dir = refactor_dir / "src"

print("Refactored project directory:")
print(refactor_dir)

print("\nSource directory:")
print(src_dir)


# ============================================================
# 1. Έλεγχος ότι υπάρχουν οι βασικοί φάκελοι
# ============================================================

required_dirs = [
    refactor_dir,
    src_dir,
    refactor_dir / "src" / "pbmcgrn",
    refactor_dir / "configs",
    refactor_dir / "scripts",
    refactor_dir / "data",
    refactor_dir / "reports",
    refactor_dir / "models",
]

print("\nChecking required folders:")

for path in required_dirs:
    print(path, "| exists:", path.exists())

missing_dirs = [path for path in required_dirs if not path.exists()]

if missing_dirs:
    raise FileNotFoundError(
        "Λείπουν βασικοί φάκελοι. Πρώτα πρέπει να τρέξει το folder setup script."
    )


# ============================================================
# 2. Προσθήκη του src/ στο Python path
# ============================================================
# Αυτό επιτρέπει imports όπως:
# from pbmcgrn.paths import ProjectPaths
#
# Χωρίς αυτό, το Python μπορεί να μη βρίσκει το package.

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

print("\nAdded src_dir to sys.path:")
print(src_dir)


# ============================================================
# 3. Import των helper modules
# ============================================================

from pbmcgrn.paths import ProjectPaths
from pbmcgrn.config import load_yaml
from pbmcgrn.io import write_json


print("\nImports successful:")
print("- ProjectPaths")
print("- load_yaml")
print("- write_json")


# ============================================================
# 4. Έλεγχος ProjectPaths
# ============================================================

paths = ProjectPaths(refactor_root=refactor_dir)
paths.ensure_output_dirs()

print("\nImportant resolved paths:")
print("refactor_root:", paths.refactor_root)
print("legacy_root:", paths.legacy_root)
print("legacy_processed:", paths.legacy_processed)
print("processed_gnn:", paths.processed_gnn)
print("reports:", paths.reports)
print("models_gnn:", paths.models_gnn)


# ============================================================
# 5. Έλεγχος config files
# ============================================================

project_config_path = paths.configs / "project.yaml"
cell_types_config_path = paths.configs / "cell_types.yaml"

print("\nChecking config files:")
print("project.yaml exists:", project_config_path.exists())
print("cell_types.yaml exists:", cell_types_config_path.exists())

if project_config_path.exists():
    project_config = load_yaml(project_config_path)
    print("\nLoaded project config keys:")
    print(project_config.keys())

if cell_types_config_path.exists():
    cell_type_config = load_yaml(cell_types_config_path)
    print("\nLoaded cell types:")
    print(cell_type_config.get("cell_types"))


# ============================================================
# 6. Γράψε μικρό test manifest
# ============================================================

manifest = {
    "step": "00E",
    "description": "Import test for refactored pbmcgrn package",
    "refactor_root": str(paths.refactor_root),
    "legacy_root": str(paths.legacy_root),
    "status": "success",
}

manifest_out = paths.manifests / "test_00_import_refactored_package_manifest.json"
write_json(manifest, manifest_out)

print("\nSaved test manifest to:")
print(manifest_out)


# ============================================================
# 7. Τελικό summary
# ============================================================

print("\nDone.")
print("Το refactored package φορτώνεται σωστά.")
print("Μπορούμε να προχωρήσουμε στον έλεγχο των παλιών GNN input files.")