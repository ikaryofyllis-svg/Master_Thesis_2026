"""Βήμα 00P — Δημιουργία project-based workspace system για καθαρά outputs ανά dataset/project"""

from pathlib import Path
import json
import shutil
import yaml


# ============================================================
# 0. Βασικές ρυθμίσεις
# ============================================================

refactor_dir = Path(
    r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project/Refactored"
)

active_project_id = "pbmc_10k_v1"

projects_dir = refactor_dir / "projects"
active_project_dir = projects_dir / active_project_id

print("Refactored root:")
print(refactor_dir)

print("\nActive project:")
print(active_project_id)

print("\nActive project directory:")
print(active_project_dir)


# ============================================================
# 1. Δημιουργία καθαρής project δομής
# ============================================================

project_dirs = [
    active_project_dir / "configs",

    active_project_dir / "data" / "raw",
    active_project_dir / "data" / "external",
    active_project_dir / "data" / "interim",
    active_project_dir / "data" / "processed",
    active_project_dir / "data" / "processed" / "rna",
    active_project_dir / "data" / "processed" / "atac",
    active_project_dir / "data" / "processed" / "integration",
    active_project_dir / "data" / "processed" / "grn",
    active_project_dir / "data" / "processed" / "gnn",

    active_project_dir / "reports",
    active_project_dir / "reports" / "tables",
    active_project_dir / "reports" / "tables" / "rna",
    active_project_dir / "reports" / "tables" / "atac",
    active_project_dir / "reports" / "tables" / "annotation",
    active_project_dir / "reports" / "tables" / "integration",
    active_project_dir / "reports" / "tables" / "grn",
    active_project_dir / "reports" / "tables" / "gnn",

    active_project_dir / "reports" / "figures",
    active_project_dir / "reports" / "figures" / "rna",
    active_project_dir / "reports" / "figures" / "atac",
    active_project_dir / "reports" / "figures" / "annotation",
    active_project_dir / "reports" / "figures" / "integration",
    active_project_dir / "reports" / "figures" / "grn",
    active_project_dir / "reports" / "figures" / "gnn",

    active_project_dir / "reports" / "manifests",
    active_project_dir / "reports" / "logs",

    active_project_dir / "models",
    active_project_dir / "models" / "annotation",
    active_project_dir / "models" / "gnn",

    active_project_dir / "notebooks",
]

for folder in project_dirs:
    folder.mkdir(parents=True, exist_ok=True)

print("\nCreated project folders.")


# ============================================================
# 2. Δημιουργία global project registry
# ============================================================
# Αυτό είναι το μόνο global config.
# Λέει ποιο project είναι ενεργό.
# Τα υπόλοιπα configs θα είναι ανά project.

registry = {
    "active_project_id": active_project_id,
    "projects": {
        active_project_id: {
            "project_name": "PBMC 10k RNA/ATAC GRN thesis project",
            "organism": "human",
            "domain": "PBMC / immune",
            "description": (
                "Refactored PBMC RNA/ATAC project for cell-type annotation, "
                "pseudo-bulk integration, GRN inference, and GNN refinement."
            ),
            "project_root": f"projects/{active_project_id}",
        }
    },
}

registry_out = refactor_dir / "project_registry.yaml"

with open(registry_out, "w", encoding="utf-8") as handle:
    yaml.safe_dump(
        registry,
        handle,
        sort_keys=False,
        allow_unicode=True,
        width=120,
    )

print("\nSaved project registry to:")
print(registry_out)


# ============================================================
# 3. Copy current configs into active project configs
# ============================================================
# Δεν σβήνουμε τα παλιά configs.
# Τα αντιγράφουμε μέσα στο active project για να είναι self-contained.

old_configs_dir = refactor_dir / "configs"
new_configs_dir = active_project_dir / "configs"

if old_configs_dir.exists():
    for config_file in old_configs_dir.glob("*.yaml"):
        dst = new_configs_dir / config_file.name

        if dst.exists():
            print("SKIP existing project config:", dst.name)
        else:
            shutil.copy2(config_file, dst)
            print("Copied config:", config_file.name)

else:
    print("\nNo old global configs directory found. Skipping config copy.")


# ============================================================
# 4. Optional migration notes
# ============================================================
# Δεν μετακινούμε αυτόματα μεγάλα h5ad files εδώ για ασφάλεια.
# Από εδώ και μετά τα scripts πρέπει να γράφουν στο project workspace.
#
# Αν θέλουμε, μετά γράφουμε ξεχωριστό migration script για να μεταφέρει:
# Refactored/data/processed/rna/*.h5ad
# Refactored/reports/tables/*.csv
# Refactored/reports/figures/rna/*.png
# στο projects/pbmc_10k_v1/.

migration_notes = {
    "important": (
        "Existing outputs were not moved automatically to avoid accidental data loss. "
        "Future scripts should write to projects/pbmc_10k_v1 using the updated ProjectPaths."
    ),
    "old_output_locations": {
        "data": str(refactor_dir / "data"),
        "reports": str(refactor_dir / "reports"),
        "models": str(refactor_dir / "models"),
    },
    "new_output_location": str(active_project_dir),
}

migration_notes_out = (
    active_project_dir
    / "reports"
    / "manifests"
    / "workspace_migration_notes.json"
)

with open(migration_notes_out, "w", encoding="utf-8") as handle:
    json.dump(
        migration_notes,
        handle,
        indent=2,
        ensure_ascii=False,
    )

print("\nSaved migration notes to:")
print(migration_notes_out)


# ============================================================
# 5. Τελικό summary
# ============================================================

print("\nDone.")
print("Η νέα project-based δομή είναι έτοιμη.")
print("Από εδώ και πέρα τα outputs πρέπει να γράφονται μέσα στο:")
print(active_project_dir)
