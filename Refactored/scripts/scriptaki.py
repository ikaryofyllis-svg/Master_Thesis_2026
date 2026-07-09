import pandas as pd
import scanpy as sc

from pathlib import Path

import sys
print(sys.executable)

project_root = Path("Refactored/projects/pbmc_10k_v1")

rna_path = project_root / "data/processed/rna/pbmc_rna_annotated_reviewed.h5ad"
atac_path = project_root / "data/processed/atac/pbmc_atac_clustered.h5ad"

rna = sc.read_h5ad(rna_path)
atac = sc.read_h5ad(atac_path)

print(rna)
print(atac)
rna_cells = pd.Index(rna.obs_names.astype(str))
atac_cells = pd.Index(atac.obs_names.astype(str))

shared = rna_cells.intersection(atac_cells)

summary = pd.DataFrame({
    "metric": [
        "n_rna_cells",
        "n_atac_cells",
        "n_shared_cells",
        "fraction_rna_shared",
        "fraction_atac_shared",
    ],
    "value": [
        len(rna_cells),
        len(atac_cells),
        len(shared),
        len(shared) / len(rna_cells),
        len(shared) / len(atac_cells),
    ],
})

print(summary)


import sys
import subprocess

subprocess.check_call([
    sys.executable,
    "-m",
    "pip",
    "install",
    "scanpy",
    "anndata",
])

import sys
from importlib.metadata import version

print("Python:", sys.executable)
print("Spyder kernels:", version("spyder-kernels"))
print("Scanpy:", version("scanpy"))

import sys
from importlib.metadata import version

print("Python:", sys.executable)
print("Spyder:", version("spyder"))
print("spyder-kernels:", version("spyder-kernels"))
print("Scanpy:", version("scanpy"))

from pathlib import Path
import pandas as pd
import scanpy as sc

refactored_root = Path(
    r"C:\Users\jkary\Desktop\bioinformatics\MSc thesis\Msc_Thesis_Project\Refactored"
)

project_root = refactored_root / "projects" / "pbmc_10k_v1"

rna_path = (
    project_root
    / "data"
    / "processed"
    / "rna"
    / "pbmc_rna_annotated_reviewed.h5ad"
)

atac_path = (
    project_root
    / "data"
    / "processed"
    / "atac"
    / "pbmc_atac_clustered.h5ad"
)

print("RNA exists:", rna_path.exists(), rna_path)
print("ATAC exists:", atac_path.exists(), atac_path)

rna = sc.read_h5ad(rna_path)
atac = sc.read_h5ad(atac_path)

rna_cells = pd.Index(rna.obs_names.astype(str))
atac_cells = pd.Index(atac.obs_names.astype(str))

shared = rna_cells.intersection(atac_cells)

summary = pd.DataFrame({
    "metric": [
        "n_rna_cells",
        "n_atac_cells",
        "n_shared_cells",
        "fraction_rna_shared",
        "fraction_atac_shared",
    ],
    "value": [
        len(rna_cells),
        len(atac_cells),
        len(shared),
        len(shared) / len(rna_cells),
        len(shared) / len(atac_cells),
    ],
})

print(summary.to_string(index=False))


import scanpy as sc

adata = sc.read_10x_h5(
    r"C:\Users\jkary\Desktop\bioinformatics\MSc thesis\Msc_Thesis_Project\data\raw\pbmc_rna_10k_v3\pbmc_10k_v3_filtered_feature_bc_matrix.h5",
    gex_only=False,
)

print(adata)
print(adata.var["feature_types"].value_counts())


import scanpy as sc

multiome = sc.read_10x_h5(
    r"C:\Users\jkary\Desktop\bioinformatics\MSc thesis\Msc_Thesis_Project\data\raw\PBMC from a Healthy Donor - Granulocytes Removed Through Cell Sorting (10k)\pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
    gex_only=False,
)

print(multiome)
print(multiome.var["feature_types"].value_counts())


multiome.var_names_make_unique()

gex_mask = multiome.var["feature_types"] == "Gene Expression"
peak_mask = multiome.var["feature_types"] == "Peaks"

rna_multiome = multiome[:, gex_mask].copy()
atac_multiome = multiome[:, peak_mask].copy()

print("RNA shape:", rna_multiome.shape)
print("ATAC shape:", atac_multiome.shape)

print(
    "Same barcodes:",
    rna_multiome.obs_names.equals(atac_multiome.obs_names)
)

from pathlib import Path


def main() -> None:
    refactored_root = Path(
        r"C:\Users\jkary\Desktop\bioinformatics\MSc thesis"
        r"\Msc_Thesis_Project\Refactored"
    )

    project_root = (
        refactored_root
        / "projects"
        / "pbmc_10k_multiome_v1"
    )

    folders = [
        project_root / "configs",
        project_root / "data" / "raw",
        project_root / "data" / "processed" / "rna",
        project_root / "data" / "processed" / "atac",
        project_root / "data" / "processed" / "multiome",
        project_root / "reports" / "tables",
        project_root / "reports" / "figures",
        project_root / "reports" / "manifests",
    ]

    print("Project root:")
    print(project_root)

    for folder in folders:
        folder.mkdir(parents=True, exist_ok=True)
        print(f"Created/exists: {folder}")

    print("\nDone.")
    print("Project structure created successfully.")


main()