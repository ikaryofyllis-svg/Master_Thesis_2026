# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 14:50:50 2026

@author: jkary
"""

import scanpy as sc
from pathlib import Path

project_root = Path(
    r"C:\Users\jkary\Desktop\bioinformatics\MSc thesis"
    r"\Msc_Thesis_Project\Refactored"
    r"\projects\pbmc_10k_multiome_v1"
)

rna = sc.read_h5ad(
    project_root
    / "data/processed/rna/pbmc_multiome_rna_raw.h5ad"
)

atac = sc.read_h5ad(
    project_root
    / "data/processed/atac/pbmc_multiome_atac_raw.h5ad"
)

print(rna.shape)
print(atac.shape)
print(rna.obs_names.equals(atac.obs_names))
print((rna.obs["cell_id"] == atac.obs["cell_id"]).all())
