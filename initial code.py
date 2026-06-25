# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 14:14:13 2026

@author: jkary
"""
import scanpy as sc
from pathlib import Path

file_path = Path("data/raw/pbmc_rna_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5")

print("Exists:", file_path.exists())
print("Size MB:", file_path.stat().st_size / 1e6)

rna = sc.read_10x_h5(file_path)
rna.var_names_make_unique()

print(rna)
print(rna.obs.head())
print(rna.var.head())
from pathlib import Path
import os

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

print("Project directory:", Path.cwd())

files = [
    "data/raw/pbmc_rna_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5",

    "data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_filtered_peak_bc_matrix.h5",
    "data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_fragments.tsv.gz",
    "data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_fragments.tsv.gz.tbi",
    "data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_peaks.bed",
    "data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_peak_annotation.tsv",
    "data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_peak_motif_mapping.bed",
]

for f in files:
    p = project_dir / f

    if p.exists():
        print("OK:", f, "| size MB:", round(p.stat().st_size / 1e6, 2))
    else:
        print("MISSING:", f)
        





from pathlib import Path
import os
import scanpy as sc

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

print("Project directory:", Path.cwd())

rna_file = project_dir / "data/raw/pbmc_rna_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
atac_file = project_dir / "data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_filtered_peak_bc_matrix.h5"

print("RNA file exists:", rna_file.exists())
print("ATAC file exists:", atac_file.exists())

rna = sc.read_10x_h5(rna_file)
rna.var_names_make_unique()

atac = sc.read_10x_h5(atac_file)
atac.var_names_make_unique()

print("\nRNA object:")
print(rna)
print(rna.var.head())

print("\nATAC object:")
print(atac)
print(atac.var.head())



from pathlib import Path
import os
import scanpy as sc
import pandas as pd

# --------------------------------------------------
# 0. Project directory
# --------------------------------------------------

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

print("Project directory:", Path.cwd())

# --------------------------------------------------
# 1. File paths
# --------------------------------------------------

rna_file = project_dir / "data/raw/pbmc_rna_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5"

atac_file = project_dir / "data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_filtered_peak_bc_matrix.h5"

processed_dir = project_dir / "data/processed"
processed_dir.mkdir(parents=True, exist_ok=True)

print("\nRNA exists:", rna_file.exists())
print("ATAC exists:", atac_file.exists())

# --------------------------------------------------
# 2. Load RNA
# --------------------------------------------------

rna = sc.read_10x_h5(rna_file)
rna.var_names_make_unique()

print("\n================ RNA ================")
print(rna)
print("\nRNA obs head:")
print(rna.obs.head())
print("\nRNA var head:")
print(rna.var.head())

if "feature_types" in rna.var.columns:
    print("\nRNA feature types:")
    print(rna.var["feature_types"].value_counts())

# --------------------------------------------------
# 3. Load ATAC
# --------------------------------------------------
# Important: for ATAC we use gex_only=False because the features are peaks, not genes.

atac = sc.read_10x_h5(atac_file, gex_only=False)
atac.var_names_make_unique()

print("\n================ ATAC ================")
print(atac)
print("\nATAC obs head:")
print(atac.obs.head())
print("\nATAC var head:")
print(atac.var.head())

if "feature_types" in atac.var.columns:
    print("\nATAC feature types:")
    print(atac.var["feature_types"].value_counts())

# --------------------------------------------------
# 4. Keep only peaks if multiple feature types exist
# --------------------------------------------------

if "feature_types" in atac.var.columns:
    print("\nFiltering ATAC object to Peaks only...")
    atac = atac[:, atac.var["feature_types"] == "Peaks"].copy()
    print("ATAC after peak filtering:")
    print(atac)

# --------------------------------------------------
# 5. Basic QC numbers
# --------------------------------------------------

print("\n================ BASIC QC ================")

print("RNA cells:", rna.n_obs)
print("RNA genes:", rna.n_vars)

print("ATAC cells:", atac.n_obs)
print("ATAC peaks:", atac.n_vars)

# RNA counts per cell
rna.obs["n_counts"] = rna.X.sum(axis=1).A1 if hasattr(rna.X, "A1") else rna.X.sum(axis=1)
rna.obs["n_genes"] = (rna.X > 0).sum(axis=1).A1 if hasattr((rna.X > 0).sum(axis=1), "A1") else (rna.X > 0).sum(axis=1)

# ATAC counts per cell
atac.obs["n_counts"] = atac.X.sum(axis=1).A1 if hasattr(atac.X, "A1") else atac.X.sum(axis=1)
atac.obs["n_peaks"] = (atac.X > 0).sum(axis=1).A1 if hasattr((atac.X > 0).sum(axis=1), "A1") else (atac.X > 0).sum(axis=1)

print("\nRNA QC summary:")
print(rna.obs[["n_counts", "n_genes"]].describe())

print("\nATAC QC summary:")
print(atac.obs[["n_counts", "n_peaks"]].describe())

# --------------------------------------------------
# 6. Save h5ad files
# --------------------------------------------------

rna_out = processed_dir / "pbmc_rna_raw.h5ad"
atac_out = processed_dir / "pbmc_atac_raw.h5ad"

rna.write_h5ad(rna_out)
atac.write_h5ad(atac_out)

print("\nSaved RNA to:", rna_out)
print("Saved ATAC to:", atac_out)



from pathlib import Path
import os
import scanpy as sc

# --------------------------------------------------
# 0. Project directory
# --------------------------------------------------

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

processed_dir = project_dir / "data/processed"

rna_in = processed_dir / "pbmc_rna_raw.h5ad"
rna_out = processed_dir / "pbmc_rna_preprocessed.h5ad"

print("Loading:", rna_in)

rna = sc.read_h5ad(rna_in)

print("Before preprocessing:")
print(rna)

# --------------------------------------------------
# 1. Basic RNA QC
# --------------------------------------------------

# Human mitochondrial genes usually start with MT-
rna.var["mt"] = rna.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    rna,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

print("\nQC columns:")
print(rna.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].describe())

# --------------------------------------------------
# 2. Filter low-quality cells and genes
# --------------------------------------------------

# Keep cells with at least 200 detected genes
sc.pp.filter_cells(rna, min_genes=200)

# Keep genes detected in at least 3 cells
sc.pp.filter_genes(rna, min_cells=3)

# Remove cells with very high mitochondrial percentage
rna = rna[rna.obs["pct_counts_mt"] < 15].copy()

# Remove likely doublets / extreme high-count cells
rna = rna[rna.obs["total_counts"] < 40000].copy()

print("\nAfter filtering:")
print(rna)
print(rna.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].describe())

# --------------------------------------------------
# 3. Store raw counts
# --------------------------------------------------

rna.layers["counts"] = rna.X.copy()

# --------------------------------------------------
# 4. Normalize and log-transform
# --------------------------------------------------

sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)

# Store normalized/log data as raw for marker gene analysis
rna.raw = rna

# --------------------------------------------------
# 5. Highly variable genes
# --------------------------------------------------

sc.pp.highly_variable_genes(
    rna,
    n_top_genes=3000,
    flavor="seurat"
)

print("\nHighly variable genes:")
print(rna.var["highly_variable"].value_counts())

# Keep HVGs for PCA/clustering
rna_hvg = rna[:, rna.var["highly_variable"]].copy()

# --------------------------------------------------
# 6. PCA / neighbors / UMAP
# --------------------------------------------------

sc.pp.scale(rna_hvg, max_value=10)
sc.tl.pca(rna_hvg, svd_solver="arpack")

sc.pp.neighbors(rna_hvg, n_neighbors=15, n_pcs=30)
sc.tl.umap(rna_hvg)

# --------------------------------------------------
# 7. Clustering
# --------------------------------------------------

try:
    sc.tl.leiden(rna_hvg, resolution=0.6)
    cluster_key = "leiden"
except Exception as e:
    print("Leiden failed:", e)
    print("Trying Louvain instead...")
    sc.tl.louvain(rna_hvg, resolution=0.6)
    cluster_key = "louvain"

print("\nCluster counts:")
print(rna_hvg.obs[cluster_key].value_counts())

# --------------------------------------------------
# 8. Save
# --------------------------------------------------

rna_hvg.write_h5ad(rna_out)

print("\nSaved preprocessed RNA to:", rna_out)
print("Final object:")
print(rna_hvg)




"""Βήμα 2A — RNA QC πριν το filtering"""

from pathlib import Path
import os
import scanpy as sc

# --------------------------------------------------
# 0. Project directory
# --------------------------------------------------

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

processed_dir = project_dir / "data/processed"

rna_in = processed_dir / "pbmc_rna_raw.h5ad"

print("Loading:", rna_in)

rna = sc.read_h5ad(rna_in)

print("\nBefore QC:")
print(rna)
# --------------------------------------------------
# 1. Calculate QC metrics
# --------------------------------------------------

rna.var["mt"] = rna.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    rna,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

print("\nRNA QC summary:")
print(rna.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].describe())

print("\nLowest n_genes_by_counts:")
print(rna.obs["n_genes_by_counts"].sort_values().head(10))

print("\nHighest total_counts:")
print(rna.obs["total_counts"].sort_values(ascending=False).head(10))

print("\nHighest mitochondrial percentage:")
print(rna.obs["pct_counts_mt"].sort_values(ascending=False).head(10))










"""Βήμα 2B — RNA filtering με βάση QC thresholds"""

from pathlib import Path
import os
import scanpy as sc

# --------------------------------------------------
# 0. Project directory
# --------------------------------------------------

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

processed_dir = project_dir / "data/processed"

rna_in = processed_dir / "pbmc_rna_raw.h5ad"
rna_out = processed_dir / "pbmc_rna_filtered.h5ad"

print("Project directory:", Path.cwd())
print("Loading:", rna_in)

# --------------------------------------------------
# 1. Load raw RNA object
# --------------------------------------------------

rna = sc.read_h5ad(rna_in)

print("\nBefore filtering:")
print(rna)

# --------------------------------------------------
# 2. Calculate QC metrics
# --------------------------------------------------

rna.var["mt"] = rna.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    rna,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

print("\nQC before filtering:")
print(rna.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].describe())

# --------------------------------------------------
# 3. Store initial cell/gene numbers
# --------------------------------------------------

n_cells_before = rna.n_obs
n_genes_before = rna.n_vars

# --------------------------------------------------
# 4. Apply filtering thresholds
# --------------------------------------------------

# Remove cells with very few detected genes
sc.pp.filter_cells(rna, min_genes=200)

# Remove genes detected in very few cells
sc.pp.filter_genes(rna, min_cells=3)

# Remove cells with high mitochondrial percentage
rna = rna[rna.obs["pct_counts_mt"] < 15].copy()

# Remove extreme high-count cells, possible doublets
rna = rna[rna.obs["total_counts"] < 40000].copy()

# --------------------------------------------------
# 5. Report filtering result
# --------------------------------------------------

n_cells_after = rna.n_obs
n_genes_after = rna.n_vars

print("\nAfter filtering:")
print(rna)

print("\nCells before:", n_cells_before)
print("Cells after:", n_cells_after)
print("Cells removed:", n_cells_before - n_cells_after)

print("\nGenes before:", n_genes_before)
print("Genes after:", n_genes_after)
print("Genes removed:", n_genes_before - n_genes_after)

print("\nQC after filtering:")
print(rna.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].describe())

# --------------------------------------------------
# 6. Save filtered RNA object
# --------------------------------------------------

rna.write_h5ad(rna_out)

print("\nSaved filtered RNA to:", rna_out)








"""Βήμα 2C — RNA normalization, log-transform και highly variable genes"""

from pathlib import Path
import os
import scanpy as sc

# --------------------------------------------------
# 0. Project directory
# --------------------------------------------------

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

processed_dir = project_dir / "data/processed"

rna_in = processed_dir / "pbmc_rna_filtered.h5ad"
rna_out = processed_dir / "pbmc_rna_hvg.h5ad"

print("Project directory:", Path.cwd())
print("Loading:", rna_in)

# --------------------------------------------------
# 1. Load filtered RNA object
# --------------------------------------------------

rna = sc.read_h5ad(rna_in)

print("\nBefore normalization:")
print(rna)

# --------------------------------------------------
# 2. Store raw counts
# --------------------------------------------------

rna.layers["counts"] = rna.X.copy()

print("\nStored raw counts in rna.layers['counts'].")

# --------------------------------------------------
# 3. Normalize total counts per cell
# --------------------------------------------------

sc.pp.normalize_total(
    rna,
    target_sum=1e4
)

print("\nAfter normalize_total:")
print("Mean total normalized counts per cell should be close to 10000.")
print(rna.X.sum(axis=1).mean())

# --------------------------------------------------
# 4. Log-transform
# --------------------------------------------------

sc.pp.log1p(rna)

print("\nApplied log1p transformation.")

# Store normalized/log-transformed full gene matrix for marker genes later
rna.raw = rna

# --------------------------------------------------
# 5. Highly variable genes
# --------------------------------------------------

sc.pp.highly_variable_genes(
    rna,
    n_top_genes=3000,
    flavor="seurat"
)

print("\nHighly variable gene counts:")
print(rna.var["highly_variable"].value_counts())

# --------------------------------------------------
# 6. Keep only HVGs for PCA/UMAP/clustering
# --------------------------------------------------

rna_hvg = rna[:, rna.var["highly_variable"]].copy()

print("\nAfter keeping HVGs:")
print(rna_hvg)

# --------------------------------------------------
# 7. Save HVG object
# --------------------------------------------------

rna_hvg.write_h5ad(rna_out)

print("\nSaved HVG RNA object to:", rna_out)





"""Βήμα 2D — PCA, neighbors, UMAP και Leiden clustering"""

from pathlib import Path
import os
import scanpy as sc

# --------------------------------------------------
# 0. Project directory
# --------------------------------------------------

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

processed_dir = project_dir / "data/processed"

rna_in = processed_dir / "pbmc_rna_hvg.h5ad"
rna_out = processed_dir / "pbmc_rna_clustered.h5ad"

print("Project directory:", Path.cwd())
print("Loading:", rna_in)

# --------------------------------------------------
# 1. Load HVG RNA object
# --------------------------------------------------

rna = sc.read_h5ad(rna_in)

print("\nBefore PCA/UMAP/clustering:")
print(rna)

# --------------------------------------------------
# 2. Scale data
# --------------------------------------------------

print("\nScaling data...")
sc.pp.scale(rna, max_value=10)

# --------------------------------------------------
# 3. PCA
# --------------------------------------------------

print("Running PCA...")
sc.tl.pca(
    rna,
    svd_solver="arpack"
)

print("\nPCA completed.")
print("PCA coordinates shape:", rna.obsm["X_pca"].shape)

# --------------------------------------------------
# 4. Neighbors graph
# --------------------------------------------------

print("\nComputing neighbors graph...")
sc.pp.neighbors(
    rna,
    n_neighbors=15,
    n_pcs=30
)

# --------------------------------------------------
# 5. UMAP
# --------------------------------------------------

print("Running UMAP...")
sc.tl.umap(rna)

print("UMAP coordinates shape:", rna.obsm["X_umap"].shape)

# --------------------------------------------------
# 6. Leiden clustering
# --------------------------------------------------

print("\nRunning Leiden clustering...")

try:
    sc.tl.leiden(
        rna,
        resolution=0.6,
        key_added="leiden_0_6"
    )
    cluster_key = "leiden_0_6"

except Exception as e:
    print("\nLeiden failed.")
    print("Error:", e)
    print("\nTrying Louvain instead...")

    sc.tl.louvain(
        rna,
        resolution=0.6,
        key_added="louvain_0_6"
    )
    cluster_key = "louvain_0_6"

print("\nCluster key:", cluster_key)

print("\nCluster counts:")
print(rna.obs[cluster_key].value_counts().sort_index())

# --------------------------------------------------
# 7. Save clustered object
# --------------------------------------------------

rna.write_h5ad(rna_out)

print("\nSaved clustered RNA object to:", rna_out)
print("Final object:")
print(rna)







"""Βήμα 2E — UMAP visualization των RNA clusters"""

from pathlib import Path
import os
import scanpy as sc
import matplotlib.pyplot as plt

# --------------------------------------------------
# 0. Project directory
# --------------------------------------------------

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

processed_dir = project_dir / "data/processed"
figures_dir = project_dir / "results/figures"
figures_dir.mkdir(parents=True, exist_ok=True)

rna_in = processed_dir / "pbmc_rna_clustered.h5ad"

print("Project directory:", Path.cwd())
print("Loading:", rna_in)

# --------------------------------------------------
# 1. Load clustered RNA object
# --------------------------------------------------

rna = sc.read_h5ad(rna_in)

print(rna)

cluster_key = "leiden_0_6"

# --------------------------------------------------
# 2. Plot UMAP clusters
# --------------------------------------------------

sc.pl.umap(
    rna,
    color=cluster_key,
    legend_loc="on data",
    frameon=False,
    show=False
)

out_file = figures_dir / "pbmc_rna_umap_leiden_0_6.png"
plt.savefig(out_file, dpi=300, bbox_inches="tight")
plt.show()

print("Saved UMAP cluster plot to:", out_file)



"""Βήμα 2F — UMAP visualization βασικών PBMC marker genes"""

from pathlib import Path
import os
import scanpy as sc
import matplotlib.pyplot as plt

# --------------------------------------------------
# 0. Project directory
# --------------------------------------------------

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

processed_dir = project_dir / "data/processed"
figures_dir = project_dir / "results/figures"
figures_dir.mkdir(parents=True, exist_ok=True)

rna_in = processed_dir / "pbmc_rna_clustered.h5ad"

print("Project directory:", Path.cwd())
print("Loading:", rna_in)

# --------------------------------------------------
# 1. Load clustered RNA object
# --------------------------------------------------

rna = sc.read_h5ad(rna_in)

# --------------------------------------------------
# 2. Define PBMC marker genes
# --------------------------------------------------

marker_genes = [
    # T cells
    "CD3D", "CD3E",

    # CD4 T cells
    "CD4", "IL7R",

    # CD8 T cells
    "CD8A", "CD8B",

    # Naive T cells
    "CCR7", "LEF1", "TCF7",

    # Cytotoxic / effector T / NK
    "NKG7", "GZMB", "PRF1", "GNLY", "CCL5",

    # Monocytes
    "CD14", "LYZ", "S100A8", "S100A9", "FCN1", "MS4A7",

    # B cells
    "MS4A1", "CD79A",

    # NK cells
    "KLRD1", "NCR1",

    # Dendritic / antigen-presenting
    "FCER1A", "CST3"
]

# --------------------------------------------------
# 3. Keep only markers present in the object
# --------------------------------------------------

if rna.raw is not None:
    available_markers = [g for g in marker_genes if g in rna.raw.var_names]
else:
    available_markers = [g for g in marker_genes if g in rna.var_names]

missing_markers = [g for g in marker_genes if g not in available_markers]

print("\nAvailable markers:")
print(available_markers)

print("\nMissing markers:")
print(missing_markers)

# --------------------------------------------------
# 4. Plot marker genes on UMAP
# --------------------------------------------------

sc.pl.umap(
    rna,
    color=available_markers,
    use_raw=True if rna.raw is not None else False,
    frameon=False,
    show=False
)

out_file = figures_dir / "pbmc_rna_umap_marker_genes.png"
plt.savefig(out_file, dpi=300, bbox_inches="tight")
plt.show()

print("Saved marker gene UMAP plot to:", out_file)



"""Βήμα 2G — Marker gene dotplot ανά Leiden cluster για PBMC annotation"""

from pathlib import Path
import os
import scanpy as sc
import matplotlib.pyplot as plt

# --------------------------------------------------
# 0. Project directory
# --------------------------------------------------

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

processed_dir = project_dir / "data/processed"
figures_dir = project_dir / "results/figures"
figures_dir.mkdir(parents=True, exist_ok=True)

rna_in = processed_dir / "pbmc_rna_clustered.h5ad"

print("Project directory:", Path.cwd())
print("Loading:", rna_in)

# --------------------------------------------------
# 1. Load clustered RNA object
# --------------------------------------------------

rna = sc.read_h5ad(rna_in)

cluster_key = "leiden_0_6"

# --------------------------------------------------
# 2. Marker gene groups
# --------------------------------------------------

marker_dict = {
    "T_cells": ["CD3D", "CD3E"],
    "CD4_T": ["CD4", "IL7R"],
    "CD8_T": ["CD8A", "CD8B"],
    "Naive_T": ["CCR7", "LEF1", "TCF7"],
    "Cytotoxic_effector": ["NKG7", "GZMB", "PRF1", "GNLY", "CCL5"],
    "CD14_Monocytes": ["CD14", "LYZ", "S100A8", "S100A9", "FCN1", "MS4A7"],
    "B_cells": ["MS4A1", "CD79A"],
    "NK_cells": ["KLRD1", "NCR1"],
    "Dendritic_APC": ["FCER1A", "CST3"]
}

# --------------------------------------------------
# 3. Dotplot
# --------------------------------------------------

sc.pl.dotplot(
    rna,
    marker_dict,
    groupby=cluster_key,
    use_raw=True if rna.raw is not None else False,
    standard_scale="var",
    show=False
)

out_file = figures_dir / "pbmc_rna_marker_dotplot_by_cluster.png"
plt.savefig(out_file, dpi=300, bbox_inches="tight")
plt.show()

print("Saved dotplot to:", out_file)

"""0  → T cell
1  → T cell
2  → unclear / likely T or lymphocyte, needs table
3  → B cell / unclear, needs table
4  → CD14+ Monocyte
5  → Cytotoxic / NK-like
6  → T cell
7  → B cell
8  → CD8 effector / cytotoxic T
9  → Naive T cell
10 → Naive T cell
11 → CD14+ Monocyte / inflammatory monocyte
12 → CD14+ Monocyte
13 → Dendritic/APC-like or unclear
14 → CD14+ Monocyte
15 → Dendritic/APC-like
16 → Monocyte/APC-like
17 → T/cytotoxic-related
18 → unclear small cluster
19 → Dendritic/APC-like / unclear"""



"""Βήμα 2H — Mean marker expression ανά Leiden cluster"""

from pathlib import Path
import os
import scanpy as sc
import pandas as pd

# --------------------------------------------------
# 0. Project directory
# --------------------------------------------------

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

processed_dir = project_dir / "data/processed"
tables_dir = project_dir / "results/tables"
tables_dir.mkdir(parents=True, exist_ok=True)

rna_in = processed_dir / "pbmc_rna_clustered.h5ad"

print("Project directory:", Path.cwd())
print("Loading:", rna_in)

# --------------------------------------------------
# 1. Load clustered RNA object
# --------------------------------------------------

rna = sc.read_h5ad(rna_in)

cluster_key = "leiden_0_6"

marker_genes = [
    "CD3D", "CD3E",
    "CD4", "IL7R",
    "CD8A", "CD8B",
    "CCR7", "LEF1", "TCF7",
    "NKG7", "GZMB", "PRF1", "GNLY", "CCL5",
    "CD14", "LYZ", "S100A8", "S100A9", "FCN1", "MS4A7",
    "MS4A1", "CD79A",
    "KLRD1", "NCR1",
    "FCER1A", "CST3"
]

# --------------------------------------------------
# 2. Use raw normalized/log expression if available
# --------------------------------------------------

if rna.raw is not None:
    expr = rna.raw[:, marker_genes].to_adata()
else:
    expr = rna[:, marker_genes].copy()

expr.obs[cluster_key] = rna.obs[cluster_key].values

# --------------------------------------------------
# 3. Convert to DataFrame
# --------------------------------------------------

df = expr.to_df()
df[cluster_key] = expr.obs[cluster_key].values

mean_expr = df.groupby(cluster_key).mean()

# --------------------------------------------------
# 4. Save and print
# --------------------------------------------------

out_file = tables_dir / "pbmc_rna_mean_marker_expression_by_cluster.csv"
mean_expr.to_csv(out_file)

print("\nMean marker expression by cluster:")
print(mean_expr.round(2))

print("\nSaved table to:", out_file)



"""0  → CD14+ Monocytes
1  → CD4 Memory / IL7R+ T cells
2  → CD4 Naive T cells
3  → B cells
4  → CD14+ Monocytes
5  → NK cells / cytotoxic lymphocytes
6  → CD8 Effector T cells
7  → B cells
8  → CD8 Effector T cells
9  → CD14+ Monocytes / intermediate monocytes
10 → CD8 Naive T cells
11 → mixed T + Monocyte signal, πιθανό doublet/ambiguous
12 → CCL5+ cytotoxic / effector-like
13 → FCER1A+ dendritic/APC-like
14 → CD14+ Monocytes
15 → GZMB+ dendritic / cytotoxic-like small cluster
16 → mixed Monocyte + B signal, πιθανό doublet/ambiguous
17 → CD4 Naive / T cell
18 → T cell, probably memory/effector-like
19 → dendritic/APC-like small cluster"""


"""Βήμα 2I — Δημιουργία RNA cell type annotations από Leiden clusters"""

from pathlib import Path
import os
import scanpy as sc
import pandas as pd

# --------------------------------------------------
# 0. Project directory
# --------------------------------------------------

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

processed_dir = project_dir / "data/processed"
tables_dir = project_dir / "results/tables"
tables_dir.mkdir(parents=True, exist_ok=True)

rna_in = processed_dir / "pbmc_rna_clustered.h5ad"
rna_out = processed_dir / "pbmc_rna_annotated.h5ad"

print("Project directory:", Path.cwd())
print("Loading:", rna_in)

# --------------------------------------------------
# 1. Load clustered RNA object
# --------------------------------------------------

rna = sc.read_h5ad(rna_in)

cluster_key = "leiden_0_6"

# --------------------------------------------------
# 2. Manual cluster-to-cell-type annotation
# --------------------------------------------------

cluster_to_celltype = {
    "0": "CD14_Monocytes",
    "1": "CD4_Memory",
    "2": "CD4_Naive",
    "3": "Other_B_cells",
    "4": "CD14_Monocytes",
    "5": "Other_NK_or_Cytotoxic",
    "6": "CD8_Effector",
    "7": "Other_B_cells",
    "8": "CD8_Effector",
    "9": "Other_Monocyte_related",
    "10": "CD8_Naive",
    "11": "Other_Ambiguous",
    "12": "Other_Cytotoxic",
    "13": "Other_Dendritic_APC",
    "14": "CD14_Monocytes",
    "15": "Other_GZMB_high",
    "16": "Other_Ambiguous",
    "17": "CD4_Naive",
    "18": "Other_T_cells",
    "19": "Other_Dendritic_APC",
}

rna.obs["cell_type"] = rna.obs[cluster_key].map(cluster_to_celltype)

# --------------------------------------------------
# 3. Report annotation counts
# --------------------------------------------------

print("\nCell type counts:")
print(rna.obs["cell_type"].value_counts())

print("\nCluster to cell type table:")
cluster_table = (
    rna.obs[[cluster_key, "cell_type"]]
    .drop_duplicates()
    .sort_values(cluster_key)
)

print(cluster_table)

# --------------------------------------------------
# 4. Save annotation table
# --------------------------------------------------

counts_out = tables_dir / "pbmc_rna_cell_type_counts.csv"
cluster_map_out = tables_dir / "pbmc_rna_cluster_to_cell_type_mapping.csv"

rna.obs["cell_type"].value_counts().to_csv(counts_out)
cluster_table.to_csv(cluster_map_out, index=False)

print("\nSaved cell type counts to:", counts_out)
print("Saved cluster mapping to:", cluster_map_out)

# --------------------------------------------------
# 5. Save annotated RNA object
# --------------------------------------------------

rna.write_h5ad(rna_out)

print("\nSaved annotated RNA object to:", rna_out)
print("Final object:")
print(rna)


"""Βήμα 2J — UMAP visualization των annotated RNA cell types"""

from pathlib import Path
import os
import scanpy as sc
import matplotlib.pyplot as plt

# --------------------------------------------------
# 0. Project directory
# --------------------------------------------------

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

processed_dir = project_dir / "data/processed"
figures_dir = project_dir / "results/figures"
figures_dir.mkdir(parents=True, exist_ok=True)

rna_in = processed_dir / "pbmc_rna_annotated.h5ad"

print("Project directory:", Path.cwd())
print("Loading:", rna_in)

# --------------------------------------------------
# 1. Load annotated RNA object
# --------------------------------------------------

rna = sc.read_h5ad(rna_in)

# --------------------------------------------------
# 2. Plot UMAP by cell type
# --------------------------------------------------

sc.pl.umap(
    rna,
    color="cell_type",
    legend_loc="right margin",
    frameon=False,
    show=False
)

out_file = figures_dir / "pbmc_rna_umap_cell_types.png"
plt.savefig(out_file, dpi=300, bbox_inches="tight")
plt.show()

print("Saved annotated UMAP to:", out_file)
"""Βήμα 2K — UMAP visualization με Leiden clusters και annotated cell types"""

from pathlib import Path
import os
import scanpy as sc
import matplotlib.pyplot as plt

# --------------------------------------------------
# 0. Project directory
# --------------------------------------------------

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

processed_dir = project_dir / "data/processed"
figures_dir = project_dir / "results/figures"
figures_dir.mkdir(parents=True, exist_ok=True)

rna_in = processed_dir / "pbmc_rna_annotated.h5ad"

print("Project directory:", Path.cwd())
print("Loading:", rna_in)

# --------------------------------------------------
# 1. Load annotated RNA object
# --------------------------------------------------

rna = sc.read_h5ad(rna_in)

print(rna)
print("\nColumns in obs:")
print(rna.obs.columns)

# --------------------------------------------------
# 2. Plot Leiden cluster numbers and cell type annotation
# --------------------------------------------------

sc.pl.umap(
    rna,
    color=["leiden_0_6", "cell_type"],
    legend_loc="on data",
    frameon=False,
    show=False
)

out_file = figures_dir / "pbmc_rna_umap_leiden_and_cell_type.png"
plt.savefig(out_file, dpi=300, bbox_inches="tight")
plt.show()

print("Saved combined UMAP to:", out_file)




"""Βήμα 2M — Cluster-aware RNA cell type annotation"""

from pathlib import Path
import os
import scanpy as sc

# --------------------------------------------------
# 0. Project directory
# --------------------------------------------------

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

processed_dir = project_dir / "data/processed"
tables_dir = project_dir / "results/tables"
tables_dir.mkdir(parents=True, exist_ok=True)

rna_in = processed_dir / "pbmc_rna_clustered.h5ad"
rna_out = processed_dir / "pbmc_rna_annotated_clusteraware.h5ad"

print("Project directory:", Path.cwd())
print("Loading:", rna_in)

# --------------------------------------------------
# 1. Load clustered RNA object
# --------------------------------------------------

rna = sc.read_h5ad(rna_in)

cluster_key = "leiden_0_6"

# --------------------------------------------------
# 2. Cluster-aware annotation
# --------------------------------------------------

cluster_to_celltype_detailed = {
    "0": "CD14_Monocytes_c0",
    "1": "CD4_Memory_c1",
    "2": "CD4_Naive_c2",
    "3": "B_cells_c3",
    "4": "CD14_Monocytes_c4",
    "5": "NK_or_Cytotoxic_c5",
    "6": "CD8_Effector_c6",
    "7": "B_cells_c7",
    "8": "CD8_Effector_c8",
    "9": "Monocyte_related_c9",
    "10": "CD8_Naive_c10",
    "11": "Ambiguous_c11",
    "12": "Cytotoxic_c12",
    "13": "Dendritic_APC_c13",
    "14": "CD14_Monocytes_c14",
    "15": "GZMB_high_c15",
    "16": "Ambiguous_c16",
    "17": "CD4_Naive_c17",
    "18": "T_cells_c18",
    "19": "Dendritic_APC_c19",
}

rna.obs["cell_type_detailed"] = rna.obs[cluster_key].map(cluster_to_celltype_detailed)

# --------------------------------------------------
# 3. Also keep broad biological labels
# --------------------------------------------------

cluster_to_celltype_broad = {
    "0": "CD14_Monocytes",
    "1": "CD4_Memory",
    "2": "CD4_Naive",
    "3": "B_cells",
    "4": "CD14_Monocytes",
    "5": "NK_or_Cytotoxic",
    "6": "CD8_Effector",
    "7": "B_cells",
    "8": "CD8_Effector",
    "9": "Monocyte_related",
    "10": "CD8_Naive",
    "11": "Ambiguous",
    "12": "Cytotoxic",
    "13": "Dendritic_APC",
    "14": "CD14_Monocytes",
    "15": "GZMB_high",
    "16": "Ambiguous",
    "17": "CD4_Naive",
    "18": "T_cells",
    "19": "Dendritic_APC",
}

rna.obs["cell_type_broad"] = rna.obs[cluster_key].map(cluster_to_celltype_broad)

# --------------------------------------------------
# 4. Print summaries
# --------------------------------------------------

print("\nDetailed cell type counts:")
print(rna.obs["cell_type_detailed"].value_counts())

print("\nBroad cell type counts:")
print(rna.obs["cell_type_broad"].value_counts())

summary = (
    rna.obs
    .groupby([cluster_key, "cell_type_detailed", "cell_type_broad"])
    .size()
    .reset_index(name="n_cells")
    .sort_values(cluster_key)
)

print("\nCluster annotation summary:")
print(summary)

summary_out = tables_dir / "pbmc_rna_clusteraware_annotation_summary.csv"
summary.to_csv(summary_out, index=False)

# --------------------------------------------------
# 5. Save object
# --------------------------------------------------

rna.write_h5ad(rna_out)

print("\nSaved cluster-aware annotated RNA object to:", rna_out)
print("Saved summary to:", summary_out)




"""Βήμα 2N — UMAP με cluster-aware RNA annotations"""

from pathlib import Path
import os
import scanpy as sc
import matplotlib.pyplot as plt

# --------------------------------------------------
# 0. Project directory
# --------------------------------------------------

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
os.chdir(project_dir)

processed_dir = project_dir / "data/processed"
figures_dir = project_dir / "results/figures"
figures_dir.mkdir(parents=True, exist_ok=True)

rna_in = processed_dir / "pbmc_rna_annotated_clusteraware.h5ad"

print("Project directory:", Path.cwd())
print("Loading:", rna_in)

# --------------------------------------------------
# 1. Load object
# --------------------------------------------------

rna = sc.read_h5ad(rna_in)

# --------------------------------------------------
# 2. Plot Leiden and detailed cell type
# --------------------------------------------------

sc.pl.umap(
    rna,
    color=["leiden_0_6", "cell_type_detailed"],
    legend_loc="on data",
    legend_fontsize=8,
   legend_fontweight="bold",
    frameon=False,
    show=False
)

out_file = figures_dir / "pbmc_rna_umap_clusteraware_annotation.png"
plt.savefig(out_file, dpi=300, bbox_inches="tight")
plt.show()

print("Saved UMAP to:", out_file)



"""Βήμα 20 — Προσθήκη final RNA annotations σε 3 επίπεδα"""

import scanpy as sc
import pandas as pd
from pathlib import Path

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

rna_path = project_dir / "data/processed/pbmc_rna_clustered.h5ad"
out_path = project_dir / "data/processed/pbmc_rna_annotated_final.h5ad"

rna = sc.read_h5ad(rna_path)

# Ensure Leiden labels are strings
rna.obs["leiden_0_6"] = rna.obs["leiden_0_6"].astype(str)

detailed_mapping = {
    "0":  "CD14_Monocytes_c0",
    "1":  "CD4_Memory_c1",
    "2":  "CD4_Naive_c2",
    "3":  "B_cells_c3",
    "4":  "CD14_Monocytes_c4",
    "5":  "NK_or_Cytotoxic_c5",
    "6":  "CD8_Effector_c6",
    "7":  "B_cells_c7",
    "8":  "CD8_Effector_c8",
    "9":  "Monocyte_related_c9",
    "10": "CD8_Naive_c10",
    "11": "Ambiguous_c11",
    "12": "Cytotoxic_c12",
    "13": "Dendritic_APC_c13",
    "14": "CD14_Monocytes_c14",
    "15": "GZMB_high_c15",
    "16": "Ambiguous_c16",
    "17": "CD4_Naive_c17",
    "18": "T_cells_c18",
    "19": "Dendritic_APC_c19",
}

broad_mapping = {
    "0":  "CD14_Monocytes",
    "1":  "CD4_Memory",
    "2":  "CD4_Naive",
    "3":  "B_cells",
    "4":  "CD14_Monocytes",
    "5":  "NK_or_Cytotoxic",
    "6":  "CD8_Effector",
    "7":  "B_cells",
    "8":  "CD8_Effector",
    "9":  "Monocyte_related",
    "10": "CD8_Naive",
    "11": "Ambiguous",
    "12": "Cytotoxic",
    "13": "Dendritic_APC",
    "14": "CD14_Monocytes",
    "15": "GZMB_high",
    "16": "Ambiguous",
    "17": "CD4_Naive",
    "18": "T_cells",
    "19": "Dendritic_APC",
}

rna.obs["cell_type_detailed"] = rna.obs["leiden_0_6"].map(detailed_mapping)
rna.obs["cell_type_broad"] = rna.obs["leiden_0_6"].map(broad_mapping)

# Sanity checks
missing_detailed = rna.obs["cell_type_detailed"].isna().sum()
missing_broad = rna.obs["cell_type_broad"].isna().sum()

print("Missing detailed annotations:", missing_detailed)
print("Missing broad annotations:", missing_broad)

print("\nCell counts by detailed annotation:")
print(rna.obs["cell_type_detailed"].value_counts())

print("\nCell counts by broad annotation:")
print(rna.obs["cell_type_broad"].value_counts())

rna.write_h5ad(out_path)

print(f"\nSaved final annotated RNA object to:\n{out_path}")
print(rna)


"""Βήμα 21 — UMAP έλεγχος RNA annotations"""

import scanpy as sc
from pathlib import Path

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
rna_path = project_dir / "data/processed/pbmc_rna_annotated_final.h5ad"

rna = sc.read_h5ad(rna_path)

sc.pl.umap(
    rna,
    color="leiden_0_6",
    legend_loc="on data",
    legend_fontsize=6,
    legend_fontweight="normal",
    frameon=False,
    title="RNA Leiden clusters"
)

sc.pl.umap(
    rna,
    color="cell_type_broad",
    legend_loc="on data",
    legend_fontsize=6,
    legend_fontweight="normal",
    frameon=False,
    title="RNA broad cell-type annotation"
)

sc.pl.umap(
    rna,
    color="cell_type_detailed",
    legend_loc="right margin",
    frameon=False,
    title="RNA detailed cluster-aware annotation"
)





"""Βήμα 23 — Δημιουργία RNA subset για τα 5 scMultiomeGRN PBMC cell types"""

import scanpy as sc
import pandas as pd
from pathlib import Path

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

rna_path = project_dir / "data/processed/pbmc_rna_annotated_final.h5ad"
out_path = project_dir / "data/processed/pbmc_rna_scmultiomegrn_5types.h5ad"
csv_out = project_dir / "data/processed/pbmc_rna_scmultiomegrn_5types_cell_counts.csv"

rna = sc.read_h5ad(rna_path)

target_mapping = {
    "CD14_Monocytes": "CD14+ Monocytes",
    "CD4_Memory": "CD4.Memory",
    "CD4_Naive": "CD4.Naive",
    "CD8_Effector": "CD8.effector",
    "CD8_Naive": "CD8.Naive",
}

target_broad = list(target_mapping.keys())

rna_5 = rna[rna.obs["cell_type_broad"].isin(target_broad)].copy()
rna_5.obs["scmultiomegrn_cell_type"] = rna_5.obs["cell_type_broad"].map(target_mapping)

print("RNA 5-type subset:")
print(rna_5)

print("\nCell counts by scMultiomeGRN label:")
counts = rna_5.obs["scmultiomegrn_cell_type"].value_counts()
print(counts)

print("\nCell counts by Leiden cluster:")
print(rna_5.obs["leiden_0_6"].value_counts().sort_index())

counts.to_csv(csv_out, header=["n_cells"])

rna_5.write_h5ad(out_path)

print(f"\nSaved RNA 5-type subset to:\n{out_path}")
print(f"Saved cell counts to:\n{csv_out}")






"""Βήμα 24 — Marker validation για τα 5 scMultiomeGRN RNA cell types"""

import scanpy as sc
from pathlib import Path

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
rna_path = project_dir / "data/processed/pbmc_rna_scmultiomegrn_5types.h5ad"

rna_5 = sc.read_h5ad(rna_path)

marker_genes = [
    # CD14+ monocytes
    "CD14", "LYZ", "S100A8", "S100A9", "FCN1", "MS4A7",
    
    # CD4 naive / memory
    "CD3D", "CD3E", "CD4", "IL7R", "CCR7", "LEF1", "TCF7",
    
    # CD8 naive / effector
    "CD8A", "CD8B", "NKG7", "GZMB", "PRF1", "GNLY", "CCL5"
]

available_markers = [g for g in marker_genes if g in rna_5.var_names]
missing_markers = [g for g in marker_genes if g not in rna_5.var_names]

print("Available markers:", available_markers)
print("Missing markers:", missing_markers)

sc.pl.dotplot(
    rna_5,
    var_names=available_markers,
    groupby="scmultiomegrn_cell_type",
    standard_scale="var",
    dendrogram=False,
    title="Marker validation — scMultiomeGRN 5 RNA cell types"
)

sc.pl.umap(
    rna_5,
    color="scmultiomegrn_cell_type",
    legend_loc="on data",
    legend_fontsize=10,
    legend_fontweight="bold",
    frameon=False,
    title="RNA subset — scMultiomeGRN 5 PBMC cell types"
)






"""Βήμα 25 — Έλεγχος barcode overlap μεταξύ RNA και ATAC"""

import scanpy as sc
from pathlib import Path

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

rna_path = project_dir / "data/processed/pbmc_rna_annotated_final.h5ad"
atac_path = project_dir / "data/processed/pbmc_atac_raw.h5ad"

rna = sc.read_h5ad(rna_path)
atac = sc.read_h5ad(atac_path)

rna_barcodes = set(rna.obs_names)
atac_barcodes = set(atac.obs_names)

overlap = rna_barcodes.intersection(atac_barcodes)

print("RNA cells:", rna.n_obs)
print("ATAC cells:", atac.n_obs)
print("Overlapping barcodes:", len(overlap))
print("RNA overlap fraction:", len(overlap) / rna.n_obs)
print("ATAC overlap fraction:", len(overlap) / atac.n_obs)

print("\nExample RNA barcodes:")
print(list(rna.obs_names[:5]))

print("\nExample ATAC barcodes:")
print(list(atac.obs_names[:5]))


""" ScATAC change"""

"""Βήμα 26 — ATAC QC inspection πριν το filtering"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

atac_path = project_dir / "data/processed/pbmc_atac_raw.h5ad"
out_csv = project_dir / "data/processed/pbmc_atac_qc_summary_before_filtering.csv"

atac = sc.read_h5ad(atac_path)

print(atac)
print("\nobs columns:")
print(atac.obs.columns.tolist())

print("\nvar columns:")
print(atac.var.columns.tolist())

# Basic ATAC QC metrics
atac.obs["total_counts"] = np.asarray(atac.X.sum(axis=1)).ravel()
atac.obs["n_peaks_by_counts"] = np.asarray((atac.X > 0).sum(axis=1)).ravel()

atac.var["n_cells_by_counts"] = np.asarray((atac.X > 0).sum(axis=0)).ravel()
atac.var["total_counts"] = np.asarray(atac.X.sum(axis=0)).ravel()

qc_summary = pd.DataFrame({
    "metric": [
        "n_cells",
        "n_peaks",
        "mean_total_counts_per_cell",
        "median_total_counts_per_cell",
        "min_total_counts_per_cell",
        "max_total_counts_per_cell",
        "mean_n_peaks_by_counts_per_cell",
        "median_n_peaks_by_counts_per_cell",
        "min_n_peaks_by_counts_per_cell",
        "max_n_peaks_by_counts_per_cell",
        "mean_cells_per_peak",
        "median_cells_per_peak",
        "min_cells_per_peak",
        "max_cells_per_peak",
    ],
    "value": [
        atac.n_obs,
        atac.n_vars,
        atac.obs["total_counts"].mean(),
        atac.obs["total_counts"].median(),
        atac.obs["total_counts"].min(),
        atac.obs["total_counts"].max(),
        atac.obs["n_peaks_by_counts"].mean(),
        atac.obs["n_peaks_by_counts"].median(),
        atac.obs["n_peaks_by_counts"].min(),
        atac.obs["n_peaks_by_counts"].max(),
        atac.var["n_cells_by_counts"].mean(),
        atac.var["n_cells_by_counts"].median(),
        atac.var["n_cells_by_counts"].min(),
        atac.var["n_cells_by_counts"].max(),
    ]
})

print("\nATAC QC summary before filtering:")
print(qc_summary)

qc_summary.to_csv(out_csv, index=False)
print(f"\nSaved QC summary to:\n{out_csv}")

# Histograms
plt.figure(figsize=(6, 4))
plt.hist(atac.obs["total_counts"], bins=100)
plt.xlabel("Total ATAC counts per cell")
plt.ylabel("Number of cells")
plt.title("ATAC total counts per cell")
plt.tight_layout()
plt.show()

plt.figure(figsize=(6, 4))
plt.hist(atac.obs["n_peaks_by_counts"], bins=100)
plt.xlabel("Accessible peaks per cell")
plt.ylabel("Number of cells")
plt.title("ATAC accessible peaks per cell")
plt.tight_layout()
plt.show()

plt.figure(figsize=(6, 4))
plt.hist(atac.var["n_cells_by_counts"], bins=100)
plt.xlabel("Cells per peak")
plt.ylabel("Number of peaks")
plt.title("ATAC cells per peak")
plt.tight_layout()
plt.show()



"""Βήμα 27 — ATAC QC quantiles για επιλογή filtering thresholds"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
atac_path = project_dir / "data/processed/pbmc_atac_raw.h5ad"

atac = sc.read_h5ad(atac_path)

atac.obs["total_counts"] = np.asarray(atac.X.sum(axis=1)).ravel()
atac.obs["n_peaks_by_counts"] = np.asarray((atac.X > 0).sum(axis=1)).ravel()
atac.var["n_cells_by_counts"] = np.asarray((atac.X > 0).sum(axis=0)).ravel()

cell_metrics = ["total_counts", "n_peaks_by_counts"]
peak_metrics = ["n_cells_by_counts"]

cell_quantiles = atac.obs[cell_metrics].quantile(
    [0.00, 0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.98, 0.99, 1.00]
)

peak_quantiles = atac.var[peak_metrics].quantile(
    [0.00, 0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.98, 0.99, 1.00]
)

print("Cell-level ATAC QC quantiles:")
print(cell_quantiles)

print("\nPeak-level ATAC QC quantiles:")
print(peak_quantiles)





"""Βήμα 28 — ATAC filtering με βάση QC quantiles"""

import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

atac_path = project_dir / "data/processed/pbmc_atac_raw.h5ad"
out_path = project_dir / "data/processed/pbmc_atac_filtered.h5ad"
summary_out = project_dir / "data/processed/pbmc_atac_filtering_summary.csv"

atac = sc.read_h5ad(atac_path)

print("ATAC before filtering:")
print(atac)

# Recompute QC metrics
atac.obs["total_counts"] = np.asarray(atac.X.sum(axis=1)).ravel()
atac.obs["n_peaks_by_counts"] = np.asarray((atac.X > 0).sum(axis=1)).ravel()
atac.var["n_cells_by_counts"] = np.asarray((atac.X > 0).sum(axis=0)).ravel()
atac.var["total_counts"] = np.asarray(atac.X.sum(axis=0)).ravel()

n_cells_before = atac.n_obs
n_peaks_before = atac.n_vars

# Filtering thresholds
min_peaks_per_cell = 1000
max_counts_per_cell = 60000
min_cells_per_peak = 20

cell_filter = (
    (atac.obs["n_peaks_by_counts"] >= min_peaks_per_cell) &
    (atac.obs["total_counts"] <= max_counts_per_cell)
)

peak_filter = atac.var["n_cells_by_counts"] >= min_cells_per_peak

atac_filtered = atac[cell_filter, peak_filter].copy()

# Recompute metrics after filtering
atac_filtered.obs["total_counts"] = np.asarray(atac_filtered.X.sum(axis=1)).ravel()
atac_filtered.obs["n_peaks_by_counts"] = np.asarray((atac_filtered.X > 0).sum(axis=1)).ravel()
atac_filtered.var["n_cells_by_counts"] = np.asarray((atac_filtered.X > 0).sum(axis=0)).ravel()
atac_filtered.var["total_counts"] = np.asarray(atac_filtered.X.sum(axis=0)).ravel()

n_cells_after = atac_filtered.n_obs
n_peaks_after = atac_filtered.n_vars

summary = pd.DataFrame({
    "metric": [
        "cells_before",
        "cells_after",
        "cells_removed",
        "peaks_before",
        "peaks_after",
        "peaks_removed",
        "min_peaks_per_cell",
        "max_counts_per_cell",
        "min_cells_per_peak",
        "mean_total_counts_after",
        "median_total_counts_after",
        "max_total_counts_after",
        "mean_n_peaks_by_counts_after",
        "median_n_peaks_by_counts_after",
        "max_n_peaks_by_counts_after",
    ],
    "value": [
        n_cells_before,
        n_cells_after,
        n_cells_before - n_cells_after,
        n_peaks_before,
        n_peaks_after,
        n_peaks_before - n_peaks_after,
        min_peaks_per_cell,
        max_counts_per_cell,
        min_cells_per_peak,
        atac_filtered.obs["total_counts"].mean(),
        atac_filtered.obs["total_counts"].median(),
        atac_filtered.obs["total_counts"].max(),
        atac_filtered.obs["n_peaks_by_counts"].mean(),
        atac_filtered.obs["n_peaks_by_counts"].median(),
        atac_filtered.obs["n_peaks_by_counts"].max(),
    ]
})

print("\nATAC after filtering:")
print(atac_filtered)

print("\nFiltering summary:")
print(summary)

summary.to_csv(summary_out, index=False)
atac_filtered.write_h5ad(out_path)

print(f"\nSaved filtered ATAC object to:\n{out_path}")
print(f"Saved filtering summary to:\n{summary_out}")




"""Βήμα 29 — ATAC TF-IDF normalization και LSI"""

import scanpy as sc
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize
from pathlib import Path

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

atac_path = project_dir / "data/processed/pbmc_atac_filtered.h5ad"
out_path = project_dir / "data/processed/pbmc_atac_lsi.h5ad"

atac = sc.read_h5ad(atac_path)

print("Loaded filtered ATAC:")
print(atac)

# Save raw binary/count matrix
atac.layers["counts"] = atac.X.copy()

# Make sure matrix is sparse CSR
X = atac.X
if not sparse.issparse(X):
    X = sparse.csr_matrix(X)
else:
    X = X.tocsr()

# Binarize accessibility matrix
X_bin = X.copy()
X_bin.data = np.ones_like(X_bin.data)

# Term frequency: counts per cell normalized by total accessible peaks per cell
cell_sums = np.asarray(X_bin.sum(axis=1)).ravel()
cell_sums[cell_sums == 0] = 1

tf = X_bin.multiply(1 / cell_sums[:, None])

# Inverse document frequency
n_cells = X_bin.shape[0]
peak_sums = np.asarray(X_bin.sum(axis=0)).ravel()
idf = np.log(1 + n_cells / (1 + peak_sums))

tfidf = tf.multiply(idf)

# Scale TF-IDF for numerical stability
tfidf = tfidf * 1e4

# LSI using truncated SVD
n_components = 50
svd = TruncatedSVD(n_components=n_components, random_state=42)
X_lsi = svd.fit_transform(tfidf)

# Usually LSI component 1 is sequencing-depth-associated.
# Keep all components for now, but later neighbors can use components 2:50.
atac.obsm["X_lsi"] = X_lsi
atac.uns["lsi_variance_ratio"] = svd.explained_variance_ratio_
atac.uns["lsi_singular_values"] = svd.singular_values_

print("LSI shape:", atac.obsm["X_lsi"].shape)
print("First 10 explained variance ratios:")
print(atac.uns["lsi_variance_ratio"][:10])

atac.write_h5ad(out_path)

print(f"\nSaved ATAC LSI object to:\n{out_path}")



"""Βήμα 30 — ATAC neighbors, UMAP και Leiden clustering με LSI components 2–50"""

import scanpy as sc
from pathlib import Path

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

atac_path = project_dir / "data/processed/pbmc_atac_lsi.h5ad"
out_path = project_dir / "data/processed/pbmc_atac_clustered.h5ad"

atac = sc.read_h5ad(atac_path)

# Use LSI components 2–50, excluding component 1
atac.obsm["X_lsi_2_50"] = atac.obsm["X_lsi"][:, 1:50]

sc.pp.neighbors(
    atac,
    n_neighbors=15,
    use_rep="X_lsi_2_50",
    metric="cosine"
)

sc.tl.umap(atac)

sc.tl.leiden(
    atac,
    resolution=0.6,
    key_added="leiden_0_6"
)

print("ATAC clustered object:")
print(atac)

print("\nATAC Leiden cluster counts:")
print(atac.obs["leiden_0_6"].value_counts().sort_index())

sc.pl.umap(
    atac,
    color="leiden_0_6",
    legend_loc="on data",
    legend_fontsize=6,
    legend_fontweight="normal",
    frameon=False,
    title="ATAC Leiden clusters"
)

atac.write_h5ad(out_path)

print(f"\nSaved clustered ATAC object to:\n{out_path}")









"""Βήμα 31 — Έλεγχος ATAC peak annotation και συμβατότητας με τα peak IDs"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

atac_path = project_dir / "data/processed/pbmc_atac_clustered.h5ad"

peak_annotation_path = (
    project_dir /
    "data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_peak_annotation.tsv"
)

# ============================================================
# Load ATAC object και peak annotation table
# ============================================================

atac = sc.read_h5ad(atac_path)
peak_anno = pd.read_csv(peak_annotation_path, sep="\t")

print("ATAC object:")
print(atac)

print("\nPeak annotation shape:")
print(peak_anno.shape)

print("\nPeak annotation columns:")
print(peak_anno.columns.tolist())

print("\nFirst rows of peak annotation:")
print(peak_anno.head())

print("\nExample ATAC peak IDs:")
print(list(atac.var_names[:10]))

# ============================================================
# Αυτόματος εντοπισμός peak column
# ============================================================

possible_peak_cols = [
    "peak", "Peak", "peak_id", "PeakID", "peak_name", "PeakName",
    "feature", "Feature"
]

peak_col = None

for col in possible_peak_cols:
    if col in peak_anno.columns:
        peak_col = col
        break

# Αν δεν υπάρχει έτοιμη peak column, δοκιμάζουμε να φτιάξουμε peak ID από chrom/start/end
if peak_col is None:
    possible_chrom_cols = ["chrom", "chr", "seqnames", "Chromosome"]
    possible_start_cols = ["start", "Start"]
    possible_end_cols = ["end", "End"]

    chrom_col = next((c for c in possible_chrom_cols if c in peak_anno.columns), None)
    start_col = next((c for c in possible_start_cols if c in peak_anno.columns), None)
    end_col = next((c for c in possible_end_cols if c in peak_anno.columns), None)

    if chrom_col is not None and start_col is not None and end_col is not None:
        peak_anno["peak_id_constructed"] = (
            peak_anno[chrom_col].astype(str) + ":" +
            peak_anno[start_col].astype(str) + "-" +
            peak_anno[end_col].astype(str)
        )
        peak_col = "peak_id_constructed"

if peak_col is None:
    raise ValueError(
        "Δεν βρέθηκε peak column. Δες τα printed columns του peak_annotation.tsv."
    )

print("\nDetected peak column:")
print(peak_col)

# ============================================================
# Αυτόματος εντοπισμός gene column
# ============================================================

possible_gene_cols = [
    "gene", "Gene", "gene_name", "GeneName", "gene_symbol",
    "symbol", "Symbol", "name", "Name"
]

gene_col = None

for col in possible_gene_cols:
    if col in peak_anno.columns:
        gene_col = col
        break

if gene_col is None:
    raise ValueError(
        "Δεν βρέθηκε gene column. Δες τα printed columns του peak_annotation.tsv."
    )

print("\nDetected gene column:")
print(gene_col)

# ============================================================
# Έλεγχος overlap μεταξύ ATAC peaks και annotated peaks
# ============================================================

atac_peaks = set(atac.var_names.astype(str))
annotated_peaks = set(peak_anno[peak_col].astype(str))

overlap = atac_peaks.intersection(annotated_peaks)

print("\nATAC peaks:", len(atac_peaks))
print("Annotated peaks:", len(annotated_peaks))
print("Overlapping peaks:", len(overlap))
print("ATAC peak overlap fraction:", len(overlap) / len(atac_peaks))

print("\nExample overlapping peaks:")
print(list(overlap)[:10])





"""Βήμα 32 — Δημιουργία ATAC gene activity matrix από peak-to-gene annotation"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from anndata import AnnData
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

atac_path = project_dir / "data/processed/pbmc_atac_clustered.h5ad"

peak_annotation_path = (
    project_dir /
    "data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_peak_annotation.tsv"
)

out_path = project_dir / "data/processed/pbmc_atac_gene_activity.h5ad"

# ============================================================
# Load ATAC object και peak annotation
# ============================================================

atac = sc.read_h5ad(atac_path)
peak_anno = pd.read_csv(peak_annotation_path, sep="\t")

print("Loaded ATAC:")
print(atac)

print("\nPeak annotation columns:")
print(peak_anno.columns.tolist())

# ============================================================
# Detect peak column
# ============================================================

possible_peak_cols = [
    "peak", "Peak", "peak_id", "PeakID", "peak_name", "PeakName",
    "feature", "Feature"
]

peak_col = None

for col in possible_peak_cols:
    if col in peak_anno.columns:
        peak_col = col
        break

if peak_col is None:
    possible_chrom_cols = ["chrom", "chr", "seqnames", "Chromosome"]
    possible_start_cols = ["start", "Start"]
    possible_end_cols = ["end", "End"]

    chrom_col = next((c for c in possible_chrom_cols if c in peak_anno.columns), None)
    start_col = next((c for c in possible_start_cols if c in peak_anno.columns), None)
    end_col = next((c for c in possible_end_cols if c in peak_anno.columns), None)

    if chrom_col is not None and start_col is not None and end_col is not None:
        peak_anno["peak_id_constructed"] = (
            peak_anno[chrom_col].astype(str) + ":" +
            peak_anno[start_col].astype(str) + "-" +
            peak_anno[end_col].astype(str)
        )
        peak_col = "peak_id_constructed"

if peak_col is None:
    raise ValueError("Δεν βρέθηκε peak column στο peak annotation file.")

# ============================================================
# Detect gene column
# ============================================================

possible_gene_cols = [
    "gene", "Gene", "gene_name", "GeneName", "gene_symbol",
    "symbol", "Symbol", "name", "Name"
]

gene_col = None

for col in possible_gene_cols:
    if col in peak_anno.columns:
        gene_col = col
        break

if gene_col is None:
    raise ValueError("Δεν βρέθηκε gene column στο peak annotation file.")

print("\nUsing peak column:", peak_col)
print("Using gene column:", gene_col)

# ============================================================
# Καθαρισμός annotation table
# ============================================================

peak_anno = peak_anno[[peak_col, gene_col]].copy()
peak_anno.columns = ["peak_id", "gene_name"]

peak_anno["peak_id"] = peak_anno["peak_id"].astype(str)
peak_anno["gene_name"] = peak_anno["gene_name"].astype(str)

# Αφαίρεση άδειων ή άκυρων gene annotations
invalid_gene_values = ["", ".", "nan", "NaN", "None", "none", "NA"]

peak_anno = peak_anno[
    ~peak_anno["gene_name"].isin(invalid_gene_values)
].copy()

# Κρατάμε μόνο peaks που υπάρχουν στο filtered/clustered ATAC object
peak_anno = peak_anno[
    peak_anno["peak_id"].isin(atac.var_names.astype(str))
].copy()

# Αφαίρεση διπλότυπων peak-gene pairs
peak_anno = peak_anno.drop_duplicates(["peak_id", "gene_name"])

print("\nUsable peak-gene links:", peak_anno.shape[0])
print("Unique annotated peaks:", peak_anno["peak_id"].nunique())
print("Unique genes:", peak_anno["gene_name"].nunique())

# ============================================================
# Δημιουργία sparse peak × gene mapping matrix
# ============================================================

peak_to_idx = pd.Series(
    np.arange(atac.n_vars),
    index=atac.var_names.astype(str)
)

gene_names = np.array(sorted(peak_anno["gene_name"].unique()))
gene_to_idx = pd.Series(
    np.arange(len(gene_names)),
    index=gene_names
)

peak_indices = peak_anno["peak_id"].map(peak_to_idx).to_numpy()
gene_indices = peak_anno["gene_name"].map(gene_to_idx).to_numpy()

data = np.ones(len(peak_anno), dtype=np.float32)

peak_gene_matrix = sparse.csr_matrix(
    (data, (peak_indices, gene_indices)),
    shape=(atac.n_vars, len(gene_names))
)

peak_gene_matrix.sum_duplicates()
peak_gene_matrix.data = np.ones_like(peak_gene_matrix.data)

print("\nPeak × gene matrix shape:")
print(peak_gene_matrix.shape)

# ============================================================
# Υπολογισμός gene activity: cells × peaks  @  peaks × genes
# ============================================================

if "counts" in atac.layers:
    X_atac = atac.layers["counts"]
else:
    X_atac = atac.X

if not sparse.issparse(X_atac):
    X_atac = sparse.csr_matrix(X_atac)
else:
    X_atac = X_atac.tocsr()

gene_activity_X = X_atac @ peak_gene_matrix

print("\nGene activity matrix shape:")
print(gene_activity_X.shape)

# ============================================================
# Δημιουργία AnnData για gene activity
# ============================================================

gene_activity = AnnData(
    X=gene_activity_X,
    obs=atac.obs.copy(),
    var=pd.DataFrame(index=gene_names)
)

# Κρατάμε UMAP και LSI coordinates από το ATAC clustering
gene_activity.obsm["X_umap"] = atac.obsm["X_umap"].copy()
gene_activity.obsm["X_lsi"] = atac.obsm["X_lsi"].copy()
gene_activity.obsm["X_lsi_2_50"] = atac.obsm["X_lsi_2_50"].copy()

# Κρατάμε raw gene activity counts
gene_activity.layers["counts"] = gene_activity.X.copy()

# ============================================================
# Normalize και log-transform gene activity
# ============================================================

sc.pp.normalize_total(gene_activity, target_sum=1e4)
sc.pp.log1p(gene_activity)

gene_activity.raw = gene_activity

print("\nFinal gene activity object:")
print(gene_activity)

print("\nGene activity obs columns:")
print(gene_activity.obs.columns.tolist())

print("\nExample genes:")
print(list(gene_activity.var_names[:20]))

# ============================================================
# Save output
# ============================================================

gene_activity.write_h5ad(out_path)

print(f"\nSaved ATAC gene activity object to:\n{out_path}")








"""Βήμα 33 — Marker gene activity dotplot για ATAC cluster annotation"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

gene_activity_path = project_dir / "data/processed/pbmc_atac_gene_activity.h5ad"

marker_table_out = (
    project_dir /
    "data/processed/pbmc_atac_gene_activity_marker_means_by_cluster.csv"
)

# ============================================================
# Load ATAC gene activity object
# ============================================================

ga = sc.read_h5ad(gene_activity_path)

print("ATAC gene activity object:")
print(ga)

# ============================================================
# Marker genes για PBMC annotation
# ============================================================

marker_genes = [
    # T cells
    "CD3D", "CD3E",

    # CD4 T cells
    "CD4", "IL7R",

    # CD8 T cells
    "CD8A", "CD8B",

    # Naive T cells
    "CCR7", "LEF1", "TCF7",

    # Cytotoxic / effector T / NK
    "NKG7", "GZMB", "PRF1", "GNLY", "CCL5",

    # Monocytes
    "CD14", "LYZ", "S100A8", "S100A9", "FCN1", "MS4A7",

    # B cells
    "MS4A1", "CD79A",

    # NK cells
    "KLRD1", "NCR1",

    # Dendritic / APC
    "FCER1A", "CST3"
]

available_markers = [g for g in marker_genes if g in ga.var_names]
missing_markers = [g for g in marker_genes if g not in ga.var_names]

print("\nAvailable markers:")
print(available_markers)

print("\nMissing markers:")
print(missing_markers)

# ============================================================
# Dotplot ανά ATAC Leiden cluster
# ============================================================

sc.pl.dotplot(
    ga,
    var_names=available_markers,
    groupby="leiden_0_6",
    standard_scale="var",
    dendrogram=False,
    title="ATAC gene activity markers by Leiden cluster"
)

# ============================================================
# UMAP με ATAC Leiden clusters
# ============================================================

sc.pl.umap(
    ga,
    color="leiden_0_6",
    legend_loc="on data",
    legend_fontsize=6,
    legend_fontweight="normal",
    frameon=False,
    title="ATAC Leiden clusters"
)

# ============================================================
# Feature plots για βασικά marker activity genes
# ============================================================

key_markers = [
    "CD14", "LYZ",
    "CD3D", "CD4",
    "CD8A", "CCR7",
    "NKG7", "GZMB",
    "MS4A1", "FCER1A"
]

key_markers = [g for g in key_markers if g in ga.var_names]

sc.pl.umap(
    ga,
    color=key_markers,
    frameon=False,
    vmax="p99",
    title=[f"{g} activity" for g in key_markers]
)

# ============================================================
# Mean marker activity table ανά cluster
# ============================================================

clusters = sorted(
    ga.obs["leiden_0_6"].astype(str).unique(),
    key=lambda x: int(x)
)

marker_means = pd.DataFrame(index=clusters, columns=available_markers)

for cluster in clusters:
    cluster_cells = ga.obs["leiden_0_6"].astype(str) == cluster
    X_cluster = ga[cluster_cells, available_markers].X

    if hasattr(X_cluster, "toarray"):
        mean_values = np.asarray(X_cluster.mean(axis=0)).ravel()
    else:
        mean_values = X_cluster.mean(axis=0)

    marker_means.loc[cluster, :] = mean_values

marker_means = marker_means.astype(float)

print("\nMean marker gene activity by ATAC Leiden cluster:")
print(marker_means)

marker_means.to_csv(marker_table_out)

print(f"\nSaved marker activity table to:\n{marker_table_out}")










"""Βήμα 34 — ATAC marker module scores για cell-type annotation"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

gene_activity_path = project_dir / "data/processed/pbmc_atac_gene_activity.h5ad"

module_score_out = (
    project_dir /
    "data/processed/pbmc_atac_gene_activity_module_scores_by_cluster.csv"
)

full_marker_table_out = (
    project_dir /
    "data/processed/pbmc_atac_gene_activity_full_marker_means_by_cluster.csv"
)

# ============================================================
# Load ATAC gene activity object
# ============================================================

ga = sc.read_h5ad(gene_activity_path)

print("Loaded ATAC gene activity object:")
print(ga)

# ============================================================
# Define marker modules
# ============================================================

marker_modules = {
    "T_cells": ["CD3D", "CD3E"],
    "CD4_T": ["CD4", "IL7R"],
    "CD4_Naive": ["CD4", "CCR7", "LEF1", "TCF7"],
    "CD4_Memory": ["CD4", "IL7R"],
    "CD8_T": ["CD8A", "CD8B"],
    "CD8_Naive": ["CD8A", "CD8B", "CCR7", "LEF1", "TCF7"],
    "CD8_Effector": ["CD8A", "CD8B", "NKG7", "GZMB", "PRF1", "GNLY", "CCL5"],
    "CD14_Monocytes": ["CD14", "LYZ", "S100A8", "S100A9", "FCN1", "MS4A7"],
    "B_cells": ["MS4A1", "CD79A"],
    "NK_cells": ["KLRD1", "NCR1", "NKG7", "GNLY"],
    "Dendritic_APC": ["FCER1A", "CST3"],
}

# ============================================================
# Keep only available markers
# ============================================================

marker_modules_available = {}

for module_name, genes in marker_modules.items():
    available_genes = [g for g in genes if g in ga.var_names]
    marker_modules_available[module_name] = available_genes

print("\nAvailable genes per marker module:")
for module_name, genes in marker_modules_available.items():
    print(module_name, ":", genes)

# ============================================================
# Compute mean module activity per cell
# ============================================================

for module_name, genes in marker_modules_available.items():
    if len(genes) == 0:
        continue

    X_module = ga[:, genes].X

    if hasattr(X_module, "toarray"):
        module_score = np.asarray(X_module.mean(axis=1)).ravel()
    else:
        module_score = np.asarray(X_module.mean(axis=1)).ravel()

    ga.obs[f"{module_name}_score"] = module_score

# ============================================================
# Summarize module scores by Leiden cluster
# ============================================================

score_cols = [f"{module_name}_score" for module_name in marker_modules_available.keys()]

module_scores_by_cluster = (
    ga.obs
    .groupby("leiden_0_6", observed=True)[score_cols]
    .mean()
)

# Remove "_score" suffix from columns
module_scores_by_cluster.columns = [
    c.replace("_score", "") for c in module_scores_by_cluster.columns
]

# Sort clusters numerically
module_scores_by_cluster = module_scores_by_cluster.loc[
    sorted(module_scores_by_cluster.index.astype(str), key=lambda x: int(x))
]

print("\nModule scores by ATAC Leiden cluster:")
print(module_scores_by_cluster.to_string())

module_scores_by_cluster.to_csv(module_score_out)

print(f"\nSaved module scores to:\n{module_score_out}")

# ============================================================
# Also print and save full marker means without truncation
# ============================================================

marker_genes = sorted(set(sum(marker_modules.values(), [])))
available_markers = [g for g in marker_genes if g in ga.var_names]

clusters = sorted(
    ga.obs["leiden_0_6"].astype(str).unique(),
    key=lambda x: int(x)
)

marker_means = pd.DataFrame(index=clusters, columns=available_markers)

for cluster in clusters:
    cluster_cells = ga.obs["leiden_0_6"].astype(str) == cluster
    X_cluster = ga[cluster_cells, available_markers].X

    if hasattr(X_cluster, "toarray"):
        mean_values = np.asarray(X_cluster.mean(axis=0)).ravel()
    else:
        mean_values = np.asarray(X_cluster.mean(axis=0)).ravel()

    marker_means.loc[cluster, :] = mean_values

marker_means = marker_means.astype(float)

print("\nFull marker means by ATAC Leiden cluster:")
print(marker_means.to_string())

marker_means.to_csv(full_marker_table_out)

print(f"\nSaved full marker means to:\n{full_marker_table_out}")




"""Βήμα 35 — Visualization των ATAC marker module scores"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

gene_activity_path = project_dir / "data/processed/pbmc_atac_gene_activity.h5ad"

# ============================================================
# Load ATAC gene activity object
# ============================================================

ga = sc.read_h5ad(gene_activity_path)

# ============================================================
# Define marker modules
# ============================================================

marker_modules = {
    "T_cells": ["CD3D", "CD3E"],
    "CD4_T": ["CD4", "IL7R"],
    "CD4_Naive": ["CD4", "CCR7", "LEF1", "TCF7"],
    "CD8_T": ["CD8A", "CD8B"],
    "CD8_Naive": ["CD8A", "CD8B", "CCR7", "LEF1", "TCF7"],
    "CD8_Effector": ["CD8A", "CD8B", "NKG7", "GZMB", "PRF1", "GNLY", "CCL5"],
    "CD14_Monocytes": ["CD14", "LYZ", "S100A8", "S100A9", "FCN1", "MS4A7"],
    "B_cells": ["MS4A1", "CD79A"],
    "NK_cells": ["KLRD1", "NCR1", "NKG7", "GNLY"],
    "Dendritic_APC": ["FCER1A", "CST3"],
}

# ============================================================
# Compute module scores
# ============================================================

for module_name, genes in marker_modules.items():
    genes = [g for g in genes if g in ga.var_names]

    if len(genes) == 0:
        continue

    X_module = ga[:, genes].X

    if hasattr(X_module, "toarray"):
        score = np.asarray(X_module.mean(axis=1)).ravel()
    else:
        score = np.asarray(X_module.mean(axis=1)).ravel()

    ga.obs[f"{module_name}_score"] = score

score_cols = [c for c in ga.obs.columns if c.endswith("_score")]

# ============================================================
# Dotplot των module scores ανά cluster
# ============================================================

sc.pl.dotplot(
    ga,
    var_names=score_cols,
    groupby="leiden_0_6",
    standard_scale="var",
    dendrogram=False,
    title="ATAC marker module scores by Leiden cluster"
)

# ============================================================
# UMAP με βασικά module scores
# ============================================================

plot_scores = [
    "CD14_Monocytes_score",
    "T_cells_score",
    "CD4_T_score",
    "CD8_T_score",
    "CD8_Effector_score",
    "B_cells_score",
    "NK_cells_score",
    "Dendritic_APC_score",
]

plot_scores = [s for s in plot_scores if s in ga.obs.columns]

sc.pl.umap(
    ga,
    color=plot_scores,
    frameon=False,
    vmax="p99",
    title=[s.replace("_score", "") for s in plot_scores]
)



"""Βήμα 36 — Προσθήκη final ATAC annotations σε peak object και gene activity object"""

import scanpy as sc
import pandas as pd
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

atac_path = project_dir / "data/processed/pbmc_atac_clustered.h5ad"
gene_activity_path = project_dir / "data/processed/pbmc_atac_gene_activity.h5ad"

atac_out = project_dir / "data/processed/pbmc_atac_annotated_final.h5ad"
gene_activity_out = project_dir / "data/processed/pbmc_atac_gene_activity_annotated_final.h5ad"

# ============================================================
# Load objects
# ============================================================

atac = sc.read_h5ad(atac_path)
ga = sc.read_h5ad(gene_activity_path)

print("Loaded ATAC peak object:")
print(atac)

print("\nLoaded ATAC gene activity object:")
print(ga)

# ============================================================
# Ensure Leiden labels are strings
# ============================================================

atac.obs["leiden_0_6"] = atac.obs["leiden_0_6"].astype(str)
ga.obs["leiden_0_6"] = ga.obs["leiden_0_6"].astype(str)

# ============================================================
# Define ATAC detailed annotation
# ============================================================

detailed_mapping = {
    "0":  "CD14_Monocytes_c0",
    "1":  "CD14_Monocytes_c1",
    "2":  "CD8_Effector_c2",
    "3":  "CD8_Naive_c3",
    "4":  "B_cells_c4",
    "5":  "CD4_Memory_c5",
    "6":  "CD4_Naive_c6",
    "7":  "NK_or_Cytotoxic_c7",
    "8":  "CD8_Naive_c8",
    "9":  "CD8_Effector_c9",
    "10": "Monocyte_related_c10",
    "11": "Dendritic_APC_c11",
    "12": "Mixed_B_Monocyte_APC_c12",
}

# ============================================================
# Define ATAC broad annotation
# ============================================================

broad_mapping = {
    "0":  "CD14_Monocytes",
    "1":  "CD14_Monocytes",
    "2":  "CD8_Effector",
    "3":  "CD8_Naive",
    "4":  "B_cells",
    "5":  "CD4_Memory",
    "6":  "CD4_Naive",
    "7":  "NK_or_Cytotoxic",
    "8":  "CD8_Naive",
    "9":  "CD8_Effector",
    "10": "Monocyte_related",
    "11": "Dendritic_APC",
    "12": "Mixed_or_Ambiguous",
}

# ============================================================
# Define scMultiomeGRN 5-type mapping
# ============================================================

scmultiomegrn_mapping = {
    "CD14_Monocytes": "CD14+ Monocytes",
    "CD4_Memory": "CD4.Memory",
    "CD4_Naive": "CD4.Naive",
    "CD8_Effector": "CD8.effector",
    "CD8_Naive": "CD8.Naive",
}

# ============================================================
# Apply annotations to both objects
# ============================================================

for obj_name, obj in {"ATAC": atac, "Gene activity": ga}.items():
    obj.obs["cell_type_detailed"] = obj.obs["leiden_0_6"].map(detailed_mapping)
    obj.obs["cell_type_broad"] = obj.obs["leiden_0_6"].map(broad_mapping)
    obj.obs["scmultiomegrn_cell_type"] = obj.obs["cell_type_broad"].map(scmultiomegrn_mapping)
    obj.obs["is_scmultiomegrn_5type"] = obj.obs["scmultiomegrn_cell_type"].notna()

    print(f"\n{obj_name} missing detailed annotations:")
    print(obj.obs["cell_type_detailed"].isna().sum())

    print(f"\n{obj_name} missing broad annotations:")
    print(obj.obs["cell_type_broad"].isna().sum())

    print(f"\n{obj_name} broad cell-type counts:")
    print(obj.obs["cell_type_broad"].value_counts())

    print(f"\n{obj_name} scMultiomeGRN 5-type counts:")
    print(obj.obs["scmultiomegrn_cell_type"].value_counts())

# ============================================================
# Save annotated objects
# ============================================================

atac.write_h5ad(atac_out)
ga.write_h5ad(gene_activity_out)

print(f"\nSaved annotated ATAC peak object to:\n{atac_out}")
print(f"Saved annotated ATAC gene activity object to:\n{gene_activity_out}")


"""Βήμα 37 — UMAP έλεγχος final ATAC annotations"""

import scanpy as sc
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

gene_activity_path = (
    project_dir /
    "data/processed/pbmc_atac_gene_activity_annotated_final.h5ad"
)

# ============================================================
# Load annotated ATAC gene activity object
# ============================================================

ga = sc.read_h5ad(gene_activity_path)

print("Loaded annotated ATAC gene activity object:")
print(ga)

# ============================================================
# UMAP με Leiden clusters
# ============================================================

sc.pl.umap(
    ga,
    color="leiden_0_6",
    legend_loc="on data",
    legend_fontsize=6,
    legend_fontweight="normal",
    frameon=False,
    title="ATAC Leiden clusters"
)

# ============================================================
# UMAP με broad cell types
# ============================================================

sc.pl.umap(
    ga,
    color="cell_type_broad",
    legend_loc="on data",
    legend_fontsize=6,
    legend_fontweight="normal",
    frameon=False,
    title="ATAC broad cell-type annotation"
)

# ============================================================
# UMAP με detailed cluster-aware labels
# ============================================================

sc.pl.umap(
    ga,
    color="cell_type_detailed",
    legend_loc="right margin",
    frameon=False,
    title="ATAC detailed cluster-aware annotation"
)

# ============================================================
# UMAP μόνο για scMultiomeGRN target labels
# ============================================================

sc.pl.umap(
    ga,
    color="scmultiomegrn_cell_type",
    legend_loc="right margin",
    frameon=False,
    title="ATAC scMultiomeGRN 5 target cell types"
)




"""Βήμα 38 — Δημιουργία ATAC subset για τα 5 scMultiomeGRN PBMC cell types"""

import scanpy as sc
import pandas as pd
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

atac_path = project_dir / "data/processed/pbmc_atac_annotated_final.h5ad"
ga_path = project_dir / "data/processed/pbmc_atac_gene_activity_annotated_final.h5ad"

atac_out = project_dir / "data/processed/pbmc_atac_scmultiomegrn_5types.h5ad"
ga_out = project_dir / "data/processed/pbmc_atac_gene_activity_scmultiomegrn_5types.h5ad"

counts_out = project_dir / "data/processed/pbmc_atac_scmultiomegrn_5types_cell_counts.csv"

# ============================================================
# Load annotated objects
# ============================================================

atac = sc.read_h5ad(atac_path)
ga = sc.read_h5ad(ga_path)

# ============================================================
# Subset to scMultiomeGRN target cell types
# ============================================================

atac_5 = atac[atac.obs["is_scmultiomegrn_5type"]].copy()
ga_5 = ga[ga.obs["is_scmultiomegrn_5type"]].copy()

# ============================================================
# Print summary
# ============================================================

print("ATAC 5-type peak object:")
print(atac_5)

print("\nATAC 5-type gene activity object:")
print(ga_5)

print("\nCell counts by scMultiomeGRN cell type:")
counts = atac_5.obs["scmultiomegrn_cell_type"].value_counts()
print(counts)

print("\nCell counts by Leiden cluster:")
print(atac_5.obs["leiden_0_6"].value_counts().sort_index())

# ============================================================
# Save outputs
# ============================================================

counts.to_csv(counts_out, header=["n_cells"])

atac_5.write_h5ad(atac_out)
ga_5.write_h5ad(ga_out)

print(f"\nSaved ATAC 5-type peak object to:\n{atac_out}")
print(f"Saved ATAC 5-type gene activity object to:\n{ga_out}")
print(f"Saved ATAC 5-type counts to:\n{counts_out}")







"""Βήμα 39 — RNA–ATAC matching summary για τα 5 scMultiomeGRN PBMC cell types"""

import scanpy as sc
import pandas as pd
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

rna_5_path = project_dir / "data/processed/pbmc_rna_scmultiomegrn_5types.h5ad"
atac_5_path = project_dir / "data/processed/pbmc_atac_scmultiomegrn_5types.h5ad"
ga_5_path = project_dir / "data/processed/pbmc_atac_gene_activity_scmultiomegrn_5types.h5ad"

summary_out = project_dir / "data/processed/pbmc_rna_atac_scmultiomegrn_5types_matching_summary.csv"

# ============================================================
# Load RNA and ATAC 5-type objects
# ============================================================

rna_5 = sc.read_h5ad(rna_5_path)
atac_5 = sc.read_h5ad(atac_5_path)
ga_5 = sc.read_h5ad(ga_5_path)

print("RNA 5-type object:")
print(rna_5)

print("\nATAC 5-type peak object:")
print(atac_5)

print("\nATAC 5-type gene activity object:")
print(ga_5)

# ============================================================
# Cell-type counts
# ============================================================

rna_counts = rna_5.obs["scmultiomegrn_cell_type"].value_counts()
atac_counts = atac_5.obs["scmultiomegrn_cell_type"].value_counts()
ga_counts = ga_5.obs["scmultiomegrn_cell_type"].value_counts()

all_cell_types = sorted(
    set(rna_counts.index)
    .union(set(atac_counts.index))
    .union(set(ga_counts.index))
)

summary = pd.DataFrame(index=all_cell_types)

summary["rna_cells"] = rna_counts.reindex(all_cell_types).fillna(0).astype(int)
summary["atac_peak_cells"] = atac_counts.reindex(all_cell_types).fillna(0).astype(int)
summary["atac_gene_activity_cells"] = ga_counts.reindex(all_cell_types).fillna(0).astype(int)

summary["rna_fraction"] = summary["rna_cells"] / summary["rna_cells"].sum()
summary["atac_fraction"] = summary["atac_peak_cells"] / summary["atac_peak_cells"].sum()

# ============================================================
# Add total row
# ============================================================

summary_with_total = summary.copy()

summary_with_total.loc["TOTAL", "rna_cells"] = summary["rna_cells"].sum()
summary_with_total.loc["TOTAL", "atac_peak_cells"] = summary["atac_peak_cells"].sum()
summary_with_total.loc["TOTAL", "atac_gene_activity_cells"] = summary["atac_gene_activity_cells"].sum()
summary_with_total.loc["TOTAL", "rna_fraction"] = 1.0
summary_with_total.loc["TOTAL", "atac_fraction"] = 1.0

# ============================================================
# Print and save summary
# ============================================================

print("\nRNA–ATAC 5-type matching summary:")
print(summary_with_total)

summary_with_total.to_csv(summary_out)

print(f"\nSaved RNA–ATAC matching summary to:\n{summary_out}")





"""Βήμα 40 — Έλεγχος RNA raw gene space και RNA–ATAC gene overlap"""

import scanpy as sc
import pandas as pd
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

rna_5_path = project_dir / "data/processed/pbmc_rna_scmultiomegrn_5types.h5ad"
ga_5_path = project_dir / "data/processed/pbmc_atac_gene_activity_scmultiomegrn_5types.h5ad"

gene_overlap_out = project_dir / "data/processed/pbmc_rna_atac_common_genes.csv"

# ============================================================
# Load objects
# ============================================================

rna_5 = sc.read_h5ad(rna_5_path)
ga_5 = sc.read_h5ad(ga_5_path)

print("RNA 5-type object:")
print(rna_5)

print("\nATAC gene activity 5-type object:")
print(ga_5)

# ============================================================
# Check RNA X, raw, and layers
# ============================================================

print("\nRNA .X gene count:")
print(rna_5.n_vars)

print("\nRNA var_names example:")
print(list(rna_5.var_names[:10]))

print("\nRNA layers:")
print(list(rna_5.layers.keys()))

print("\nRNA raw exists:")
print(rna_5.raw is not None)

if rna_5.raw is not None:
    print("\nRNA raw shape:")
    print(rna_5.raw.shape)

    print("\nRNA raw var_names example:")
    print(list(rna_5.raw.var_names[:10]))

# ============================================================
# Define RNA gene universe
# ============================================================

if rna_5.raw is not None:
    rna_genes = set(rna_5.raw.var_names)
    rna_gene_source = "rna_5.raw.var_names"
else:
    rna_genes = set(rna_5.var_names)
    rna_gene_source = "rna_5.var_names"

atac_ga_genes = set(ga_5.var_names)

common_genes = sorted(rna_genes.intersection(atac_ga_genes))

print("\nRNA gene source used:")
print(rna_gene_source)

print("\nRNA genes:")
print(len(rna_genes))

print("\nATAC gene activity genes:")
print(len(atac_ga_genes))

print("\nCommon RNA–ATAC genes:")
print(len(common_genes))

print("\nExample common genes:")
print(common_genes[:30])

# ============================================================
# Check important marker genes in common space
# ============================================================

important_genes = [
    "CD14", "LYZ", "S100A8", "S100A9", "FCN1", "MS4A7",
    "CD3D", "CD3E", "CD4", "IL7R", "CCR7", "LEF1", "TCF7",
    "CD8A", "CD8B", "NKG7", "GZMB", "PRF1", "GNLY", "CCL5"
]

important_status = pd.DataFrame({
    "gene": important_genes,
    "in_rna": [g in rna_genes for g in important_genes],
    "in_atac_gene_activity": [g in atac_ga_genes for g in important_genes],
    "in_common": [g in common_genes for g in important_genes],
})

print("\nImportant marker gene overlap status:")
print(important_status)

# ============================================================
# Save common genes
# ============================================================

common_genes_df = pd.DataFrame({"gene": common_genes})
common_genes_df.to_csv(gene_overlap_out, index=False)

print(f"\nSaved common RNA–ATAC genes to:\n{gene_overlap_out}")


"""Βήμα 41 — Εξαγωγή RNA και ATAC metadata για τα 5 scMultiomeGRN cell types"""

import scanpy as sc
import pandas as pd
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

rna_5_path = project_dir / "data/processed/pbmc_rna_scmultiomegrn_5types.h5ad"
atac_5_path = project_dir / "data/processed/pbmc_atac_scmultiomegrn_5types.h5ad"
ga_5_path = project_dir / "data/processed/pbmc_atac_gene_activity_scmultiomegrn_5types.h5ad"

rna_meta_out = project_dir / "data/processed/pbmc_rna_scmultiomegrn_5types_metadata.csv"
atac_meta_out = project_dir / "data/processed/pbmc_atac_scmultiomegrn_5types_metadata.csv"
ga_meta_out = project_dir / "data/processed/pbmc_atac_gene_activity_scmultiomegrn_5types_metadata.csv"

# ============================================================
# Load objects
# ============================================================

rna_5 = sc.read_h5ad(rna_5_path)
atac_5 = sc.read_h5ad(atac_5_path)
ga_5 = sc.read_h5ad(ga_5_path)

# ============================================================
# Select metadata columns
# ============================================================

rna_cols = [
    "leiden_0_6",
    "cell_type_detailed",
    "cell_type_broad",
    "scmultiomegrn_cell_type",
]

atac_cols = [
    "leiden_0_6",
    "cell_type_detailed",
    "cell_type_broad",
    "scmultiomegrn_cell_type",
    "total_counts",
    "n_peaks_by_counts",
]

# ============================================================
# Export metadata
# ============================================================

rna_meta = rna_5.obs[rna_cols].copy()
rna_meta.index.name = "cell_barcode"

atac_meta = atac_5.obs[atac_cols].copy()
atac_meta.index.name = "cell_barcode"

ga_meta = ga_5.obs[atac_cols].copy()
ga_meta.index.name = "cell_barcode"

rna_meta.to_csv(rna_meta_out)
atac_meta.to_csv(atac_meta_out)
ga_meta.to_csv(ga_meta_out)

print("Saved RNA metadata to:")
print(rna_meta_out)

print("\nSaved ATAC peak metadata to:")
print(atac_meta_out)

print("\nSaved ATAC gene activity metadata to:")
print(ga_meta_out)

print("\nRNA metadata preview:")
print(rna_meta.head())

print("\nATAC metadata preview:")
print(atac_meta.head())






"""Βήμα 42 — Δημιουργία pseudo-bulk RNA expression και ATAC gene activity ανά cell type"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

rna_5_path = project_dir / "data/processed/pbmc_rna_scmultiomegrn_5types.h5ad"
ga_5_path = project_dir / "data/processed/pbmc_atac_gene_activity_scmultiomegrn_5types.h5ad"
common_genes_path = project_dir / "data/processed/pbmc_rna_atac_common_genes.csv"

rna_pseudobulk_out = project_dir / "data/processed/pbmc_rna_pseudobulk_5types_common_genes.csv"
atac_ga_pseudobulk_out = project_dir / "data/processed/pbmc_atac_gene_activity_pseudobulk_5types_common_genes.csv"
combined_summary_out = project_dir / "data/processed/pbmc_rna_atac_pseudobulk_5types_summary.csv"

# ============================================================
# Load objects
# ============================================================

rna_5 = sc.read_h5ad(rna_5_path)
ga_5 = sc.read_h5ad(ga_5_path)

common_genes = pd.read_csv(common_genes_path)["gene"].tolist()

print("Loaded RNA 5-type object:")
print(rna_5)

print("\nLoaded ATAC gene activity 5-type object:")
print(ga_5)

print("\nNumber of common genes:")
print(len(common_genes))

# ============================================================
# Use RNA raw gene space, not scaled HVG-only .X
# ============================================================

if rna_5.raw is None:
    raise ValueError("rna_5.raw is missing. Δεν μπορούμε να χρησιμοποιήσουμε full RNA gene space.")

rna_full = rna_5.raw.to_adata()
rna_full.obs = rna_5.obs.copy()

rna_common_genes = [g for g in common_genes if g in rna_full.var_names]
ga_common_genes = [g for g in common_genes if g in ga_5.var_names]

final_common_genes = sorted(set(rna_common_genes).intersection(ga_common_genes))

print("\nRNA common genes available:")
print(len(rna_common_genes))

print("\nATAC gene activity common genes available:")
print(len(ga_common_genes))

print("\nFinal common genes used:")
print(len(final_common_genes))

rna_full = rna_full[:, final_common_genes].copy()
ga_common = ga_5[:, final_common_genes].copy()

# ============================================================
# Define cell-type order
# ============================================================

cell_type_order = [
    "CD14+ Monocytes",
    "CD4.Memory",
    "CD4.Naive",
    "CD8.effector",
    "CD8.Naive",
]

# ============================================================
# Helper function: mean pseudo-bulk by cell type
# ============================================================

def compute_mean_pseudobulk(adata, group_col, groups, genes):
    pseudobulk = pd.DataFrame(index=groups, columns=genes, dtype=float)

    for group in groups:
        group_mask = adata.obs[group_col].astype(str) == group
        n_cells = int(group_mask.sum())

        print(f"Computing pseudo-bulk for {group}: {n_cells} cells")

        if n_cells == 0:
            raise ValueError(f"No cells found for group: {group}")

        X_group = adata[group_mask, genes].X

        if sparse.issparse(X_group):
            mean_values = np.asarray(X_group.mean(axis=0)).ravel()
        else:
            mean_values = np.asarray(X_group.mean(axis=0)).ravel()

        pseudobulk.loc[group, :] = mean_values

    return pseudobulk

# ============================================================
# Compute RNA pseudo-bulk expression
# ============================================================

rna_pseudobulk = compute_mean_pseudobulk(
    adata=rna_full,
    group_col="scmultiomegrn_cell_type",
    groups=cell_type_order,
    genes=final_common_genes,
)

print("\nRNA pseudo-bulk shape:")
print(rna_pseudobulk.shape)

# ============================================================
# Compute ATAC gene activity pseudo-bulk
# ============================================================

atac_ga_pseudobulk = compute_mean_pseudobulk(
    adata=ga_common,
    group_col="scmultiomegrn_cell_type",
    groups=cell_type_order,
    genes=final_common_genes,
)

print("\nATAC gene activity pseudo-bulk shape:")
print(atac_ga_pseudobulk.shape)

# ============================================================
# Save pseudo-bulk matrices
# ============================================================

rna_pseudobulk.to_csv(rna_pseudobulk_out)
atac_ga_pseudobulk.to_csv(atac_ga_pseudobulk_out)

print("\nSaved RNA pseudo-bulk to:")
print(rna_pseudobulk_out)

print("\nSaved ATAC gene activity pseudo-bulk to:")
print(atac_ga_pseudobulk_out)

# ============================================================
# Create summary table
# ============================================================

summary = pd.DataFrame({
    "cell_type": cell_type_order,
    "rna_cells": [
        int((rna_5.obs["scmultiomegrn_cell_type"].astype(str) == ct).sum())
        for ct in cell_type_order
    ],
    "atac_gene_activity_cells": [
        int((ga_5.obs["scmultiomegrn_cell_type"].astype(str) == ct).sum())
        for ct in cell_type_order
    ],
})

summary["n_common_genes"] = len(final_common_genes)

summary.to_csv(combined_summary_out, index=False)

print("\nPseudo-bulk summary:")
print(summary)

print("\nSaved pseudo-bulk summary to:")
print(combined_summary_out)









"""Βήμα 43 — Sanity check των pseudo-bulk RNA και ATAC gene activity marker profiles"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

rna_pseudobulk_path = project_dir / "data/processed/pbmc_rna_pseudobulk_5types_common_genes.csv"
atac_pseudobulk_path = project_dir / "data/processed/pbmc_atac_gene_activity_pseudobulk_5types_common_genes.csv"

marker_summary_out = project_dir / "data/processed/pbmc_pseudobulk_marker_sanity_check.csv"
rna_marker_zscore_out = project_dir / "data/processed/pbmc_rna_pseudobulk_marker_zscores.csv"
atac_marker_zscore_out = project_dir / "data/processed/pbmc_atac_gene_activity_pseudobulk_marker_zscores.csv"

# ============================================================
# Load pseudo-bulk matrices
# ============================================================

rna_pb = pd.read_csv(rna_pseudobulk_path, index_col=0)
atac_pb = pd.read_csv(atac_pseudobulk_path, index_col=0)

print("RNA pseudo-bulk shape:")
print(rna_pb.shape)

print("\nATAC gene activity pseudo-bulk shape:")
print(atac_pb.shape)

print("\nRNA cell types:")
print(rna_pb.index.tolist())

print("\nATAC cell types:")
print(atac_pb.index.tolist())

# ============================================================
# Define marker genes
# ============================================================

marker_sets = {
    "CD14_Monocytes": ["CD14", "LYZ", "S100A8", "S100A9", "FCN1", "MS4A7"],
    "CD4_Naive": ["CD4", "CCR7", "LEF1", "TCF7"],
    "CD4_Memory": ["CD4", "IL7R"],
    "CD8_Naive": ["CD8A", "CD8B", "CCR7", "LEF1", "TCF7"],
    "CD8_Effector": ["CD8A", "CD8B", "NKG7", "GZMB", "PRF1", "GNLY", "CCL5"],
}

all_markers = sorted(set(sum(marker_sets.values(), [])))

available_markers = [
    gene for gene in all_markers
    if gene in rna_pb.columns and gene in atac_pb.columns
]

missing_markers = [
    gene for gene in all_markers
    if gene not in rna_pb.columns or gene not in atac_pb.columns
]

print("\nAvailable markers:")
print(available_markers)

print("\nMissing markers:")
print(missing_markers)

# ============================================================
# Extract marker pseudo-bulk values
# ============================================================

rna_markers = rna_pb[available_markers].copy()
atac_markers = atac_pb[available_markers].copy()

print("\nRNA marker pseudo-bulk values:")
print(rna_markers.to_string())

print("\nATAC gene activity marker pseudo-bulk values:")
print(atac_markers.to_string())

# ============================================================
# Z-score ανά gene μέσα σε κάθε modality
# ============================================================

def zscore_by_gene(df):
    gene_means = df.mean(axis=0)
    gene_stds = df.std(axis=0, ddof=0)
    gene_stds = gene_stds.replace(0, np.nan)
    return (df - gene_means) / gene_stds

rna_marker_z = zscore_by_gene(rna_markers)
atac_marker_z = zscore_by_gene(atac_markers)

print("\nRNA marker z-scores:")
print(rna_marker_z.to_string())

print("\nATAC gene activity marker z-scores:")
print(atac_marker_z.to_string())

# ============================================================
# Save marker tables
# ============================================================

rna_marker_z.to_csv(rna_marker_zscore_out)
atac_marker_z.to_csv(atac_marker_zscore_out)

combined_marker_summary = []

for gene in available_markers:
    for cell_type in rna_pb.index:
        combined_marker_summary.append({
            "gene": gene,
            "cell_type": cell_type,
            "rna_value": rna_markers.loc[cell_type, gene],
            "atac_gene_activity_value": atac_markers.loc[cell_type, gene],
            "rna_zscore": rna_marker_z.loc[cell_type, gene],
            "atac_gene_activity_zscore": atac_marker_z.loc[cell_type, gene],
        })

combined_marker_summary = pd.DataFrame(combined_marker_summary)
combined_marker_summary.to_csv(marker_summary_out, index=False)

print("\nSaved combined marker sanity check table to:")
print(marker_summary_out)

# ============================================================
# Plot RNA marker z-score heatmap
# ============================================================

plt.figure(figsize=(12, 4))
plt.imshow(rna_marker_z, aspect="auto")
plt.xticks(range(len(available_markers)), available_markers, rotation=90)
plt.yticks(range(len(rna_marker_z.index)), rna_marker_z.index)
plt.colorbar(label="RNA marker z-score")
plt.title("RNA pseudo-bulk marker z-scores")
plt.tight_layout()
plt.show()

# ============================================================
# Plot ATAC gene activity marker z-score heatmap
# ============================================================

plt.figure(figsize=(12, 4))
plt.imshow(atac_marker_z, aspect="auto")
plt.xticks(range(len(available_markers)), available_markers, rotation=90)
plt.yticks(range(len(atac_marker_z.index)), atac_marker_z.index)
plt.colorbar(label="ATAC gene activity marker z-score")
plt.title("ATAC gene activity pseudo-bulk marker z-scores")
plt.tight_layout()
plt.show()



"""Βήμα 44B — Καθάρισμα και diagnosis των RNA–ATAC pseudo-bulk gene correlations"""

import pandas as pd
import numpy as np
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

correlation_path = project_dir / "data/processed/pbmc_rna_atac_pseudobulk_gene_correlations.csv"

clean_correlation_out = project_dir / "data/processed/pbmc_rna_atac_pseudobulk_gene_correlations_finite.csv"
nan_correlation_out = project_dir / "data/processed/pbmc_rna_atac_pseudobulk_gene_correlations_nan_genes.csv"
positive_out = project_dir / "data/processed/pbmc_rna_atac_pseudobulk_gene_correlations_top_positive.csv"
negative_out = project_dir / "data/processed/pbmc_rna_atac_pseudobulk_gene_correlations_top_negative.csv"
concordant_out = project_dir / "data/processed/pbmc_rna_atac_pseudobulk_gene_correlations_concordant_genes.csv"

# ============================================================
# Load correlation table
# ============================================================

corr = pd.read_csv(correlation_path)

print("Original correlation table shape:")
print(corr.shape)

# ============================================================
# Split finite and non-finite correlations
# ============================================================

finite_corr = corr[
    np.isfinite(corr["rna_atac_pseudobulk_correlation"])
].copy()

nan_corr = corr[
    ~np.isfinite(corr["rna_atac_pseudobulk_correlation"])
].copy()

print("\nFinite correlations:")
print(finite_corr.shape[0])

print("\nNaN correlations:")
print(nan_corr.shape[0])

print("\nNaN genes with n_cell_types counts:")
print(nan_corr["n_cell_types"].value_counts(dropna=False).sort_index())

# ============================================================
# Sort finite correlations
# ============================================================

finite_corr = finite_corr.sort_values(
    "rna_atac_pseudobulk_correlation",
    ascending=False
)

top_positive = finite_corr.head(100).copy()

top_negative = finite_corr.sort_values(
    "rna_atac_pseudobulk_correlation",
    ascending=True
).head(100).copy()

# ============================================================
# Define concordant genes
# ============================================================

concordant = finite_corr[
    finite_corr["rna_atac_pseudobulk_correlation"] >= 0.40
].copy()

strong_concordant = finite_corr[
    finite_corr["rna_atac_pseudobulk_correlation"] >= 0.70
].copy()

print("\nGenes with correlation >= 0.40:")
print(concordant.shape[0])

print("\nGenes with correlation >= 0.70:")
print(strong_concordant.shape[0])

# ============================================================
# Save cleaned outputs
# ============================================================

finite_corr.to_csv(clean_correlation_out, index=False)
nan_corr.to_csv(nan_correlation_out, index=False)
top_positive.to_csv(positive_out, index=False)
top_negative.to_csv(negative_out, index=False)
concordant.to_csv(concordant_out, index=False)

print("\nSaved finite correlations to:")
print(clean_correlation_out)

print("\nSaved NaN genes to:")
print(nan_correlation_out)

print("\nSaved top positive genes to:")
print(positive_out)

print("\nSaved top negative genes to:")
print(negative_out)

print("\nSaved concordant genes to:")
print(concordant_out)

# ============================================================
# Inspect actual top positive and actual top negative genes
# ============================================================

print("\nTop 30 finite positively correlated genes:")
print(
    finite_corr.head(30).to_string(index=False)
)

print("\nTop 30 finite negatively correlated genes:")
print(
    finite_corr.sort_values(
        "rna_atac_pseudobulk_correlation",
        ascending=True
    ).head(30).to_string(index=False)
)

# ============================================================
# Marker gene inspection
# ============================================================

important_markers = [
    "CD14", "LYZ", "S100A8", "S100A9", "FCN1", "MS4A7",
    "CD4", "IL7R", "CCR7", "LEF1", "TCF7",
    "CD8A", "CD8B", "NKG7", "GZMB", "PRF1", "GNLY", "CCL5"
]

marker_corr = finite_corr[
    finite_corr["gene"].isin(important_markers)
].copy()

marker_corr = marker_corr.sort_values(
    "rna_atac_pseudobulk_correlation",
    ascending=False
)

print("\nMarker gene finite RNA–ATAC correlations:")
print(marker_corr.to_string(index=False))


"""Βήμα 44C — Έλεγχος raw pseudo-bulk profiles για markers με χαμηλότερη RNA–ATAC concordance"""

import pandas as pd
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

rna_pseudobulk_path = project_dir / "data/processed/pbmc_rna_pseudobulk_5types_common_genes.csv"
atac_pseudobulk_path = project_dir / "data/processed/pbmc_atac_gene_activity_pseudobulk_5types_common_genes.csv"

marker_profile_out = project_dir / "data/processed/pbmc_rna_atac_selected_marker_raw_profiles.csv"

# ============================================================
# Load pseudo-bulk matrices
# ============================================================

rna_pb = pd.read_csv(rna_pseudobulk_path, index_col=0)
atac_pb = pd.read_csv(atac_pseudobulk_path, index_col=0)

# ============================================================
# Select markers to inspect
# ============================================================

markers_to_check = [
    "CD8A", "CD8B", "GNLY",
    "NKG7", "GZMB", "PRF1", "CCL5",
    "CD14", "FCN1", "LYZ",
    "LEF1", "TCF7", "CCR7", "IL7R"
]

markers_to_check = [
    g for g in markers_to_check
    if g in rna_pb.columns and g in atac_pb.columns
]

# ============================================================
# Build long-format comparison table
# ============================================================

records = []

for gene in markers_to_check:
    for cell_type in rna_pb.index:
        if cell_type not in atac_pb.index:
            continue

        records.append({
            "gene": gene,
            "cell_type": cell_type,
            "rna_pseudobulk": rna_pb.loc[cell_type, gene],
            "atac_gene_activity_pseudobulk": atac_pb.loc[cell_type, gene],
        })

marker_profiles = pd.DataFrame(records)

marker_profiles.to_csv(marker_profile_out, index=False)

print("Saved selected marker raw profiles to:")
print(marker_profile_out)

print("\nSelected marker raw profiles:")
print(marker_profiles.to_string(index=False))


"""Βήμα 45 — Peak-to-gene annotation table για τα 5 ATAC cell types"""

import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path
from scipy import sparse

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

atac_5types_path = project_dir / "data/processed/pbmc_atac_scmultiomegrn_5types.h5ad"

peak_annotation_path = project_dir / "data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_peak_annotation.tsv"

gene_corr_path = project_dir / "data/processed/pbmc_rna_atac_pseudobulk_gene_correlations_finite.csv"
concordant_gene_path = project_dir / "data/processed/pbmc_rna_atac_pseudobulk_gene_correlations_concordant_genes.csv"

peak_accessibility_out = project_dir / "data/processed/pbmc_atac_peak_accessibility_pseudobulk_5types.csv"
peak_detection_out = project_dir / "data/processed/pbmc_atac_peak_detection_fraction_5types.csv"
peak_to_gene_out = project_dir / "data/processed/pbmc_atac_peak_to_gene_annotation_5types.csv"
peak_to_gene_summary_out = project_dir / "data/processed/pbmc_atac_peak_to_gene_annotation_5types_summary.csv"

# ============================================================
# Load ATAC 5-type object
# ============================================================

atac = sc.read_h5ad(atac_5types_path)

print("ATAC 5-type object:")
print(atac)

print("\nATAC obs columns:")
print(list(atac.obs.columns))

print("\nATAC var names example:")
print(atac.var_names[:5].tolist())

# ============================================================
# Define expected cell types
# ============================================================

cell_type_order = [
    "CD14+ Monocytes",
    "CD4.Memory",
    "CD4.Naive",
    "CD8.effector",
    "CD8.Naive",
]

# ============================================================
# Detect cell type column
# ============================================================

candidate_cell_type_columns = [
    "scmultiomegrn_cell_type",
    "cell_type_scmultiomegrn",
    "cell_type_5types",
    "cell_type",
    "cell_type_broad",
    "cell_type_detailed",
]

cell_type_col = None

for col in candidate_cell_type_columns:
    if col in atac.obs.columns:
        values = set(atac.obs[col].astype(str).unique())
        if set(cell_type_order).issubset(values):
            cell_type_col = col
            break

if cell_type_col is None:
    print("\nCould not automatically detect the 5-type cell type column.")
    print("Available obs columns and example values:")
    for col in atac.obs.columns:
        vals = atac.obs[col].astype(str).value_counts().head(10)
        print(f"\n{col}:")
        print(vals)
    raise ValueError("Please check which obs column contains the 5 scMultiomeGRN cell type labels.")

print("\nUsing cell type column:")
print(cell_type_col)

print("\nCell type counts:")
print(atac.obs[cell_type_col].value_counts())

# ============================================================
# Select matrix for peak accessibility
# ============================================================

if "counts" in atac.layers:
    X = atac.layers["counts"]
    matrix_source = "atac.layers['counts']"
else:
    X = atac.X
    matrix_source = "atac.X"

if sparse.issparse(X):
    X = X.tocsr()

print("\nUsing matrix source for peak accessibility:")
print(matrix_source)

print("\nMatrix shape:")
print(X.shape)

# ============================================================
# Compute ATAC peak pseudo-bulk mean accessibility and detection fraction
# ============================================================

peak_ids = atac.var_names.astype(str).tolist()

mean_accessibility = pd.DataFrame(index=peak_ids)
detection_fraction = pd.DataFrame(index=peak_ids)

cell_counts = {}

for cell_type in cell_type_order:
    mask = atac.obs[cell_type_col].astype(str).values == cell_type
    n_cells = int(mask.sum())
    cell_counts[cell_type] = n_cells

    if n_cells == 0:
        raise ValueError(f"No cells found for cell type: {cell_type}")

    X_ct = X[mask, :]

    mean_values = np.asarray(X_ct.mean(axis=0)).ravel()

    if sparse.issparse(X_ct):
        detected_values = np.asarray((X_ct > 0).mean(axis=0)).ravel()
    else:
        detected_values = (X_ct > 0).mean(axis=0)

    safe_name = (
        cell_type
        .replace("+", "plus")
        .replace(" ", "_")
        .replace(".", "_")
        .replace("-", "_")
    )

    mean_accessibility[f"mean_accessibility_{safe_name}"] = mean_values
    detection_fraction[f"detection_fraction_{safe_name}"] = detected_values

mean_accessibility.index.name = "peak_id"
detection_fraction.index.name = "peak_id"

mean_accessibility.to_csv(peak_accessibility_out)
detection_fraction.to_csv(peak_detection_out)

print("\nSaved peak mean accessibility pseudo-bulk to:")
print(peak_accessibility_out)

print("\nSaved peak detection fraction pseudo-bulk to:")
print(peak_detection_out)

# ============================================================
# Load 10x peak annotation table
# ============================================================

peak_annot_raw = pd.read_csv(peak_annotation_path, sep="\t")

peak_annot_raw.columns = [
    str(c).strip() for c in peak_annot_raw.columns
]

print("\nRaw peak annotation shape:")
print(peak_annot_raw.shape)

print("\nRaw peak annotation columns:")
print(list(peak_annot_raw.columns))

print("\nRaw peak annotation head:")
print(peak_annot_raw.head())

# ============================================================
# Helper for flexible column detection
# ============================================================

def find_column(columns, candidates):
    lower_to_original = {c.lower(): c for c in columns}

    for candidate in candidates:
        if candidate.lower() in lower_to_original:
            return lower_to_original[candidate.lower()]

    return None

columns = list(peak_annot_raw.columns)

peak_col = find_column(
    columns,
    ["peak", "peak_id", "interval", "region"]
)

gene_col = find_column(
    columns,
    ["gene", "gene_name", "target_gene", "symbol"]
)

distance_col = find_column(
    columns,
    ["distance", "dist", "distance_to_tss", "tss_distance"]
)

peak_type_col = find_column(
    columns,
    ["peak_type", "annotation", "feature_type", "type"]
)

chrom_col = find_column(
    columns,
    ["chrom", "chr", "chromosome", "#chrom"]
)

start_col = find_column(
    columns,
    ["start", "peak_start", "chromStart"]
)

end_col = find_column(
    columns,
    ["end", "peak_end", "chromEnd"]
)

if gene_col is None:
    raise ValueError("Could not find gene column in peak annotation table.")

# ============================================================
# Build peak_id column
# ============================================================

peak_annot = peak_annot_raw.copy()

if peak_col is not None:
    peak_annot["peak_id"] = peak_annot[peak_col].astype(str).str.strip()
else:
    if chrom_col is None or start_col is None or end_col is None:
        raise ValueError(
            "Could not find a peak column or chrom/start/end columns in peak annotation table."
        )

    peak_annot["peak_id"] = (
        peak_annot[chrom_col].astype(str).str.strip()
        + ":"
        + peak_annot[start_col].astype(int).astype(str)
        + "-"
        + peak_annot[end_col].astype(int).astype(str)
    )

peak_annot["target_gene"] = peak_annot[gene_col].astype(str).str.strip()

# ============================================================
# Clean target genes
# ============================================================

peak_annot = peak_annot[
    peak_annot["target_gene"].notna()
].copy()

peak_annot = peak_annot[
    ~peak_annot["target_gene"].isin(["", "nan", "None", "NA"])
].copy()

# Some annotation tables may contain multiple genes per peak
peak_annot["target_gene"] = peak_annot["target_gene"].str.split(";")
peak_annot = peak_annot.explode("target_gene")
peak_annot["target_gene"] = peak_annot["target_gene"].astype(str).str.strip()

peak_annot = peak_annot[
    ~peak_annot["target_gene"].isin(["", "nan", "None", "NA"])
].copy()

# ============================================================
# Keep annotation metadata columns
# ============================================================

keep_cols = ["peak_id", "target_gene"]

if distance_col is not None:
    keep_cols.append(distance_col)

if peak_type_col is not None:
    keep_cols.append(peak_type_col)

for optional_col in [chrom_col, start_col, end_col]:
    if optional_col is not None and optional_col not in keep_cols:
        keep_cols.append(optional_col)

peak_to_gene = peak_annot[keep_cols].copy()

rename_dict = {}

if distance_col is not None:
    rename_dict[distance_col] = "distance_to_gene"

if peak_type_col is not None:
    rename_dict[peak_type_col] = "peak_annotation_type"

if chrom_col is not None:
    rename_dict[chrom_col] = "chrom"

if start_col is not None:
    rename_dict[start_col] = "start"

if end_col is not None:
    rename_dict[end_col] = "end"

peak_to_gene = peak_to_gene.rename(columns=rename_dict)

# ============================================================
# Filter to peaks present in the 5-type ATAC object
# ============================================================

atac_peak_set = set(peak_ids)

peak_to_gene = peak_to_gene[
    peak_to_gene["peak_id"].isin(atac_peak_set)
].copy()

peak_to_gene = peak_to_gene.drop_duplicates()

print("\nPeak-to-gene rows after filtering to ATAC 5-type peaks:")
print(peak_to_gene.shape[0])

print("\nUnique annotated peaks:")
print(peak_to_gene["peak_id"].nunique())

print("\nUnique target genes:")
print(peak_to_gene["target_gene"].nunique())

# ============================================================
# Add ATAC peak accessibility pseudo-bulk values
# ============================================================

peak_features = mean_accessibility.join(detection_fraction, how="left")
peak_features = peak_features.reset_index()

peak_to_gene = peak_to_gene.merge(
    peak_features,
    on="peak_id",
    how="left"
)

# ============================================================
# Add top accessible cell type per peak
# ============================================================

mean_cols = [
    col for col in peak_to_gene.columns
    if col.startswith("mean_accessibility_")
]

detect_cols = [
    col for col in peak_to_gene.columns
    if col.startswith("detection_fraction_")
]

def extract_cell_type_from_safe_col(col, prefix):
    safe = col.replace(prefix, "")

    safe_to_cell_type = {
        (
            ct.replace("+", "plus")
            .replace(" ", "_")
            .replace(".", "_")
            .replace("-", "_")
        ): ct
        for ct in cell_type_order
    }

    return safe_to_cell_type.get(safe, safe)

peak_to_gene["peak_top_accessibility_cell_type"] = (
    peak_to_gene[mean_cols]
    .idxmax(axis=1)
    .apply(lambda x: extract_cell_type_from_safe_col(x, "mean_accessibility_"))
)

peak_to_gene["peak_max_mean_accessibility"] = peak_to_gene[mean_cols].max(axis=1)

peak_to_gene["peak_top_detection_cell_type"] = (
    peak_to_gene[detect_cols]
    .idxmax(axis=1)
    .apply(lambda x: extract_cell_type_from_safe_col(x, "detection_fraction_"))
)

peak_to_gene["peak_max_detection_fraction"] = peak_to_gene[detect_cols].max(axis=1)

# ============================================================
# Add RNA–ATAC gene-level concordance from Step 44
# ============================================================

if gene_corr_path.exists():
    gene_corr = pd.read_csv(gene_corr_path)

    gene_corr = gene_corr.rename(columns={
        "gene": "target_gene",
        "rna_atac_pseudobulk_correlation": "target_gene_rna_atac_pseudobulk_correlation",
        "n_cell_types": "target_gene_correlation_n_cell_types",
    })

    peak_to_gene = peak_to_gene.merge(
        gene_corr[
            [
                "target_gene",
                "target_gene_rna_atac_pseudobulk_correlation",
                "target_gene_correlation_n_cell_types",
            ]
        ],
        on="target_gene",
        how="left"
    )
else:
    peak_to_gene["target_gene_rna_atac_pseudobulk_correlation"] = np.nan
    peak_to_gene["target_gene_correlation_n_cell_types"] = np.nan

if concordant_gene_path.exists():
    concordant_genes = pd.read_csv(concordant_gene_path)
    concordant_gene_set = set(concordant_genes["gene"].astype(str))

    peak_to_gene["target_gene_concordant_ge_0_40"] = (
        peak_to_gene["target_gene"].astype(str).isin(concordant_gene_set)
    )
else:
    peak_to_gene["target_gene_concordant_ge_0_40"] = (
        peak_to_gene["target_gene_rna_atac_pseudobulk_correlation"] >= 0.40
    )

# ============================================================
# Sort table for downstream GRN construction
# ============================================================

sort_cols = [
    "target_gene_concordant_ge_0_40",
    "target_gene_rna_atac_pseudobulk_correlation",
    "peak_max_mean_accessibility",
    "peak_max_detection_fraction",
]

peak_to_gene = peak_to_gene.sort_values(
    sort_cols,
    ascending=[False, False, False, False]
)

# ============================================================
# Save peak-to-gene annotation table
# ============================================================

peak_to_gene.to_csv(peak_to_gene_out, index=False)

print("\nSaved peak-to-gene annotation table to:")
print(peak_to_gene_out)

print("\nPeak-to-gene table shape:")
print(peak_to_gene.shape)

print("\nPeak-to-gene table columns:")
print(list(peak_to_gene.columns))

print("\nTop 20 peak-to-gene rows:")
print(peak_to_gene.head(20).to_string(index=False))

# ============================================================
# Summary table
# ============================================================

summary_records = []

summary_records.append({
    "metric": "atac_cells",
    "value": atac.n_obs,
})

summary_records.append({
    "metric": "atac_peaks",
    "value": atac.n_vars,
})

summary_records.append({
    "metric": "peak_to_gene_rows",
    "value": peak_to_gene.shape[0],
})

summary_records.append({
    "metric": "unique_annotated_peaks",
    "value": peak_to_gene["peak_id"].nunique(),
})

summary_records.append({
    "metric": "unique_target_genes",
    "value": peak_to_gene["target_gene"].nunique(),
})

summary_records.append({
    "metric": "targets_with_rna_atac_correlation",
    "value": peak_to_gene["target_gene_rna_atac_pseudobulk_correlation"].notna().sum(),
})

summary_records.append({
    "metric": "peak_gene_links_with_concordant_target_gene_ge_0_40",
    "value": int(peak_to_gene["target_gene_concordant_ge_0_40"].sum()),
})

for cell_type, n_cells in cell_counts.items():
    summary_records.append({
        "metric": f"cells_{cell_type}",
        "value": n_cells,
    })

summary = pd.DataFrame(summary_records)

summary.to_csv(peak_to_gene_summary_out, index=False)

print("\nSaved Step 45 summary to:")
print(peak_to_gene_summary_out)

print("\nStep 45 summary:")
print(summary.to_string(index=False))





"""Βήμα 46 — Motif-to-TF table από το peak motif mapping για τα 5 ATAC cell types"""

import pandas as pd
import numpy as np
import re
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

motif_mapping_path = project_dir / "data/raw/pbmc_atac_10k/atac_pbmc_10k_nextgem_peak_motif_mapping.bed"

peak_to_gene_path = project_dir / "data/processed/pbmc_atac_peak_to_gene_annotation_5types.csv"

peak_motif_out = project_dir / "data/processed/pbmc_atac_peak_motif_mapping_5types.csv"
motif_to_tf_out = project_dir / "data/processed/pbmc_atac_motif_to_tf_table.csv"
step46_summary_out = project_dir / "data/processed/pbmc_atac_motif_to_tf_step46_summary.csv"

# ============================================================
# Load Step 45 peak-to-gene table to know valid 5-type peaks
# ============================================================

peak_to_gene = pd.read_csv(
    peak_to_gene_path,
    usecols=["peak_id"]
)

valid_peak_set = set(peak_to_gene["peak_id"].astype(str))

print("Valid annotated peaks from Step 45:")
print(len(valid_peak_set))

# ============================================================
# Inspect motif mapping BED structure
# ============================================================

preview = pd.read_csv(
    motif_mapping_path,
    sep="\t",
    header=None,
    comment="#",
    nrows=10
)

print("\nMotif mapping preview:")
print(preview)

print("\nNumber of columns in motif mapping BED:")
print(preview.shape[1])

if preview.shape[1] < 4:
    raise ValueError(
        "Expected at least 4 BED columns: chrom, start, end, motif."
    )

# ============================================================
# Column names
# ============================================================

n_cols = preview.shape[1]

col_names = ["chrom", "start", "end", "motif_raw"]

for i in range(4, n_cols):
    col_names.append(f"extra_col_{i}")

print("\nAssigned column names:")
print(col_names)

# ============================================================
# Motif parsing helper
# ============================================================

def parse_motif_record(motif_value):
    """
    Converts raw motif labels into:
    - motif_id
    - motif_name
    - list of candidate TF symbols

    This is intentionally flexible because 10x/JASPAR motif labels can appear as:
    MA0139.1_CTCF
    CTCF_MA0139.1
    FOS::JUN
    MA0477.1_FOSL1::JUN
    """
    s = str(motif_value).strip()

    if s == "" or s.lower() in {"nan", "none", "na"}:
        return np.nan, np.nan, []

    # Find JASPAR-like motif ID if present
    motif_id_match = re.search(r"(MA\d+(?:\.\d+)?)", s)

    if motif_id_match:
        motif_id = motif_id_match.group(1)

        # Remove motif ID from label to infer TF name
        tf_label = (
            s[:motif_id_match.start()]
            + s[motif_id_match.end():]
        )
        tf_label = re.sub(r"^[\s_\-:|]+|[\s_\-:|]+$", "", tf_label)
    else:
        motif_id = s
        tf_label = s

    if tf_label == "":
        tf_label = s

    motif_name = tf_label

    # Clean common decorations
    tf_label = re.sub(r"\(.*?\)", "", tf_label)
    tf_label = tf_label.replace(" ", "")

    # Split TF complexes such as FOS::JUN
    tf_candidates = re.split(r"::", tf_label)

    cleaned_tfs = []

    for tf in tf_candidates:
        tf = str(tf).strip()
        tf = re.sub(r"^[^A-Za-z0-9]+|[^A-Za-z0-9]+$", "", tf)

        # Remove remaining motif IDs if present
        tf = re.sub(r"MA\d+(?:\.\d+)?", "", tf)
        tf = re.sub(r"^[\s_\-:|]+|[\s_\-:|]+$", "", tf)

        # Human TF symbols are usually uppercase
        tf = tf.upper()

        if tf != "" and tf.lower() not in {"nan", "none", "na"}:
            cleaned_tfs.append(tf)

    cleaned_tfs = list(dict.fromkeys(cleaned_tfs))

    return motif_id, motif_name, cleaned_tfs

# ============================================================
# Process motif mapping in chunks
# ============================================================

chunk_size = 500_000

first_write = True

total_rows = 0
kept_rows = 0

unique_peaks_with_motifs = set()
unique_motifs = set()
unique_tfs = set()

for chunk_idx, chunk in enumerate(
    pd.read_csv(
        motif_mapping_path,
        sep="\t",
        header=None,
        comment="#",
        names=col_names,
        chunksize=chunk_size
    )
):
    total_rows += chunk.shape[0]

    # ------------------------------------------------------------
    # Clean BED coordinates
    # ------------------------------------------------------------

    chunk = chunk.copy()

    chunk["chrom"] = chunk["chrom"].astype(str).str.strip()
    chunk["start"] = pd.to_numeric(chunk["start"], errors="coerce")
    chunk["end"] = pd.to_numeric(chunk["end"], errors="coerce")

    chunk = chunk[
        chunk["chrom"].notna() &
        chunk["start"].notna() &
        chunk["end"].notna()
    ].copy()

    chunk["start"] = chunk["start"].astype(int)
    chunk["end"] = chunk["end"].astype(int)

    chunk["peak_id"] = (
        chunk["chrom"].astype(str)
        + ":"
        + chunk["start"].astype(str)
        + "-"
        + chunk["end"].astype(str)
    )

    # ------------------------------------------------------------
    # Keep only peaks used in Step 45
    # ------------------------------------------------------------

    chunk = chunk[
        chunk["peak_id"].isin(valid_peak_set)
    ].copy()

    if chunk.empty:
        print(f"Chunk {chunk_idx}: no valid Step 45 peaks.")
        continue

    # ------------------------------------------------------------
    # Handle possible comma-separated motif labels
    # ------------------------------------------------------------

    chunk["motif_raw"] = chunk["motif_raw"].astype(str).str.strip()
    chunk["motif_raw"] = chunk["motif_raw"].str.split(",")

    chunk = chunk.explode("motif_raw")
    chunk["motif_raw"] = chunk["motif_raw"].astype(str).str.strip()

    # ------------------------------------------------------------
    # Parse motif labels into motif ID and candidate TF
    # ------------------------------------------------------------

    parsed = chunk["motif_raw"].map(parse_motif_record)

    parsed_df = pd.DataFrame(
        parsed.tolist(),
        columns=["motif_id", "motif_name", "tf_list"],
        index=chunk.index
    )

    chunk = pd.concat([chunk, parsed_df], axis=1)

    chunk = chunk.explode("tf_list")
    chunk = chunk.rename(columns={"tf_list": "tf"})

    chunk["tf"] = chunk["tf"].astype(str).str.strip().str.upper()

    chunk = chunk[
        chunk["tf"].notna() &
        ~chunk["tf"].isin(["", "NAN", "NONE", "NA"])
    ].copy()

    # ------------------------------------------------------------
    # Keep essential columns
    # ------------------------------------------------------------

    keep_cols = [
        "peak_id",
        "chrom",
        "start",
        "end",
        "motif_raw",
        "motif_id",
        "motif_name",
        "tf"
    ]

    chunk_out = chunk[keep_cols].drop_duplicates().copy()

    # ------------------------------------------------------------
    # Update summary sets
    # ------------------------------------------------------------

    kept_rows += chunk_out.shape[0]

    unique_peaks_with_motifs.update(chunk_out["peak_id"].astype(str).unique())
    unique_motifs.update(chunk_out["motif_id"].astype(str).unique())
    unique_tfs.update(chunk_out["tf"].astype(str).unique())

    # ------------------------------------------------------------
    # Append to output file
    # ------------------------------------------------------------

    chunk_out.to_csv(
        peak_motif_out,
        index=False,
        mode="w" if first_write else "a",
        header=first_write
    )

    first_write = False

    print(
        f"Chunk {chunk_idx}: kept {chunk_out.shape[0]} peak-motif-TF rows"
    )

print("\nFinished processing motif mapping.")

# ============================================================
# Load final peak-motif table and remove global duplicates
# ============================================================

peak_motif = pd.read_csv(peak_motif_out)

print("\nPeak-motif-TF table before global deduplication:")
print(peak_motif.shape)

peak_motif = peak_motif.drop_duplicates()

peak_motif.to_csv(peak_motif_out, index=False)

print("\nPeak-motif-TF table after global deduplication:")
print(peak_motif.shape)

# ============================================================
# Build motif-to-TF table
# ============================================================

motif_to_tf = (
    peak_motif[
        ["motif_raw", "motif_id", "motif_name", "tf"]
    ]
    .drop_duplicates()
    .sort_values(["tf", "motif_id", "motif_raw"])
)

motif_to_tf.to_csv(motif_to_tf_out, index=False)

print("\nSaved peak-motif-TF table to:")
print(peak_motif_out)

print("\nSaved motif-to-TF table to:")
print(motif_to_tf_out)

print("\nMotif-to-TF table preview:")
print(motif_to_tf.head(30).to_string(index=False))

# ============================================================
# Save Step 46 summary
# ============================================================

summary = pd.DataFrame([
    {
        "metric": "motif_mapping_total_rows_processed",
        "value": total_rows,
    },
    {
        "metric": "peak_motif_tf_rows_kept_before_global_deduplication",
        "value": kept_rows,
    },
    {
        "metric": "peak_motif_tf_rows_after_global_deduplication",
        "value": peak_motif.shape[0],
    },
    {
        "metric": "unique_peaks_with_motifs",
        "value": peak_motif["peak_id"].nunique(),
    },
    {
        "metric": "unique_motifs",
        "value": peak_motif["motif_id"].nunique(),
    },
    {
        "metric": "unique_candidate_tfs",
        "value": peak_motif["tf"].nunique(),
    },
])

summary.to_csv(step46_summary_out, index=False)

print("\nSaved Step 46 summary to:")
print(step46_summary_out)

print("\nStep 46 summary:")
print(summary.to_string(index=False))

"""Βήμα 47 — Candidate TF–target edge construction για GRN inference"""

import pandas as pd
import numpy as np
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

peak_to_gene_path = project_dir / "data/processed/pbmc_atac_peak_to_gene_annotation_5types.csv"
peak_motif_path = project_dir / "data/processed/pbmc_atac_peak_motif_mapping_5types.csv"

rna_pseudobulk_path = project_dir / "data/processed/pbmc_rna_pseudobulk_5types_common_genes.csv"

peak_level_edges_out = project_dir / "data/processed/pbmc_grn_candidate_peak_level_edges_filtered.csv"
collapsed_edges_out = project_dir / "data/processed/pbmc_grn_candidate_tf_target_edges_collapsed.csv"
step47_summary_out = project_dir / "data/processed/pbmc_grn_candidate_edges_step47_summary.csv"

# ============================================================
# Parameters
# ============================================================

MIN_TARGET_RNA_ATAC_CORR = 0.40
MIN_PEAK_MAX_DETECTION_FRACTION = 0.01
MIN_PEAK_MAX_MEAN_ACCESSIBILITY = 0.005
MIN_TF_RNA_IN_PEAK_TOP_CELL_TYPE = 0.05
MIN_TARGET_RNA_IN_PEAK_TOP_CELL_TYPE = 0.05

EXCLUDE_SELF_LOOPS = True

CHUNK_SIZE = 500_000

# ============================================================
# Cell type order
# ============================================================

cell_type_order = [
    "CD14+ Monocytes",
    "CD4.Memory",
    "CD4.Naive",
    "CD8.effector",
    "CD8.Naive",
]

cell_type_to_safe = {
    "CD14+ Monocytes": "CD14plus_Monocytes",
    "CD4.Memory": "CD4_Memory",
    "CD4.Naive": "CD4_Naive",
    "CD8.effector": "CD8_effector",
    "CD8.Naive": "CD8_Naive",
}

# ============================================================
# Load Step 45 peak-to-gene table
# ============================================================

peak_to_gene = pd.read_csv(peak_to_gene_path)

print("Peak-to-gene table:")
print(peak_to_gene.shape)

required_peak_to_gene_cols = [
    "peak_id",
    "target_gene",
    "peak_annotation_type",
    "distance_to_gene",
    "peak_top_accessibility_cell_type",
    "peak_max_mean_accessibility",
    "peak_top_detection_cell_type",
    "peak_max_detection_fraction",
    "target_gene_rna_atac_pseudobulk_correlation",
    "target_gene_concordant_ge_0_40",
]

missing_cols = [
    col for col in required_peak_to_gene_cols
    if col not in peak_to_gene.columns
]

if missing_cols:
    raise ValueError(f"Missing required columns from Step 45 table: {missing_cols}")

# ============================================================
# Pre-filter Step 45 peak-to-gene links
# ============================================================

peak_to_gene["target_gene"] = peak_to_gene["target_gene"].astype(str).str.strip()

peak_to_gene["target_gene_rna_atac_pseudobulk_correlation"] = pd.to_numeric(
    peak_to_gene["target_gene_rna_atac_pseudobulk_correlation"],
    errors="coerce"
)

peak_to_gene["peak_max_detection_fraction"] = pd.to_numeric(
    peak_to_gene["peak_max_detection_fraction"],
    errors="coerce"
)

peak_to_gene["peak_max_mean_accessibility"] = pd.to_numeric(
    peak_to_gene["peak_max_mean_accessibility"],
    errors="coerce"
)

peak_to_gene_filtered = peak_to_gene[
    (peak_to_gene["target_gene_rna_atac_pseudobulk_correlation"] >= MIN_TARGET_RNA_ATAC_CORR) &
    (peak_to_gene["peak_max_detection_fraction"] >= MIN_PEAK_MAX_DETECTION_FRACTION) &
    (peak_to_gene["peak_max_mean_accessibility"] >= MIN_PEAK_MAX_MEAN_ACCESSIBILITY)
].copy()

print("\nPeak-to-gene links before filtering:")
print(peak_to_gene.shape[0])

print("\nPeak-to-gene links after Step 47 pre-filtering:")
print(peak_to_gene_filtered.shape[0])

# ============================================================
# Load RNA pseudo-bulk for TF and target expression evidence
# ============================================================

rna_pb = pd.read_csv(rna_pseudobulk_path, index_col=0)

rna_pb.index = rna_pb.index.astype(str).str.strip()
rna_pb.columns = rna_pb.columns.astype(str).str.strip().str.upper()

# Make target gene and TF comparisons uppercase
rna_gene_set = set(rna_pb.columns)

print("\nRNA pseudo-bulk shape:")
print(rna_pb.shape)

# ============================================================
# Helper dictionaries for fast expression lookup
# ============================================================

rna_expr_by_cell_type = {}

for cell_type in cell_type_order:
    if cell_type not in rna_pb.index:
        raise ValueError(f"Missing RNA pseudo-bulk cell type: {cell_type}")

    rna_expr_by_cell_type[cell_type] = rna_pb.loc[cell_type].to_dict()

# ============================================================
# Add target-gene uppercase column to peak-to-gene
# ============================================================

peak_to_gene_filtered["target_gene_original"] = peak_to_gene_filtered["target_gene"]
peak_to_gene_filtered["target_gene"] = (
    peak_to_gene_filtered["target_gene"]
    .astype(str)
    .str.strip()
    .str.upper()
)

# ============================================================
# Keep useful columns from Step 45
# ============================================================

mean_accessibility_cols = [
    col for col in peak_to_gene_filtered.columns
    if col.startswith("mean_accessibility_")
]

detection_fraction_cols = [
    col for col in peak_to_gene_filtered.columns
    if col.startswith("detection_fraction_")
]

step45_keep_cols = [
    "peak_id",
    "target_gene",
    "target_gene_original",
    "distance_to_gene",
    "peak_annotation_type",
    "peak_top_accessibility_cell_type",
    "peak_max_mean_accessibility",
    "peak_top_detection_cell_type",
    "peak_max_detection_fraction",
    "target_gene_rna_atac_pseudobulk_correlation",
    "target_gene_concordant_ge_0_40",
] + mean_accessibility_cols + detection_fraction_cols

peak_to_gene_filtered = peak_to_gene_filtered[step45_keep_cols].copy()

# ============================================================
# Function to add RNA expression evidence
# ============================================================

def add_expression_evidence(edge_chunk):
    """
    Adds:
    - TF RNA expression in the cell type where the peak is most accessible
    - target gene RNA expression in the cell type where the peak is most accessible
    """
    edge_chunk = edge_chunk.copy()

    edge_chunk["tf"] = edge_chunk["tf"].astype(str).str.strip().str.upper()
    edge_chunk["target_gene"] = edge_chunk["target_gene"].astype(str).str.strip().str.upper()

    edge_chunk["tf_in_rna_pseudobulk"] = edge_chunk["tf"].isin(rna_gene_set)
    edge_chunk["target_gene_in_rna_pseudobulk"] = edge_chunk["target_gene"].isin(rna_gene_set)

    edge_chunk["tf_rna_in_peak_top_cell_type"] = np.nan
    edge_chunk["target_gene_rna_in_peak_top_cell_type"] = np.nan

    for cell_type in cell_type_order:
        mask = edge_chunk["peak_top_accessibility_cell_type"].astype(str) == cell_type

        if mask.sum() == 0:
            continue

        expr_dict = rna_expr_by_cell_type[cell_type]

        edge_chunk.loc[mask, "tf_rna_in_peak_top_cell_type"] = (
            edge_chunk.loc[mask, "tf"].map(expr_dict)
        )

        edge_chunk.loc[mask, "target_gene_rna_in_peak_top_cell_type"] = (
            edge_chunk.loc[mask, "target_gene"].map(expr_dict)
        )

    return edge_chunk

# ============================================================
# Function to score candidate peak-level edges
# ============================================================

def score_candidate_edges(edge_chunk):
    """
    Heuristic score for ranking candidate regulatory evidence.

    This is not a statistical probability.
    It combines:
    - target gene RNA–ATAC concordance
    - peak accessibility
    - peak detection fraction
    - TF RNA expression in the peak-top cell type
    - target gene RNA expression in the peak-top cell type
    """
    edge_chunk = edge_chunk.copy()

    corr_score = edge_chunk["target_gene_rna_atac_pseudobulk_correlation"].clip(lower=0, upper=1)

    accessibility_score = np.log1p(
        edge_chunk["peak_max_mean_accessibility"].fillna(0)
    )

    detection_score = edge_chunk["peak_max_detection_fraction"].fillna(0)

    tf_expr_score = np.log1p(
        edge_chunk["tf_rna_in_peak_top_cell_type"].fillna(0)
    )

    target_expr_score = np.log1p(
        edge_chunk["target_gene_rna_in_peak_top_cell_type"].fillna(0)
    )

    edge_chunk["edge_evidence_score"] = (
        corr_score *
        accessibility_score *
        detection_score *
        tf_expr_score *
        target_expr_score
    )

    edge_chunk["is_promoter_peak"] = (
        edge_chunk["peak_annotation_type"]
        .astype(str)
        .str.lower()
        .str.contains("promoter", na=False)
    )

    edge_chunk["is_distal_peak"] = (
        edge_chunk["peak_annotation_type"]
        .astype(str)
        .str.lower()
        .str.contains("distal", na=False)
    )

    edge_chunk["is_self_loop"] = (
        edge_chunk["tf"].astype(str).str.upper()
        ==
        edge_chunk["target_gene"].astype(str).str.upper()
    )

    return edge_chunk

# ============================================================
# Merge Step 45 and Step 46 in chunks
# ============================================================

first_write = True

total_peak_motif_rows = 0
total_merged_rows = 0
total_filtered_edge_rows = 0

for chunk_idx, peak_motif_chunk in enumerate(
    pd.read_csv(
        peak_motif_path,
        chunksize=CHUNK_SIZE
    )
):
    total_peak_motif_rows += peak_motif_chunk.shape[0]

    peak_motif_chunk["peak_id"] = peak_motif_chunk["peak_id"].astype(str)
    peak_motif_chunk["tf"] = peak_motif_chunk["tf"].astype(str).str.strip().str.upper()

    # ------------------------------------------------------------
    # Merge: peak → gene with peak → motif → TF
    # ------------------------------------------------------------

    edge_chunk = peak_to_gene_filtered.merge(
        peak_motif_chunk,
        on="peak_id",
        how="inner"
    )

    total_merged_rows += edge_chunk.shape[0]

    if edge_chunk.empty:
        print(f"Chunk {chunk_idx}: no merged edge rows.")
        continue

    # ------------------------------------------------------------
    # Add RNA expression evidence and score
    # ------------------------------------------------------------

    edge_chunk = add_expression_evidence(edge_chunk)
    edge_chunk = score_candidate_edges(edge_chunk)

    # ------------------------------------------------------------
    # Final filtering for candidate regulatory evidence
    # ------------------------------------------------------------

    edge_chunk = edge_chunk[
        (edge_chunk["tf_in_rna_pseudobulk"]) &
        (edge_chunk["target_gene_in_rna_pseudobulk"]) &
        (edge_chunk["tf_rna_in_peak_top_cell_type"] >= MIN_TF_RNA_IN_PEAK_TOP_CELL_TYPE) &
        (edge_chunk["target_gene_rna_in_peak_top_cell_type"] >= MIN_TARGET_RNA_IN_PEAK_TOP_CELL_TYPE)
    ].copy()

    if EXCLUDE_SELF_LOOPS:
        edge_chunk = edge_chunk[
            ~edge_chunk["is_self_loop"]
        ].copy()

    edge_chunk = edge_chunk.sort_values(
        "edge_evidence_score",
        ascending=False
    )

    total_filtered_edge_rows += edge_chunk.shape[0]

    # ------------------------------------------------------------
    # Save filtered peak-level candidate edges
    # ------------------------------------------------------------

    edge_chunk.to_csv(
        peak_level_edges_out,
        index=False,
        mode="w" if first_write else "a",
        header=first_write
    )

    first_write = False

    print(
        f"Chunk {chunk_idx}: merged {edge_chunk.shape[0]} filtered candidate peak-level edges"
    )

print("\nFinished constructing filtered peak-level candidate edges.")

print("\nTotal peak-motif rows processed:")
print(total_peak_motif_rows)

print("\nTotal merged peak-level rows before final TF/target expression filtering:")
print(total_merged_rows)

print("\nTotal filtered peak-level candidate edges:")
print(total_filtered_edge_rows)

# ============================================================
# Load filtered peak-level edges
# ============================================================

if total_filtered_edge_rows == 0:
    raise ValueError(
        "No candidate edges passed the filters. "
        "Try lowering MIN_TF_RNA_IN_PEAK_TOP_CELL_TYPE or MIN_PEAK_MAX_DETECTION_FRACTION."
    )

peak_level_edges = pd.read_csv(peak_level_edges_out)

print("\nFiltered peak-level candidate edges:")
print(peak_level_edges.shape)

# ============================================================
# Collapse peak-level evidence to TF-target edges
# ============================================================

peak_level_edges["tf"] = peak_level_edges["tf"].astype(str).str.upper()
peak_level_edges["target_gene"] = peak_level_edges["target_gene"].astype(str).str.upper()

# Top evidence row per TF-target edge
top_evidence = (
    peak_level_edges
    .sort_values("edge_evidence_score", ascending=False)
    .drop_duplicates(["tf", "target_gene"])
    [
        [
            "tf",
            "target_gene",
            "peak_id",
            "motif_id",
            "motif_name",
            "peak_annotation_type",
            "distance_to_gene",
            "peak_top_accessibility_cell_type",
            "peak_max_mean_accessibility",
            "peak_max_detection_fraction",
            "target_gene_rna_atac_pseudobulk_correlation",
            "tf_rna_in_peak_top_cell_type",
            "target_gene_rna_in_peak_top_cell_type",
            "edge_evidence_score",
        ]
    ]
    .rename(columns={
        "peak_id": "top_evidence_peak_id",
        "motif_id": "top_evidence_motif_id",
        "motif_name": "top_evidence_motif_name",
        "peak_annotation_type": "top_evidence_peak_annotation_type",
        "distance_to_gene": "top_evidence_distance_to_gene",
        "peak_top_accessibility_cell_type": "top_evidence_cell_type",
        "peak_max_mean_accessibility": "top_evidence_peak_max_mean_accessibility",
        "peak_max_detection_fraction": "top_evidence_peak_max_detection_fraction",
        "target_gene_rna_atac_pseudobulk_correlation": "top_evidence_target_gene_rna_atac_corr",
        "tf_rna_in_peak_top_cell_type": "top_evidence_tf_rna",
        "target_gene_rna_in_peak_top_cell_type": "top_evidence_target_gene_rna",
        "edge_evidence_score": "top_evidence_score",
    })
)

# Group-level edge evidence
collapsed = (
    peak_level_edges
    .groupby(["tf", "target_gene"])
    .agg(
        n_supporting_peak_motif_links=("peak_id", "count"),
        n_unique_supporting_peaks=("peak_id", "nunique"),
        n_unique_motifs=("motif_id", "nunique"),
        n_promoter_peak_links=("is_promoter_peak", "sum"),
        n_distal_peak_links=("is_distal_peak", "sum"),
        max_edge_evidence_score=("edge_evidence_score", "max"),
        mean_edge_evidence_score=("edge_evidence_score", "mean"),
        max_target_gene_rna_atac_corr=("target_gene_rna_atac_pseudobulk_correlation", "max"),
        max_peak_mean_accessibility=("peak_max_mean_accessibility", "max"),
        max_peak_detection_fraction=("peak_max_detection_fraction", "max"),
        max_tf_rna_in_peak_top_cell_type=("tf_rna_in_peak_top_cell_type", "max"),
        max_target_gene_rna_in_peak_top_cell_type=("target_gene_rna_in_peak_top_cell_type", "max"),
        supported_cell_types=(
            "peak_top_accessibility_cell_type",
            lambda x: ";".join(sorted(set(x.astype(str))))
        ),
    )
    .reset_index()
)

collapsed = collapsed.merge(
    top_evidence,
    on=["tf", "target_gene"],
    how="left"
)

collapsed = collapsed.sort_values(
    [
        "max_edge_evidence_score",
        "n_unique_supporting_peaks",
        "max_target_gene_rna_atac_corr",
    ],
    ascending=[False, False, False]
)

collapsed.to_csv(collapsed_edges_out, index=False)

print("\nSaved collapsed TF-target candidate edges to:")
print(collapsed_edges_out)

print("\nCollapsed candidate TF-target edges:")
print(collapsed.shape)

print("\nTop 30 candidate TF-target edges:")
print(collapsed.head(30).to_string(index=False))

# ============================================================
# Save Step 47 summary
# ============================================================

summary = pd.DataFrame([
    {
        "metric": "peak_to_gene_rows_input",
        "value": peak_to_gene.shape[0],
    },
    {
        "metric": "peak_to_gene_rows_after_prefiltering",
        "value": peak_to_gene_filtered.shape[0],
    },
    {
        "metric": "peak_motif_rows_processed",
        "value": total_peak_motif_rows,
    },
    {
        "metric": "merged_peak_level_rows_before_expression_filtering",
        "value": total_merged_rows,
    },
    {
        "metric": "filtered_peak_level_candidate_edges",
        "value": peak_level_edges.shape[0],
    },
    {
        "metric": "collapsed_tf_target_edges",
        "value": collapsed.shape[0],
    },
    {
        "metric": "unique_tfs_in_candidate_edges",
        "value": collapsed["tf"].nunique(),
    },
    {
        "metric": "unique_target_genes_in_candidate_edges",
        "value": collapsed["target_gene"].nunique(),
    },
    {
        "metric": "min_target_rna_atac_corr",
        "value": MIN_TARGET_RNA_ATAC_CORR,
    },
    {
        "metric": "min_peak_max_detection_fraction",
        "value": MIN_PEAK_MAX_DETECTION_FRACTION,
    },
    {
        "metric": "min_peak_max_mean_accessibility",
        "value": MIN_PEAK_MAX_MEAN_ACCESSIBILITY,
    },
    {
        "metric": "min_tf_rna_in_peak_top_cell_type",
        "value": MIN_TF_RNA_IN_PEAK_TOP_CELL_TYPE,
    },
    {
        "metric": "min_target_rna_in_peak_top_cell_type",
        "value": MIN_TARGET_RNA_IN_PEAK_TOP_CELL_TYPE,
    },
    {
        "metric": "exclude_self_loops",
        "value": EXCLUDE_SELF_LOOPS,
    },
])

summary.to_csv(step47_summary_out, index=False)

print("\nSaved Step 47 summary to:")
print(step47_summary_out)

print("\nStep 47 summary:")
print(summary.to_string(index=False))





"""Βήμα 48 — Quality control και biological inspection των candidate GRN edges"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

collapsed_edges_path = project_dir / "data/processed/pbmc_grn_candidate_tf_target_edges_collapsed.csv"
peak_level_edges_path = project_dir / "data/processed/pbmc_grn_candidate_peak_level_edges_filtered.csv"

qc_table_dir = project_dir / "data/processed/grn_qc"
qc_figure_dir = project_dir / "figures/grn_qc"

qc_table_dir.mkdir(parents=True, exist_ok=True)
qc_figure_dir.mkdir(parents=True, exist_ok=True)

# ============================================================
# Load candidate GRN edge tables
# ============================================================

collapsed = pd.read_csv(collapsed_edges_path)

print("Collapsed TF-target candidate edges:")
print(collapsed.shape)

print("\nCollapsed edge columns:")
print(list(collapsed.columns))

if peak_level_edges_path.exists():
    peak_level = pd.read_csv(peak_level_edges_path)
    print("\nPeak-level candidate edges:")
    print(peak_level.shape)
else:
    peak_level = None
    print("\nPeak-level edge table not found. Continuing with collapsed table only.")

# ============================================================
# Basic cleaning
# ============================================================

collapsed["tf"] = collapsed["tf"].astype(str).str.strip().str.upper()
collapsed["target_gene"] = collapsed["target_gene"].astype(str).str.strip().str.upper()

numeric_cols = [
    "n_supporting_peak_motif_links",
    "n_unique_supporting_peaks",
    "n_unique_motifs",
    "n_promoter_peak_links",
    "n_distal_peak_links",
    "max_edge_evidence_score",
    "mean_edge_evidence_score",
    "max_target_gene_rna_atac_corr",
    "max_peak_mean_accessibility",
    "max_peak_detection_fraction",
    "max_tf_rna_in_peak_top_cell_type",
    "max_target_gene_rna_in_peak_top_cell_type",
    "top_evidence_score",
]

for col in numeric_cols:
    if col in collapsed.columns:
        collapsed[col] = pd.to_numeric(collapsed[col], errors="coerce")

# ============================================================
# Overall GRN summary
# ============================================================

score_q50 = collapsed["max_edge_evidence_score"].quantile(0.50)
score_q75 = collapsed["max_edge_evidence_score"].quantile(0.75)
score_q90 = collapsed["max_edge_evidence_score"].quantile(0.90)
score_q95 = collapsed["max_edge_evidence_score"].quantile(0.95)
score_q99 = collapsed["max_edge_evidence_score"].quantile(0.99)

overall_summary = pd.DataFrame([
    {"metric": "collapsed_tf_target_edges", "value": collapsed.shape[0]},
    {"metric": "unique_tfs", "value": collapsed["tf"].nunique()},
    {"metric": "unique_target_genes", "value": collapsed["target_gene"].nunique()},
    {"metric": "median_max_edge_evidence_score", "value": score_q50},
    {"metric": "q75_max_edge_evidence_score", "value": score_q75},
    {"metric": "q90_max_edge_evidence_score", "value": score_q90},
    {"metric": "q95_max_edge_evidence_score", "value": score_q95},
    {"metric": "q99_max_edge_evidence_score", "value": score_q99},
    {"metric": "edges_with_at_least_2_supporting_peaks", "value": int((collapsed["n_unique_supporting_peaks"] >= 2).sum())},
    {"metric": "edges_with_promoter_support", "value": int((collapsed["n_promoter_peak_links"] > 0).sum())},
    {"metric": "edges_with_distal_support", "value": int((collapsed["n_distal_peak_links"] > 0).sum())},
    {"metric": "edges_with_target_corr_ge_0_70", "value": int((collapsed["max_target_gene_rna_atac_corr"] >= 0.70).sum())},
])

overall_summary.to_csv(qc_table_dir / "step48_grn_overall_summary.csv", index=False)

print("\nOverall GRN summary:")
print(overall_summary.to_string(index=False))

# ============================================================
# Top TFs by number of target genes
# ============================================================

tf_summary = (
    collapsed
    .groupby("tf")
    .agg(
        n_target_genes=("target_gene", "nunique"),
        n_edges=("target_gene", "count"),
        total_supporting_peak_motif_links=("n_supporting_peak_motif_links", "sum"),
        total_unique_peak_support_across_edges=("n_unique_supporting_peaks", "sum"),
        total_unique_motif_support_across_edges=("n_unique_motifs", "sum"),
        n_promoter_peak_links=("n_promoter_peak_links", "sum"),
        n_distal_peak_links=("n_distal_peak_links", "sum"),
        max_edge_evidence_score=("max_edge_evidence_score", "max"),
        mean_edge_evidence_score=("max_edge_evidence_score", "mean"),
        median_edge_evidence_score=("max_edge_evidence_score", "median"),
        max_target_gene_rna_atac_corr=("max_target_gene_rna_atac_corr", "max"),
        mean_target_gene_rna_atac_corr=("max_target_gene_rna_atac_corr", "mean"),
    )
    .reset_index()
    .sort_values(
        ["n_target_genes", "max_edge_evidence_score"],
        ascending=[False, False]
    )
)

tf_summary.to_csv(qc_table_dir / "step48_top_tfs_by_targets.csv", index=False)

print("\nTop 30 TFs by number of target genes:")
print(tf_summary.head(30).to_string(index=False))

# ============================================================
# Top target genes by number of regulators
# ============================================================

target_summary = (
    collapsed
    .groupby("target_gene")
    .agg(
        n_regulator_tfs=("tf", "nunique"),
        n_edges=("tf", "count"),
        total_supporting_peak_motif_links=("n_supporting_peak_motif_links", "sum"),
        total_unique_peak_support_across_edges=("n_unique_supporting_peaks", "sum"),
        n_promoter_peak_links=("n_promoter_peak_links", "sum"),
        n_distal_peak_links=("n_distal_peak_links", "sum"),
        max_edge_evidence_score=("max_edge_evidence_score", "max"),
        mean_edge_evidence_score=("max_edge_evidence_score", "mean"),
        target_gene_rna_atac_corr=("max_target_gene_rna_atac_corr", "max"),
        max_target_gene_rna=("max_target_gene_rna_in_peak_top_cell_type", "max"),
    )
    .reset_index()
    .sort_values(
        ["n_regulator_tfs", "max_edge_evidence_score"],
        ascending=[False, False]
    )
)

target_summary.to_csv(qc_table_dir / "step48_top_target_genes_by_regulators.csv", index=False)

print("\nTop 30 target genes by number of regulator TFs:")
print(target_summary.head(30).to_string(index=False))

# ============================================================
# Top candidate edges by evidence score
# ============================================================

top_edges = (
    collapsed
    .sort_values(
        ["max_edge_evidence_score", "n_unique_supporting_peaks", "max_target_gene_rna_atac_corr"],
        ascending=[False, False, False]
    )
    .copy()
)

top_edges.to_csv(qc_table_dir / "step48_candidate_edges_ranked_by_evidence_score.csv", index=False)

print("\nTop 30 candidate TF-target edges by evidence score:")
cols_to_show = [
    "tf",
    "target_gene",
    "max_edge_evidence_score",
    "n_unique_supporting_peaks",
    "n_unique_motifs",
    "n_promoter_peak_links",
    "n_distal_peak_links",
    "max_target_gene_rna_atac_corr",
    "supported_cell_types",
    "top_evidence_peak_id",
    "top_evidence_cell_type",
]

available_cols_to_show = [col for col in cols_to_show if col in top_edges.columns]
print(top_edges[available_cols_to_show].head(30).to_string(index=False))

# ============================================================
# High-confidence edge subsets
# ============================================================

high_confidence_q90 = collapsed[
    (collapsed["max_edge_evidence_score"] >= score_q90) &
    (collapsed["max_target_gene_rna_atac_corr"] >= 0.70)
].copy()

high_confidence_q95 = collapsed[
    (collapsed["max_edge_evidence_score"] >= score_q95) &
    (collapsed["max_target_gene_rna_atac_corr"] >= 0.70)
].copy()

multi_peak_edges = collapsed[
    (collapsed["n_unique_supporting_peaks"] >= 2) &
    (collapsed["max_target_gene_rna_atac_corr"] >= 0.70)
].copy()

high_confidence_q90 = high_confidence_q90.sort_values("max_edge_evidence_score", ascending=False)
high_confidence_q95 = high_confidence_q95.sort_values("max_edge_evidence_score", ascending=False)
multi_peak_edges = multi_peak_edges.sort_values(
    ["n_unique_supporting_peaks", "max_edge_evidence_score"],
    ascending=[False, False]
)

high_confidence_q90.to_csv(qc_table_dir / "step48_high_confidence_edges_q90_score_corr_ge_0_70.csv", index=False)
high_confidence_q95.to_csv(qc_table_dir / "step48_high_confidence_edges_q95_score_corr_ge_0_70.csv", index=False)
multi_peak_edges.to_csv(qc_table_dir / "step48_multi_peak_supported_edges_corr_ge_0_70.csv", index=False)

print("\nHigh-confidence edge counts:")
print("Q90 score + target corr >= 0.70:", high_confidence_q90.shape[0])
print("Q95 score + target corr >= 0.70:", high_confidence_q95.shape[0])
print("Multi-peak supported + target corr >= 0.70:", multi_peak_edges.shape[0])

# ============================================================
# Promoter versus distal support
# ============================================================

promoter_distal_summary = pd.DataFrame([
    {
        "category": "edges_with_promoter_support",
        "n_edges": int((collapsed["n_promoter_peak_links"] > 0).sum()),
    },
    {
        "category": "edges_with_distal_support",
        "n_edges": int((collapsed["n_distal_peak_links"] > 0).sum()),
    },
    {
        "category": "edges_with_both_promoter_and_distal_support",
        "n_edges": int(
            ((collapsed["n_promoter_peak_links"] > 0) &
             (collapsed["n_distal_peak_links"] > 0)).sum()
        ),
    },
    {
        "category": "edges_with_only_promoter_support",
        "n_edges": int(
            ((collapsed["n_promoter_peak_links"] > 0) &
             (collapsed["n_distal_peak_links"] == 0)).sum()
        ),
    },
    {
        "category": "edges_with_only_distal_support",
        "n_edges": int(
            ((collapsed["n_promoter_peak_links"] == 0) &
             (collapsed["n_distal_peak_links"] > 0)).sum()
        ),
    },
])

promoter_distal_summary.to_csv(qc_table_dir / "step48_promoter_distal_support_summary.csv", index=False)

print("\nPromoter/distal support summary:")
print(promoter_distal_summary.to_string(index=False))

# ============================================================
# Supported cell type distribution
# ============================================================

celltype_edges = collapsed.copy()
celltype_edges["supported_cell_types"] = celltype_edges["supported_cell_types"].fillna("unknown")
celltype_edges["supported_cell_type"] = celltype_edges["supported_cell_types"].astype(str).str.split(";")
celltype_edges = celltype_edges.explode("supported_cell_type")
celltype_edges["supported_cell_type"] = celltype_edges["supported_cell_type"].astype(str).str.strip()

celltype_summary = (
    celltype_edges
    .groupby("supported_cell_type")
    .agg(
        n_edges=("target_gene", "count"),
        n_tfs=("tf", "nunique"),
        n_target_genes=("target_gene", "nunique"),
        max_edge_evidence_score=("max_edge_evidence_score", "max"),
        mean_edge_evidence_score=("max_edge_evidence_score", "mean"),
    )
    .reset_index()
    .sort_values("n_edges", ascending=False)
)

celltype_summary.to_csv(qc_table_dir / "step48_supported_cell_type_summary.csv", index=False)

print("\nSupported cell type summary:")
print(celltype_summary.to_string(index=False))

# ============================================================
# Known immune TF inspection
# ============================================================

known_immune_tfs = [
    "FOS", "FOSB", "FOSL1", "FOSL2",
    "JUN", "JUNB", "JUND",
    "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9",
    "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
    "NFKB1", "NFKB2", "RELA", "RELB", "REL",
    "SPI1", "CEBPA", "CEBPB", "CEBPD",
    "KLF2", "KLF3", "KLF4", "KLF6",
    "IKZF1", "IKZF2", "IKZF3",
    "ETS1", "ELF1", "ELK1",
    "RUNX1", "RUNX2", "RUNX3",
    "TCF7", "LEF1",
    "GATA3", "TBX21", "EOMES",
    "PRDM1", "BCL6", "FOXP1", "FOXO1",
    "MEF2A", "MEF2C",
    "NR4A1", "NR4A2", "NR4A3",
    "AHR", "BATF", "MAF", "MAFB",
]

immune_edges = collapsed[
    collapsed["tf"].isin(known_immune_tfs)
].copy()

immune_tf_summary = (
    immune_edges
    .groupby("tf")
    .agg(
        n_target_genes=("target_gene", "nunique"),
        n_edges=("target_gene", "count"),
        n_unique_supporting_peaks_total=("n_unique_supporting_peaks", "sum"),
        n_promoter_peak_links=("n_promoter_peak_links", "sum"),
        n_distal_peak_links=("n_distal_peak_links", "sum"),
        max_edge_evidence_score=("max_edge_evidence_score", "max"),
        mean_edge_evidence_score=("max_edge_evidence_score", "mean"),
        max_target_gene_rna_atac_corr=("max_target_gene_rna_atac_corr", "max"),
    )
    .reset_index()
    .sort_values(
        ["n_target_genes", "max_edge_evidence_score"],
        ascending=[False, False]
    )
)

immune_edges_ranked = immune_edges.sort_values(
    ["max_edge_evidence_score", "n_unique_supporting_peaks"],
    ascending=[False, False]
)

immune_tf_summary.to_csv(qc_table_dir / "step48_known_immune_tf_summary.csv", index=False)
immune_edges_ranked.to_csv(qc_table_dir / "step48_known_immune_tf_edges_ranked.csv", index=False)

print("\nKnown immune TFs found in candidate GRN:")
print(sorted(immune_edges["tf"].unique()))

print("\nKnown immune TF summary:")
print(immune_tf_summary.to_string(index=False))

print("\nTop 50 known immune TF candidate edges:")
print(immune_edges_ranked[available_cols_to_show].head(50).to_string(index=False))

# ============================================================
# Housekeeping/ribosomal/mitochondrial target gene inspection
# ============================================================

def classify_target_gene(gene):
    gene = str(gene).upper()

    if gene.startswith("RPL") or gene.startswith("RPS"):
        return "ribosomal"
    if gene.startswith("MT-"):
        return "mitochondrial"
    if gene.startswith("HSP") or gene.startswith("HSPA") or gene.startswith("HSPB") or gene.startswith("HSPD") or gene.startswith("HSPH"):
        return "heat_shock"
    if gene.startswith("HIST") or gene.startswith("H1-") or gene.startswith("H2A") or gene.startswith("H2B") or gene.startswith("H3") or gene.startswith("H4"):
        return "histone"
    if gene in {"MALAT1", "NEAT1", "B2M", "ACTB", "GAPDH", "TUBB", "TUBA1B"}:
        return "common_high_expression"
    return "other"

collapsed["target_gene_category"] = collapsed["target_gene"].map(classify_target_gene)

target_category_summary = (
    collapsed
    .groupby("target_gene_category")
    .agg(
        n_edges=("target_gene", "count"),
        n_target_genes=("target_gene", "nunique"),
        n_tfs=("tf", "nunique"),
        max_edge_evidence_score=("max_edge_evidence_score", "max"),
        mean_edge_evidence_score=("max_edge_evidence_score", "mean"),
    )
    .reset_index()
    .sort_values("n_edges", ascending=False)
)

target_category_summary.to_csv(qc_table_dir / "step48_target_gene_category_summary.csv", index=False)

non_housekeeping_edges = collapsed[
    collapsed["target_gene_category"] == "other"
].copy()

non_housekeeping_edges = non_housekeeping_edges.sort_values(
    "max_edge_evidence_score",
    ascending=False
)

non_housekeeping_edges.to_csv(qc_table_dir / "step48_non_housekeeping_candidate_edges.csv", index=False)

print("\nTarget gene category summary:")
print(target_category_summary.to_string(index=False))

print("\nTop 30 non-housekeeping-like candidate edges:")
print(non_housekeeping_edges[available_cols_to_show].head(30).to_string(index=False))

# ============================================================
# Optional exact peak-level TF summary if peak-level table is available
# ============================================================

if peak_level is not None:
    peak_level["tf"] = peak_level["tf"].astype(str).str.strip().str.upper()
    peak_level["target_gene"] = peak_level["target_gene"].astype(str).str.strip().str.upper()

    peak_level_tf_summary = (
        peak_level
        .groupby("tf")
        .agg(
            n_peak_level_rows=("target_gene", "count"),
            n_target_genes=("target_gene", "nunique"),
            n_unique_peaks=("peak_id", "nunique"),
            n_unique_motifs=("motif_id", "nunique"),
            max_edge_evidence_score=("edge_evidence_score", "max"),
            mean_edge_evidence_score=("edge_evidence_score", "mean"),
        )
        .reset_index()
        .sort_values(
            ["n_target_genes", "n_unique_peaks"],
            ascending=[False, False]
        )
    )

    peak_level_tf_summary.to_csv(
        qc_table_dir / "step48_peak_level_tf_summary.csv",
        index=False
    )

    print("\nPeak-level TF summary top 30:")
    print(peak_level_tf_summary.head(30).to_string(index=False))

# ============================================================
# Plot helper
# ============================================================

def save_current_plot(output_path):
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()

# ============================================================
# Plot 1 — Evidence score distribution
# ============================================================

plt.figure(figsize=(7, 5))
plt.hist(collapsed["max_edge_evidence_score"].dropna(), bins=60)
plt.xlabel("Max edge evidence score")
plt.ylabel("Number of TF-target edges")
plt.title("Distribution of candidate GRN edge evidence scores")
save_current_plot(qc_figure_dir / "step48_edge_evidence_score_distribution.png")

# ============================================================
# Plot 2 — Top TFs by number of targets
# ============================================================

top_n = 30
plot_df = tf_summary.head(top_n).sort_values("n_target_genes", ascending=True)

plt.figure(figsize=(8, 8))
plt.barh(plot_df["tf"], plot_df["n_target_genes"])
plt.xlabel("Number of target genes")
plt.ylabel("TF")
plt.title(f"Top {top_n} TFs by number of candidate target genes")
save_current_plot(qc_figure_dir / "step48_top_tfs_by_targets.png")

# ============================================================
# Plot 3 — Top target genes by number of regulators
# ============================================================

plot_df = target_summary.head(top_n).sort_values("n_regulator_tfs", ascending=True)

plt.figure(figsize=(8, 8))
plt.barh(plot_df["target_gene"], plot_df["n_regulator_tfs"])
plt.xlabel("Number of candidate regulator TFs")
plt.ylabel("Target gene")
plt.title(f"Top {top_n} target genes by number of candidate regulators")
save_current_plot(qc_figure_dir / "step48_top_targets_by_regulators.png")

# ============================================================
# Plot 4 — Promoter/distal support
# ============================================================

plot_df = promoter_distal_summary.sort_values("n_edges", ascending=True)

plt.figure(figsize=(8, 5))
plt.barh(plot_df["category"], plot_df["n_edges"])
plt.xlabel("Number of edges")
plt.ylabel("Support category")
plt.title("Promoter and distal support among candidate GRN edges")
save_current_plot(qc_figure_dir / "step48_promoter_distal_support.png")

# ============================================================
# Plot 5 — Supported cell type distribution
# ============================================================

plot_df = celltype_summary.sort_values("n_edges", ascending=True)

plt.figure(figsize=(8, 5))
plt.barh(plot_df["supported_cell_type"], plot_df["n_edges"])
plt.xlabel("Number of candidate edges")
plt.ylabel("Supported cell type")
plt.title("Cell-type support distribution of candidate GRN edges")
save_current_plot(qc_figure_dir / "step48_supported_cell_type_distribution.png")

# ============================================================
# Plot 6 — Known immune TFs by number of targets
# ============================================================

if immune_tf_summary.shape[0] > 0:
    plot_df = immune_tf_summary.head(30).sort_values("n_target_genes", ascending=True)

    plt.figure(figsize=(8, 8))
    plt.barh(plot_df["tf"], plot_df["n_target_genes"])
    plt.xlabel("Number of candidate target genes")
    plt.ylabel("Known immune TF")
    plt.title("Known immune TFs in candidate GRN")
    save_current_plot(qc_figure_dir / "step48_known_immune_tfs_by_targets.png")

# ============================================================
# Plot 7 — Small TF-target evidence heatmap
# ============================================================

top_heatmap_tfs = tf_summary.head(15)["tf"].tolist()
top_heatmap_targets = target_summary.head(20)["target_gene"].tolist()

heatmap_edges = collapsed[
    collapsed["tf"].isin(top_heatmap_tfs) &
    collapsed["target_gene"].isin(top_heatmap_targets)
].copy()

if heatmap_edges.shape[0] > 0:
    heatmap_matrix = (
        heatmap_edges
        .pivot_table(
            index="tf",
            columns="target_gene",
            values="max_edge_evidence_score",
            aggfunc="max",
            fill_value=0
        )
        .reindex(index=top_heatmap_tfs, columns=top_heatmap_targets)
        .fillna(0)
    )

    plt.figure(figsize=(12, 7))
    plt.imshow(heatmap_matrix.values, aspect="auto")
    plt.xticks(
        ticks=np.arange(heatmap_matrix.shape[1]),
        labels=heatmap_matrix.columns,
        rotation=90
    )
    plt.yticks(
        ticks=np.arange(heatmap_matrix.shape[0]),
        labels=heatmap_matrix.index
    )
    plt.colorbar(label="Max edge evidence score")
    plt.title("Candidate TF-target evidence heatmap")
    save_current_plot(qc_figure_dir / "step48_tf_target_evidence_heatmap.png")

# ============================================================
# Save annotated collapsed table with target categories
# ============================================================

collapsed.to_csv(qc_table_dir / "step48_candidate_edges_with_qc_annotations.csv", index=False)

print("\nSaved Step 48 QC tables to:")
print(qc_table_dir)

print("\nSaved Step 48 QC figures to:")
print(qc_figure_dir)

print("\nMain files to inspect next:")
print(qc_table_dir / "step48_grn_overall_summary.csv")
print(qc_table_dir / "step48_top_tfs_by_targets.csv")
print(qc_table_dir / "step48_known_immune_tf_summary.csv")
print(qc_table_dir / "step48_high_confidence_edges_q95_score_corr_ge_0_70.csv")
print(qc_table_dir / "step48_non_housekeeping_candidate_edges.csv")



"""Βήμα 49 — Thesis-ready high-confidence GRN export και visualization"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

grn_qc_dir = project_dir / "data/processed/grn_qc"

high_confidence_path = grn_qc_dir / "step48_high_confidence_edges_q95_score_corr_ge_0_70.csv"
all_qc_edges_path = grn_qc_dir / "step48_candidate_edges_with_qc_annotations.csv"

final_grn_dir = project_dir / "data/processed/final_grn"
final_fig_dir = project_dir / "figures/final_grn"

final_grn_dir.mkdir(parents=True, exist_ok=True)
final_fig_dir.mkdir(parents=True, exist_ok=True)

# ============================================================
# Output paths
# ============================================================

high_conf_all_out = final_grn_dir / "step49_high_confidence_grn_all_edges.csv"
high_conf_clean_out = final_grn_dir / "step49_high_confidence_grn_non_housekeeping_edges.csv"
immune_grn_out = final_grn_dir / "step49_high_confidence_known_immune_tf_grn.csv"

top_edges_out = final_grn_dir / "step49_top_100_high_confidence_edges.csv"
top_immune_edges_out = final_grn_dir / "step49_top_100_known_immune_tf_edges.csv"

tf_summary_out = final_grn_dir / "step49_high_confidence_tf_summary.csv"
target_summary_out = final_grn_dir / "step49_high_confidence_target_summary.csv"
immune_tf_summary_out = final_grn_dir / "step49_known_immune_tf_summary.csv"

celltype_edge_long_out = final_grn_dir / "step49_cell_type_edge_long_table.csv"
celltype_summary_out = final_grn_dir / "step49_cell_type_specific_grn_summary.csv"

cytoscape_edges_out = final_grn_dir / "step49_cytoscape_edges_high_confidence_immune_grn.csv"
cytoscape_nodes_out = final_grn_dir / "step49_cytoscape_nodes_high_confidence_immune_grn.csv"

thesis_summary_out = final_grn_dir / "step49_thesis_ready_grn_summary.csv"

# ============================================================
# Load high-confidence GRN
# ============================================================

high_conf = pd.read_csv(high_confidence_path)

print("Loaded high-confidence GRN:")
print(high_conf.shape)

print("Columns:")
print(list(high_conf.columns))

# ============================================================
# Basic cleaning
# ============================================================

high_conf["tf"] = high_conf["tf"].astype(str).str.strip().str.upper()
high_conf["target_gene"] = high_conf["target_gene"].astype(str).str.strip().str.upper()

numeric_cols = [
    "n_supporting_peak_motif_links",
    "n_unique_supporting_peaks",
    "n_unique_motifs",
    "n_promoter_peak_links",
    "n_distal_peak_links",
    "max_edge_evidence_score",
    "mean_edge_evidence_score",
    "max_target_gene_rna_atac_corr",
    "max_peak_mean_accessibility",
    "max_peak_detection_fraction",
    "max_tf_rna_in_peak_top_cell_type",
    "max_target_gene_rna_in_peak_top_cell_type",
    "top_evidence_score",
    "top_evidence_distance_to_gene",
    "top_evidence_peak_max_mean_accessibility",
    "top_evidence_peak_max_detection_fraction",
    "top_evidence_target_gene_rna_atac_corr",
    "top_evidence_tf_rna",
    "top_evidence_target_gene_rna",
]

for col in numeric_cols:
    if col in high_conf.columns:
        high_conf[col] = pd.to_numeric(high_conf[col], errors="coerce")

# ============================================================
# Target-gene category classification
# ============================================================

def classify_target_gene(gene):
    gene = str(gene).upper()

    if gene.startswith("RPL") or gene.startswith("RPS"):
        return "ribosomal"

    if gene.startswith("MT-"):
        return "mitochondrial"

    if gene.startswith("HSP"):
        return "heat_shock"

    if (
        gene.startswith("HIST") or
        gene.startswith("H1-") or
        gene.startswith("H2A") or
        gene.startswith("H2B") or
        gene.startswith("H3") or
        gene.startswith("H4")
    ):
        return "histone"

    if gene in {"MALAT1", "NEAT1", "B2M", "ACTB", "GAPDH", "TUBB", "TUBA1B"}:
        return "common_high_expression"

    return "other"

high_conf["target_gene_category"] = high_conf["target_gene"].map(classify_target_gene)

# ============================================================
# Clean non-housekeeping GRN
# ============================================================

clean_grn = high_conf[
    high_conf["target_gene_category"] == "other"
].copy()

clean_grn = clean_grn.sort_values(
    [
        "max_edge_evidence_score",
        "n_unique_supporting_peaks",
        "max_target_gene_rna_atac_corr",
    ],
    ascending=[False, False, False]
)

high_conf.to_csv(high_conf_all_out, index=False)
clean_grn.to_csv(high_conf_clean_out, index=False)

print()
print("High-confidence all edges:")
print(high_conf.shape)

print()
print("High-confidence non-housekeeping edges:")
print(clean_grn.shape)

# ============================================================
# Known immune TF list
# ============================================================

known_immune_tfs = [
    "FOS", "FOSB", "FOSL1", "FOSL2",
    "JUN", "JUNB", "JUND",
    "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9",
    "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
    "NFKB1", "NFKB2", "RELA", "RELB", "REL",
    "SPI1", "CEBPA", "CEBPB", "CEBPD",
    "KLF2", "KLF3", "KLF4", "KLF6",
    "IKZF1", "IKZF2", "IKZF3",
    "ETS1", "ELF1", "ELK1",
    "RUNX1", "RUNX2", "RUNX3",
    "TCF7", "LEF1",
    "GATA3", "TBX21", "EOMES",
    "PRDM1", "BCL6", "FOXP1", "FOXO1",
    "MEF2A", "MEF2C",
    "NR4A1", "NR4A2", "NR4A3",
    "AHR", "BATF", "MAF", "MAFB",
]

immune_grn = clean_grn[
    clean_grn["tf"].isin(known_immune_tfs)
].copy()

immune_grn = immune_grn.sort_values(
    [
        "max_edge_evidence_score",
        "n_unique_supporting_peaks",
        "max_target_gene_rna_atac_corr",
    ],
    ascending=[False, False, False]
)

immune_grn.to_csv(immune_grn_out, index=False)

print()
print("High-confidence known immune TF GRN:")
print(immune_grn.shape)

print()
print("Known immune TFs retained:")
print(sorted(immune_grn["tf"].unique()))

# ============================================================
# Top edge exports
# ============================================================

top_edges = clean_grn.head(100).copy()
top_immune_edges = immune_grn.head(100).copy()

top_edges.to_csv(top_edges_out, index=False)
top_immune_edges.to_csv(top_immune_edges_out, index=False)

print()
print("Top 30 high-confidence non-housekeeping edges:")
show_cols = [
    "tf",
    "target_gene",
    "max_edge_evidence_score",
    "n_unique_supporting_peaks",
    "n_unique_motifs",
    "n_promoter_peak_links",
    "n_distal_peak_links",
    "max_target_gene_rna_atac_corr",
    "supported_cell_types",
    "top_evidence_peak_id",
    "top_evidence_cell_type",
]

show_cols = [col for col in show_cols if col in clean_grn.columns]
print(clean_grn[show_cols].head(30).to_string(index=False))

print()
print("Top 30 known immune TF high-confidence edges:")
print(immune_grn[show_cols].head(30).to_string(index=False))

# ============================================================
# TF summary
# ============================================================

tf_summary = (
    clean_grn
    .groupby("tf")
    .agg(
        n_target_genes=("target_gene", "nunique"),
        n_edges=("target_gene", "count"),
        n_supporting_peak_motif_links=("n_supporting_peak_motif_links", "sum"),
        n_unique_supporting_peaks_total=("n_unique_supporting_peaks", "sum"),
        n_unique_motifs_total=("n_unique_motifs", "sum"),
        n_promoter_peak_links=("n_promoter_peak_links", "sum"),
        n_distal_peak_links=("n_distal_peak_links", "sum"),
        max_edge_evidence_score=("max_edge_evidence_score", "max"),
        mean_edge_evidence_score=("max_edge_evidence_score", "mean"),
        median_edge_evidence_score=("max_edge_evidence_score", "median"),
        mean_target_gene_rna_atac_corr=("max_target_gene_rna_atac_corr", "mean"),
        max_target_gene_rna_atac_corr=("max_target_gene_rna_atac_corr", "max"),
    )
    .reset_index()
    .sort_values(
        ["n_target_genes", "max_edge_evidence_score"],
        ascending=[False, False]
    )
)

tf_summary.to_csv(tf_summary_out, index=False)

print()
print("Top 30 TFs in clean high-confidence GRN:")
print(tf_summary.head(30).to_string(index=False))

# ============================================================
# Target summary
# ============================================================

target_summary = (
    clean_grn
    .groupby("target_gene")
    .agg(
        n_regulator_tfs=("tf", "nunique"),
        n_edges=("tf", "count"),
        n_supporting_peak_motif_links=("n_supporting_peak_motif_links", "sum"),
        n_unique_supporting_peaks_total=("n_unique_supporting_peaks", "sum"),
        n_promoter_peak_links=("n_promoter_peak_links", "sum"),
        n_distal_peak_links=("n_distal_peak_links", "sum"),
        max_edge_evidence_score=("max_edge_evidence_score", "max"),
        mean_edge_evidence_score=("max_edge_evidence_score", "mean"),
        target_gene_rna_atac_corr=("max_target_gene_rna_atac_corr", "max"),
    )
    .reset_index()
    .sort_values(
        ["n_regulator_tfs", "max_edge_evidence_score"],
        ascending=[False, False]
    )
)

target_summary.to_csv(target_summary_out, index=False)

print()
print("Top 30 target genes by number of regulator TFs:")
print(target_summary.head(30).to_string(index=False))

# ============================================================
# Known immune TF summary
# ============================================================

immune_tf_summary = (
    immune_grn
    .groupby("tf")
    .agg(
        n_target_genes=("target_gene", "nunique"),
        n_edges=("target_gene", "count"),
        n_supporting_peak_motif_links=("n_supporting_peak_motif_links", "sum"),
        n_unique_supporting_peaks_total=("n_unique_supporting_peaks", "sum"),
        n_unique_motifs_total=("n_unique_motifs", "sum"),
        n_promoter_peak_links=("n_promoter_peak_links", "sum"),
        n_distal_peak_links=("n_distal_peak_links", "sum"),
        max_edge_evidence_score=("max_edge_evidence_score", "max"),
        mean_edge_evidence_score=("max_edge_evidence_score", "mean"),
        mean_target_gene_rna_atac_corr=("max_target_gene_rna_atac_corr", "mean"),
        max_target_gene_rna_atac_corr=("max_target_gene_rna_atac_corr", "max"),
    )
    .reset_index()
    .sort_values(
        ["n_target_genes", "max_edge_evidence_score"],
        ascending=[False, False]
    )
)

immune_tf_summary.to_csv(immune_tf_summary_out, index=False)

print()
print("Known immune TF summary:")
print(immune_tf_summary.to_string(index=False))

# ============================================================
# Cell-type-specific long table
# ============================================================

celltype_long = clean_grn.copy()

if "supported_cell_types" not in celltype_long.columns:
    if "top_evidence_cell_type" in celltype_long.columns:
        celltype_long["supported_cell_types"] = celltype_long["top_evidence_cell_type"]
    else:
        celltype_long["supported_cell_types"] = "unknown"

celltype_long["supported_cell_types"] = (
    celltype_long["supported_cell_types"]
    .fillna("unknown")
    .astype(str)
)

celltype_long["supported_cell_type"] = (
    celltype_long["supported_cell_types"]
    .str.split(";")
)

celltype_long = celltype_long.explode("supported_cell_type")
celltype_long["supported_cell_type"] = (
    celltype_long["supported_cell_type"]
    .astype(str)
    .str.strip()
)

celltype_long.to_csv(celltype_edge_long_out, index=False)

# ============================================================
# Cell-type-specific summary
# ============================================================

celltype_summary = (
    celltype_long
    .groupby("supported_cell_type")
    .agg(
        n_edges=("target_gene", "count"),
        n_tfs=("tf", "nunique"),
        n_target_genes=("target_gene", "nunique"),
        max_edge_evidence_score=("max_edge_evidence_score", "max"),
        mean_edge_evidence_score=("max_edge_evidence_score", "mean"),
        n_promoter_peak_links=("n_promoter_peak_links", "sum"),
        n_distal_peak_links=("n_distal_peak_links", "sum"),
    )
    .reset_index()
    .sort_values("n_edges", ascending=False)
)

celltype_summary.to_csv(celltype_summary_out, index=False)

print()
print("Cell-type-specific GRN summary:")
print(celltype_summary.to_string(index=False))

# ============================================================
# Export per-cell-type edge tables
# ============================================================

def safe_name(x):
    return (
        str(x)
        .replace("+", "plus")
        .replace(" ", "_")
        .replace(".", "_")
        .replace("-", "_")
        .replace("/", "_")
    )

for cell_type in sorted(celltype_long["supported_cell_type"].dropna().unique()):
    sub = celltype_long[
        celltype_long["supported_cell_type"] == cell_type
    ].copy()

    sub = sub.sort_values(
        ["max_edge_evidence_score", "n_unique_supporting_peaks"],
        ascending=[False, False]
    )

    out_path = final_grn_dir / f"step49_cell_type_specific_edges_{safe_name(cell_type)}.csv"
    sub.to_csv(out_path, index=False)

# ============================================================
# Promoter/distal summary
# ============================================================

promoter_distal_summary = pd.DataFrame([
    {
        "category": "edges_with_promoter_support",
        "n_edges": int((clean_grn["n_promoter_peak_links"] > 0).sum()),
    },
    {
        "category": "edges_with_distal_support",
        "n_edges": int((clean_grn["n_distal_peak_links"] > 0).sum()),
    },
    {
        "category": "edges_with_both_promoter_and_distal_support",
        "n_edges": int(
            (
                (clean_grn["n_promoter_peak_links"] > 0) &
                (clean_grn["n_distal_peak_links"] > 0)
            ).sum()
        ),
    },
    {
        "category": "edges_with_only_promoter_support",
        "n_edges": int(
            (
                (clean_grn["n_promoter_peak_links"] > 0) &
                (clean_grn["n_distal_peak_links"] == 0)
            ).sum()
        ),
    },
    {
        "category": "edges_with_only_distal_support",
        "n_edges": int(
            (
                (clean_grn["n_promoter_peak_links"] == 0) &
                (clean_grn["n_distal_peak_links"] > 0)
            ).sum()
        ),
    },
])

promoter_distal_summary.to_csv(
    final_grn_dir / "step49_promoter_distal_summary_clean_grn.csv",
    index=False
)

print()
print("Promoter/distal summary for clean high-confidence GRN:")
print(promoter_distal_summary.to_string(index=False))

# ============================================================
# Cytoscape-ready immune GRN export
# ============================================================

cytoscape_edges = immune_grn.copy()

cytoscape_edges = cytoscape_edges.rename(
    columns={
        "tf": "source",
        "target_gene": "target",
        "max_edge_evidence_score": "score",
    }
)

cytoscape_edges["interaction"] = "candidate_regulates"
cytoscape_edges["method"] = "multiome_informed_candidate_grn"
cytoscape_edges["source_type"] = "TF"
cytoscape_edges["target_type"] = "gene"

cytoscape_keep_cols = [
    "source",
    "target",
    "interaction",
    "score",
    "method",
    "source_type",
    "target_type",
    "n_unique_supporting_peaks",
    "n_unique_motifs",
    "n_promoter_peak_links",
    "n_distal_peak_links",
    "max_target_gene_rna_atac_corr",
    "supported_cell_types",
    "top_evidence_peak_id",
    "top_evidence_cell_type",
]

cytoscape_keep_cols = [
    col for col in cytoscape_keep_cols
    if col in cytoscape_edges.columns
]

cytoscape_edges = cytoscape_edges[cytoscape_keep_cols].copy()
cytoscape_edges.to_csv(cytoscape_edges_out, index=False)

tf_nodes = pd.DataFrame({
    "id": sorted(immune_grn["tf"].unique()),
    "node_type": "TF",
})

target_nodes = pd.DataFrame({
    "id": sorted(immune_grn["target_gene"].unique()),
    "node_type": "target_gene",
})

cytoscape_nodes = pd.concat(
    [tf_nodes, target_nodes],
    axis=0,
    ignore_index=True
).drop_duplicates("id")

cytoscape_nodes.to_csv(cytoscape_nodes_out, index=False)

print()
print("Cytoscape edges:")
print(cytoscape_edges.shape)

print()
print("Cytoscape nodes:")
print(cytoscape_nodes.shape)

# ============================================================
# Thesis-ready summary table
# ============================================================

thesis_summary = pd.DataFrame([
    {
        "metric": "high_confidence_edges_all",
        "value": high_conf.shape[0],
    },
    {
        "metric": "high_confidence_unique_tfs_all",
        "value": high_conf["tf"].nunique(),
    },
    {
        "metric": "high_confidence_unique_targets_all",
        "value": high_conf["target_gene"].nunique(),
    },
    {
        "metric": "clean_non_housekeeping_edges",
        "value": clean_grn.shape[0],
    },
    {
        "metric": "clean_non_housekeeping_unique_tfs",
        "value": clean_grn["tf"].nunique(),
    },
    {
        "metric": "clean_non_housekeeping_unique_targets",
        "value": clean_grn["target_gene"].nunique(),
    },
    {
        "metric": "known_immune_tf_edges",
        "value": immune_grn.shape[0],
    },
    {
        "metric": "known_immune_tfs",
        "value": immune_grn["tf"].nunique(),
    },
    {
        "metric": "known_immune_tf_targets",
        "value": immune_grn["target_gene"].nunique(),
    },
    {
        "metric": "clean_edges_with_at_least_2_supporting_peaks",
        "value": int((clean_grn["n_unique_supporting_peaks"] >= 2).sum()),
    },
    {
        "metric": "clean_edges_with_promoter_support",
        "value": int((clean_grn["n_promoter_peak_links"] > 0).sum()),
    },
    {
        "metric": "clean_edges_with_distal_support",
        "value": int((clean_grn["n_distal_peak_links"] > 0).sum()),
    },
    {
        "metric": "median_clean_edge_score",
        "value": clean_grn["max_edge_evidence_score"].median(),
    },
    {
        "metric": "mean_clean_edge_score",
        "value": clean_grn["max_edge_evidence_score"].mean(),
    },
    {
        "metric": "median_clean_target_rna_atac_corr",
        "value": clean_grn["max_target_gene_rna_atac_corr"].median(),
    },
    {
        "metric": "mean_clean_target_rna_atac_corr",
        "value": clean_grn["max_target_gene_rna_atac_corr"].mean(),
    },
])

thesis_summary.to_csv(thesis_summary_out, index=False)

print()
print("Thesis-ready GRN summary:")
print(thesis_summary.to_string(index=False))

# ============================================================
# Plot helper
# ============================================================

def save_current_plot(output_path):
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()

# ============================================================
# Plot 1 — Edge evidence score distribution
# ============================================================

plt.figure(figsize=(7, 5))
plt.hist(clean_grn["max_edge_evidence_score"].dropna(), bins=40)
plt.xlabel("Max edge evidence score")
plt.ylabel("Number of edges")
plt.title("High-confidence non-housekeeping GRN edge scores")
save_current_plot(final_fig_dir / "step49_clean_grn_edge_score_distribution.png")

# ============================================================
# Plot 2 — Top TFs by number of target genes
# ============================================================

top_n = 25
plot_df = tf_summary.head(top_n).sort_values("n_target_genes", ascending=True)

plt.figure(figsize=(8, 7))
plt.barh(plot_df["tf"], plot_df["n_target_genes"])
plt.xlabel("Number of target genes")
plt.ylabel("TF")
plt.title("Top TFs in high-confidence non-housekeeping GRN")
save_current_plot(final_fig_dir / "step49_top_tfs_clean_grn.png")

# ============================================================
# Plot 3 — Known immune TFs by number of targets
# ============================================================

if immune_tf_summary.shape[0] > 0:
    plot_df = immune_tf_summary.head(top_n).sort_values("n_target_genes", ascending=True)

    plt.figure(figsize=(8, 7))
    plt.barh(plot_df["tf"], plot_df["n_target_genes"])
    plt.xlabel("Number of target genes")
    plt.ylabel("Known immune TF")
    plt.title("Known immune TFs in high-confidence GRN")
    save_current_plot(final_fig_dir / "step49_known_immune_tfs.png")

# ============================================================
# Plot 4 — Cell-type edge counts
# ============================================================

plot_df = celltype_summary.sort_values("n_edges", ascending=True)

plt.figure(figsize=(8, 5))
plt.barh(plot_df["supported_cell_type"], plot_df["n_edges"])
plt.xlabel("Number of high-confidence edges")
plt.ylabel("Supported cell type")
plt.title("Cell-type support in high-confidence GRN")
save_current_plot(final_fig_dir / "step49_cell_type_edge_counts.png")

# ============================================================
# Plot 5 — Promoter/distal support
# ============================================================

plot_df = promoter_distal_summary.sort_values("n_edges", ascending=True)

plt.figure(figsize=(8, 5))
plt.barh(plot_df["category"], plot_df["n_edges"])
plt.xlabel("Number of edges")
plt.ylabel("Support category")
plt.title("Promoter and distal support in clean high-confidence GRN")
save_current_plot(final_fig_dir / "step49_promoter_distal_support.png")

# ============================================================
# Plot 6 — Top target genes by number of candidate regulators
# ============================================================

plot_df = target_summary.head(top_n).sort_values("n_regulator_tfs", ascending=True)

plt.figure(figsize=(8, 7))
plt.barh(plot_df["target_gene"], plot_df["n_regulator_tfs"])
plt.xlabel("Number of regulator TFs")
plt.ylabel("Target gene")
plt.title("Top target genes by number of candidate regulators")
save_current_plot(final_fig_dir / "step49_top_targets_by_regulators.png")

# ============================================================
# Optional network plot using networkx
# ============================================================

try:
    import networkx as nx

    network_edges = immune_grn.head(60).copy()

    if network_edges.shape[0] > 0:
        G = nx.DiGraph()

        for _, row in network_edges.iterrows():
            tf = row["tf"]
            target = row["target_gene"]
            score = row["max_edge_evidence_score"]

            G.add_node(tf, node_type="TF")
            G.add_node(target, node_type="target_gene")
            G.add_edge(tf, target, weight=score)

        pos = nx.spring_layout(G, seed=42, k=0.8)

        tf_nodes = [
            node for node, data in G.nodes(data=True)
            if data.get("node_type") == "TF"
        ]

        target_nodes = [
            node for node, data in G.nodes(data=True)
            if data.get("node_type") == "target_gene"
        ]

        edge_weights = [
            G[u][v].get("weight", 0.01)
            for u, v in G.edges()
        ]

        max_weight = max(edge_weights) if len(edge_weights) > 0 else 1.0

        edge_widths = [
            0.5 + 3.0 * (w / max_weight)
            for w in edge_weights
        ]

        plt.figure(figsize=(13, 10))

        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=tf_nodes,
            node_shape="s",
            node_size=700,
        )

        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=target_nodes,
            node_shape="o",
            node_size=350,
        )

        nx.draw_networkx_edges(
            G,
            pos,
            arrows=True,
            width=edge_widths,
            alpha=0.5,
        )

        nx.draw_networkx_labels(
            G,
            pos,
            font_size=7,
        )

        plt.title("Top known immune TF high-confidence candidate GRN")
        plt.axis("off")
        save_current_plot(final_fig_dir / "step49_top_known_immune_tf_network.png")

        print()
        print("Saved network plot using networkx.")

except ImportError:
    print()
    print("networkx is not installed. Skipping network plot.")
    print("You can install it with: pip install networkx")

# ============================================================
# Final messages
# ============================================================

print()
print("Saved Step 49 final GRN tables to:")
print(final_grn_dir)

print()
print("Saved Step 49 final GRN figures to:")
print(final_fig_dir)

print()
print("Main thesis-ready files:")
print(high_conf_clean_out)
print(immune_grn_out)
print(tf_summary_out)
print(immune_tf_summary_out)
print(celltype_summary_out)
print(cytoscape_edges_out)
print(cytoscape_nodes_out)
print(thesis_summary_out)










"""Βήμα 47R — Per-cell-type candidate GRN construction για τα 5 PBMC cell types"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

peak_to_gene_path = project_dir / "data/processed/pbmc_atac_peak_to_gene_annotation_5types.csv"
peak_motif_path = project_dir / "data/processed/pbmc_atac_peak_motif_mapping_5types.csv"

rna_pseudobulk_path = project_dir / "data/processed/pbmc_rna_pseudobulk_5types_common_genes.csv"
atac_gene_activity_pseudobulk_path = project_dir / "data/processed/pbmc_atac_gene_activity_pseudobulk_5types_common_genes.csv"

out_dir = project_dir / "data/processed/grn_per_cell_type"
fig_dir = project_dir / "figures/grn_per_cell_type"

out_dir.mkdir(parents=True, exist_ok=True)
fig_dir.mkdir(parents=True, exist_ok=True)

# ============================================================
# Output paths
# ============================================================

all_peak_level_out = out_dir / "step47R_all_cell_type_peak_level_candidate_edges.csv"
all_collapsed_out = out_dir / "step47R_all_cell_type_tf_target_edges_collapsed.csv"
summary_out = out_dir / "step47R_per_cell_type_grn_summary.csv"
tf_degree_matrix_out = out_dir / "step47R_tf_degree_matrix_by_cell_type.csv"
tf_degree_similarity_out = out_dir / "step47R_tf_degree_similarity_by_cell_type.csv"

# ============================================================
# Parameters
# ============================================================

MIN_PEAK_DETECTION_FRACTION = 0.01
MIN_PEAK_MEAN_ACCESSIBILITY = 0.005

MIN_TF_RNA_EXPRESSION = 0.05
MIN_TARGET_RNA_EXPRESSION = 0.05

EXCLUDE_SELF_LOOPS = True

# Keep this False for the rewritten per-cell-type biological logic.
# The cross-cell-type RNA-ATAC correlation is useful as QC, but not as a core GRN filter.
USE_GLOBAL_RNA_ATAC_CORR_AS_FILTER = False
MIN_GLOBAL_RNA_ATAC_CORR = 0.40

# Optional: keep ATAC gene activity as annotation, not required.
USE_TARGET_ATAC_GENE_ACTIVITY_IN_SCORE = False

# ============================================================
# Cell type definitions
# ============================================================

cell_type_order = [
    "CD14+ Monocytes",
    "CD4.Memory",
    "CD4.Naive",
    "CD8.effector",
    "CD8.Naive",
]

cell_type_to_safe = {
    "CD14+ Monocytes": "CD14plus_Monocytes",
    "CD4.Memory": "CD4_Memory",
    "CD4.Naive": "CD4_Naive",
    "CD8.effector": "CD8_effector",
    "CD8.Naive": "CD8_Naive",
}

# ============================================================
# Load input tables
# ============================================================

peak_to_gene = pd.read_csv(peak_to_gene_path)
peak_motif = pd.read_csv(peak_motif_path)

rna_pb = pd.read_csv(rna_pseudobulk_path, index_col=0)
rna_pb.index = rna_pb.index.astype(str).str.strip()
rna_pb.columns = rna_pb.columns.astype(str).str.strip().str.upper()

if atac_gene_activity_pseudobulk_path.exists():
    atac_ga_pb = pd.read_csv(atac_gene_activity_pseudobulk_path, index_col=0)
    atac_ga_pb.index = atac_ga_pb.index.astype(str).str.strip()
    atac_ga_pb.columns = atac_ga_pb.columns.astype(str).str.strip().str.upper()
else:
    atac_ga_pb = None

print("Peak-to-gene table:")
print(peak_to_gene.shape)

print()
print("Peak-motif-TF table:")
print(peak_motif.shape)

print()
print("RNA pseudo-bulk:")
print(rna_pb.shape)

if atac_ga_pb is not None:
    print()
    print("ATAC gene activity pseudo-bulk:")
    print(atac_ga_pb.shape)

# ============================================================
# Basic cleaning
# ============================================================

peak_to_gene["peak_id"] = peak_to_gene["peak_id"].astype(str).str.strip()
peak_to_gene["target_gene_original"] = peak_to_gene["target_gene"].astype(str).str.strip()
peak_to_gene["target_gene"] = peak_to_gene["target_gene_original"].str.upper()

peak_motif["peak_id"] = peak_motif["peak_id"].astype(str).str.strip()
peak_motif["tf"] = peak_motif["tf"].astype(str).str.strip().str.upper()

if "motif_id" not in peak_motif.columns:
    peak_motif["motif_id"] = "unknown"

if "motif_name" not in peak_motif.columns:
    peak_motif["motif_name"] = peak_motif["tf"]

rna_gene_set = set(rna_pb.columns)

# ============================================================
# Helper functions
# ============================================================

def get_expression_dict(matrix_df, cell_type):
    if matrix_df is None:
        return {}

    if cell_type not in matrix_df.index:
        raise ValueError(f"Missing cell type in matrix: {cell_type}")

    return matrix_df.loc[cell_type].to_dict()


def add_gene_expression_columns(df, cell_type):
    df = df.copy()

    rna_expr = get_expression_dict(rna_pb, cell_type)

    df["tf_in_rna"] = df["tf"].isin(rna_gene_set)
    df["target_gene_in_rna"] = df["target_gene"].isin(rna_gene_set)

    df["tf_rna_expression"] = df["tf"].map(rna_expr)
    df["target_gene_rna_expression"] = df["target_gene"].map(rna_expr)

    if atac_ga_pb is not None:
        atac_ga_expr = get_expression_dict(atac_ga_pb, cell_type)
        df["target_gene_atac_gene_activity"] = df["target_gene"].map(atac_ga_expr)
    else:
        df["target_gene_atac_gene_activity"] = np.nan

    return df


def score_edges(df):
    df = df.copy()

    accessibility_score = np.log1p(df["cell_type_peak_mean_accessibility"].fillna(0))
    detection_score = df["cell_type_peak_detection_fraction"].fillna(0)

    tf_expr_score = np.log1p(df["tf_rna_expression"].fillna(0))
    target_expr_score = np.log1p(df["target_gene_rna_expression"].fillna(0))

    score = (
        accessibility_score *
        detection_score *
        tf_expr_score *
        target_expr_score
    )

    if USE_TARGET_ATAC_GENE_ACTIVITY_IN_SCORE:
        target_atac_score = np.log1p(df["target_gene_atac_gene_activity"].fillna(0))
        score = score * target_atac_score

    df["cell_type_edge_evidence_score"] = score

    df["is_promoter_peak"] = (
        df["peak_annotation_type"]
        .astype(str)
        .str.lower()
        .str.contains("promoter", na=False)
    )

    df["is_distal_peak"] = (
        df["peak_annotation_type"]
        .astype(str)
        .str.lower()
        .str.contains("distal", na=False)
    )

    df["is_self_loop"] = (
        df["tf"].astype(str).str.upper()
        ==
        df["target_gene"].astype(str).str.upper()
    )

    return df


def safe_cell_type_name(cell_type):
    return (
        str(cell_type)
        .replace("+", "plus")
        .replace(" ", "_")
        .replace(".", "_")
        .replace("-", "_")
        .replace("/", "_")
    )


def collapse_cell_type_edges(edge_df, cell_type):
    if edge_df.empty:
        return pd.DataFrame()

    edge_df = edge_df.copy()

    top_evidence = (
        edge_df
        .sort_values("cell_type_edge_evidence_score", ascending=False)
        .drop_duplicates(["cell_type", "tf", "target_gene"])
        [
            [
                "cell_type",
                "tf",
                "target_gene",
                "peak_id",
                "motif_id",
                "motif_name",
                "peak_annotation_type",
                "distance_to_gene",
                "cell_type_peak_mean_accessibility",
                "cell_type_peak_detection_fraction",
                "tf_rna_expression",
                "target_gene_rna_expression",
                "target_gene_atac_gene_activity",
                "cell_type_edge_evidence_score",
            ]
        ]
        .rename(columns={
            "peak_id": "top_evidence_peak_id",
            "motif_id": "top_evidence_motif_id",
            "motif_name": "top_evidence_motif_name",
            "peak_annotation_type": "top_evidence_peak_annotation_type",
            "distance_to_gene": "top_evidence_distance_to_gene",
            "cell_type_peak_mean_accessibility": "top_evidence_peak_mean_accessibility",
            "cell_type_peak_detection_fraction": "top_evidence_peak_detection_fraction",
            "tf_rna_expression": "top_evidence_tf_rna_expression",
            "target_gene_rna_expression": "top_evidence_target_gene_rna_expression",
            "target_gene_atac_gene_activity": "top_evidence_target_gene_atac_gene_activity",
            "cell_type_edge_evidence_score": "top_evidence_score",
        })
    )

    agg_dict = {
        "n_supporting_peak_motif_links": ("peak_id", "count"),
        "n_unique_supporting_peaks": ("peak_id", "nunique"),
        "n_unique_motifs": ("motif_id", "nunique"),
        "n_promoter_peak_links": ("is_promoter_peak", "sum"),
        "n_distal_peak_links": ("is_distal_peak", "sum"),
        "max_cell_type_edge_evidence_score": ("cell_type_edge_evidence_score", "max"),
        "mean_cell_type_edge_evidence_score": ("cell_type_edge_evidence_score", "mean"),
        "max_peak_mean_accessibility": ("cell_type_peak_mean_accessibility", "max"),
        "max_peak_detection_fraction": ("cell_type_peak_detection_fraction", "max"),
        "max_tf_rna_expression": ("tf_rna_expression", "max"),
        "max_target_gene_rna_expression": ("target_gene_rna_expression", "max"),
        "max_target_gene_atac_gene_activity": ("target_gene_atac_gene_activity", "max"),
        "min_abs_distance_to_gene": ("abs_distance_to_gene", "min"),
    }

    if "target_gene_rna_atac_pseudobulk_correlation" in edge_df.columns:
        agg_dict["global_target_gene_rna_atac_corr"] = (
            "target_gene_rna_atac_pseudobulk_correlation",
            "max",
        )

    collapsed = (
        edge_df
        .groupby(["cell_type", "tf", "target_gene"])
        .agg(**agg_dict)
        .reset_index()
    )

    collapsed = collapsed.merge(
        top_evidence,
        on=["cell_type", "tf", "target_gene"],
        how="left"
    )

    collapsed = collapsed.sort_values(
        [
            "max_cell_type_edge_evidence_score",
            "n_unique_supporting_peaks",
            "max_peak_detection_fraction",
        ],
        ascending=[False, False, False]
    )

    return collapsed

# ============================================================
# Per-cell-type GRN construction
# ============================================================

all_peak_level_edges = []
all_collapsed_edges = []
summary_records = []

for cell_type in cell_type_order:
    print()
    print("============================================================")
    print("Constructing GRN for cell type:")
    print(cell_type)
    print("============================================================")

    safe_ct = cell_type_to_safe[cell_type]

    mean_col = f"mean_accessibility_{safe_ct}"
    detection_col = f"detection_fraction_{safe_ct}"

    if mean_col not in peak_to_gene.columns:
        raise ValueError(f"Missing column: {mean_col}")

    if detection_col not in peak_to_gene.columns:
        raise ValueError(f"Missing column: {detection_col}")

    # ------------------------------------------------------------
    # Cell-type-specific peak-to-gene filtering
    # ------------------------------------------------------------

    ct_peak_to_gene = peak_to_gene.copy()

    ct_peak_to_gene["cell_type"] = cell_type

    ct_peak_to_gene["cell_type_peak_mean_accessibility"] = pd.to_numeric(
        ct_peak_to_gene[mean_col],
        errors="coerce"
    )

    ct_peak_to_gene["cell_type_peak_detection_fraction"] = pd.to_numeric(
        ct_peak_to_gene[detection_col],
        errors="coerce"
    )

    if "distance_to_gene" in ct_peak_to_gene.columns:
        ct_peak_to_gene["distance_to_gene"] = pd.to_numeric(
            ct_peak_to_gene["distance_to_gene"],
            errors="coerce"
        )
        ct_peak_to_gene["abs_distance_to_gene"] = ct_peak_to_gene["distance_to_gene"].abs()
    else:
        ct_peak_to_gene["distance_to_gene"] = np.nan
        ct_peak_to_gene["abs_distance_to_gene"] = np.nan

    if "target_gene_rna_atac_pseudobulk_correlation" in ct_peak_to_gene.columns:
        ct_peak_to_gene["target_gene_rna_atac_pseudobulk_correlation"] = pd.to_numeric(
            ct_peak_to_gene["target_gene_rna_atac_pseudobulk_correlation"],
            errors="coerce"
        )

    ct_peak_to_gene_prefiltered = ct_peak_to_gene[
        (ct_peak_to_gene["cell_type_peak_mean_accessibility"] >= MIN_PEAK_MEAN_ACCESSIBILITY) &
        (ct_peak_to_gene["cell_type_peak_detection_fraction"] >= MIN_PEAK_DETECTION_FRACTION)
    ].copy()

    if USE_GLOBAL_RNA_ATAC_CORR_AS_FILTER:
        ct_peak_to_gene_prefiltered = ct_peak_to_gene_prefiltered[
            ct_peak_to_gene_prefiltered["target_gene_rna_atac_pseudobulk_correlation"] >= MIN_GLOBAL_RNA_ATAC_CORR
        ].copy()

    print("Peak-to-gene rows before cell-type filtering:")
    print(ct_peak_to_gene.shape[0])

    print("Peak-to-gene rows after cell-type peak accessibility filtering:")
    print(ct_peak_to_gene_prefiltered.shape[0])

    # ------------------------------------------------------------
    # Merge with peak-motif-TF evidence
    # ------------------------------------------------------------

    ct_edges = ct_peak_to_gene_prefiltered.merge(
        peak_motif,
        on="peak_id",
        how="inner"
    )

    print("Merged peak-level rows before expression filtering:")
    print(ct_edges.shape[0])

    if ct_edges.empty:
        summary_records.append({
            "cell_type": cell_type,
            "peak_to_gene_rows_input": ct_peak_to_gene.shape[0],
            "peak_to_gene_rows_after_cell_type_peak_filtering": ct_peak_to_gene_prefiltered.shape[0],
            "merged_peak_level_rows_before_expression_filtering": 0,
            "filtered_peak_level_edges": 0,
            "collapsed_tf_target_edges": 0,
            "unique_tfs": 0,
            "unique_target_genes": 0,
        })
        continue

    # ------------------------------------------------------------
    # Add cell-type-specific RNA expression evidence
    # ------------------------------------------------------------

    ct_edges = add_gene_expression_columns(ct_edges, cell_type)

    ct_edges = ct_edges[
        (ct_edges["tf_in_rna"]) &
        (ct_edges["target_gene_in_rna"]) &
        (ct_edges["tf_rna_expression"] >= MIN_TF_RNA_EXPRESSION) &
        (ct_edges["target_gene_rna_expression"] >= MIN_TARGET_RNA_EXPRESSION)
    ].copy()

    if EXCLUDE_SELF_LOOPS:
        ct_edges["is_self_loop"] = (
            ct_edges["tf"].astype(str).str.upper()
            ==
            ct_edges["target_gene"].astype(str).str.upper()
        )
        ct_edges = ct_edges[~ct_edges["is_self_loop"]].copy()

    # ------------------------------------------------------------
    # Score cell-type-specific candidate edges
    # ------------------------------------------------------------

    ct_edges = score_edges(ct_edges)

    ct_edges = ct_edges[
        ct_edges["cell_type_edge_evidence_score"] > 0
    ].copy()

    ct_edges = ct_edges.sort_values(
        "cell_type_edge_evidence_score",
        ascending=False
    )

    print("Filtered peak-level edges after TF/target expression filtering:")
    print(ct_edges.shape[0])

    # ------------------------------------------------------------
    # Save per-cell-type peak-level table
    # ------------------------------------------------------------

    ct_safe_file = safe_cell_type_name(cell_type)

    ct_peak_level_out = out_dir / f"step47R_peak_level_edges_{ct_safe_file}.csv"
    ct_edges.to_csv(ct_peak_level_out, index=False)

    # ------------------------------------------------------------
    # Collapse to TF-target edges for this cell type
    # ------------------------------------------------------------

    ct_collapsed = collapse_cell_type_edges(ct_edges, cell_type)

    ct_collapsed_out = out_dir / f"step47R_tf_target_edges_collapsed_{ct_safe_file}.csv"
    ct_collapsed.to_csv(ct_collapsed_out, index=False)

    print("Collapsed TF-target edges:")
    print(ct_collapsed.shape[0])

    print("Unique TFs:")
    print(ct_collapsed["tf"].nunique() if not ct_collapsed.empty else 0)

    print("Unique target genes:")
    print(ct_collapsed["target_gene"].nunique() if not ct_collapsed.empty else 0)

    print("Top 10 TF-target edges:")
    if not ct_collapsed.empty:
        print(
            ct_collapsed[
                [
                    "cell_type",
                    "tf",
                    "target_gene",
                    "max_cell_type_edge_evidence_score",
                    "n_unique_supporting_peaks",
                    "n_unique_motifs",
                    "max_peak_detection_fraction",
                    "max_tf_rna_expression",
                    "max_target_gene_rna_expression",
                    "top_evidence_peak_id",
                ]
            ].head(10).to_string(index=False)
        )

    # ------------------------------------------------------------
    # Store for all-cell-type outputs
    # ------------------------------------------------------------

    all_peak_level_edges.append(ct_edges)
    all_collapsed_edges.append(ct_collapsed)

    summary_records.append({
        "cell_type": cell_type,
        "peak_to_gene_rows_input": ct_peak_to_gene.shape[0],
        "peak_to_gene_rows_after_cell_type_peak_filtering": ct_peak_to_gene_prefiltered.shape[0],
        "merged_peak_level_rows_before_expression_filtering": ct_edges.shape[0],
        "filtered_peak_level_edges": ct_edges.shape[0],
        "collapsed_tf_target_edges": ct_collapsed.shape[0],
        "unique_tfs": ct_collapsed["tf"].nunique() if not ct_collapsed.empty else 0,
        "unique_target_genes": ct_collapsed["target_gene"].nunique() if not ct_collapsed.empty else 0,
        "min_peak_detection_fraction": MIN_PEAK_DETECTION_FRACTION,
        "min_peak_mean_accessibility": MIN_PEAK_MEAN_ACCESSIBILITY,
        "min_tf_rna_expression": MIN_TF_RNA_EXPRESSION,
        "min_target_rna_expression": MIN_TARGET_RNA_EXPRESSION,
        "used_global_rna_atac_corr_as_filter": USE_GLOBAL_RNA_ATAC_CORR_AS_FILTER,
        "exclude_self_loops": EXCLUDE_SELF_LOOPS,
    })

# ============================================================
# Save combined outputs
# ============================================================

if len(all_peak_level_edges) > 0:
    all_peak_level = pd.concat(all_peak_level_edges, axis=0, ignore_index=True)
else:
    all_peak_level = pd.DataFrame()

if len(all_collapsed_edges) > 0:
    all_collapsed = pd.concat(all_collapsed_edges, axis=0, ignore_index=True)
else:
    all_collapsed = pd.DataFrame()

all_peak_level.to_csv(all_peak_level_out, index=False)
all_collapsed.to_csv(all_collapsed_out, index=False)

summary = pd.DataFrame(summary_records)
summary.to_csv(summary_out, index=False)

print()
print("============================================================")
print("Saved combined per-cell-type peak-level edges to:")
print(all_peak_level_out)

print()
print("Saved combined per-cell-type collapsed TF-target edges to:")
print(all_collapsed_out)

print()
print("Saved per-cell-type GRN summary to:")
print(summary_out)

print()
print("Per-cell-type GRN summary:")
print(summary.to_string(index=False))

# ============================================================
# TF degree matrix by cell type
# ============================================================

if not all_collapsed.empty:
    tf_degree = (
        all_collapsed
        .groupby(["cell_type", "tf"])
        .agg(
            n_target_genes=("target_gene", "nunique"),
            n_edges=("target_gene", "count"),
            max_edge_score=("max_cell_type_edge_evidence_score", "max"),
        )
        .reset_index()
    )

    tf_degree_matrix = (
        tf_degree
        .pivot_table(
            index="tf",
            columns="cell_type",
            values="n_target_genes",
            fill_value=0
        )
        .reindex(columns=cell_type_order)
        .fillna(0)
    )

    tf_degree_matrix.to_csv(tf_degree_matrix_out)

    tf_similarity = tf_degree_matrix.corr(method="pearson")
    tf_similarity.to_csv(tf_degree_similarity_out)

    print()
    print("Saved TF degree matrix to:")
    print(tf_degree_matrix_out)

    print()
    print("Saved TF degree similarity matrix to:")
    print(tf_degree_similarity_out)

    print()
    print("TF degree similarity matrix:")
    print(tf_similarity.to_string())

    # ------------------------------------------------------------
    # Plot TF degree similarity heatmap
    # ------------------------------------------------------------

    plt.figure(figsize=(6, 5))
    plt.imshow(tf_similarity.values, vmin=-1, vmax=1)
    plt.xticks(
        ticks=np.arange(tf_similarity.shape[1]),
        labels=tf_similarity.columns,
        rotation=45,
        ha="right"
    )
    plt.yticks(
        ticks=np.arange(tf_similarity.shape[0]),
        labels=tf_similarity.index
    )
    plt.colorbar(label="Pearson correlation")
    plt.title("Similarity between per-cell-type GRNs based on TF degree")
    plt.tight_layout()
    plt.savefig(
        fig_dir / "step47R_tf_degree_similarity_heatmap.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.show()

# ============================================================
# Top TFs per cell type
# ============================================================

if not all_collapsed.empty:
    top_tfs_per_cell_type = []

    for cell_type in cell_type_order:
        sub = all_collapsed[
            all_collapsed["cell_type"] == cell_type
        ].copy()

        if sub.empty:
            continue

        tf_summary = (
            sub
            .groupby("tf")
            .agg(
                n_target_genes=("target_gene", "nunique"),
                n_edges=("target_gene", "count"),
                n_unique_supporting_peaks_total=("n_unique_supporting_peaks", "sum"),
                max_edge_score=("max_cell_type_edge_evidence_score", "max"),
                mean_edge_score=("max_cell_type_edge_evidence_score", "mean"),
            )
            .reset_index()
            .sort_values(
                ["n_target_genes", "max_edge_score"],
                ascending=[False, False]
            )
        )

        tf_summary["cell_type"] = cell_type

        top_tfs_per_cell_type.append(tf_summary)

        tf_summary.to_csv(
            out_dir / f"step47R_top_tfs_{safe_cell_type_name(cell_type)}.csv",
            index=False
        )

    top_tfs_all = pd.concat(top_tfs_per_cell_type, axis=0, ignore_index=True)
    top_tfs_all.to_csv(
        out_dir / "step47R_top_tfs_all_cell_types.csv",
        index=False
    )

    print()
    print("Top 10 TFs per cell type:")
    for cell_type in cell_type_order:
        sub = top_tfs_all[top_tfs_all["cell_type"] == cell_type]
        print()
        print(cell_type)
        print(sub[["tf", "n_target_genes", "n_edges", "max_edge_score"]].head(10).to_string(index=False))

print()
print("Βήμα 47R complete.")



"""Βήμα 48R — High-confidence filtering και QC των per-cell-type candidate GRNs"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

in_dir = project_dir / "data/processed/grn_per_cell_type"
out_dir = project_dir / "data/processed/grn_per_cell_type/high_confidence"
fig_dir = project_dir / "figures/grn_per_cell_type/high_confidence"

out_dir.mkdir(parents=True, exist_ok=True)
fig_dir.mkdir(parents=True, exist_ok=True)

collapsed_path = in_dir / "step47R_all_cell_type_tf_target_edges_collapsed.csv"

high_conf_all_out = out_dir / "step48R_high_confidence_per_cell_type_grn_all_edges.csv"
high_conf_clean_out = out_dir / "step48R_high_confidence_per_cell_type_grn_non_housekeeping_edges.csv"
summary_out = out_dir / "step48R_high_confidence_per_cell_type_summary.csv"
tf_summary_out = out_dir / "step48R_high_confidence_top_tfs_all_cell_types.csv"
target_summary_out = out_dir / "step48R_high_confidence_top_targets_all_cell_types.csv"
tf_degree_matrix_out = out_dir / "step48R_high_confidence_tf_degree_matrix_by_cell_type.csv"
tf_degree_similarity_out = out_dir / "step48R_high_confidence_tf_degree_similarity_by_cell_type.csv"

# ============================================================
# Parameters
# ============================================================

SCORE_QUANTILE_PER_CELL_TYPE = 0.95

MIN_UNIQUE_SUPPORTING_PEAKS = 1
MIN_UNIQUE_MOTIFS = 1
MIN_PEAK_DETECTION_FRACTION = 0.01
MIN_TF_RNA_EXPRESSION = 0.05
MIN_TARGET_RNA_EXPRESSION = 0.05

REMOVE_HOUSEKEEPING_LIKE_TARGETS = True

# Optional additional filter. Keep False initially.
USE_GLOBAL_TARGET_RNA_ATAC_CORR_FILTER = False
MIN_GLOBAL_TARGET_RNA_ATAC_CORR = 0.40

# ============================================================
# Cell type order
# ============================================================

cell_type_order = [
    "CD14+ Monocytes",
    "CD4.Memory",
    "CD4.Naive",
    "CD8.effector",
    "CD8.Naive",
]

# ============================================================
# Load per-cell-type collapsed GRN
# ============================================================

grn = pd.read_csv(collapsed_path)

grn["cell_type"] = grn["cell_type"].astype(str).str.strip()
grn["tf"] = grn["tf"].astype(str).str.strip().str.upper()
grn["target_gene"] = grn["target_gene"].astype(str).str.strip().str.upper()

print("Loaded per-cell-type collapsed GRN:")
print(grn.shape)

print()
print("Cell types:")
print(grn["cell_type"].value_counts())

# ============================================================
# Ensure numeric columns
# ============================================================

numeric_cols = [
    "n_supporting_peak_motif_links",
    "n_unique_supporting_peaks",
    "n_unique_motifs",
    "n_promoter_peak_links",
    "n_distal_peak_links",
    "max_cell_type_edge_evidence_score",
    "mean_cell_type_edge_evidence_score",
    "max_peak_mean_accessibility",
    "max_peak_detection_fraction",
    "max_tf_rna_expression",
    "max_target_gene_rna_expression",
    "max_target_gene_atac_gene_activity",
    "min_abs_distance_to_gene",
    "global_target_gene_rna_atac_corr",
    "top_evidence_score",
]

for col in numeric_cols:
    if col in grn.columns:
        grn[col] = pd.to_numeric(grn[col], errors="coerce")

# ============================================================
# Housekeeping-like target gene classification
# ============================================================

def classify_target_gene(gene):
    gene = str(gene).upper()

    if gene.startswith("RPL") or gene.startswith("RPS"):
        return "ribosomal"

    if gene.startswith("MT-"):
        return "mitochondrial"

    if gene.startswith("HSP"):
        return "heat_shock"

    if (
        gene.startswith("HIST") or
        gene.startswith("H1-") or
        gene.startswith("H2A") or
        gene.startswith("H2B") or
        gene.startswith("H3") or
        gene.startswith("H4")
    ):
        return "histone"

    if gene in {
        "MALAT1", "NEAT1", "B2M", "ACTB", "GAPDH",
        "TUBB", "TUBA1B", "TMSB4X", "TMSB10"
    }:
        return "common_high_expression"

    return "other"

grn["target_gene_category"] = grn["target_gene"].map(classify_target_gene)

# ============================================================
# Per-cell-type high-confidence filtering
# ============================================================

high_conf_parts = []
summary_records = []

for cell_type in cell_type_order:
    sub = grn[grn["cell_type"] == cell_type].copy()

    if sub.empty:
        summary_records.append({
            "cell_type": cell_type,
            "raw_edges": 0,
            "score_threshold_q95": np.nan,
            "high_confidence_edges_all_targets": 0,
            "high_confidence_edges_non_housekeeping": 0,
            "unique_tfs_non_housekeeping": 0,
            "unique_targets_non_housekeeping": 0,
        })
        continue

    score_threshold = sub["max_cell_type_edge_evidence_score"].quantile(
        SCORE_QUANTILE_PER_CELL_TYPE
    )

    filtered = sub[
        (sub["max_cell_type_edge_evidence_score"] >= score_threshold) &
        (sub["n_unique_supporting_peaks"] >= MIN_UNIQUE_SUPPORTING_PEAKS) &
        (sub["n_unique_motifs"] >= MIN_UNIQUE_MOTIFS) &
        (sub["max_peak_detection_fraction"] >= MIN_PEAK_DETECTION_FRACTION) &
        (sub["max_tf_rna_expression"] >= MIN_TF_RNA_EXPRESSION) &
        (sub["max_target_gene_rna_expression"] >= MIN_TARGET_RNA_EXPRESSION)
    ].copy()

    if USE_GLOBAL_TARGET_RNA_ATAC_CORR_FILTER and "global_target_gene_rna_atac_corr" in filtered.columns:
        filtered = filtered[
            filtered["global_target_gene_rna_atac_corr"] >= MIN_GLOBAL_TARGET_RNA_ATAC_CORR
        ].copy()

    filtered["score_threshold_used"] = score_threshold

    high_conf_parts.append(filtered)

    filtered_clean = filtered[
        filtered["target_gene_category"] == "other"
    ].copy()

    summary_records.append({
        "cell_type": cell_type,
        "raw_edges": sub.shape[0],
        "score_threshold_q95": score_threshold,
        "high_confidence_edges_all_targets": filtered.shape[0],
        "high_confidence_edges_non_housekeeping": filtered_clean.shape[0],
        "unique_tfs_non_housekeeping": filtered_clean["tf"].nunique(),
        "unique_targets_non_housekeeping": filtered_clean["target_gene"].nunique(),
        "edges_with_at_least_2_supporting_peaks": int((filtered_clean["n_unique_supporting_peaks"] >= 2).sum()),
        "edges_with_promoter_support": int((filtered_clean["n_promoter_peak_links"] > 0).sum()),
        "edges_with_distal_support": int((filtered_clean["n_distal_peak_links"] > 0).sum()),
        "mean_edge_score_non_housekeeping": filtered_clean["max_cell_type_edge_evidence_score"].mean(),
        "median_edge_score_non_housekeeping": filtered_clean["max_cell_type_edge_evidence_score"].median(),
    })

# ============================================================
# Combine high-confidence outputs
# ============================================================

high_conf = pd.concat(high_conf_parts, axis=0, ignore_index=True)

if REMOVE_HOUSEKEEPING_LIKE_TARGETS:
    high_conf_clean = high_conf[
        high_conf["target_gene_category"] == "other"
    ].copy()
else:
    high_conf_clean = high_conf.copy()

high_conf = high_conf.sort_values(
    ["cell_type", "max_cell_type_edge_evidence_score"],
    ascending=[True, False]
)

high_conf_clean = high_conf_clean.sort_values(
    ["cell_type", "max_cell_type_edge_evidence_score"],
    ascending=[True, False]
)

high_conf.to_csv(high_conf_all_out, index=False)
high_conf_clean.to_csv(high_conf_clean_out, index=False)

summary = pd.DataFrame(summary_records)
summary.to_csv(summary_out, index=False)

print()
print("High-confidence all-target edges:")
print(high_conf.shape)

print()
print("High-confidence non-housekeeping edges:")
print(high_conf_clean.shape)

print()
print("High-confidence per-cell-type summary:")
print(summary.to_string(index=False))

# ============================================================
# Save separate high-confidence GRN per cell type
# ============================================================

def safe_name(x):
    return (
        str(x)
        .replace("+", "plus")
        .replace(" ", "_")
        .replace(".", "_")
        .replace("-", "_")
        .replace("/", "_")
    )

for cell_type in cell_type_order:
    sub = high_conf_clean[
        high_conf_clean["cell_type"] == cell_type
    ].copy()

    sub = sub.sort_values(
        ["max_cell_type_edge_evidence_score", "n_unique_supporting_peaks"],
        ascending=[False, False]
    )

    out_path = out_dir / f"step48R_high_confidence_grn_{safe_name(cell_type)}.csv"
    sub.to_csv(out_path, index=False)

# ============================================================
# TF summary per cell type
# ============================================================

tf_summary = (
    high_conf_clean
    .groupby(["cell_type", "tf"])
    .agg(
        n_target_genes=("target_gene", "nunique"),
        n_edges=("target_gene", "count"),
        n_supporting_peak_motif_links=("n_supporting_peak_motif_links", "sum"),
        n_unique_supporting_peaks_total=("n_unique_supporting_peaks", "sum"),
        n_unique_motifs_total=("n_unique_motifs", "sum"),
        n_promoter_peak_links=("n_promoter_peak_links", "sum"),
        n_distal_peak_links=("n_distal_peak_links", "sum"),
        max_edge_score=("max_cell_type_edge_evidence_score", "max"),
        mean_edge_score=("max_cell_type_edge_evidence_score", "mean"),
        median_edge_score=("max_cell_type_edge_evidence_score", "median"),
        mean_tf_rna_expression=("max_tf_rna_expression", "mean"),
        mean_target_gene_rna_expression=("max_target_gene_rna_expression", "mean"),
    )
    .reset_index()
    .sort_values(
        ["cell_type", "n_target_genes", "max_edge_score"],
        ascending=[True, False, False]
    )
)

tf_summary.to_csv(tf_summary_out, index=False)

print()
print("Top 10 TFs per cell type after high-confidence filtering:")
for cell_type in cell_type_order:
    sub = tf_summary[tf_summary["cell_type"] == cell_type].copy()
    print()
    print(cell_type)
    print(
        sub[
            [
                "tf",
                "n_target_genes",
                "n_edges",
                "max_edge_score",
                "mean_edge_score",
                "mean_tf_rna_expression",
            ]
        ].head(10).to_string(index=False)
    )

# ============================================================
# Target summary per cell type
# ============================================================

target_summary = (
    high_conf_clean
    .groupby(["cell_type", "target_gene"])
    .agg(
        n_regulator_tfs=("tf", "nunique"),
        n_edges=("tf", "count"),
        n_supporting_peak_motif_links=("n_supporting_peak_motif_links", "sum"),
        n_unique_supporting_peaks_total=("n_unique_supporting_peaks", "sum"),
        max_edge_score=("max_cell_type_edge_evidence_score", "max"),
        mean_edge_score=("max_cell_type_edge_evidence_score", "mean"),
        max_target_gene_rna_expression=("max_target_gene_rna_expression", "max"),
    )
    .reset_index()
    .sort_values(
        ["cell_type", "n_regulator_tfs", "max_edge_score"],
        ascending=[True, False, False]
    )
)

target_summary.to_csv(target_summary_out, index=False)

# ============================================================
# TF degree matrix and similarity after high-confidence filtering
# ============================================================

tf_degree_matrix = (
    tf_summary
    .pivot_table(
        index="tf",
        columns="cell_type",
        values="n_target_genes",
        fill_value=0
    )
    .reindex(columns=cell_type_order)
    .fillna(0)
)

tf_degree_matrix.to_csv(tf_degree_matrix_out)

tf_similarity = tf_degree_matrix.corr(method="pearson")
tf_similarity.to_csv(tf_degree_similarity_out)

print()
print("High-confidence TF degree similarity matrix:")
print(tf_similarity.to_string())

# ============================================================
# Known immune TF-focused summary
# ============================================================

known_immune_tfs = [
    "SPI1", "CEBPA", "CEBPB", "CEBPD",
    "FOS", "FOSB", "FOSL1", "FOSL2",
    "JUN", "JUNB", "JUND",
    "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9",
    "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
    "NFKB1", "NFKB2", "RELA", "RELB", "REL",
    "KLF2", "KLF3", "KLF4", "KLF6", "KLF12",
    "IKZF1", "IKZF2", "IKZF3",
    "ETS1", "ELF1", "ELK1",
    "RUNX1", "RUNX2", "RUNX3",
    "TCF7", "LEF1",
    "GATA3", "TBX21", "EOMES",
    "PRDM1", "BCL6", "FOXP1", "FOXO1",
    "MEF2A", "MEF2C",
    "NR4A1", "NR4A2", "NR4A3",
    "AHR", "ARNT", "BATF", "MAF", "MAFB",
    "NFATC1", "NFATC2", "NFATC3",
]

immune_tf_summary = tf_summary[
    tf_summary["tf"].isin(known_immune_tfs)
].copy()

immune_tf_summary = immune_tf_summary.sort_values(
    ["cell_type", "n_target_genes", "max_edge_score"],
    ascending=[True, False, False]
)

immune_tf_summary.to_csv(
    out_dir / "step48R_known_immune_tf_summary_per_cell_type.csv",
    index=False
)

print()
print("Known immune TFs after high-confidence filtering:")
for cell_type in cell_type_order:
    sub = immune_tf_summary[immune_tf_summary["cell_type"] == cell_type]
    print()
    print(cell_type)
    print(
        sub[
            [
                "tf",
                "n_target_genes",
                "n_edges",
                "max_edge_score",
                "mean_tf_rna_expression",
            ]
        ].head(15).to_string(index=False)
    )

# ============================================================
# Plot helper
# ============================================================

def save_current_plot(output_path):
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()

# ============================================================
# Plot 1 — High-confidence edge counts per cell type
# ============================================================

plot_df = summary.sort_values("high_confidence_edges_non_housekeeping", ascending=True)

plt.figure(figsize=(8, 5))
plt.barh(plot_df["cell_type"], plot_df["high_confidence_edges_non_housekeeping"])
plt.xlabel("Number of high-confidence non-housekeeping edges")
plt.ylabel("Cell type")
plt.title("High-confidence per-cell-type GRN sizes")
save_current_plot(fig_dir / "step48R_high_confidence_grn_sizes.png")

# ============================================================
# Plot 2 — TF degree similarity heatmap
# ============================================================

plt.figure(figsize=(6, 5))
plt.imshow(tf_similarity.values, vmin=-1, vmax=1)
plt.xticks(
    ticks=np.arange(tf_similarity.shape[1]),
    labels=tf_similarity.columns,
    rotation=45,
    ha="right"
)
plt.yticks(
    ticks=np.arange(tf_similarity.shape[0]),
    labels=tf_similarity.index
)
plt.colorbar(label="Pearson correlation")
plt.title("High-confidence GRN similarity based on TF degree")
save_current_plot(fig_dir / "step48R_high_confidence_tf_degree_similarity_heatmap.png")

# ============================================================
# Plot 3 — Top TFs per cell type
# ============================================================

for cell_type in cell_type_order:
    sub = tf_summary[tf_summary["cell_type"] == cell_type].head(20).copy()

    if sub.empty:
        continue

    sub = sub.sort_values("n_target_genes", ascending=True)

    plt.figure(figsize=(8, 6))
    plt.barh(sub["tf"], sub["n_target_genes"])
    plt.xlabel("Number of target genes")
    plt.ylabel("TF")
    plt.title(f"Top TFs in {cell_type} high-confidence GRN")
    save_current_plot(
        fig_dir / f"step48R_top_tfs_{safe_name(cell_type)}.png"
    )

# ============================================================
# Final messages
# ============================================================

print()
print("Saved Step 48R high-confidence GRN tables to:")
print(out_dir)

print()
print("Saved Step 48R figures to:")
print(fig_dir)

print()
print("Main files:")
print(high_conf_clean_out)
print(summary_out)
print(tf_summary_out)
print(tf_degree_similarity_out)

print()
print("Βήμα 48R complete.")



"""Βήμα 49S — TF-centered network visualization των final per-cell-type GRNs"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# Optional networkx import
# ============================================================

try:
    import networkx as nx
except ImportError:
    raise ImportError(
        "Το package networkx δεν είναι εγκατεστημένο. "
        "Τρέξε πρώτα: pip install networkx"
    )

# ============================================================
# Paths
# ============================================================

project_dir = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

final_grn_path = (
    project_dir
    / "data/processed/final_grn_per_cell_type/step49R_final_per_cell_type_grn_edges.csv"
)

network_fig_dir = (
    project_dir
    / "figures/final_grn_per_cell_type/tf_centered_networks"
)

network_table_dir = (
    project_dir
    / "data/processed/final_grn_per_cell_type/tf_centered_networks"
)

network_fig_dir.mkdir(parents=True, exist_ok=True)
network_table_dir.mkdir(parents=True, exist_ok=True)

# ============================================================
# Load final GRN
# ============================================================

grn = pd.read_csv(final_grn_path)

grn["cell_type"] = grn["cell_type"].astype(str).str.strip()
grn["tf"] = grn["tf"].astype(str).str.strip().str.upper()
grn["target_gene"] = grn["target_gene"].astype(str).str.strip().str.upper()

numeric_cols = [
    "max_cell_type_edge_evidence_score",
    "n_unique_supporting_peaks",
    "n_unique_motifs",
    "n_promoter_peak_links",
    "n_distal_peak_links",
    "max_tf_rna_expression",
    "max_target_gene_rna_expression",
    "max_peak_detection_fraction",
]

for col in numeric_cols:
    if col in grn.columns:
        grn[col] = pd.to_numeric(grn[col], errors="coerce")

print("Loaded final per-cell-type GRN:")
print(grn.shape)

print()
print("Cell types:")
print(grn["cell_type"].value_counts())

# ============================================================
# Helper functions
# ============================================================

def safe_name(x):
    return (
        str(x)
        .replace("+", "plus")
        .replace(" ", "_")
        .replace(".", "_")
        .replace("-", "_")
        .replace("/", "_")
    )


def get_tf_edges(cell_type, tf, max_targets=25):
    """
    Select top target edges for a given TF in a given cell type.
    """
    sub = grn[
        (grn["cell_type"] == cell_type) &
        (grn["tf"] == tf.upper())
    ].copy()

    if sub.empty:
        return sub

    sub = sub.sort_values(
        [
            "max_cell_type_edge_evidence_score",
            "n_unique_supporting_peaks",
            "max_peak_detection_fraction",
        ],
        ascending=[False, False, False]
    )

    return sub.head(max_targets).copy()


def plot_tf_network(cell_type, tf, max_targets=25):
    """
    Draw a TF-centered directed network:
    TF -> target genes.
    """
    edges = get_tf_edges(
        cell_type=cell_type,
        tf=tf,
        max_targets=max_targets
    )

    if edges.empty:
        print()
        print(f"No edges found for {tf} in {cell_type}. Skipping.")
        return None

    # ------------------------------------------------------------
    # Build graph
    # ------------------------------------------------------------

    G = nx.DiGraph()

    tf_node = tf.upper()
    G.add_node(tf_node, node_type="TF")

    for _, row in edges.iterrows():
        target = row["target_gene"]
        score = row["max_cell_type_edge_evidence_score"]

        G.add_node(target, node_type="target_gene")
        G.add_edge(tf_node, target, weight=score)

    # ------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------

    pos = nx.spring_layout(
        G,
        seed=42,
        k=1.2,
        iterations=200
    )

    # Force TF approximately to center
    pos[tf_node] = np.array([0.0, 0.0])

    # ------------------------------------------------------------
    # Node/edge style
    # ------------------------------------------------------------

    target_nodes = [
        node for node, data in G.nodes(data=True)
        if data.get("node_type") == "target_gene"
    ]

    tf_nodes = [tf_node]

    edge_weights = [
        G[u][v].get("weight", 0.0)
        for u, v in G.edges()
    ]

    max_weight = max(edge_weights) if len(edge_weights) > 0 else 1.0

    edge_widths = [
        0.8 + 4.0 * (w / max_weight)
        for w in edge_weights
    ]

    target_scores = {
        row["target_gene"]: row["max_cell_type_edge_evidence_score"]
        for _, row in edges.iterrows()
    }

    target_sizes = [
        250 + 900 * (target_scores.get(node, 0) / max_weight)
        for node in target_nodes
    ]

    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------

    plt.figure(figsize=(11, 9))

    nx.draw_networkx_edges(
        G,
        pos,
        arrows=True,
        arrowstyle="-|>",
        arrowsize=12,
        width=edge_widths,
        alpha=0.55,
        connectionstyle="arc3,rad=0.08",
    )

    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=target_nodes,
        node_size=target_sizes,
        node_shape="o",
        alpha=0.85,
    )

    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=tf_nodes,
        node_size=1800,
        node_shape="s",
        alpha=0.95,
    )

    nx.draw_networkx_labels(
        G,
        pos,
        font_size=8,
        font_weight="normal",
    )

    plt.title(
        f"{tf_node}-centered candidate GRN\n{cell_type} | top {edges.shape[0]} targets",
        fontsize=13,
    )

    plt.axis("off")
    plt.tight_layout()

    # ------------------------------------------------------------
    # Save files
    # ------------------------------------------------------------

    fig_base = f"step49S_network_{safe_name(cell_type)}_{tf_node}"

    png_path = network_fig_dir / f"{fig_base}.png"
    pdf_path = network_fig_dir / f"{fig_base}.pdf"
    edge_table_path = network_table_dir / f"{fig_base}_edges.csv"

    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.show()

    edges.to_csv(edge_table_path, index=False)

    print()
    print(f"Saved TF-centered network for {tf_node} in {cell_type}:")
    print(png_path)
    print(pdf_path)
    print(edge_table_path)

    return edges


# ============================================================
# TFs to visualize per cell type
# ============================================================

tf_panels = {
    "CD14+ Monocytes": ["SPI1", "FOS", "JUN", "KLF4", "IRF1", "AHR"],
    "CD4.Memory": ["GATA3", "KLF2", "JUN", "FOS", "IKZF1", "IRF1"],
    "CD4.Naive": ["TCF7", "LEF1", "IKZF1", "JUN", "IRF1", "KLF2"],
    "CD8.effector": ["IRF1", "JUN", "JUNB", "FOSB", "GATA3", "IKZF1"],
    "CD8.Naive": ["TCF7", "LEF1", "IRF1", "JUN", "IKZF1", "FOS"],
}

# ============================================================
# Generate TF-centered network plots
# ============================================================

all_visualized_edges = []

for cell_type, tf_list in tf_panels.items():
    for tf in tf_list:
        edges = plot_tf_network(
            cell_type=cell_type,
            tf=tf,
            max_targets=25
        )

        if edges is not None and not edges.empty:
            all_visualized_edges.append(edges)

if len(all_visualized_edges) > 0:
    all_visualized_edges = pd.concat(
        all_visualized_edges,
        axis=0,
        ignore_index=True
    )

    all_visualized_edges.to_csv(
        network_table_dir / "step49S_all_visualized_tf_centered_edges.csv",
        index=False
    )

# ============================================================
# Combined immune TF network per cell type
# ============================================================

def plot_cell_type_top_tf_network(cell_type, top_tfs=6, top_targets_per_tf=8):
    """
    Draw a combined network for the top TFs in one cell type.
    """
    sub = grn[grn["cell_type"] == cell_type].copy()

    if sub.empty:
        print(f"No edges found for {cell_type}. Skipping combined plot.")
        return None

    tf_summary = (
        sub
        .groupby("tf")
        .agg(
            n_targets=("target_gene", "nunique"),
            max_score=("max_cell_type_edge_evidence_score", "max"),
        )
        .reset_index()
        .sort_values(["n_targets", "max_score"], ascending=[False, False])
    )

    selected_tfs = tf_summary.head(top_tfs)["tf"].tolist()

    selected_edges = []

    for tf in selected_tfs:
        tf_edges = (
            sub[sub["tf"] == tf]
            .sort_values(
                [
                    "max_cell_type_edge_evidence_score",
                    "n_unique_supporting_peaks",
                ],
                ascending=[False, False]
            )
            .head(top_targets_per_tf)
        )

        selected_edges.append(tf_edges)

    selected_edges = pd.concat(selected_edges, axis=0, ignore_index=True)

    # ------------------------------------------------------------
    # Build graph
    # ------------------------------------------------------------

    G = nx.DiGraph()

    for _, row in selected_edges.iterrows():
        tf = row["tf"]
        target = row["target_gene"]
        score = row["max_cell_type_edge_evidence_score"]

        G.add_node(tf, node_type="TF")
        G.add_node(target, node_type="target_gene")
        G.add_edge(tf, target, weight=score)

    pos = nx.spring_layout(
        G,
        seed=42,
        k=0.9,
        iterations=250
    )

    tf_nodes = [
        node for node, data in G.nodes(data=True)
        if data.get("node_type") == "TF"
    ]

    target_nodes = [
        node for node, data in G.nodes(data=True)
        if data.get("node_type") == "target_gene"
    ]

    edge_weights = [
        G[u][v].get("weight", 0.0)
        for u, v in G.edges()
    ]

    max_weight = max(edge_weights) if len(edge_weights) > 0 else 1.0

    edge_widths = [
        0.5 + 3.5 * (w / max_weight)
        for w in edge_weights
    ]

    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------

    plt.figure(figsize=(14, 11))

    nx.draw_networkx_edges(
        G,
        pos,
        arrows=True,
        arrowstyle="-|>",
        arrowsize=10,
        width=edge_widths,
        alpha=0.45,
        connectionstyle="arc3,rad=0.05",
    )

    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=target_nodes,
        node_size=300,
        node_shape="o",
        alpha=0.80,
    )

    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=tf_nodes,
        node_size=1200,
        node_shape="s",
        alpha=0.95,
    )

    nx.draw_networkx_labels(
        G,
        pos,
        font_size=7,
        font_weight="normal",
    )

    plt.title(
        f"Top TF-centered candidate GRN\n{cell_type}",
        fontsize=14,
    )

    plt.axis("off")
    plt.tight_layout()

    fig_base = f"step49S_combined_top_TF_network_{safe_name(cell_type)}"

    png_path = network_fig_dir / f"{fig_base}.png"
    pdf_path = network_fig_dir / f"{fig_base}.pdf"
    edge_table_path = network_table_dir / f"{fig_base}_edges.csv"

    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.show()

    selected_edges.to_csv(edge_table_path, index=False)

    print()
    print(f"Saved combined top-TF network for {cell_type}:")
    print(png_path)
    print(pdf_path)
    print(edge_table_path)

    return selected_edges


# ============================================================
# Generate combined cell-type network plots
# ============================================================

for cell_type in tf_panels.keys():
    plot_cell_type_top_tf_network(
        cell_type=cell_type,
        top_tfs=6,
        top_targets_per_tf=8
    )

# ============================================================
# Final messages
# ============================================================

print()
print("Saved TF-centered network figures to:")
print(network_fig_dir)

print()
print("Saved TF-centered network edge tables to:")
print(network_table_dir)

print()
print("Βήμα 49S complete.")








"""Βήμα 50G — Build pseudo-supervised GNN input tables"""

# ============================================================
# Imports
# ============================================================

from pathlib import Path
import re
import json
import random
import warnings

import numpy as np
import pandas as pd


# ============================================================
# Paths and configuration
# ============================================================

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

PROJECT_DIR = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
PROCESSED_DIR = PROJECT_DIR / "data" / "processed"
OUT_DIR = PROCESSED_DIR / "gnn"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Optional manual overrides.
# Αν ξέρεις ακριβώς το permissive per-cell-type GRN table, βάλε path εδώ.
MANUAL_CANDIDATE_EDGE_FILE = None

# Αν έχεις ξεχωριστό final/high-confidence per-cell-type GRN table, βάλε path εδώ.
MANUAL_HIGH_CONFIDENCE_EDGE_FILE = None

POSITIVE_QUANTILE_IF_NO_HC_FILE = 0.80
WEAK_NEGATIVE_QUANTILE = 0.20
SHUFFLED_NEGATIVES_PER_POSITIVE = 1.0
MAX_SHUFFLED_NEGATIVES_PER_CELL_TYPE = 100_000

REQUIRED_CELL_TYPES = [
    "CD14+ Monocytes",
    "CD4.Memory",
    "CD4.Naive",
    "CD8.effector",
    "CD8.Naive",
]


# ============================================================
# Helper functions
# ============================================================

def normalize_colname(x):
    return re.sub(r"[^a-z0-9]+", "_", str(x).strip().lower()).strip("_")


def safe_name(x):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x))


def find_column(df, candidates, required=True):
    norm_map = {normalize_colname(c): c for c in df.columns}
    candidate_norms = [normalize_colname(c) for c in candidates]

    for cn in candidate_norms:
        if cn in norm_map:
            return norm_map[cn]

    for cn in candidate_norms:
        for nc, original in norm_map.items():
            if cn in nc:
                return original

    if required:
        raise ValueError(f"Could not find required column. Tried: {candidates}")
    return None


def preview_csv_columns(path):
    try:
        return pd.read_csv(path, nrows=50).columns.tolist()
    except Exception:
        return []


def score_edge_file(path):
    cols = preview_csv_columns(path)
    if not cols:
        return -999

    norm_cols = {normalize_colname(c) for c in cols}
    name = normalize_colname(path.name)

    has_tf = bool(norm_cols & {"tf", "transcription_factor", "source_tf", "regulator", "source"})
    has_target = bool(norm_cols & {"target_gene", "target", "gene", "target_symbol", "target_name"})
    has_cell = bool(norm_cols & {"cell_type", "celltype", "cluster", "annotation"})
    has_score = any("score" in c or "evidence" in c for c in norm_cols)

    score = 0
    score += 10 if has_tf else 0
    score += 10 if has_target else 0
    score += 10 if has_cell else 0
    score += 5 if has_score else 0

    for token in ["per_cell", "cell_type", "candidate", "grn", "edge"]:
        if token in name:
            score += 2

    for token in ["cytoscape", "node", "summary", "tf_summary", "correlation"]:
        if token in name:
            score -= 4

    return score


def auto_find_edge_file(prefer_high_confidence=False):
    csvs = list(PROCESSED_DIR.rglob("*.csv"))
    if not csvs:
        raise FileNotFoundError(f"No CSV files found under: {PROCESSED_DIR}")

    scored = []
    for p in csvs:
        s = score_edge_file(p)
        name = normalize_colname(p.name)

        if prefer_high_confidence:
            if any(t in name for t in ["high_conf", "high_confidence", "final", "filtered"]):
                s += 8
            if "candidate" in name and not any(t in name for t in ["high", "final", "filtered"]):
                s -= 3
        else:
            if "candidate" in name:
                s += 8
            if any(t in name for t in ["high_conf", "high_confidence", "final", "filtered"]):
                s -= 3

        scored.append((s, p))

    scored = sorted(scored, key=lambda x: x[0], reverse=True)
    best_score, best_path = scored[0]

    if best_score < 25:
        raise ValueError(
            "Could not confidently auto-detect a per-cell-type TF-target edge table. "
            "Set MANUAL_CANDIDATE_EDGE_FILE manually."
        )

    return best_path


def canonicalize_edge_table(df):
    df = df.copy()

    tf_col = find_column(
        df,
        ["tf", "TF", "transcription_factor", "source_tf", "regulator", "source"],
        required=True,
    )
    target_col = find_column(
        df,
        ["target_gene", "target", "gene", "target_symbol", "target_name"],
        required=True,
    )
    cell_col = find_column(
        df,
        ["cell_type", "celltype", "cluster", "annotation"],
        required=True,
    )

    rename = {
        tf_col: "tf",
        target_col: "target_gene",
        cell_col: "cell_type",
    }
    df = df.rename(columns=rename)

    df["tf"] = df["tf"].astype(str).str.strip()
    df["target_gene"] = df["target_gene"].astype(str).str.strip()
    df["cell_type"] = df["cell_type"].astype(str).str.strip()

    df = df[
        (df["tf"] != "")
        & (df["target_gene"] != "")
        & (df["tf"].str.lower() != "nan")
        & (df["target_gene"].str.lower() != "nan")
    ].copy()

    return df


def choose_score_column(df):
    priority = [
        "edge_evidence_score",
        "final_edge_score",
        "combined_edge_score",
        "evidence_score",
        "grn_edge_score",
        "score",
    ]

    norm_to_original = {normalize_colname(c): c for c in df.columns}
    for p in priority:
        if p in norm_to_original:
            return norm_to_original[p]

    numeric_cols = [
        c for c in df.columns
        if c not in {"tf", "target_gene", "cell_type"}
        and pd.api.types.is_numeric_dtype(df[c])
    ]

    numeric_score_like = [
        c for c in numeric_cols
        if any(t in normalize_colname(c) for t in ["score", "evidence", "access", "detect", "expr", "support"])
    ]

    if numeric_score_like:
        synthetic = "__synthetic_evidence_score"
        z = df[numeric_score_like].copy()
        z = z.apply(pd.to_numeric, errors="coerce")
        z = z.rank(pct=True)
        df[synthetic] = z.mean(axis=1)
        return synthetic

    synthetic = "__synthetic_evidence_score"
    df[synthetic] = 1.0
    return synthetic


def collapse_to_tf_target_edges(df, score_col):
    df = df.copy()
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce").fillna(0)

    # Keep strongest row if duplicate TF-target-cell_type rows exist.
    df = df.sort_values(["cell_type", "tf", "target_gene", score_col])
    df = df.drop_duplicates(["cell_type", "tf", "target_gene"], keep="last").copy()
    return df


def make_shuffled_negatives(ct_df, n_requested):
    ct = ct_df["cell_type"].iloc[0]
    tfs = sorted(ct_df["tf"].dropna().astype(str).unique())
    targets = sorted(ct_df["target_gene"].dropna().astype(str).unique())
    existing = set(zip(ct_df["tf"].astype(str), ct_df["target_gene"].astype(str)))

    if len(tfs) == 0 or len(targets) == 0:
        return pd.DataFrame()

    possible_upper = len(tfs) * len(targets) - len(existing)
    n_to_sample = int(min(n_requested, possible_upper, MAX_SHUFFLED_NEGATIVES_PER_CELL_TYPE))

    rows = []
    seen = set()
    attempts = 0
    max_attempts = max(10_000, n_to_sample * 30)

    while len(rows) < n_to_sample and attempts < max_attempts:
        attempts += 1
        tf = random.choice(tfs)
        target = random.choice(targets)

        if tf == target:
            continue
        if (tf, target) in existing:
            continue
        if (tf, target) in seen:
            continue

        seen.add((tf, target))
        rows.append({
            "cell_type": ct,
            "tf": tf,
            "target_gene": target,
            "pseudo_label": 0,
            "edge_source_type": "shuffled_negative",
            "is_pseudo_positive": 0,
            "is_weak_negative": 0,
            "is_shuffled_negative": 1,
            "edge_evidence_score_for_labeling": 0.0,
            "edge_evidence_rank_for_labeling": 0.0,
        })

    return pd.DataFrame(rows)


# ============================================================
# Load candidate and optional high-confidence GRN tables
# ============================================================

candidate_path = Path(MANUAL_CANDIDATE_EDGE_FILE) if MANUAL_CANDIDATE_EDGE_FILE else auto_find_edge_file(False)
print(f"[50G] Candidate edge table: {candidate_path}")

candidate_df = pd.read_csv(candidate_path)
candidate_df = canonicalize_edge_table(candidate_df)

score_col = choose_score_column(candidate_df)
candidate_df = collapse_to_tf_target_edges(candidate_df, score_col)

candidate_df["edge_evidence_score_for_labeling"] = pd.to_numeric(
    candidate_df[score_col], errors="coerce"
).fillna(0)

candidate_df["edge_evidence_rank_for_labeling"] = (
    candidate_df
    .groupby("cell_type")["edge_evidence_score_for_labeling"]
    .rank(method="average", pct=True)
)

candidate_df["edge_source_type"] = "candidate"
candidate_df["is_shuffled_negative"] = 0


# ============================================================
# Mark pseudo-positive edges
# ============================================================

candidate_df["is_high_confidence_input"] = 0

if MANUAL_HIGH_CONFIDENCE_EDGE_FILE is not None:
    high_conf_path = Path(MANUAL_HIGH_CONFIDENCE_EDGE_FILE)
else:
    try:
        high_conf_path = auto_find_edge_file(prefer_high_confidence=True)
        if high_conf_path == candidate_path:
            high_conf_path = None
    except Exception:
        high_conf_path = None

if high_conf_path is not None and high_conf_path.exists():
    print(f"[50G] High-confidence edge table: {high_conf_path}")
    high_df = pd.read_csv(high_conf_path)
    high_df = canonicalize_edge_table(high_df)

    high_keys = high_df[["cell_type", "tf", "target_gene"]].drop_duplicates()
    high_keys["is_high_confidence_input"] = 1

    candidate_df = candidate_df.drop(columns=["is_high_confidence_input"]).merge(
        high_keys,
        on=["cell_type", "tf", "target_gene"],
        how="left",
    )
    candidate_df["is_high_confidence_input"] = candidate_df["is_high_confidence_input"].fillna(0).astype(int)
else:
    print("[50G] No separate high-confidence table detected. Using top quantile as pseudo-positive labels.")

candidate_df["is_pseudo_positive"] = (
    (candidate_df["is_high_confidence_input"] == 1)
    | (candidate_df["edge_evidence_rank_for_labeling"] >= POSITIVE_QUANTILE_IF_NO_HC_FILE)
).astype(int)

candidate_df["is_weak_negative"] = (
    (candidate_df["is_pseudo_positive"] == 0)
    & (candidate_df["edge_evidence_rank_for_labeling"] <= WEAK_NEGATIVE_QUANTILE)
).astype(int)

candidate_df["pseudo_label"] = -1
candidate_df.loc[candidate_df["is_pseudo_positive"] == 1, "pseudo_label"] = 1
candidate_df.loc[candidate_df["is_weak_negative"] == 1, "pseudo_label"] = 0


# ============================================================
# Add shuffled/background negatives
# ============================================================

shuffle_tables = []

for ct, ct_df in candidate_df.groupby("cell_type", sort=False):
    n_pos = int((ct_df["is_pseudo_positive"] == 1).sum())
    n_shuffle = int(np.ceil(n_pos * SHUFFLED_NEGATIVES_PER_POSITIVE))

    if n_pos == 0:
        warnings.warn(f"[50G] No pseudo-positives for {ct}; skipping shuffled negatives.")
        continue

    shuffled = make_shuffled_negatives(ct_df, n_shuffle)
    shuffle_tables.append(shuffled)
    print(f"[50G] {ct}: positives={n_pos:,}, shuffled_negatives={len(shuffled):,}")

shuffled_df = pd.concat(shuffle_tables, ignore_index=True) if shuffle_tables else pd.DataFrame()


# ============================================================
# Build candidate scoring table and training examples table
# ============================================================

candidate_df["edge_id"] = (
    candidate_df["cell_type"].map(safe_name)
    + "|"
    + candidate_df["tf"].astype(str)
    + "|"
    + candidate_df["target_gene"].astype(str)
)

labeled_candidate_df = candidate_df[candidate_df["pseudo_label"].isin([0, 1])].copy()

if not shuffled_df.empty:
    for col in labeled_candidate_df.columns:
        if col not in shuffled_df.columns:
            shuffled_df[col] = np.nan

    shuffled_df["edge_id"] = (
        shuffled_df["cell_type"].map(safe_name)
        + "|"
        + shuffled_df["tf"].astype(str)
        + "|"
        + shuffled_df["target_gene"].astype(str)
    )

    edge_examples_df = pd.concat(
        [labeled_candidate_df, shuffled_df[labeled_candidate_df.columns]],
        ignore_index=True,
    )
else:
    edge_examples_df = labeled_candidate_df.copy()

candidate_out = OUT_DIR / "pbmc_gnn_candidate_edges_for_scoring.csv"
examples_out = OUT_DIR / "pbmc_gnn_edge_examples.csv"
manifest_out = OUT_DIR / "pbmc_gnn_step50_manifest.json"

candidate_df.to_csv(candidate_out, index=False)
edge_examples_df.to_csv(examples_out, index=False)

manifest = {
    "candidate_edge_file_detected": str(candidate_path),
    "high_confidence_edge_file_detected": str(high_conf_path) if high_conf_path is not None else None,
    "score_column_used_for_labeling": score_col,
    "positive_quantile_if_no_hc_file": POSITIVE_QUANTILE_IF_NO_HC_FILE,
    "weak_negative_quantile": WEAK_NEGATIVE_QUANTILE,
    "shuffled_negatives_per_positive": SHUFFLED_NEGATIVES_PER_POSITIVE,
    "candidate_edges_for_scoring": str(candidate_out),
    "edge_examples": str(examples_out),
}

with open(manifest_out, "w", encoding="utf-8") as f:
    json.dump(manifest, f, indent=2)

print("\n[50G] Done.")
print(f"[50G] Candidate edges for scoring: {candidate_out}")
print(f"[50G] Training examples: {examples_out}")
print("\n[50G] Label summary:")
print(edge_examples_df.groupby(["cell_type", "pseudo_label", "edge_source_type"]).size())



"""Βήμα 51G.0 — Inspect Step 50G outputs before feature matrix construction"""

# ============================================================
# Imports
# ============================================================

from pathlib import Path
import pandas as pd
import numpy as np
import re


# ============================================================
# Paths
# ============================================================

PROJECT_DIR = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
GNN_DIR = PROJECT_DIR / "data" / "processed" / "gnn"

CANDIDATE_FILE = GNN_DIR / "pbmc_gnn_candidate_edges_for_scoring.csv"
EXAMPLES_FILE = GNN_DIR / "pbmc_gnn_edge_examples.csv"

REPORT_OUT = GNN_DIR / "pbmc_gnn_step51_column_inspection.csv"


# ============================================================
# Helper functions
# ============================================================

def normalize_colname(x):
    return re.sub(r"[^a-z0-9]+", "_", str(x).strip().lower()).strip("_")


def feature_column_category(col):
    nc = normalize_colname(col)

    if nc in {"cell_type", "tf", "target_gene", "edge_id"}:
        return "core_id_column"

    if nc in {
        "pseudo_label",
        "is_pseudo_positive",
        "is_weak_negative",
        "is_shuffled_negative",
        "is_high_confidence_input",
        "edge_source_type",
    }:
        return "label_or_source_column"

    if any(x in nc for x in ["access", "accessible", "accessibility"]):
        return "accessibility_feature_candidate"

    if any(x in nc for x in ["detect", "detection"]):
        return "detection_feature_candidate"

    if any(x in nc for x in ["tf_expr", "tf_rna", "tf_expression"]):
        return "tf_expression_feature_candidate"

    if any(x in nc for x in ["target_expr", "target_rna", "target_expression"]):
        return "target_expression_feature_candidate"

    if any(x in nc for x in ["gene_activity", "activity", "atac"]):
        return "atac_gene_activity_feature_candidate"

    if any(x in nc for x in ["motif"]):
        return "motif_feature_candidate"

    if any(x in nc for x in ["peak"]):
        return "peak_feature_candidate"

    if any(x in nc for x in ["distance", "dist_to", "tss"]):
        return "distance_feature_candidate"

    if any(x in nc for x in ["promoter", "distal"]):
        return "regulatory_region_feature_candidate"

    if any(x in nc for x in ["support"]):
        return "support_feature_candidate"

    if any(x in nc for x in ["score", "evidence", "rank"]):
        return "score_or_evidence_feature_candidate"

    return "other"


def inspect_table(df, table_name):
    rows = []

    for col in df.columns:
        rows.append({
            "table": table_name,
            "column": col,
            "normalized_column": normalize_colname(col),
            "dtype": str(df[col].dtype),
            "n_missing": int(df[col].isna().sum()),
            "n_unique": int(df[col].nunique(dropna=True)),
            "category": feature_column_category(col),
            "example_values": ", ".join(
                df[col].dropna().astype(str).head(5).tolist()
            )
        })

    return pd.DataFrame(rows)


# ============================================================
# Load Step 50G outputs
# ============================================================

if not CANDIDATE_FILE.exists():
    raise FileNotFoundError(f"Missing candidate file: {CANDIDATE_FILE}")

if not EXAMPLES_FILE.exists():
    raise FileNotFoundError(f"Missing examples file: {EXAMPLES_FILE}")

candidate_df = pd.read_csv(CANDIDATE_FILE)
examples_df = pd.read_csv(EXAMPLES_FILE)


# ============================================================
# Basic table summaries
# ============================================================

print("Candidate edges for scoring")
print("Path:", CANDIDATE_FILE)
print("Shape:", candidate_df.shape)
print()

print("Training edge examples")
print("Path:", EXAMPLES_FILE)
print("Shape:", examples_df.shape)
print()


# ============================================================
# Required core columns check
# ============================================================

required_cols = ["cell_type", "tf", "target_gene", "edge_id"]

for table_name, df in {
    "candidate_df": candidate_df,
    "examples_df": examples_df,
}.items():
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f"[WARNING] {table_name} is missing required columns: {missing}")
    else:
        print(f"[OK] {table_name} has all required core columns.")

print()


# ============================================================
# Label summary
# ============================================================

if "pseudo_label" in examples_df.columns:
    print("Training label summary:")
    print(
        examples_df
        .groupby(["cell_type", "pseudo_label", "edge_source_type"])
        .size()
        .rename("n_edges")
        .reset_index()
    )
else:
    print("[WARNING] examples_df has no pseudo_label column.")

print()


# ============================================================
# Duplicate edge checks
# ============================================================

for table_name, df in {
    "candidate_df": candidate_df,
    "examples_df": examples_df,
}.items():
    duplicated_edges = df.duplicated(["cell_type", "tf", "target_gene"]).sum()
    duplicated_edge_ids = df.duplicated(["edge_id"]).sum() if "edge_id" in df.columns else np.nan

    print(f"{table_name} duplicate check:")
    print("Duplicated cell_type/tf/target_gene rows:", duplicated_edges)
    print("Duplicated edge_id rows:", duplicated_edge_ids)
    print()


# ============================================================
# Column inspection report
# ============================================================

candidate_report = inspect_table(candidate_df, "candidate_edges_for_scoring")
examples_report = inspect_table(examples_df, "training_edge_examples")

column_report = pd.concat([candidate_report, examples_report], ignore_index=True)
column_report.to_csv(REPORT_OUT, index=False)

print("Column category counts:")
print(
    column_report
    .groupby(["table", "category"])
    .size()
    .rename("n_columns")
    .reset_index()
    .sort_values(["table", "category"])
)

print()
print("Feature-like columns detected:")
feature_like = column_report[
    ~column_report["category"].isin(["core_id_column", "label_or_source_column", "other"])
].copy()

print(feature_like[["table", "column", "dtype", "category", "example_values"]].to_string(index=False))

print()
print(f"[51G.0] Column inspection saved to: {REPORT_OUT}")





"""Βήμα 51G — Build node and edge feature matrices for pseudo-supervised GNN"""

# ============================================================
# Imports
# ============================================================

from pathlib import Path
import json
import re
import numpy as np
import pandas as pd


# ============================================================
# Paths
# ============================================================

PROJECT_DIR = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
GNN_DIR = PROJECT_DIR / "data" / "processed" / "gnn"
GNN_DIR.mkdir(parents=True, exist_ok=True)

CANDIDATE_FILE = GNN_DIR / "pbmc_gnn_candidate_edges_for_scoring.csv"
EXAMPLES_FILE = GNN_DIR / "pbmc_gnn_edge_examples.csv"

NODE_FEATURE_OUT = GNN_DIR / "pbmc_gnn_node_features.csv"
CANDIDATE_EDGE_FEATURE_OUT = GNN_DIR / "pbmc_gnn_edge_features_candidates.csv"
TRAINING_EDGE_FEATURE_OUT = GNN_DIR / "pbmc_gnn_edge_features_training_examples.csv"
FEATURE_MANIFEST_OUT = GNN_DIR / "pbmc_gnn_step51_feature_manifest.json"


# ============================================================
# Biological reference lists
# ============================================================

KNOWN_IMMUNE_TFS = {
    "SPI1", "RUNX1", "RUNX2",
    "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF7", "IRF8",
    "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
    "JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2",
    "NFKB1", "NFKB2", "RELA", "REL", "RELB",
    "IKZF1", "IKZF2", "IKZF3",
    "ETS1", "ELF1", "ELK1",
    "KLF2", "KLF3", "KLF6", "KLF12",
    "TCF7", "LEF1", "GATA3", "BCL11B",
    "BATF", "BATF3",
    "NFATC1", "NFATC2", "NFATC3",
    "EOMES", "TBX21", "PRDM1", "BACH2",
    "EGR1", "EGR2", "EGR3", "EGR4",
    "NR4A1", "NR4A2", "NR4A3",
    "CEBPA", "CEBPB", "CEBPD",
}

BROAD_MOTIF_TFS = {
    "MAZ", "THAP11", "ARNT", "MGA", "ZNF384", "CTCF",
    "SP1", "SP2", "SP3", "SP4",
    "YY1", "MAX", "MXI1", "USF1", "USF2",
    "E2F1", "E2F2", "E2F3", "E2F4",
}

MARKER_GENES_BY_CELL_TYPE = {
    "CD14+ Monocytes": {
        "FCN1", "CD14", "LYZ", "S100A8", "S100A9", "MS4A7", "CTSS", "LST1", "LGALS3",
    },
    "CD4.Memory": {
        "IL7R", "CD4", "LTB", "CCR7", "ANXA1", "S100A4", "CD40LG", "TRAC",
    },
    "CD4.Naive": {
        "CCR7", "LEF1", "TCF7", "IL7R", "CD4", "SELL", "LTB", "MAL",
    },
    "CD8.effector": {
        "NKG7", "PRF1", "GZMB", "CCL5", "GNLY", "CD8A", "CD8B", "GZMH", "CST7",
    },
    "CD8.Naive": {
        "CD8A", "CD8B", "CCR7", "LEF1", "TCF7", "SELL", "IL7R", "LTB",
    },
}


# ============================================================
# Helper functions
# ============================================================

def normalize_colname(x):
    return re.sub(r"[^a-z0-9]+", "_", str(x).strip().lower()).strip("_")


def safe_numeric(df, col, default=0.0):
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def log1p_abs(series):
    return np.log1p(pd.to_numeric(series, errors="coerce").fillna(0).abs())


def rank01_by_cell_type(df, col):
    return (
        df
        .groupby("cell_type")[col]
        .rank(method="average", pct=True)
        .fillna(0)
    )


def compute_gene_specificity(long_df, gene_col, value_col):
    mat = long_df.pivot_table(
        index=gene_col,
        columns="cell_type",
        values=value_col,
        aggfunc="mean",
        fill_value=0,
    )

    mat = mat.clip(lower=0)
    denom = mat.sum(axis=1).replace(0, np.nan)
    specificity = (mat.max(axis=1) / denom).fillna(0)

    return specificity.rename("specificity").reset_index()


# ============================================================
# Load Step 50G outputs
# ============================================================

candidate_df = pd.read_csv(CANDIDATE_FILE)
examples_df = pd.read_csv(EXAMPLES_FILE)

for df_name, df in {"candidate_df": candidate_df, "examples_df": examples_df}.items():
    required = ["cell_type", "tf", "target_gene", "edge_id"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{df_name} is missing required columns: {missing}")

candidate_df["tf"] = candidate_df["tf"].astype(str)
candidate_df["target_gene"] = candidate_df["target_gene"].astype(str)
candidate_df["cell_type"] = candidate_df["cell_type"].astype(str)

examples_df["tf"] = examples_df["tf"].astype(str)
examples_df["target_gene"] = examples_df["target_gene"].astype(str)
examples_df["cell_type"] = examples_df["cell_type"].astype(str)


# ============================================================
# Edge feature construction
# ============================================================

def build_edge_features(df, table_name):
    out_cols = ["edge_id", "cell_type", "tf", "target_gene"]

    if "pseudo_label" in df.columns:
        out_cols.append("pseudo_label")

    out = df[out_cols].copy()

    # Direct biological evidence features
    out["ef_cell_type_peak_mean_accessibility"] = safe_numeric(df, "cell_type_peak_mean_accessibility")
    out["ef_cell_type_peak_detection_fraction"] = safe_numeric(df, "cell_type_peak_detection_fraction")
    out["ef_peak_max_mean_accessibility"] = safe_numeric(df, "peak_max_mean_accessibility")
    out["ef_peak_max_detection_fraction"] = safe_numeric(df, "peak_max_detection_fraction")

    out["ef_tf_rna_expression"] = safe_numeric(df, "tf_rna_expression")
    out["ef_target_gene_atac_gene_activity"] = safe_numeric(df, "target_gene_atac_gene_activity")
    out["ef_target_rna_atac_pseudobulk_correlation"] = safe_numeric(df, "target_gene_rna_atac_pseudobulk_correlation")

    out["ef_cell_type_edge_evidence_score"] = safe_numeric(df, "cell_type_edge_evidence_score")
    out["ef_edge_evidence_score_for_labeling"] = safe_numeric(df, "edge_evidence_score_for_labeling")
    out["ef_edge_evidence_rank_for_labeling"] = safe_numeric(df, "edge_evidence_rank_for_labeling")

    # Distance / regulatory region features
    out["ef_log1p_abs_distance_to_gene"] = log1p_abs(safe_numeric(df, "abs_distance_to_gene"))
    out["ef_is_promoter_peak"] = safe_numeric(df, "is_promoter_peak")
    out["ef_is_distal_peak"] = safe_numeric(df, "is_distal_peak")

    # Source / label-derived indicator, useful mainly to keep shuffled negatives recognizable
    out["ef_is_shuffled_negative"] = safe_numeric(df, "is_shuffled_negative")

    # TF and target biological priors
    out["ef_tf_is_known_immune_tf"] = out["tf"].isin(KNOWN_IMMUNE_TFS).astype(float)
    out["ef_tf_is_broad_motif_candidate"] = out["tf"].isin(BROAD_MOTIF_TFS).astype(float)

    out["ef_target_is_marker_gene"] = out.apply(
        lambda r: float(r["target_gene"] in MARKER_GENES_BY_CELL_TYPE.get(r["cell_type"], set())),
        axis=1,
    )

    # Peak annotation one-hot features
    peak_annotation = df["peak_annotation_type"].fillna("unknown").astype(str).str.lower() if "peak_annotation_type" in df.columns else pd.Series("unknown", index=df.index)

    out["ef_peak_annotation_promoter"] = peak_annotation.eq("promoter").astype(float)
    out["ef_peak_annotation_distal"] = peak_annotation.eq("distal").astype(float)
    out["ef_peak_annotation_other"] = (~peak_annotation.isin(["promoter", "distal"])).astype(float)

    # Cell-type match features: whether top peak accessibility/detection matches current cell type
    if "peak_top_accessibility_cell_type" in df.columns:
        out["ef_peak_top_accessibility_matches_cell_type"] = (
            df["peak_top_accessibility_cell_type"].astype(str) == df["cell_type"].astype(str)
        ).astype(float)
    else:
        out["ef_peak_top_accessibility_matches_cell_type"] = 0.0

    if "peak_top_detection_cell_type" in df.columns:
        out["ef_peak_top_detection_matches_cell_type"] = (
            df["peak_top_detection_cell_type"].astype(str) == df["cell_type"].astype(str)
        ).astype(float)
    else:
        out["ef_peak_top_detection_matches_cell_type"] = 0.0

    # Rank-normalized versions within cell type
    for col in [
        "ef_cell_type_peak_mean_accessibility",
        "ef_cell_type_peak_detection_fraction",
        "ef_tf_rna_expression",
        "ef_target_gene_atac_gene_activity",
        "ef_cell_type_edge_evidence_score",
        "ef_edge_evidence_score_for_labeling",
    ]:
        out[f"{col}_rank01"] = rank01_by_cell_type(out, col)

    # Fill numeric feature columns
    ef_cols = [c for c in out.columns if c.startswith("ef_")]
    out[ef_cols] = out[ef_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    print(f"[51G] {table_name}: {out.shape[0]:,} rows, {len(ef_cols)} edge features")
    return out


candidate_edge_features = build_edge_features(candidate_df, "candidate edge features")
training_edge_features = build_edge_features(examples_df, "training edge features")


# ============================================================
# TF-level graph statistics
# ============================================================

tf_stats = (
    candidate_df
    .groupby(["cell_type", "tf"])
    .agg(
        nf_tf_candidate_out_degree=("target_gene", "nunique"),
        nf_tf_candidate_edge_count=("target_gene", "size"),
        nf_tf_mean_edge_evidence=("edge_evidence_score_for_labeling", "mean"),
        nf_tf_mean_accessibility=("cell_type_peak_mean_accessibility", "mean"),
        nf_tf_mean_detection=("cell_type_peak_detection_fraction", "mean"),
        nf_tf_mean_expression=("tf_rna_expression", "mean"),
    )
    .reset_index()
    .rename(columns={"tf": "node_id"})
)

target_stats = (
    candidate_df
    .groupby(["cell_type", "target_gene"])
    .agg(
        nf_target_candidate_in_degree=("tf", "nunique"),
        nf_target_candidate_edge_count=("tf", "size"),
        nf_target_mean_edge_evidence=("edge_evidence_score_for_labeling", "mean"),
        nf_target_mean_atac_gene_activity=("target_gene_atac_gene_activity", "mean"),
        nf_target_mean_rna_atac_correlation=("target_gene_rna_atac_pseudobulk_correlation", "mean"),
    )
    .reset_index()
    .rename(columns={"target_gene": "node_id"})
)

total_targets_by_cell_type = candidate_df.groupby("cell_type")["target_gene"].nunique().to_dict()

tf_stats["nf_tf_motif_broadness"] = tf_stats.apply(
    lambda r: r["nf_tf_candidate_out_degree"] / max(total_targets_by_cell_type.get(r["cell_type"], 1), 1),
    axis=1,
)

tf_celltype_breadth = (
    candidate_df
    .groupby("tf")["cell_type"]
    .nunique()
    .rename("nf_tf_celltype_breadth")
    .reset_index()
    .rename(columns={"tf": "node_id"})
)

n_cell_types = max(candidate_df["cell_type"].nunique(), 1)
tf_celltype_breadth["nf_tf_celltype_breadth"] = tf_celltype_breadth["nf_tf_celltype_breadth"] / n_cell_types


# ============================================================
# Node feature construction
# ============================================================

tf_nodes = candidate_df[["cell_type", "tf"]].rename(columns={"tf": "node_id"}).drop_duplicates()
tf_nodes["nf_is_tf"] = 1.0
tf_nodes["nf_is_target_gene"] = 0.0

target_nodes = candidate_df[["cell_type", "target_gene"]].rename(columns={"target_gene": "node_id"}).drop_duplicates()
target_nodes["nf_is_tf"] = 0.0
target_nodes["nf_is_target_gene"] = 1.0

nodes = pd.concat([tf_nodes, target_nodes], ignore_index=True)

nodes = (
    nodes
    .groupby(["cell_type", "node_id"], as_index=False)
    .agg(
        nf_is_tf=("nf_is_tf", "max"),
        nf_is_target_gene=("nf_is_target_gene", "max"),
    )
)

nodes = nodes.merge(tf_stats, on=["cell_type", "node_id"], how="left")
nodes = nodes.merge(target_stats, on=["cell_type", "node_id"], how="left")
nodes = nodes.merge(tf_celltype_breadth, on="node_id", how="left")

nodes["nf_is_known_immune_tf"] = nodes["node_id"].isin(KNOWN_IMMUNE_TFS).astype(float)
nodes["nf_is_broad_motif_candidate"] = nodes["node_id"].isin(BROAD_MOTIF_TFS).astype(float)

nodes["nf_is_marker_gene"] = nodes.apply(
    lambda r: float(r["node_id"] in MARKER_GENES_BY_CELL_TYPE.get(r["cell_type"], set())),
    axis=1,
)

# Combined expression/activity proxy
nodes["nf_expression_or_activity_proxy"] = nodes[
    [
        "nf_tf_mean_expression",
        "nf_target_mean_atac_gene_activity",
    ]
].max(axis=1, skipna=True)

nodes["nf_total_degree"] = (
    nodes["nf_tf_candidate_out_degree"].fillna(0)
    + nodes["nf_target_candidate_in_degree"].fillna(0)
)

nodes["nf_log1p_total_degree"] = np.log1p(nodes["nf_total_degree"])
nodes["nf_log1p_tf_out_degree"] = np.log1p(nodes["nf_tf_candidate_out_degree"].fillna(0))
nodes["nf_log1p_target_in_degree"] = np.log1p(nodes["nf_target_candidate_in_degree"].fillna(0))

# Expression/activity specificity across cell types
specificity_df = compute_gene_specificity(
    nodes[["cell_type", "node_id", "nf_expression_or_activity_proxy"]].copy(),
    gene_col="node_id",
    value_col="nf_expression_or_activity_proxy",
).rename(columns={"specificity": "nf_expression_or_activity_specificity"})

nodes = nodes.merge(specificity_df, on="node_id", how="left")

nf_cols = [c for c in nodes.columns if c.startswith("nf_")]
nodes[nf_cols] = nodes[nf_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

print(f"[51G] node features: {nodes.shape[0]:,} rows, {len(nf_cols)} node features")


# ============================================================
# Save outputs
# ============================================================

candidate_edge_features.to_csv(CANDIDATE_EDGE_FEATURE_OUT, index=False)
training_edge_features.to_csv(TRAINING_EDGE_FEATURE_OUT, index=False)
nodes.to_csv(NODE_FEATURE_OUT, index=False)

edge_feature_cols = [c for c in candidate_edge_features.columns if c.startswith("ef_")]
node_feature_cols = [c for c in nodes.columns if c.startswith("nf_")]

manifest = {
    "node_feature_file": str(NODE_FEATURE_OUT),
    "candidate_edge_feature_file": str(CANDIDATE_EDGE_FEATURE_OUT),
    "training_edge_feature_file": str(TRAINING_EDGE_FEATURE_OUT),
    "node_feature_cols": node_feature_cols,
    "edge_feature_cols": edge_feature_cols,
    "strategy": {
        "strategy_1": "evidence_based_edges_from_step49R",
        "strategy_2": "pseudo_supervised_tabular_edge_features",
        "strategy_3": "gnn_refined_edge_classifier_inputs",
    },
}

with open(FEATURE_MANIFEST_OUT, "w", encoding="utf-8") as f:
    json.dump(manifest, f, indent=2)

print()
print("[51G] Done.")
print(f"[51G] Node features: {NODE_FEATURE_OUT}")
print(f"[51G] Candidate edge features: {CANDIDATE_EDGE_FEATURE_OUT}")
print(f"[51G] Training edge features: {TRAINING_EDGE_FEATURE_OUT}")
print(f"[51G] Feature manifest: {FEATURE_MANIFEST_OUT}")

print()
print("[51G] Candidate edge feature shape:", candidate_edge_features.shape)
print("[51G] Training edge feature shape:", training_edge_features.shape)
print("[51G] Node feature shape:", nodes.shape)

print()
print("[51G] Training label check:")
print(
    training_edge_features
    .groupby(["cell_type", "pseudo_label"])
    .size()
    .rename("n_edges")
    .reset_index()
)


"""Βήμα 52G — Train per-cell-type pseudo-supervised GNN edge classifier"""

# ============================================================
# Imports
# ============================================================

from pathlib import Path
import json
import re
import random
import warnings

import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score
from sklearn.linear_model import LogisticRegression
import joblib

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
except ImportError as e:
    raise ImportError(
        "PyTorch is required for Step 52G. Install PyTorch before running this step."
    ) from e


# ============================================================
# Paths and configuration
# ============================================================

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)

PROJECT_DIR = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
GNN_DIR = PROJECT_DIR / "data" / "processed" / "gnn"
MODEL_DIR = PROJECT_DIR / "models" / "gnn"
MODEL_DIR.mkdir(parents=True, exist_ok=True)

NODE_FEATURE_FILE = GNN_DIR / "pbmc_gnn_node_features.csv"
CANDIDATE_EDGE_FEATURE_FILE = GNN_DIR / "pbmc_gnn_edge_features_candidates.csv"
TRAINING_EDGE_FEATURE_FILE = GNN_DIR / "pbmc_gnn_edge_features_training_examples.csv"
FEATURE_MANIFEST_FILE = GNN_DIR / "pbmc_gnn_step51_feature_manifest.json"

METRICS_OUT = GNN_DIR / "pbmc_gnn_step52_training_metrics.csv"
USED_FEATURES_OUT = GNN_DIR / "pbmc_gnn_step52_used_features.json"

MAX_EPOCHS = 250
PATIENCE = 30
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 1e-4
HIDDEN_DIM = 64
DROPOUT = 0.25

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

print(f"[52G] Device: {DEVICE}")


# ============================================================
# Leakage-aware feature filtering
# ============================================================

EDGE_FEATURES_TO_EXCLUDE = {
    # direct label/rank leakage from Step 50G labeling
    "ef_edge_evidence_score_for_labeling",
    "ef_edge_evidence_rank_for_labeling",
    "ef_edge_evidence_score_for_labeling_rank01",

    # source leakage: this tells model which rows are artificial shuffled negatives
    "ef_is_shuffled_negative",
}

NODE_FEATURES_TO_EXCLUDE = set()


# ============================================================
# Model definitions
# ============================================================

class GraphSAGELayer(nn.Module):
    def __init__(self, in_dim, out_dim):
        super().__init__()
        self.linear = nn.Linear(in_dim * 2, out_dim)

    def forward(self, x, edge_index):
        # edge_index[0] = destination, edge_index[1] = source
        dst, src = edge_index

        agg = torch.zeros_like(x)
        deg = torch.zeros((x.shape[0], 1), device=x.device)

        agg.index_add_(0, dst, x[src])
        deg.index_add_(0, dst, torch.ones((src.shape[0], 1), device=x.device))

        agg = agg / deg.clamp(min=1.0)

        return self.linear(torch.cat([x, agg], dim=1))


class EdgeGNN(nn.Module):
    def __init__(self, node_dim, edge_dim, hidden_dim=64, dropout=0.25):
        super().__init__()

        self.layer1 = GraphSAGELayer(node_dim, hidden_dim)
        self.layer2 = GraphSAGELayer(hidden_dim, hidden_dim)
        self.dropout = nn.Dropout(dropout)

        decoder_in = hidden_dim * 4 + edge_dim

        self.decoder = nn.Sequential(
            nn.Linear(decoder_in, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1),
        )

    def encode(self, x, edge_index):
        h = self.layer1(x, edge_index)
        h = F.relu(h)
        h = self.dropout(h)

        h = self.layer2(h, edge_index)
        h = F.relu(h)
        h = self.dropout(h)

        return h

    def decode(self, h, pairs, edge_attr):
        src = pairs[:, 0]
        dst = pairs[:, 1]

        h_src = h[src]
        h_dst = h[dst]

        z = torch.cat(
            [
                h_src,
                h_dst,
                torch.abs(h_src - h_dst),
                h_src * h_dst,
                edge_attr,
            ],
            dim=1,
        )

        return self.decoder(z).squeeze(1)

    def forward(self, x, edge_index, pairs, edge_attr):
        h = self.encode(x, edge_index)
        logits = self.decode(h, pairs, edge_attr)
        return logits


# ============================================================
# Helper functions
# ============================================================

def safe_cell_name(x):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x))


def safe_metric(metric_func, y_true, y_score):
    try:
        if len(np.unique(y_true)) < 2:
            return np.nan
        return metric_func(y_true, y_score)
    except Exception:
        return np.nan


def build_node_index(nodes_ct, candidate_ct, examples_ct):
    all_nodes = sorted(
        set(nodes_ct["node_id"].astype(str))
        | set(candidate_ct["tf"].astype(str))
        | set(candidate_ct["target_gene"].astype(str))
        | set(examples_ct["tf"].astype(str))
        | set(examples_ct["target_gene"].astype(str))
    )

    return {node: idx for idx, node in enumerate(all_nodes)}


def matrix_for_nodes(nodes_ct, node_index, node_feature_cols):
    ordered_nodes = [node for node, idx in sorted(node_index.items(), key=lambda x: x[1])]

    base = pd.DataFrame({"node_id": ordered_nodes})
    merged = base.merge(nodes_ct, on="node_id", how="left")

    for col in node_feature_cols:
        if col not in merged.columns:
            merged[col] = 0.0

    x = (
        merged[node_feature_cols]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0)
        .values
    )

    return x


def pairs_from_edges(df, node_index):
    src = df["tf"].astype(str).map(node_index)
    dst = df["target_gene"].astype(str).map(node_index)

    valid = src.notna() & dst.notna()

    pairs = np.vstack(
        [
            src[valid].astype(int).values,
            dst[valid].astype(int).values,
        ]
    ).T

    return pairs, valid.values


def make_message_passing_edge_index(candidate_ct, node_index):
    pairs, valid = pairs_from_edges(candidate_ct, node_index)

    if pairs.shape[0] == 0:
        raise ValueError("No valid candidate edges for message passing.")

    src = pairs[:, 0]
    dst = pairs[:, 1]

    # Use undirected message passing over directed biological edges.
    # Edge prediction remains directional via TF -> target pairs.
    mp_dst = np.concatenate([dst, src])
    mp_src = np.concatenate([src, dst])

    edge_index = np.vstack([mp_dst, mp_src])

    return edge_index


def evaluate_gnn(model, x_tensor, edge_index_tensor, pairs_tensor, edge_attr_tensor, y_np, indices):
    model.eval()

    with torch.no_grad():
        logits = model(
            x_tensor,
            edge_index_tensor,
            pairs_tensor[indices],
            edge_attr_tensor[indices],
        )
        probs = torch.sigmoid(logits).detach().cpu().numpy()

    y_true = y_np[indices]
    preds = (probs >= 0.5).astype(int)

    return {
        "auroc": safe_metric(roc_auc_score, y_true, probs),
        "aupr": safe_metric(average_precision_score, y_true, probs),
        "accuracy": accuracy_score(y_true, preds),
    }


def evaluate_tabular(model, x_np, y_np, indices):
    probs = model.predict_proba(x_np[indices])[:, 1]
    y_true = y_np[indices]
    preds = (probs >= 0.5).astype(int)

    return {
        "auroc": safe_metric(roc_auc_score, y_true, probs),
        "aupr": safe_metric(average_precision_score, y_true, probs),
        "accuracy": accuracy_score(y_true, preds),
    }


# ============================================================
# Load data
# ============================================================

with open(FEATURE_MANIFEST_FILE, "r", encoding="utf-8") as f:
    feature_manifest = json.load(f)

node_feature_cols_all = feature_manifest["node_feature_cols"]
edge_feature_cols_all = feature_manifest["edge_feature_cols"]

node_feature_cols = [
    c for c in node_feature_cols_all
    if c not in NODE_FEATURES_TO_EXCLUDE
]

edge_feature_cols = [
    c for c in edge_feature_cols_all
    if c not in EDGE_FEATURES_TO_EXCLUDE
]

nodes = pd.read_csv(NODE_FEATURE_FILE)
candidate_edges = pd.read_csv(CANDIDATE_EDGE_FEATURE_FILE)
training_edges = pd.read_csv(TRAINING_EDGE_FEATURE_FILE)

print(f"[52G] Node features loaded: {nodes.shape}")
print(f"[52G] Candidate edge features loaded: {candidate_edges.shape}")
print(f"[52G] Training edge features loaded: {training_edges.shape}")

print(f"[52G] Node features used: {len(node_feature_cols)}")
print(f"[52G] Edge features used: {len(edge_feature_cols)}")

print("[52G] Excluded edge features:")
print(sorted(set(edge_feature_cols_all) - set(edge_feature_cols)))


# ============================================================
# Save used feature list
# ============================================================

with open(USED_FEATURES_OUT, "w", encoding="utf-8") as f:
    json.dump(
        {
            "node_feature_cols": node_feature_cols,
            "edge_feature_cols": edge_feature_cols,
            "excluded_edge_features": sorted(set(edge_feature_cols_all) - set(edge_feature_cols)),
            "excluded_node_features": sorted(set(node_feature_cols_all) - set(node_feature_cols)),
        },
        f,
        indent=2,
    )


# ============================================================
# Train models per cell type
# ============================================================

metrics = []

for cell_type in sorted(training_edges["cell_type"].unique()):
    print()
    print("=" * 80)
    print(f"[52G] Training cell type: {cell_type}")
    print("=" * 80)

    nodes_ct = nodes[nodes["cell_type"] == cell_type].copy()
    candidate_ct = candidate_edges[candidate_edges["cell_type"] == cell_type].copy()
    examples_ct = training_edges[training_edges["cell_type"] == cell_type].copy()

    examples_ct = examples_ct[examples_ct["pseudo_label"].isin([0, 1])].copy()
    examples_ct["pseudo_label"] = examples_ct["pseudo_label"].astype(int)

    print("[52G] Raw class counts:")
    print(examples_ct["pseudo_label"].value_counts().sort_index())

    if examples_ct.shape[0] < 100 or examples_ct["pseudo_label"].nunique() < 2:
        warnings.warn(f"[52G] Skipping {cell_type}: insufficient labeled examples.")
        continue

    # ============================================================
    # Node and edge matrices
    # ============================================================

    node_index = build_node_index(nodes_ct, candidate_ct, examples_ct)

    x_raw = matrix_for_nodes(
        nodes_ct[["node_id"] + node_feature_cols],
        node_index,
        node_feature_cols,
    )

    node_scaler = StandardScaler()
    x_scaled = node_scaler.fit_transform(x_raw)

    for col in edge_feature_cols:
        if col not in examples_ct.columns:
            examples_ct[col] = 0.0
        if col not in candidate_ct.columns:
            candidate_ct[col] = 0.0

    edge_raw = (
        examples_ct[edge_feature_cols]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0)
        .values
    )

    edge_scaler = StandardScaler()
    edge_scaled = edge_scaler.fit_transform(edge_raw)

    pairs, valid = pairs_from_edges(examples_ct, node_index)

    examples_ct = examples_ct.loc[valid].reset_index(drop=True)
    edge_scaled = edge_scaled[valid]
    y = examples_ct["pseudo_label"].values.astype(int)

    if len(np.unique(y)) < 2:
        warnings.warn(f"[52G] Skipping {cell_type}: only one class after valid-pair filtering.")
        continue

    mp_edge_index = make_message_passing_edge_index(candidate_ct, node_index)

    print(f"[52G] Nodes: {len(node_index):,}")
    print(f"[52G] Message-passing edges, directed both ways: {mp_edge_index.shape[1]:,}")
    print(f"[52G] Labeled edge examples: {len(y):,}")
    print("[52G] Final class counts:")
    print(pd.Series(y).value_counts().sort_index())

    # ============================================================
    # Train / validation / test split
    # ============================================================

    indices = np.arange(len(y))

    train_idx, temp_idx = train_test_split(
        indices,
        test_size=0.30,
        random_state=SEED,
        stratify=y,
    )

    val_idx, test_idx = train_test_split(
        temp_idx,
        test_size=0.50,
        random_state=SEED,
        stratify=y[temp_idx],
    )

    # ============================================================
    # Tensor conversion
    # ============================================================

    x_tensor = torch.tensor(x_scaled, dtype=torch.float32, device=DEVICE)
    edge_index_tensor = torch.tensor(mp_edge_index, dtype=torch.long, device=DEVICE)
    pairs_tensor = torch.tensor(pairs, dtype=torch.long, device=DEVICE)
    edge_attr_tensor = torch.tensor(edge_scaled, dtype=torch.float32, device=DEVICE)
    y_tensor = torch.tensor(y, dtype=torch.float32, device=DEVICE)

    # ============================================================
    # GNN training
    # ============================================================

    n_pos = max(float((y[train_idx] == 1).sum()), 1.0)
    n_neg = max(float((y[train_idx] == 0).sum()), 1.0)

    pos_weight = torch.tensor([n_neg / n_pos], dtype=torch.float32, device=DEVICE)

    model = EdgeGNN(
        node_dim=len(node_feature_cols),
        edge_dim=len(edge_feature_cols),
        hidden_dim=HIDDEN_DIM,
        dropout=DROPOUT,
    ).to(DEVICE)

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=LEARNING_RATE,
        weight_decay=WEIGHT_DECAY,
    )

    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)

    best_val_aupr = -np.inf
    best_state = None
    patience_counter = 0

    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        optimizer.zero_grad()

        logits = model(
            x_tensor,
            edge_index_tensor,
            pairs_tensor[train_idx],
            edge_attr_tensor[train_idx],
        )

        loss = criterion(logits, y_tensor[train_idx])
        loss.backward()
        optimizer.step()

        val_metrics = evaluate_gnn(
            model,
            x_tensor,
            edge_index_tensor,
            pairs_tensor,
            edge_attr_tensor,
            y,
            val_idx,
        )

        val_aupr = val_metrics["aupr"]
        monitor_score = val_aupr if not np.isnan(val_aupr) else -np.inf

        if monitor_score > best_val_aupr:
            best_val_aupr = monitor_score
            best_state = {
                k: v.detach().cpu().clone()
                for k, v in model.state_dict().items()
            }
            patience_counter = 0
        else:
            patience_counter += 1

        if epoch == 1 or epoch % 25 == 0:
            print(
                f"[52G] epoch={epoch:03d} "
                f"loss={loss.item():.4f} "
                f"val_AUPR={val_metrics['aupr']:.4f} "
                f"val_AUROC={val_metrics['auroc']:.4f} "
                f"val_ACC={val_metrics['accuracy']:.4f}"
            )

        if patience_counter >= PATIENCE:
            print(f"[52G] Early stopping at epoch {epoch}")
            break

    if best_state is not None:
        model.load_state_dict(best_state)

    gnn_val_metrics = evaluate_gnn(
        model,
        x_tensor,
        edge_index_tensor,
        pairs_tensor,
        edge_attr_tensor,
        y,
        val_idx,
    )

    gnn_test_metrics = evaluate_gnn(
        model,
        x_tensor,
        edge_index_tensor,
        pairs_tensor,
        edge_attr_tensor,
        y,
        test_idx,
    )

    # ============================================================
    # Strategy 2 baseline: pseudo-supervised tabular classifier
    # ============================================================

    tabular_model = LogisticRegression(
        max_iter=3000,
        class_weight="balanced",
        random_state=SEED,
        n_jobs=None,
    )

    tabular_model.fit(edge_scaled[train_idx], y[train_idx])

    tabular_val_metrics = evaluate_tabular(
        tabular_model,
        edge_scaled,
        y,
        val_idx,
    )

    tabular_test_metrics = evaluate_tabular(
        tabular_model,
        edge_scaled,
        y,
        test_idx,
    )

    # ============================================================
    # Save model and preprocessing bundle
    # ============================================================

    cell_safe = safe_cell_name(cell_type)

    checkpoint_path = MODEL_DIR / f"edge_gnn_{cell_safe}.pt"
    bundle_path = MODEL_DIR / f"edge_gnn_{cell_safe}_preprocessing.joblib"

    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "node_dim": len(node_feature_cols),
            "edge_dim": len(edge_feature_cols),
            "hidden_dim": HIDDEN_DIM,
            "dropout": DROPOUT,
            "cell_type": cell_type,
            "node_feature_cols": node_feature_cols,
            "edge_feature_cols": edge_feature_cols,
        },
        checkpoint_path,
    )

    joblib.dump(
        {
            "node_scaler": node_scaler,
            "edge_scaler": edge_scaler,
            "tabular_model": tabular_model,
            "node_feature_cols": node_feature_cols,
            "edge_feature_cols": edge_feature_cols,
            "node_index": node_index,
            "checkpoint_path": str(checkpoint_path),
            "cell_type": cell_type,
        },
        bundle_path,
    )

    # ============================================================
    # Collect metrics
    # ============================================================

    metrics.extend([
        {
            "cell_type": cell_type,
            "model": "gnn_refined",
            "split": "validation",
            **gnn_val_metrics,
            "n_examples": len(y),
            "n_train": len(train_idx),
            "n_validation": len(val_idx),
            "n_test": len(test_idx),
            "checkpoint": str(checkpoint_path),
        },
        {
            "cell_type": cell_type,
            "model": "gnn_refined",
            "split": "test",
            **gnn_test_metrics,
            "n_examples": len(y),
            "n_train": len(train_idx),
            "n_validation": len(val_idx),
            "n_test": len(test_idx),
            "checkpoint": str(checkpoint_path),
        },
        {
            "cell_type": cell_type,
            "model": "pseudo_supervised_tabular",
            "split": "validation",
            **tabular_val_metrics,
            "n_examples": len(y),
            "n_train": len(train_idx),
            "n_validation": len(val_idx),
            "n_test": len(test_idx),
            "checkpoint": str(bundle_path),
        },
        {
            "cell_type": cell_type,
            "model": "pseudo_supervised_tabular",
            "split": "test",
            **tabular_test_metrics,
            "n_examples": len(y),
            "n_train": len(train_idx),
            "n_validation": len(val_idx),
            "n_test": len(test_idx),
            "checkpoint": str(bundle_path),
        },
    ])

    print()
    print(f"[52G] Saved GNN checkpoint: {checkpoint_path}")
    print(f"[52G] Saved preprocessing bundle: {bundle_path}")
    print(f"[52G] GNN validation metrics: {gnn_val_metrics}")
    print(f"[52G] GNN test metrics: {gnn_test_metrics}")
    print(f"[52G] Tabular validation metrics: {tabular_val_metrics}")
    print(f"[52G] Tabular test metrics: {tabular_test_metrics}")


# ============================================================
# Save final metrics
# ============================================================

metrics_df = pd.DataFrame(metrics)
metrics_df.to_csv(METRICS_OUT, index=False)

print()
print("[52G] Done.")
print(f"[52G] Training metrics saved to: {METRICS_OUT}")
print(f"[52G] Used features saved to: {USED_FEATURES_OUT}")

print()
print("[52G] Metrics summary:")
print(metrics_df)





"""Βήμα 52G.1 — Inspect pseudo-supervised GNN training metrics and interpretation checks"""

# ============================================================
# Imports
# ============================================================

from pathlib import Path
import pandas as pd
import numpy as np


# ============================================================
# Paths
# ============================================================

PROJECT_DIR = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
GNN_DIR = PROJECT_DIR / "data" / "processed" / "gnn"

METRICS_FILE = GNN_DIR / "pbmc_gnn_step52_training_metrics.csv"
TRAINING_EDGE_FEATURE_FILE = GNN_DIR / "pbmc_gnn_edge_features_training_examples.csv"


# ============================================================
# Load outputs
# ============================================================

metrics_df = pd.read_csv(METRICS_FILE)
training_edges = pd.read_csv(TRAINING_EDGE_FEATURE_FILE)


# ============================================================
# Display full training metrics
# ============================================================

print("[52G.1] Full metrics table:")
print(
    metrics_df[
        ["cell_type", "model", "split", "auroc", "aupr", "accuracy", "n_examples", "n_train", "n_validation", "n_test"]
    ]
    .sort_values(["cell_type", "model", "split"])
    .to_string(index=False)
)


# ============================================================
# Test metrics only
# ============================================================

test_metrics = metrics_df[metrics_df["split"] == "test"].copy()

print()
print("[52G.1] Test metrics summary:")
print(
    test_metrics[
        ["cell_type", "model", "auroc", "aupr", "accuracy"]
    ]
    .sort_values(["cell_type", "model"])
    .to_string(index=False)
)


# ============================================================
# Label and source composition
# ============================================================

print()
print("[52G.1] Training example composition by label and source:")
print(
    training_edges
    .groupby(["cell_type", "pseudo_label", "ef_is_shuffled_negative"])
    .size()
    .rename("n_edges")
    .reset_index()
    .sort_values(["cell_type", "pseudo_label", "ef_is_shuffled_negative"])
    .to_string(index=False)
)


# ============================================================
# Accuracy interpretation helper
# ============================================================

print()
print("[52G.1] Interpretation:")
print(
    "Accuracy here measures agreement with pseudo-labels, not external biological ground truth. "
    "High AUROC/AUPR/accuracy means the model separates high-confidence evidence-based edges "
    "from weak/shuffled examples very well."
)




"""Βήμα 52G.2 — Strict feature ablation to test pseudo-label leakage"""

# ============================================================
# Imports
# ============================================================

from pathlib import Path
import json
import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score


# ============================================================
# Paths
# ============================================================

PROJECT_DIR = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")
GNN_DIR = PROJECT_DIR / "data" / "processed" / "gnn"

TRAINING_EDGE_FEATURE_FILE = GNN_DIR / "pbmc_gnn_edge_features_training_examples.csv"
FEATURE_MANIFEST_FILE = GNN_DIR / "pbmc_gnn_step51_feature_manifest.json"

ABLATION_OUT = GNN_DIR / "pbmc_gnn_step52G2_strict_ablation_metrics.csv"


# ============================================================
# Load data
# ============================================================

training_edges = pd.read_csv(TRAINING_EDGE_FEATURE_FILE)

with open(FEATURE_MANIFEST_FILE, "r", encoding="utf-8") as f:
    manifest = json.load(f)

edge_feature_cols_all = manifest["edge_feature_cols"]


# ============================================================
# Define feature sets
# ============================================================

label_leakage_or_prior_features = {
    "ef_edge_evidence_score_for_labeling",
    "ef_edge_evidence_rank_for_labeling",
    "ef_edge_evidence_score_for_labeling_rank01",
    "ef_cell_type_edge_evidence_score",
    "ef_cell_type_edge_evidence_score_rank01",
    "ef_is_shuffled_negative",
    "ef_tf_is_known_immune_tf",
    "ef_tf_is_broad_motif_candidate",
    "ef_target_is_marker_gene",
}

strict_edge_features = [
    c for c in edge_feature_cols_all
    if c not in label_leakage_or_prior_features
]

print("[52G.2] Number of all edge features:", len(edge_feature_cols_all))
print("[52G.2] Number of strict edge features:", len(strict_edge_features))

print()
print("[52G.2] Removed features:")
for c in sorted(set(edge_feature_cols_all) - set(strict_edge_features)):
    print(" -", c)

print()
print("[52G.2] Strict features used:")
for c in strict_edge_features:
    print(" -", c)


# ============================================================
# Helper functions
# ============================================================

def safe_metric(metric_func, y_true, y_score):
    if len(np.unique(y_true)) < 2:
        return np.nan
    return metric_func(y_true, y_score)


def evaluate_model(model, x, y, idx):
    probs = model.predict_proba(x[idx])[:, 1]
    preds = (probs >= 0.5).astype(int)

    return {
        "auroc": safe_metric(roc_auc_score, y[idx], probs),
        "aupr": safe_metric(average_precision_score, y[idx], probs),
        "accuracy": accuracy_score(y[idx], preds),
    }


# ============================================================
# Train strict tabular model per cell type
# ============================================================

rows = []

for cell_type in sorted(training_edges["cell_type"].unique()):
    ct_df = training_edges[training_edges["cell_type"] == cell_type].copy()
    ct_df = ct_df[ct_df["pseudo_label"].isin([0, 1])].copy()
    ct_df["pseudo_label"] = ct_df["pseudo_label"].astype(int)

    x = (
        ct_df[strict_edge_features]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0)
        .values
    )

    y = ct_df["pseudo_label"].values.astype(int)

    scaler = StandardScaler()
    x_scaled = scaler.fit_transform(x)

    indices = np.arange(len(y))

    train_idx, temp_idx = train_test_split(
        indices,
        test_size=0.30,
        random_state=42,
        stratify=y,
    )

    val_idx, test_idx = train_test_split(
        temp_idx,
        test_size=0.50,
        random_state=42,
        stratify=y[temp_idx],
    )

    model = LogisticRegression(
        max_iter=3000,
        class_weight="balanced",
        random_state=42,
    )

    model.fit(x_scaled[train_idx], y[train_idx])

    val_metrics = evaluate_model(model, x_scaled, y, val_idx)
    test_metrics = evaluate_model(model, x_scaled, y, test_idx)

    rows.append({
        "cell_type": cell_type,
        "model": "strict_tabular_no_evidence_score_no_priors",
        "split": "validation",
        **val_metrics,
        "n_examples": len(y),
        "n_features": len(strict_edge_features),
    })

    rows.append({
        "cell_type": cell_type,
        "model": "strict_tabular_no_evidence_score_no_priors",
        "split": "test",
        **test_metrics,
        "n_examples": len(y),
        "n_features": len(strict_edge_features),
    })

    print()
    print(f"[52G.2] {cell_type}")
    print("Validation:", val_metrics)
    print("Test:", test_metrics)


# ============================================================
# Save metrics
# ============================================================

ablation_df = pd.DataFrame(rows)
ablation_df.to_csv(ABLATION_OUT, index=False)

print()
print("[52G.2] Done.")
print(f"[52G.2] Strict ablation metrics saved to: {ABLATION_OUT}")

print()
print("[52G.2] Test metrics:")
print(
    ablation_df[ablation_df["split"] == "test"]
    .sort_values("cell_type")
    .to_string(index=False)
)





"""Βήμα 53G — Score candidate TF-target edges with per-cell-type GNN models"""

# ============================================================
# Imports
# ============================================================

from pathlib import Path
import re
import numpy as np
import pandas as pd
import joblib

import torch
import torch.nn as nn
import torch.nn.functional as F


# ============================================================
# Paths and configuration
# ============================================================

PROJECT_DIR = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

GNN_DIR = PROJECT_DIR / "data" / "processed" / "gnn"
MODEL_DIR = PROJECT_DIR / "models" / "gnn"

NODE_FEATURE_FILE = GNN_DIR / "pbmc_gnn_node_features.csv"
CANDIDATE_EDGE_FEATURE_FILE = GNN_DIR / "pbmc_gnn_edge_features_candidates.csv"
RAW_CANDIDATE_FILE = GNN_DIR / "pbmc_gnn_candidate_edges_for_scoring.csv"

SCORED_OUT = GNN_DIR / "pbmc_grn_gnn_refined_edges_scored.csv"
HIGH_CONF_GNN_OUT = GNN_DIR / "pbmc_grn_gnn_refined_edges_high_confidence.csv"
TF_SUMMARY_OUT = GNN_DIR / "pbmc_grn_gnn_refined_tf_summary.csv"

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# Keep approximately the same number of selected GNN-refined edges
# as the Strategy 1 evidence-based high-confidence edges per cell type.
FALLBACK_TOP_FRACTION = 0.10

print(f"[53G] Device: {DEVICE}")


# ============================================================
# Model definitions — must match Step 52G
# ============================================================

class GraphSAGELayer(nn.Module):
    def __init__(self, in_dim, out_dim):
        super().__init__()
        self.linear = nn.Linear(in_dim * 2, out_dim)

    def forward(self, x, edge_index):
        # edge_index[0] = destination, edge_index[1] = source
        dst, src = edge_index

        agg = torch.zeros_like(x)
        deg = torch.zeros((x.shape[0], 1), device=x.device)

        agg.index_add_(0, dst, x[src])
        deg.index_add_(0, dst, torch.ones((src.shape[0], 1), device=x.device))

        agg = agg / deg.clamp(min=1.0)

        return self.linear(torch.cat([x, agg], dim=1))


class EdgeGNN(nn.Module):
    def __init__(self, node_dim, edge_dim, hidden_dim=64, dropout=0.25):
        super().__init__()

        self.layer1 = GraphSAGELayer(node_dim, hidden_dim)
        self.layer2 = GraphSAGELayer(hidden_dim, hidden_dim)
        self.dropout = nn.Dropout(dropout)

        decoder_in = hidden_dim * 4 + edge_dim

        self.decoder = nn.Sequential(
            nn.Linear(decoder_in, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1),
        )

    def encode(self, x, edge_index):
        h = self.layer1(x, edge_index)
        h = F.relu(h)
        h = self.dropout(h)

        h = self.layer2(h, edge_index)
        h = F.relu(h)
        h = self.dropout(h)

        return h

    def decode(self, h, pairs, edge_attr):
        src = pairs[:, 0]
        dst = pairs[:, 1]

        h_src = h[src]
        h_dst = h[dst]

        z = torch.cat(
            [
                h_src,
                h_dst,
                torch.abs(h_src - h_dst),
                h_src * h_dst,
                edge_attr,
            ],
            dim=1,
        )

        return self.decoder(z).squeeze(1)

    def forward(self, x, edge_index, pairs, edge_attr):
        h = self.encode(x, edge_index)
        logits = self.decode(h, pairs, edge_attr)
        return logits


# ============================================================
# Helper functions
# ============================================================

def safe_cell_name(x):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x))


def rank01(series):
    return (
        pd.to_numeric(series, errors="coerce")
        .fillna(0)
        .rank(method="average", pct=True)
        .fillna(0)
    )


def pairs_from_edges(df, node_index):
    src = df["tf"].astype(str).map(node_index)
    dst = df["target_gene"].astype(str).map(node_index)

    valid = src.notna() & dst.notna()

    pairs = np.vstack(
        [
            src[valid].astype(int).values,
            dst[valid].astype(int).values,
        ]
    ).T

    return pairs, valid.values


def make_message_passing_edge_index(candidate_ct, node_index):
    pairs, valid = pairs_from_edges(candidate_ct, node_index)

    if pairs.shape[0] == 0:
        raise ValueError("No valid candidate edges for message passing.")

    src = pairs[:, 0]
    dst = pairs[:, 1]

    # Undirected message passing over directed biological candidate edges.
    mp_dst = np.concatenate([dst, src])
    mp_src = np.concatenate([src, dst])

    edge_index = np.vstack([mp_dst, mp_src])

    return edge_index


def matrix_for_saved_node_order(nodes_ct, node_index, node_feature_cols):
    ordered_nodes = [node for node, idx in sorted(node_index.items(), key=lambda x: x[1])]

    base = pd.DataFrame({"node_id": ordered_nodes})
    merged = base.merge(nodes_ct, on="node_id", how="left")

    for col in node_feature_cols:
        if col not in merged.columns:
            merged[col] = 0.0

    x = (
        merged[node_feature_cols]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0)
        .values
    )

    return x


def safe_torch_load(path, device):
    try:
        return torch.load(path, map_location=device, weights_only=False)
    except TypeError:
        return torch.load(path, map_location=device)


# ============================================================
# Load feature tables
# ============================================================

nodes = pd.read_csv(NODE_FEATURE_FILE)
candidate_features = pd.read_csv(CANDIDATE_EDGE_FEATURE_FILE)
raw_candidates = pd.read_csv(RAW_CANDIDATE_FILE)

print("[53G] Loaded:")
print("Node features:", nodes.shape)
print("Candidate edge features:", candidate_features.shape)
print("Raw candidates:", raw_candidates.shape)


# ============================================================
# Score candidate edges per cell type
# ============================================================

scored_tables = []

key_cols = ["edge_id", "cell_type", "tf", "target_gene"]

for cell_type in sorted(candidate_features["cell_type"].unique()):
    print()
    print("=" * 80)
    print(f"[53G] Scoring cell type: {cell_type}")
    print("=" * 80)

    cell_safe = safe_cell_name(cell_type)

    checkpoint_path = MODEL_DIR / f"edge_gnn_{cell_safe}.pt"
    bundle_path = MODEL_DIR / f"edge_gnn_{cell_safe}_preprocessing.joblib"

    if not checkpoint_path.exists():
        print(f"[53G] Missing checkpoint, skipping: {checkpoint_path}")
        continue

    if not bundle_path.exists():
        print(f"[53G] Missing preprocessing bundle, skipping: {bundle_path}")
        continue

    checkpoint = safe_torch_load(checkpoint_path, DEVICE)
    bundle = joblib.load(bundle_path)

    node_feature_cols = bundle["node_feature_cols"]
    edge_feature_cols = bundle["edge_feature_cols"]
    node_index = bundle["node_index"]
    node_scaler = bundle["node_scaler"]
    edge_scaler = bundle["edge_scaler"]
    tabular_model = bundle["tabular_model"]

    nodes_ct = nodes[nodes["cell_type"] == cell_type].copy()
    candidate_ct = candidate_features[candidate_features["cell_type"] == cell_type].copy()
    raw_ct = raw_candidates[raw_candidates["cell_type"] == cell_type].copy()

    print(f"[53G] Candidate edges before node-index filtering: {candidate_ct.shape[0]:,}")

    # ============================================================
    # Matrix construction
    # ============================================================

    for col in edge_feature_cols:
        if col not in candidate_ct.columns:
            candidate_ct[col] = 0.0

    x_raw = matrix_for_saved_node_order(
        nodes_ct=nodes_ct,
        node_index=node_index,
        node_feature_cols=node_feature_cols,
    )

    x_scaled = node_scaler.transform(x_raw)

    edge_raw = (
        candidate_ct[edge_feature_cols]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0)
        .values
    )

    edge_scaled = edge_scaler.transform(edge_raw)

    pairs, valid = pairs_from_edges(candidate_ct, node_index)

    candidate_ct_valid = candidate_ct.loc[valid].reset_index(drop=True)
    edge_scaled_valid = edge_scaled[valid]

    print(f"[53G] Candidate edges after node-index filtering: {candidate_ct_valid.shape[0]:,}")

    if candidate_ct_valid.empty:
        print(f"[53G] No valid candidate edges for {cell_type}, skipping.")
        continue

    mp_edge_index = make_message_passing_edge_index(candidate_ct_valid, node_index)

    # ============================================================
    # Tensor conversion
    # ============================================================

    x_tensor = torch.tensor(x_scaled, dtype=torch.float32, device=DEVICE)
    edge_index_tensor = torch.tensor(mp_edge_index, dtype=torch.long, device=DEVICE)
    pairs_tensor = torch.tensor(pairs, dtype=torch.long, device=DEVICE)
    edge_attr_tensor = torch.tensor(edge_scaled_valid, dtype=torch.float32, device=DEVICE)

    # ============================================================
    # Load and run GNN
    # ============================================================

    model = EdgeGNN(
        node_dim=checkpoint["node_dim"],
        edge_dim=checkpoint["edge_dim"],
        hidden_dim=checkpoint["hidden_dim"],
        dropout=checkpoint["dropout"],
    ).to(DEVICE)

    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    with torch.no_grad():
        logits = model(
            x_tensor,
            edge_index_tensor,
            pairs_tensor,
            edge_attr_tensor,
        )
        gnn_probability = torch.sigmoid(logits).detach().cpu().numpy()

    # ============================================================
    # Strategy 2 tabular pseudo-supervised probability
    # ============================================================

    pseudo_supervised_probability = tabular_model.predict_proba(edge_scaled_valid)[:, 1]

    # ============================================================
    # Merge scores with original candidate annotations
    # ============================================================

    scored_ct = candidate_ct_valid[key_cols].copy()

    scored_ct["gnn_probability"] = gnn_probability
    scored_ct["pseudo_supervised_probability"] = pseudo_supervised_probability

    raw_extra_cols = [c for c in raw_ct.columns if c not in key_cols]

    scored_ct = scored_ct.merge(
        raw_ct[key_cols + raw_extra_cols],
        on=key_cols,
        how="left",
    )

    # Evidence rank from Strategy 1, used only as a small stabilizing component.
    if "edge_evidence_score_for_labeling" in scored_ct.columns:
        scored_ct["evidence_rank_score"] = rank01(scored_ct["edge_evidence_score_for_labeling"])
    elif "cell_type_edge_evidence_score" in scored_ct.columns:
        scored_ct["evidence_rank_score"] = rank01(scored_ct["cell_type_edge_evidence_score"])
    else:
        scored_ct["evidence_rank_score"] = 0.0

    # Composite score:
    # GNN is primary; tabular pseudo-supervised score and original evidence score stabilize ranking.
    scored_ct["gnn_refined_composite_score"] = (
        0.70 * scored_ct["gnn_probability"]
        + 0.20 * scored_ct["pseudo_supervised_probability"]
        + 0.10 * scored_ct["evidence_rank_score"]
    )

    # ============================================================
    # Select GNN high-confidence edges
    # ============================================================

    if "is_pseudo_positive" in scored_ct.columns:
        n_strategy1_edges = int(
            pd.to_numeric(scored_ct["is_pseudo_positive"], errors="coerce")
            .fillna(0)
            .sum()
        )
    else:
        n_strategy1_edges = 0

    if n_strategy1_edges <= 0:
        n_keep = max(1, int(np.ceil(scored_ct.shape[0] * FALLBACK_TOP_FRACTION)))
    else:
        n_keep = min(n_strategy1_edges, scored_ct.shape[0])

    scored_ct["gnn_rank_within_cell_type"] = (
        scored_ct["gnn_refined_composite_score"]
        .rank(ascending=False, method="first")
    )

    scored_ct["is_gnn_high_confidence"] = (
        scored_ct["gnn_rank_within_cell_type"] <= n_keep
    ).astype(int)

    scored_tables.append(scored_ct)

    print(f"[53G] Strategy 1 evidence-based high-confidence edges: {n_strategy1_edges:,}")
    print(f"[53G] GNN high-confidence edges selected: {n_keep:,}")
    print(f"[53G] Mean GNN probability: {scored_ct['gnn_probability'].mean():.4f}")
    print(f"[53G] Mean pseudo-supervised probability: {scored_ct['pseudo_supervised_probability'].mean():.4f}")


# ============================================================
# Combine and save scored edge tables
# ============================================================

if not scored_tables:
    raise RuntimeError("[53G] No scored tables were produced. Check model files and input tables.")

scored = pd.concat(scored_tables, ignore_index=True)

high_conf = scored[scored["is_gnn_high_confidence"] == 1].copy()

scored.to_csv(SCORED_OUT, index=False)
high_conf.to_csv(HIGH_CONF_GNN_OUT, index=False)

print()
print("[53G] Saved scored candidate edges:")
print(SCORED_OUT)
print("Shape:", scored.shape)

print()
print("[53G] Saved GNN high-confidence edges:")
print(HIGH_CONF_GNN_OUT)
print("Shape:", high_conf.shape)


# ============================================================
# TF summary
# ============================================================

tf_summary = (
    scored
    .groupby(["cell_type", "tf"])
    .agg(
        n_candidate_targets=("target_gene", "nunique"),
        n_gnn_high_confidence_targets=("is_gnn_high_confidence", "sum"),
        mean_gnn_probability=("gnn_probability", "mean"),
        max_gnn_probability=("gnn_probability", "max"),
        mean_pseudo_supervised_probability=("pseudo_supervised_probability", "mean"),
        mean_evidence_rank_score=("evidence_rank_score", "mean"),
        mean_gnn_refined_composite_score=("gnn_refined_composite_score", "mean"),
    )
    .reset_index()
)

tf_summary["tf_rank_by_gnn_high_confidence_targets"] = (
    tf_summary
    .groupby("cell_type")["n_gnn_high_confidence_targets"]
    .rank(ascending=False, method="dense")
)

tf_summary["tf_rank_by_mean_gnn_composite_score"] = (
    tf_summary
    .groupby("cell_type")["mean_gnn_refined_composite_score"]
    .rank(ascending=False, method="dense")
)

tf_summary = tf_summary.sort_values(
    [
        "cell_type",
        "tf_rank_by_gnn_high_confidence_targets",
        "tf_rank_by_mean_gnn_composite_score",
    ]
)

tf_summary.to_csv(TF_SUMMARY_OUT, index=False)

print()
print("[53G] Saved GNN-refined TF summary:")
print(TF_SUMMARY_OUT)
print("Shape:", tf_summary.shape)


# ============================================================
# Quick sanity summaries
# ============================================================

print()
print("[53G] Edge count summary:")
print(
    scored
    .groupby("cell_type")
    .agg(
        n_candidate_edges=("edge_id", "nunique"),
        n_strategy1_pseudo_positive_edges=("is_pseudo_positive", "sum"),
        n_gnn_high_confidence_edges=("is_gnn_high_confidence", "sum"),
        mean_gnn_probability=("gnn_probability", "mean"),
        median_gnn_probability=("gnn_probability", "median"),
        mean_composite_score=("gnn_refined_composite_score", "mean"),
    )
    .reset_index()
    .to_string(index=False)
)

print()
print("[53G] Top 15 TFs per cell type by GNN high-confidence targets:")

top_tfs = (
    tf_summary
    .sort_values(
        [
            "cell_type",
            "n_gnn_high_confidence_targets",
            "mean_gnn_refined_composite_score",
        ],
        ascending=[True, False, False],
    )
    .groupby("cell_type")
    .head(15)
)

print(
    top_tfs[
        [
            "cell_type",
            "tf",
            "n_candidate_targets",
            "n_gnn_high_confidence_targets",
            "mean_gnn_probability",
            "mean_gnn_refined_composite_score",
            "tf_rank_by_gnn_high_confidence_targets",
        ]
    ].to_string(index=False)
)

print()
print("[53G] Done.")



"""Βήμα 54G — Compare evidence-based, pseudo-supervised, and GNN-refined GRNs"""

# ============================================================
# Imports
# ============================================================

from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from scipy.stats import fisher_exact
except Exception:
    fisher_exact = None


# ============================================================
# Paths
# ============================================================

PROJECT_DIR = Path(r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project")

GNN_DIR = PROJECT_DIR / "data" / "processed" / "gnn"
REPORT_TABLE_DIR = PROJECT_DIR / "reports" / "tables"
REPORT_FIGURE_DIR = PROJECT_DIR / "reports" / "figures" / "gnn_method_comparison"

REPORT_TABLE_DIR.mkdir(parents=True, exist_ok=True)
REPORT_FIGURE_DIR.mkdir(parents=True, exist_ok=True)

SCORED_FILE = GNN_DIR / "pbmc_grn_gnn_refined_edges_scored.csv"

TF_COMPARISON_OUT = REPORT_TABLE_DIR / "pbmc_grn_strategy_comparison_tf_summary.csv"
EDGE_JACCARD_OUT = REPORT_TABLE_DIR / "pbmc_grn_strategy_edge_jaccard.csv"
TF_TARGET_JACCARD_OUT = REPORT_TABLE_DIR / "pbmc_grn_strategy_tf_target_jaccard.csv"
IMMUNE_TF_ENRICHMENT_OUT = REPORT_TABLE_DIR / "pbmc_grn_strategy_known_immune_tf_enrichment.csv"
MARKER_TARGET_ENRICHMENT_OUT = REPORT_TABLE_DIR / "pbmc_grn_strategy_marker_target_enrichment.csv"
BROAD_TF_RANK_CHANGE_OUT = REPORT_TABLE_DIR / "pbmc_grn_strategy_broad_tf_rank_change.csv"
KEY_TF_RANK_TABLE_OUT = REPORT_TABLE_DIR / "pbmc_grn_strategy_key_tf_rank_table.csv"
METHOD_SIZE_OUT = REPORT_TABLE_DIR / "pbmc_grn_strategy_size_summary.csv"


# ============================================================
# Biological reference lists
# ============================================================

KNOWN_IMMUNE_TFS = {
    "SPI1", "RUNX1", "RUNX2",
    "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF7", "IRF8",
    "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
    "JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2",
    "NFKB1", "NFKB2", "RELA", "REL", "RELB",
    "IKZF1", "IKZF2", "IKZF3",
    "ETS1", "ELF1", "ELK1",
    "KLF2", "KLF3", "KLF6", "KLF12",
    "TCF7", "LEF1", "GATA3", "BCL11B",
    "BATF", "BATF3",
    "NFATC1", "NFATC2", "NFATC3",
    "EOMES", "TBX21", "PRDM1", "BACH2",
    "EGR1", "EGR2", "EGR3", "EGR4",
    "NR4A1", "NR4A2", "NR4A3",
    "CEBPA", "CEBPB", "CEBPD",
}

BROAD_MOTIF_TFS = {
    "ARNT", "AHR", "ARID3A", "MAZ", "THAP11", "MGA", "ZNF384", "CTCF",
    "SP1", "SP2", "SP3", "SP4",
    "YY1", "MAX", "MXI1", "USF1", "USF2",
    "E2F1", "E2F2", "E2F3", "E2F4",
    "NRF1",
}

KEY_TFS_TO_TRACK = [
    "SPI1", "IRF1", "STAT1", "STAT2",
    "JUN", "JUNB", "JUND", "FOS", "FOSB",
    "IKZF1", "ETS1", "KLF2", "TCF7", "LEF1",
    "GATA3", "NFATC2", "NR4A2",
    "ARNT", "AHR", "ARID3A", "THAP11", "ZNF384", "MAZ", "MGA", "CTCF", "NRF1",
]

MARKER_GENES_BY_CELL_TYPE = {
    "CD14+ Monocytes": {
        "FCN1", "CD14", "LYZ", "S100A8", "S100A9", "MS4A7", "CTSS", "LST1", "LGALS3",
    },
    "CD4.Memory": {
        "IL7R", "CD4", "LTB", "CCR7", "ANXA1", "S100A4", "CD40LG", "TRAC",
    },
    "CD4.Naive": {
        "CCR7", "LEF1", "TCF7", "IL7R", "CD4", "SELL", "LTB", "MAL",
    },
    "CD8.effector": {
        "NKG7", "PRF1", "GZMB", "CCL5", "GNLY", "CD8A", "CD8B", "GZMH", "CST7",
    },
    "CD8.Naive": {
        "CD8A", "CD8B", "CCR7", "LEF1", "TCF7", "SELL", "IL7R", "LTB",
    },
}


# ============================================================
# Helper functions
# ============================================================

def edge_key_df(df):
    out = df[["cell_type", "tf", "target_gene"]].copy()
    out["edge_key"] = (
        out["cell_type"].astype(str)
        + "|"
        + out["tf"].astype(str)
        + "|"
        + out["target_gene"].astype(str)
    )
    return out


def jaccard(a, b):
    a = set(a)
    b = set(b)
    if len(a | b) == 0:
        return np.nan
    return len(a & b) / len(a | b)


def select_top_n_by_cell_type(df, score_col, n_by_ct):
    selected = []

    for cell_type, ct_df in df.groupby("cell_type", sort=False):
        n_keep = int(n_by_ct.get(cell_type, max(1, int(np.ceil(ct_df.shape[0] * 0.10)))))
        selected.append(
            ct_df
            .sort_values(score_col, ascending=False)
            .head(n_keep)
            .copy()
        )

    return pd.concat(selected, ignore_index=True)


def fisher_enrichment(selected_items, background_items, positive_set):
    selected_items = set(selected_items)
    background_items = set(background_items)
    positive_set = set(positive_set)

    a = len(selected_items & positive_set)
    b = len(selected_items - positive_set)
    c = len((background_items - selected_items) & positive_set)
    d = len((background_items - selected_items) - positive_set)

    odds_ratio = np.nan
    p_value = np.nan

    if fisher_exact is not None:
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")

    return {
        "selected_positive": a,
        "selected_negative": b,
        "background_positive_not_selected": c,
        "background_negative_not_selected": d,
        "odds_ratio": odds_ratio,
        "p_value": p_value,
    }


def safe_rank(df, group_cols, value_col, rank_col):
    df[rank_col] = (
        df
        .groupby(group_cols)[value_col]
        .rank(ascending=False, method="dense")
    )
    return df


# ============================================================
# Load scored edge table
# ============================================================

scored = pd.read_csv(SCORED_FILE)

required_cols = [
    "cell_type", "tf", "target_gene", "edge_id",
    "is_pseudo_positive",
    "is_gnn_high_confidence",
    "pseudo_supervised_probability",
    "gnn_refined_composite_score",
    "gnn_probability",
]

missing = [c for c in required_cols if c not in scored.columns]
if missing:
    raise ValueError(f"Scored edge table is missing required columns: {missing}")

scored["is_pseudo_positive"] = pd.to_numeric(
    scored["is_pseudo_positive"], errors="coerce"
).fillna(0).astype(int)

scored["is_gnn_high_confidence"] = pd.to_numeric(
    scored["is_gnn_high_confidence"], errors="coerce"
).fillna(0).astype(int)


# ============================================================
# Define three strategy-specific selected edge sets
# ============================================================

evidence_edges = scored[scored["is_pseudo_positive"] == 1].copy()
gnn_edges = scored[scored["is_gnn_high_confidence"] == 1].copy()

n_evidence_by_cell_type = evidence_edges.groupby("cell_type").size().to_dict()

pseudo_supervised_edges = select_top_n_by_cell_type(
    scored,
    score_col="pseudo_supervised_probability",
    n_by_ct=n_evidence_by_cell_type,
)

method_edges = {
    "strategy1_evidence_based": evidence_edges,
    "strategy2_pseudo_supervised": pseudo_supervised_edges,
    "strategy3_gnn_refined": gnn_edges,
}


# ============================================================
# Method size summary
# ============================================================

size_rows = []

for method, df in method_edges.items():
    for cell_type, ct_df in df.groupby("cell_type"):
        size_rows.append({
            "method": method,
            "cell_type": cell_type,
            "n_edges": ct_df[["cell_type", "tf", "target_gene"]].drop_duplicates().shape[0],
            "n_tfs": ct_df["tf"].nunique(),
            "n_targets": ct_df["target_gene"].nunique(),
            "mean_gnn_probability": ct_df["gnn_probability"].mean(),
            "mean_pseudo_supervised_probability": ct_df["pseudo_supervised_probability"].mean(),
            "mean_composite_score": ct_df["gnn_refined_composite_score"].mean(),
        })

size_df = pd.DataFrame(size_rows)
size_df.to_csv(METHOD_SIZE_OUT, index=False)


# ============================================================
# TF summaries per strategy
# ============================================================

tf_summary_tables = []

for method, df in method_edges.items():
    summary = (
        df
        .groupby(["cell_type", "tf"])
        .agg(
            n_targets=("target_gene", "nunique"),
            n_edges=("target_gene", "size"),
            mean_gnn_probability=("gnn_probability", "mean"),
            mean_pseudo_supervised_probability=("pseudo_supervised_probability", "mean"),
            mean_gnn_refined_composite_score=("gnn_refined_composite_score", "mean"),
        )
        .reset_index()
    )

    summary["method"] = method
    summary["is_known_immune_tf"] = summary["tf"].isin(KNOWN_IMMUNE_TFS).astype(int)
    summary["is_broad_motif_tf"] = summary["tf"].isin(BROAD_MOTIF_TFS).astype(int)

    summary = safe_rank(
        summary,
        group_cols=["method", "cell_type"],
        value_col="n_targets",
        rank_col="rank_by_n_targets",
    )

    summary = safe_rank(
        summary,
        group_cols=["method", "cell_type"],
        value_col="mean_gnn_refined_composite_score",
        rank_col="rank_by_mean_score",
    )

    tf_summary_tables.append(summary)

tf_comparison = pd.concat(tf_summary_tables, ignore_index=True)
tf_comparison.to_csv(TF_COMPARISON_OUT, index=False)


# ============================================================
# Edge-level Jaccard overlap between strategies
# ============================================================

edge_jaccard_rows = []

for cell_type in sorted(scored["cell_type"].unique()):
    ct_sets = {
        method: set(edge_key_df(df[df["cell_type"] == cell_type])["edge_key"])
        for method, df in method_edges.items()
    }

    for method_1, method_2 in combinations(ct_sets.keys(), 2):
        edge_jaccard_rows.append({
            "cell_type": cell_type,
            "method_1": method_1,
            "method_2": method_2,
            "n_edges_method_1": len(ct_sets[method_1]),
            "n_edges_method_2": len(ct_sets[method_2]),
            "n_intersection": len(ct_sets[method_1] & ct_sets[method_2]),
            "jaccard": jaccard(ct_sets[method_1], ct_sets[method_2]),
        })

edge_jaccard_df = pd.DataFrame(edge_jaccard_rows)
edge_jaccard_df.to_csv(EDGE_JACCARD_OUT, index=False)


# ============================================================
# TF target-set Jaccard overlap between strategies
# ============================================================

tf_target_jaccard_rows = []

for cell_type in sorted(scored["cell_type"].unique()):
    all_tfs = sorted(scored.loc[scored["cell_type"] == cell_type, "tf"].dropna().unique())

    for tf in all_tfs:
        target_sets = {
            method: set(df[(df["cell_type"] == cell_type) & (df["tf"] == tf)]["target_gene"])
            for method, df in method_edges.items()
        }

        for method_1, method_2 in combinations(target_sets.keys(), 2):
            if len(target_sets[method_1]) == 0 and len(target_sets[method_2]) == 0:
                continue

            tf_target_jaccard_rows.append({
                "cell_type": cell_type,
                "tf": tf,
                "method_1": method_1,
                "method_2": method_2,
                "n_targets_method_1": len(target_sets[method_1]),
                "n_targets_method_2": len(target_sets[method_2]),
                "n_intersection": len(target_sets[method_1] & target_sets[method_2]),
                "jaccard": jaccard(target_sets[method_1], target_sets[method_2]),
                "is_known_immune_tf": int(tf in KNOWN_IMMUNE_TFS),
                "is_broad_motif_tf": int(tf in BROAD_MOTIF_TFS),
            })

tf_target_jaccard_df = pd.DataFrame(tf_target_jaccard_rows)
tf_target_jaccard_df.to_csv(TF_TARGET_JACCARD_OUT, index=False)


# ============================================================
# Known immune TF enrichment per strategy
# ============================================================

immune_enrichment_rows = []

for cell_type in sorted(scored["cell_type"].unique()):
    background_tfs = set(scored.loc[scored["cell_type"] == cell_type, "tf"].dropna())

    for method, df in method_edges.items():
        selected_tfs = set(df.loc[df["cell_type"] == cell_type, "tf"].dropna())

        enrichment = fisher_enrichment(
            selected_items=selected_tfs,
            background_items=background_tfs,
            positive_set=KNOWN_IMMUNE_TFS,
        )

        immune_enrichment_rows.append({
            "cell_type": cell_type,
            "method": method,
            "n_selected_tfs": len(selected_tfs),
            "n_background_tfs": len(background_tfs),
            **enrichment,
        })

immune_enrichment_df = pd.DataFrame(immune_enrichment_rows)
immune_enrichment_df.to_csv(IMMUNE_TF_ENRICHMENT_OUT, index=False)


# ============================================================
# Marker target enrichment per strategy
# ============================================================

marker_enrichment_rows = []

for cell_type in sorted(scored["cell_type"].unique()):
    background_targets = set(scored.loc[scored["cell_type"] == cell_type, "target_gene"].dropna())
    marker_targets = MARKER_GENES_BY_CELL_TYPE.get(cell_type, set())

    for method, df in method_edges.items():
        selected_targets = set(df.loc[df["cell_type"] == cell_type, "target_gene"].dropna())

        enrichment = fisher_enrichment(
            selected_items=selected_targets,
            background_items=background_targets,
            positive_set=marker_targets,
        )

        marker_enrichment_rows.append({
            "cell_type": cell_type,
            "method": method,
            "n_selected_targets": len(selected_targets),
            "n_background_targets": len(background_targets),
            "n_marker_targets_reference": len(marker_targets),
            **enrichment,
        })

marker_enrichment_df = pd.DataFrame(marker_enrichment_rows)
marker_enrichment_df.to_csv(MARKER_TARGET_ENRICHMENT_OUT, index=False)


# ============================================================
# Broad TF rank-change analysis
# ============================================================

rank_table = tf_comparison.pivot_table(
    index=["cell_type", "tf"],
    columns="method",
    values="rank_by_n_targets",
    aggfunc="min",
).reset_index()

for method in method_edges.keys():
    if method not in rank_table.columns:
        rank_table[method] = np.nan

rank_table["is_known_immune_tf"] = rank_table["tf"].isin(KNOWN_IMMUNE_TFS).astype(int)
rank_table["is_broad_motif_tf"] = rank_table["tf"].isin(BROAD_MOTIF_TFS).astype(int)

rank_table["rank_change_evidence_to_pseudo"] = (
    rank_table["strategy2_pseudo_supervised"]
    - rank_table["strategy1_evidence_based"]
)

rank_table["rank_change_evidence_to_gnn"] = (
    rank_table["strategy3_gnn_refined"]
    - rank_table["strategy1_evidence_based"]
)

broad_rank_change = rank_table[rank_table["is_broad_motif_tf"] == 1].copy()
broad_rank_change.to_csv(BROAD_TF_RANK_CHANGE_OUT, index=False)


# ============================================================
# Key TF rank table
# ============================================================

key_tf_rank_table = rank_table[rank_table["tf"].isin(KEY_TFS_TO_TRACK)].copy()
key_tf_rank_table.to_csv(KEY_TF_RANK_TABLE_OUT, index=False)


# ============================================================
# Figures — top TF target counts by strategy
# ============================================================

for cell_type in sorted(tf_comparison["cell_type"].unique()):
    ct_df = tf_comparison[tf_comparison["cell_type"] == cell_type].copy()

    top_tfs = (
        ct_df
        .sort_values(["method", "n_targets"], ascending=[True, False])
        .groupby("method")
        .head(12)["tf"]
        .unique()
    )

    plot_df = (
        ct_df[ct_df["tf"].isin(top_tfs)]
        .pivot_table(
            index="tf",
            columns="method",
            values="n_targets",
            aggfunc="max",
            fill_value=0,
        )
    )

    plot_df = plot_df.loc[plot_df.sum(axis=1).sort_values(ascending=False).index]

    plt.figure(figsize=(12, 5))
    plot_df.plot(kind="bar", ax=plt.gca())
    plt.title(f"Top TF target counts by GRN strategy — {cell_type}")
    plt.xlabel("TF")
    plt.ylabel("Number of selected targets")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    fig_path = REPORT_FIGURE_DIR / f"top_tf_targets_by_strategy_{cell_type.replace('+', 'plus').replace(' ', '_')}.png"
    plt.savefig(fig_path, dpi=300)
    plt.close()


# ============================================================
# Figures — edge Jaccard heatmaps
# ============================================================

for cell_type in sorted(edge_jaccard_df["cell_type"].unique()):
    ct_df = edge_jaccard_df[edge_jaccard_df["cell_type"] == cell_type].copy()

    labels = sorted(set(ct_df["method_1"]) | set(ct_df["method_2"]))
    mat = pd.DataFrame(np.eye(len(labels)), index=labels, columns=labels)

    for _, row in ct_df.iterrows():
        mat.loc[row["method_1"], row["method_2"]] = row["jaccard"]
        mat.loc[row["method_2"], row["method_1"]] = row["jaccard"]

    plt.figure(figsize=(6, 5))
    plt.imshow(mat.values)
    plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
    plt.yticks(range(len(labels)), labels)
    plt.title(f"Edge-set Jaccard overlap — {cell_type}")
    plt.colorbar(label="Jaccard")

    for i in range(len(labels)):
        for j in range(len(labels)):
            plt.text(j, i, f"{mat.values[i, j]:.2f}", ha="center", va="center")

    plt.tight_layout()

    fig_path = REPORT_FIGURE_DIR / f"edge_jaccard_by_strategy_{cell_type.replace('+', 'plus').replace(' ', '_')}.png"
    plt.savefig(fig_path, dpi=300)
    plt.close()


# ============================================================
# Print concise results for interpretation
# ============================================================

print("[54G] Done.")
print()
print("[54G] Saved tables:")
print(f"Method size summary: {METHOD_SIZE_OUT}")
print(f"TF comparison: {TF_COMPARISON_OUT}")
print(f"Edge Jaccard: {EDGE_JACCARD_OUT}")
print(f"TF target Jaccard: {TF_TARGET_JACCARD_OUT}")
print(f"Known immune TF enrichment: {IMMUNE_TF_ENRICHMENT_OUT}")
print(f"Marker target enrichment: {MARKER_TARGET_ENRICHMENT_OUT}")
print(f"Broad TF rank changes: {BROAD_TF_RANK_CHANGE_OUT}")
print(f"Key TF rank table: {KEY_TF_RANK_TABLE_OUT}")
print(f"Figures: {REPORT_FIGURE_DIR}")

print()
print("[54G] Method size summary:")
print(size_df.sort_values(["cell_type", "method"]).to_string(index=False))

print()
print("[54G] Edge Jaccard overlap:")
print(edge_jaccard_df.sort_values(["cell_type", "method_1", "method_2"]).to_string(index=False))

print()
print("[54G] Known immune TF enrichment:")
print(immune_enrichment_df.sort_values(["cell_type", "method"]).to_string(index=False))

print()
print("[54G] Marker target enrichment:")
print(marker_enrichment_df.sort_values(["cell_type", "method"]).to_string(index=False))

print()
print("[54G] Broad TF rank changes:")
print(
    broad_rank_change
    .sort_values(["cell_type", "strategy3_gnn_refined"])
    .to_string(index=False)
)

print()
print("[54G] Key TF ranks:")
print(
    key_tf_rank_table
    .sort_values(["cell_type", "strategy3_gnn_refined"])
    .to_string(index=False)
)






