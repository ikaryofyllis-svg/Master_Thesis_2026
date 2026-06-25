"""Βήμα 02A — Reusable RNA preprocessing functions για το refactored PBMC GRN project"""

from __future__ import annotations  # επιτρέπει καθαρότερα type hints

from typing import Iterable  # για type hints σε λίστες/στήλες

import pandas as pd  # για summary tables
import scanpy as sc  # βασικό scRNA-seq preprocessing API


# ============================================================
# 0. Internal validation helpers
# ============================================================

def _require_obs_columns(adata, required_columns: Iterable[str], step_name: str) -> None:
    """Check ότι υπάρχουν απαραίτητα columns στο adata.obs."""

    missing_columns = [col for col in required_columns if col not in adata.obs.columns]  # βρίσκει missing obs columns

    if missing_columns:  # αν λείπει κάτι, σταματάμε καθαρά
        raise ValueError(  # δίνει χρήσιμο error αντί για silent failure
            f"Missing required adata.obs columns before {step_name}: {missing_columns}. "
            "Run the previous preprocessing step first."
        )


def _require_var_columns(adata, required_columns: Iterable[str], step_name: str) -> None:
    """Check ότι υπάρχουν απαραίτητα columns στο adata.var."""

    missing_columns = [col for col in required_columns if col not in adata.var.columns]  # βρίσκει missing var columns

    if missing_columns:  # αν λείπει κάτι, σταματάμε καθαρά
        raise ValueError(  # καθαρό diagnostic error
            f"Missing required adata.var columns before {step_name}: {missing_columns}. "
            "Run the previous preprocessing step first."
        )


# ============================================================
# 1. RNA QC metrics
# ============================================================

def calculate_rna_qc_metrics(
    adata,
    mito_prefix: str = "MT-",
):
    """
    Calculate basic RNA QC metrics.

    Adds:
    - adata.var["mt"]
    - adata.obs["total_counts"]
    - adata.obs["n_genes_by_counts"]
    - adata.obs["pct_counts_mt"]
    """

    adata = adata.copy()  # κρατάμε τις συναρτήσεις pure, χωρίς να αλλάζουν το original object

    adata.var["mt"] = adata.var_names.astype(str).str.startswith(mito_prefix)  # flag για mitochondrial genes

    sc.pp.calculate_qc_metrics(  # υπολογίζει βασικά RNA QC metrics
        adata,  # AnnData object
        qc_vars=["mt"],  # χρησιμοποιεί το mitochondrial gene flag
        percent_top=None,  # δεν χρειάζονται top-N gene percentages εδώ
        log1p=False,  # κρατάμε raw-scale QC metrics
        inplace=True,  # γράφει τα metrics στο adata.obs και adata.var
    )

    return adata  # επιστρέφει AnnData με QC columns


# ============================================================
# 2. RNA QC summary table
# ============================================================

def summarize_rna_qc(
    adata,
    label: str,
):
    """
    Create compact QC summary table for an RNA AnnData object.
    Useful before/after filtering comparisons.
    """

    required_columns = [  # βασικά QC metrics που περιμένουμε
        "total_counts",  # library size per cell
        "n_genes_by_counts",  # detected genes per cell
        "pct_counts_mt",  # mitochondrial percentage per cell
    ]

    _require_obs_columns(adata, required_columns, "summarize_rna_qc")  # προστασία από λάθος σειρά βημάτων

    summary = adata.obs[required_columns].describe().T  # descriptive statistics ανά QC metric

    summary.insert(0, "stage", label)  # π.χ. raw, filtered
    summary.insert(1, "n_cells", adata.n_obs)  # αριθμός cells στο συγκεκριμένο stage
    summary.insert(2, "n_genes", adata.n_vars)  # αριθμός genes στο συγκεκριμένο stage

    summary = summary.reset_index(names="metric")  # μετατρέπει το metric index σε κανονική στήλη

    return summary  # επιστρέφει dataframe για reports/tables


# ============================================================
# 3. RNA filtering
# ============================================================

def filter_rna_cells_and_genes(
    adata,
    min_genes_per_cell: int,
    min_cells_per_gene: int,
    max_pct_mito: float,
    max_total_counts: float,
):
    """
    Apply basic RNA filtering.

    Filtering logic:
    - remove cells with few detected genes
    - remove genes detected in very few cells
    - remove cells with high mitochondrial percentage
    - remove cells with extreme total counts, likely doublets
    """

    required_columns = [  # columns που πρέπει να υπάρχουν από QC
        "total_counts",  # χρειάζεται για high-count filtering
        "n_genes_by_counts",  # χρήσιμο για diagnostics
        "pct_counts_mt",  # χρειάζεται για mito filtering
    ]

    _require_obs_columns(adata, required_columns, "filter_rna_cells_and_genes")  # εξασφαλίζει ότι έχει τρέξει QC

    adata = adata.copy()  # δουλεύουμε σε copy για predictable behavior

    n_cells_before = adata.n_obs  # cells πριν το filtering
    n_genes_before = adata.n_vars  # genes πριν το filtering

    sc.pp.filter_cells(  # αφαιρεί cells με πολύ λίγα detected genes
        adata,  # AnnData object
        min_genes=min_genes_per_cell,  # threshold από YAML/config
    )

    sc.pp.filter_genes(  # αφαιρεί genes που εμφανίζονται σε πολύ λίγα cells
        adata,  # AnnData object
        min_cells=min_cells_per_gene,  # threshold από YAML/config
    )

    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mito].copy()  # αφαιρεί high-mito low-quality cells

    adata = adata[adata.obs["total_counts"] < max_total_counts].copy()  # αφαιρεί extreme-count πιθανές doublets

    filtering_summary = pd.DataFrame(  # compact table για thesis/report
        [
            {
                "cells_before": n_cells_before,  # cells πριν
                "cells_after": adata.n_obs,  # cells μετά
                "cells_removed": n_cells_before - adata.n_obs,  # πόσα cells αφαιρέθηκαν
                "cell_retention_fraction": adata.n_obs / n_cells_before if n_cells_before else 0,  # ποσοστό cells που έμειναν
                "genes_before": n_genes_before,  # genes πριν
                "genes_after": adata.n_vars,  # genes μετά
                "genes_removed": n_genes_before - adata.n_vars,  # πόσα genes αφαιρέθηκαν
                "gene_retention_fraction": adata.n_vars / n_genes_before if n_genes_before else 0,  # ποσοστό genes που έμειναν
                "min_genes_per_cell": min_genes_per_cell,  # recorded parameter
                "min_cells_per_gene": min_cells_per_gene,  # recorded parameter
                "max_pct_mito": max_pct_mito,  # recorded parameter
                "max_total_counts": max_total_counts,  # recorded parameter
            }
        ]
    )

    return adata, filtering_summary  # επιστρέφει filtered AnnData και report table


# ============================================================
# 4. RNA normalization and HVG selection
# ============================================================

def normalize_log_hvg(
    adata,
    target_sum: float,
    apply_log1p: bool,
    n_top_genes: int,
    hvg_flavor: str,
):
    """
    Store counts, normalize, log-transform, and select highly variable genes.

    Important:
    - raw counts are stored in adata.layers["counts"]
    - normalized/log expression for all genes is stored in adata.raw
    - returned object keeps only HVGs for PCA/UMAP/clustering
    """

    adata = adata.copy()  # κρατάμε το input filtered object ανέπαφο

    adata.layers["counts"] = adata.X.copy()  # αποθηκεύει raw counts πριν normalization

    sc.pp.normalize_total(  # library-size normalization
        adata,  # AnnData object
        target_sum=target_sum,  # συνήθως 10000
    )

    sc.pp.log1p(adata)  # log-transform για πιο σταθερή έκφραση

    adata.raw = adata.copy()  # κρατά normalized/log expression όλων των genes για markers/plots

    if hvg_flavor in {"seurat_v3", "seurat_v3_paper"}:  # αυτά τα flavors προτιμούν counts
        sc.pp.highly_variable_genes(  # HVG selection από raw counts layer
            adata,  # AnnData object
            n_top_genes=n_top_genes,  # πόσα HVGs να κρατήσουμε
            flavor=hvg_flavor,  # method από YAML/config
            layer="counts",  # χρησιμοποιεί raw counts αντί για log data
        )
    else:  # π.χ. seurat ή cell_ranger πάνω σε log-normalized data
        sc.pp.highly_variable_genes(  # HVG selection από το current log-normalized X
            adata,  # AnnData object
            n_top_genes=n_top_genes,  # πόσα HVGs να κρατήσουμε
            flavor=hvg_flavor,  # method από YAML/config
        )

    _require_var_columns(adata, ["highly_variable"], "normalize_log_hvg")  # επιβεβαιώνει ότι δημιουργήθηκε HVG flag

    n_hvg = int(adata.var["highly_variable"].sum())  # μετρά selected HVGs

    hvg_summary = pd.DataFrame(  # compact report table
        [
            {
                "n_genes_total": adata.n_vars,  # genes πριν το HVG subset
                "n_hvg": n_hvg,  # selected highly variable genes
                "n_non_hvg": adata.n_vars - n_hvg,  # genes που δεν κρατήθηκαν για dimensionality reduction
                "target_sum": target_sum,  # recorded normalization parameter
                "n_top_genes": n_top_genes,  # recorded HVG parameter
                "hvg_flavor": hvg_flavor,  # recorded HVG method
            }
        ]
    )

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()  # κρατά μόνο HVGs για PCA/UMAP/clustering

    return adata_hvg, hvg_summary  # επιστρέφει HVG AnnData και summary table


# ============================================================
# 5. PCA, neighbors, UMAP, Leiden clustering
# ============================================================

def run_rna_dimensionality_reduction_and_clustering(
    adata,
    scale_max_value: float,
    pca_solver: str,
    n_comps: int,
    n_neighbors: int,
    n_pcs: int,
    leiden_resolution: float,
    leiden_key: str,
    random_state: int,
):
    """
    Run scaling, PCA, neighbors, UMAP and Leiden clustering.
    This function assumes the input AnnData is already HVG-filtered.
    """

    adata = adata.copy()  # δεν αλλάζουμε το input object απευθείας

    sc.pp.scale(  # standardization ανά gene
        adata,  # AnnData object
        max_value=scale_max_value,  # clipping για extreme scaled values
    )

    sc.tl.pca(
        adata,  # HVG AnnData
        n_comps=n_comps,  # από YAML, π.χ. 50
        svd_solver=pca_solver,  # από YAML, π.χ. arpack
        random_state=random_state,  # reproducibility
    )
    sc.pp.neighbors(  # φτιάχνει kNN graph
        adata,  # AnnData object
        n_neighbors=n_neighbors,  # πόσοι γείτονες ανά cell
        n_pcs=n_pcs,  # πόσα PCs χρησιμοποιούνται
    )

    sc.tl.umap(  # UMAP embedding για visualization
        adata,  # AnnData object
        random_state=random_state,  # reproducible layout
    )

    sc.tl.leiden(  # Leiden clustering πάνω στο neighbor graph
        adata,  # AnnData object
        resolution=leiden_resolution,  # granularity των clusters
        key_added=leiden_key,  # π.χ. leiden_0_6
        random_state=random_state,  # reproducible clustering
    )

    return adata  # επιστρέφει clustered AnnData


# ============================================================
# 6. Cluster summary
# ============================================================

def summarize_clusters(
    adata,
    cluster_key: str,
):
    """
    Summarize number and fraction of cells per RNA cluster.
    """

    _require_obs_columns(adata, [cluster_key], "summarize_clusters")  # ελέγχει ότι υπάρχει το cluster column

    cluster_summary = (  # group-level table
        adata.obs[cluster_key]  # παίρνει cluster labels
        .value_counts()  # μετρά cells ανά cluster
        .sort_index()  # ταξινομεί τα cluster ids
        .rename_axis(cluster_key)  # ονομάζει τη στήλη των clusters
        .reset_index(name="n_cells")  # κάνει tidy dataframe
    )

    cluster_summary["fraction_cells"] = cluster_summary["n_cells"] / adata.n_obs  # ποσοστό dataset ανά cluster

    return cluster_summary  # επιστρέφει report table

# ============================================================
# 7. RNA PCA-only diagnostics
# ============================================================

def run_rna_pca_only(
    adata,
    scale_max_value: float,
    pca_solver: str,
    n_comps: int,
    random_state: int,
):
    """
    Scale RNA HVG data and run PCA only.
    Useful before choosing n_pcs for neighbors/UMAP/Leiden.
    """

    adata = adata.copy()  # keep input object unchanged

    sc.pp.scale(
        adata,
        max_value=scale_max_value,
    )  # standardize genes and clip extreme values

    sc.tl.pca(
        adata,
        n_comps=n_comps,
        svd_solver=pca_solver,
        random_state=random_state,
    )  # compute PCA for elbow diagnostics

    return adata  # returns AnnData with PCA stored in .obsm and .uns

# ============================================================
# 8. PCA variance table
# ============================================================

def summarize_pca_variance(adata):
    """
    Create a PCA variance table for choosing how many PCs to use.
    """

    variance_ratio = adata.uns["pca"]["variance_ratio"]  # variance explained per PC

    pca_table = pd.DataFrame(
        {
            "pc": range(1, len(variance_ratio) + 1),
            "variance_ratio": variance_ratio,
            "cumulative_variance_ratio": variance_ratio.cumsum(),
        }
    )  # tidy table for reports and parameter justification

    return pca_table  # return PCA variance summary



# ============================================================
# 9. UMAP and Leiden from existing PCA
# ============================================================

def run_rna_neighbors_umap_leiden(
    adata,
    n_neighbors: int,
    n_pcs: int,
    leiden_resolution: float,
    leiden_key: str,
    random_state: int,
):
    """
    Run neighbors, UMAP and Leiden using an already PCA-computed AnnData.
    """

    adata = adata.copy()  # avoid modifying PCA diagnostic object directly

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
    )  # build kNN graph using selected PCs

    sc.tl.umap(
        adata,
        random_state=random_state,
    )  # compute UMAP layout

    sc.tl.leiden(
        adata,
        resolution=leiden_resolution,
        key_added=leiden_key,
        random_state=random_state,
    )  # add Leiden cluster labels to adata.obs

    return adata  # return UMAP + Leiden object