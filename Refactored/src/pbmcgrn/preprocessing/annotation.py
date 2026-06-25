"""Βήμα 03B — Reusable annotation functions για RNA marker validation και cell-type annotation"""

import pandas as pd


# ============================================================
# 1. Flatten marker dictionary
# ============================================================

def flatten_marker_dict(marker_dict):
    """
    Convert marker dictionary into a unique marker gene list.

    Example:
    {"T_cells": ["CD3D", "CD3E"], "Monocytes": ["CD14"]}
    -> ["CD3D", "CD3E", "CD14"]
    """

    markers = []

    for genes in marker_dict.values():
        markers.extend(genes)

    markers = list(dict.fromkeys(markers))

    return markers


# ============================================================
# 2. Find available and missing markers
# ============================================================

def get_available_markers(adata, marker_genes, use_raw=True):
    """
    Split markers into available and missing markers.

    If use_raw=True and adata.raw exists, marker presence is checked in adata.raw.var_names.
    Otherwise, marker presence is checked in adata.var_names.
    """

    if use_raw and adata.raw is not None:
        var_names = adata.raw.var_names
    else:
        var_names = adata.var_names

    available = [gene for gene in marker_genes if gene in var_names]
    missing = [gene for gene in marker_genes if gene not in var_names]

    return available, missing


# ============================================================
# 3. Mean marker expression by cluster
# ============================================================

def compute_mean_marker_expression(adata, marker_genes, cluster_key, use_raw=True):
    """
    Compute mean marker expression per cluster.
    """

    available_markers, missing_markers = get_available_markers(
        adata=adata,
        marker_genes=marker_genes,
        use_raw=use_raw,
    )

    if len(available_markers) == 0:
        raise ValueError("No marker genes found in AnnData object.")

    if use_raw and adata.raw is not None:
        marker_adata = adata.raw[:, available_markers].to_adata()
    else:
        marker_adata = adata[:, available_markers].copy()

    marker_df = marker_adata.to_df()
    marker_df[cluster_key] = adata.obs[cluster_key].values

    mean_expr = marker_df.groupby(cluster_key).mean()

    return mean_expr, available_markers, missing_markers


# ============================================================
# 4. Apply cluster annotation
# ============================================================

def apply_cluster_annotation(
    adata,
    cluster_key,
    detailed_mapping,
    broad_mapping,
):
    """
    Add detailed and broad cell-type annotations to adata.obs.

    Adds:
    - cell_type_detailed
    - cell_type_broad
    """

    adata = adata.copy()

    adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)

    adata.obs["cell_type_detailed"] = (
        adata.obs[cluster_key]
        .map(detailed_mapping)
    )

    adata.obs["cell_type_broad"] = (
        adata.obs[cluster_key]
        .map(broad_mapping)
    )

    missing_detailed = int(adata.obs["cell_type_detailed"].isna().sum())
    missing_broad = int(adata.obs["cell_type_broad"].isna().sum())

    if missing_detailed > 0 or missing_broad > 0:
        raise ValueError(
            "Some clusters were not annotated. "
            f"Missing detailed: {missing_detailed}; missing broad: {missing_broad}"
        )

    return adata


# ============================================================
# 5. Annotation summary
# ============================================================

def summarize_annotation(adata, cluster_key):
    """
    Create summary table for cluster-level annotation.
    """

    summary = (
        adata.obs
        .groupby([cluster_key, "cell_type_detailed", "cell_type_broad"])
        .size()
        .reset_index(name="n_cells")
        .sort_values(cluster_key)
    )

    return summary


# ============================================================
# 6. Cell type counts
# ============================================================

def summarize_cell_type_counts(adata):
    """
    Create broad and detailed cell-type count table.
    """

    broad_counts = (
        adata.obs["cell_type_broad"]
        .value_counts()
        .rename_axis("cell_type")
        .reset_index(name="n_cells")
    )

    broad_counts.insert(0, "annotation_level", "broad")

    detailed_counts = (
        adata.obs["cell_type_detailed"]
        .value_counts()
        .rename_axis("cell_type")
        .reset_index(name="n_cells")
    )

    detailed_counts.insert(0, "annotation_level", "detailed")

    counts = pd.concat(
        [broad_counts, detailed_counts],
        axis=0,
        ignore_index=True,
    )

    return counts


# ============================================================
# 7. Subset to scMultiomeGRN 5 cell types
# ============================================================

def subset_to_scmultiomegrn_5types(adata, mapping):
    """
    Subset RNA object to the 5 PBMC cell types used in the thesis GRN pipeline.

    Input mapping example:
    {
        "CD14_Monocytes": "CD14+ Monocytes",
        "CD4_Memory": "CD4.Memory",
        ...
    }
    """

    adata = adata.copy()

    target_broad_cell_types = list(mapping.keys())

    adata_5 = adata[
        adata.obs["cell_type_broad"].isin(target_broad_cell_types)
    ].copy()

    adata_5.obs["scmultiomegrn_cell_type"] = (
        adata_5.obs["cell_type_broad"]
        .map(mapping)
    )

    counts = (
        adata_5.obs["scmultiomegrn_cell_type"]
        .value_counts()
        .rename_axis("scmultiomegrn_cell_type")
        .reset_index(name="n_cells")
    )

    return adata_5, counts