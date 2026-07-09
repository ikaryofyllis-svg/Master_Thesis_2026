# -*- coding: utf-8 -*-
"""
Step 03A — RNA-informed marker-based annotation evidence.

Purpose
-------
1. Load the final selected clustered RNA object from Step 02B.
2. Read the selected Leiden cluster key dynamically from
   configs/preprocessing_rna.yaml.
3. Read marker groups from configs/rna_marker_validation.yaml.
4. Generate marker-based annotation evidence:
   - cluster sizes
   - marker availability
   - marker mean expression per cluster
   - marker percent-expressing per cluster
   - marker-group scores
   - marker dotplot
   - marker matrixplot
   - marker UMAP batches
   - manual annotation worksheet
5. Do NOT:
   - rerun QC
   - rerun normalization
   - rerun HVG selection
   - rerun PCA
   - rerun neighbors
   - rerun UMAP
   - rerun Leiden
   - assign final cell-type labels
   - write a final annotated h5ad

The selected clustering source of truth is:

configs/preprocessing_rna.yaml
    dimensionality_reduction:
        leiden_key
        n_pcs
        n_neighbors
        leiden_resolution
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any
import re
import sys

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from anndata import AnnData


# ============================================================
# 0. Resolve Refactored root and import path
# ============================================================

try:
    # Expected location:
    # Refactored/scripts/run_03A_rna_informed_annotation.py
    refactor_dir = Path(__file__).resolve().parents[2]
except NameError:
    # Fallback for Spyder interactive execution
    refactor_dir = Path.cwd().resolve()

src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))


# ============================================================
# 1. Import project helpers
# ============================================================

from pbmcgrn.config import load_yaml
from pbmcgrn.io import write_csv, write_json
from pbmcgrn.paths import ProjectPaths


# ============================================================
# 2. Generic helpers
# ============================================================

def assert_path_inside_project(
    path: Path,
    project_root: Path,
    label: str,
) -> None:
    """
    Fail if a path is outside the active project workspace.
    """

    resolved_path = Path(path).resolve()
    resolved_project_root = Path(project_root).resolve()

    try:
        resolved_path.relative_to(resolved_project_root)
    except ValueError as exc:
        raise ValueError(
            f"{label} is outside the active project workspace.\n"
            f"Path: {resolved_path}\n"
            f"Active project root: {resolved_project_root}"
        ) from exc


def require_config_key(
    config: dict[str, Any],
    key: str,
    config_name: str,
) -> Any:
    """
    Return a required config key or fail clearly.
    """

    if key not in config:
        raise KeyError(
            f"Missing required key '{key}' in {config_name}."
        )

    return config[key]


def ensure_parent_dir(path: Path) -> None:
    """
    Create the parent directory of an output file.
    """

    Path(path).parent.mkdir(
        parents=True,
        exist_ok=True,
    )


def safe_filename_token(text: str) -> str:
    """
    Convert text into a filename-safe token.
    """

    safe = re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        str(text),
    )

    return safe.strip("_")


def sort_cluster_table(
    df: pd.DataFrame,
    cluster_col: str = "cluster",
) -> pd.DataFrame:
    """
    Sort Leiden cluster labels numerically when possible.
    """

    sorted_df = df.copy()

    sorted_df["_cluster_numeric"] = pd.to_numeric(
        sorted_df[cluster_col],
        errors="coerce",
    )

    sorted_df = sorted_df.sort_values(
        by=["_cluster_numeric", cluster_col],
        na_position="last",
    )

    sorted_df = sorted_df.drop(
        columns=["_cluster_numeric"],
    )

    return sorted_df.reset_index(drop=True)


def save_current_figure(
    output_path: Path,
    dpi: int,
    display_figures: bool,
) -> None:
    """
    Save the current Matplotlib figure and optionally display it.
    """

    ensure_parent_dir(output_path)

    plt.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
    )

    if display_figures:
        plt.show()

    plt.close("all")


# ============================================================
# 3. Marker helpers
# ============================================================

def flatten_marker_dict(
    marker_dict: dict[str, list[str]],
) -> list[str]:
    """
    Flatten marker groups into one unique ordered list.
    """

    markers: list[str] = []

    for genes in marker_dict.values():
        markers.extend(
            [str(gene) for gene in genes]
        )

    # Remove duplicates while preserving order
    return list(dict.fromkeys(markers))


def get_marker_var_names(
    adata: AnnData,
    use_raw: bool,
) -> pd.Index:
    """
    Return the gene universe used for marker lookup.
    """

    if use_raw and adata.raw is not None:
        return adata.raw.var_names

    return adata.var_names


def split_available_missing_markers(
    adata: AnnData,
    marker_genes: list[str],
    use_raw: bool,
) -> tuple[list[str], list[str]]:
    """
    Split requested markers into available and missing genes.
    """

    var_names = set(
        get_marker_var_names(
            adata=adata,
            use_raw=use_raw,
        ).astype(str)
    )

    available = [
        gene
        for gene in marker_genes
        if gene in var_names
    ]

    missing = [
        gene
        for gene in marker_genes
        if gene not in var_names
    ]

    return available, missing


def filter_marker_dict_to_available(
    marker_dict: dict[str, list[str]],
    available_markers: list[str],
) -> dict[str, list[str]]:
    """
    Keep only available genes in each marker group.
    """

    available_set = set(available_markers)

    filtered: dict[str, list[str]] = {}

    for marker_group, genes in marker_dict.items():
        kept_genes = [
            str(gene)
            for gene in genes
            if str(gene) in available_set
        ]

        # Remove duplicates inside each biological group
        kept_genes = list(dict.fromkeys(kept_genes))

        if kept_genes:
            filtered[str(marker_group)] = kept_genes

    return filtered


def build_marker_availability_table(
    marker_dict: dict[str, list[str]],
    available_markers: list[str],
) -> pd.DataFrame:
    """
    Create one row per marker-group and marker-gene combination.
    """

    available_set = set(available_markers)

    rows: list[dict[str, Any]] = []

    for marker_group, genes in marker_dict.items():
        for gene in genes:
            gene = str(gene)

            rows.append(
                {
                    "marker_group": str(marker_group),
                    "marker_gene": gene,
                    "available": gene in available_set,
                }
            )

    availability = pd.DataFrame(rows)

    availability = availability.sort_values(
        by=["marker_group", "marker_gene"],
    ).reset_index(drop=True)

    return availability


def get_expression_dataframe(
    adata: AnnData,
    marker_genes: list[str],
    use_raw: bool,
) -> pd.DataFrame:
    """
    Extract cells × marker-gene expression as a DataFrame.

    This is safe here because the marker panel is small, even if the
    full RNA matrix is sparse.
    """

    if not marker_genes:
        raise ValueError(
            "No available marker genes were provided."
        )

    if use_raw and adata.raw is not None:
        expression = adata.raw[:, marker_genes].to_adata()
    else:
        expression = adata[:, marker_genes].copy()

    return expression.to_df()


def compute_cluster_sizes(
    adata: AnnData,
    cluster_key: str,
) -> pd.DataFrame:
    """
    Compute number and fraction of cells per selected cluster.
    """

    cluster_sizes = (
        adata.obs[cluster_key]
        .astype(str)
        .value_counts(dropna=False)
        .rename_axis("cluster")
        .reset_index(name="n_cells")
    )

    cluster_sizes["fraction_cells"] = (
        cluster_sizes["n_cells"] / adata.n_obs
    )

    cluster_sizes["percent_cells"] = (
        cluster_sizes["fraction_cells"] * 100.0
    )

    return sort_cluster_table(
        cluster_sizes,
        cluster_col="cluster",
    )


def compute_mean_marker_expression(
    adata: AnnData,
    marker_genes: list[str],
    cluster_key: str,
    use_raw: bool,
) -> pd.DataFrame:
    """
    Compute mean marker expression per selected Leiden cluster.
    """

    expression_df = get_expression_dataframe(
        adata=adata,
        marker_genes=marker_genes,
        use_raw=use_raw,
    )

    expression_df["cluster"] = (
        adata.obs[cluster_key]
        .astype(str)
        .values
    )

    mean_expression = (
        expression_df
        .groupby("cluster", observed=False)
        .mean(numeric_only=True)
        .reset_index()
    )

    return sort_cluster_table(
        mean_expression,
        cluster_col="cluster",
    )


def compute_pct_marker_expressing(
    adata: AnnData,
    marker_genes: list[str],
    cluster_key: str,
    use_raw: bool,
) -> pd.DataFrame:
    """
    Compute percentage of cells expressing each marker in each cluster.

    Expression is defined as expression > 0 in the selected expression
    representation.
    """

    expression_df = get_expression_dataframe(
        adata=adata,
        marker_genes=marker_genes,
        use_raw=use_raw,
    )

    binary_expression = (
        expression_df > 0
    ).astype(float)

    binary_expression["cluster"] = (
        adata.obs[cluster_key]
        .astype(str)
        .values
    )

    pct_expression = (
        binary_expression
        .groupby("cluster", observed=False)
        .mean(numeric_only=True)
        .mul(100.0)
        .reset_index()
    )

    return sort_cluster_table(
        pct_expression,
        cluster_col="cluster",
    )


def compute_marker_group_scores(
    mean_expression: pd.DataFrame,
    marker_dict_available: dict[str, list[str]],
) -> pd.DataFrame:
    """
    Compute a simple mean-expression score for each marker group.

    These scores are annotation evidence only. They are not a trained
    classifier and must not be interpreted as final labels.
    """

    score_table = mean_expression[["cluster"]].copy()

    for marker_group, genes in marker_dict_available.items():
        genes_in_table = [
            gene
            for gene in genes
            if gene in mean_expression.columns
        ]

        if not genes_in_table:
            continue

        score_table[marker_group] = (
            mean_expression[genes_in_table]
            .mean(axis=1)
        )

    return score_table


def compute_marker_group_pct_scores(
    pct_expression: pd.DataFrame,
    marker_dict_available: dict[str, list[str]],
) -> pd.DataFrame:
    """
    Compute average percent-expressing score for each marker group.
    """

    score_table = pct_expression[["cluster"]].copy()

    for marker_group, genes in marker_dict_available.items():
        genes_in_table = [
            gene
            for gene in genes
            if gene in pct_expression.columns
        ]

        if not genes_in_table:
            continue

        score_table[marker_group] = (
            pct_expression[genes_in_table]
            .mean(axis=1)
        )

    return score_table


def build_annotation_template(
    cluster_sizes: pd.DataFrame,
    marker_group_scores: pd.DataFrame,
    marker_group_pct_scores: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build the Step 03A manual annotation worksheet.

    The strongest marker signal is included only as a review aid.
    Final annotation fields remain empty.
    """

    template = cluster_sizes.copy()

    expression_score_columns = [
        column
        for column in marker_group_scores.columns
        if column != "cluster"
    ]

    pct_score_columns = [
        column
        for column in marker_group_pct_scores.columns
        if column != "cluster"
    ]

    if expression_score_columns:
        expression_top = marker_group_scores.copy()

        expression_top["top_mean_expression_marker_group"] = (
            expression_top[expression_score_columns]
            .idxmax(axis=1)
        )

        expression_top["top_mean_expression_score"] = (
            expression_top[expression_score_columns]
            .max(axis=1)
        )

        template = template.merge(
            expression_top[
                [
                    "cluster",
                    "top_mean_expression_marker_group",
                    "top_mean_expression_score",
                ]
            ],
            on="cluster",
            how="left",
        )
    else:
        template["top_mean_expression_marker_group"] = ""
        template["top_mean_expression_score"] = pd.NA

    if pct_score_columns:
        pct_top = marker_group_pct_scores.copy()

        pct_top["top_pct_expressing_marker_group"] = (
            pct_top[pct_score_columns]
            .idxmax(axis=1)
        )

        pct_top["top_pct_expressing_score"] = (
            pct_top[pct_score_columns]
            .max(axis=1)
        )

        template = template.merge(
            pct_top[
                [
                    "cluster",
                    "top_pct_expressing_marker_group",
                    "top_pct_expressing_score",
                ]
            ],
            on="cluster",
            how="left",
        )
    else:
        template["top_pct_expressing_marker_group"] = ""
        template["top_pct_expressing_score"] = pd.NA

    # Manual-review fields: intentionally empty
    template["candidate_cell_type_broad"] = ""
    template["candidate_cell_type_detailed"] = ""
    template["annotation_confidence"] = ""
    template["positive_marker_evidence"] = ""
    template["negative_marker_evidence"] = ""
    template["possible_mixed_cluster"] = ""
    template["needs_subclustering"] = ""
    template["include_in_downstream_provisional"] = ""
    template["review_notes"] = ""

    column_order = [
        "cluster",
        "n_cells",
        "fraction_cells",
        "percent_cells",
        "top_mean_expression_marker_group",
        "top_mean_expression_score",
        "top_pct_expressing_marker_group",
        "top_pct_expressing_score",
        "candidate_cell_type_broad",
        "candidate_cell_type_detailed",
        "annotation_confidence",
        "positive_marker_evidence",
        "negative_marker_evidence",
        "possible_mixed_cluster",
        "needs_subclustering",
        "include_in_downstream_provisional",
        "review_notes",
    ]

    return template[column_order]


# ============================================================
# 4. Plotting helpers
# ============================================================

def save_marker_dotplot(
    adata: AnnData,
    marker_dict_available: dict[str, list[str]],
    cluster_key: str,
    use_raw: bool,
    standard_scale: str | None,
    output_path: Path,
    dpi: int,
    display_figures: bool,
) -> None:
    """
    Create and save grouped marker dotplot.
    """

    print("\nCreating marker dotplot...")

    dotplot = sc.pl.dotplot(
        adata,
        var_names=marker_dict_available,
        groupby=cluster_key,
        use_raw=use_raw,
        standard_scale=standard_scale,
        dendrogram=False,
        swap_axes=False,
        show=False,
        return_fig=True,
    )

    ensure_parent_dir(output_path)

    dotplot.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
    )

    if display_figures:
        dotplot.show()

    plt.close("all")


def save_marker_matrixplot(
    adata: AnnData,
    marker_dict_available: dict[str, list[str]],
    cluster_key: str,
    use_raw: bool,
    standard_scale: str | None,
    output_path: Path,
    dpi: int,
    display_figures: bool,
) -> None:
    """
    Create and save grouped marker matrixplot.
    """

    print("\nCreating marker matrixplot...")

    matrixplot = sc.pl.matrixplot(
        adata,
        var_names=marker_dict_available,
        groupby=cluster_key,
        use_raw=use_raw,
        standard_scale=standard_scale,
        dendrogram=False,
        swap_axes=False,
        show=False,
        return_fig=True,
    )

    ensure_parent_dir(output_path)

    matrixplot.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
    )

    if display_figures:
        matrixplot.show()

    plt.close("all")


def save_cluster_umap(
    adata: AnnData,
    cluster_key: str,
    output_path: Path,
    dpi: int,
    display_figures: bool,
) -> None:
    """
    Save the existing Step 02B cluster UMAP for annotation reference.
    """

    sc.pl.umap(
        adata,
        color=cluster_key,
        legend_loc="on data",
        legend_fontsize=7,
        legend_fontweight="normal",
        frameon=False,
        title=f"Selected RNA clusters: {cluster_key}",
        show=False,
    )

    save_current_figure(
        output_path=output_path,
        dpi=dpi,
        display_figures=display_figures,
    )


def save_marker_umap_batches(
    adata: AnnData,
    available_markers: list[str],
    marker_umaps_dir: Path,
    batch_size: int,
    use_raw: bool,
    dpi: int,
    cluster_key_token: str,
    display_figures: bool,
) -> list[str]:
    """
    Save marker-expression UMAP panels in manageable batches.
    """

    if "X_umap" not in adata.obsm:
        raise KeyError(
            "UMAP coordinates were not found in adata.obsm['X_umap'].\n"
            "Step 03A must not recompute UMAP. Rerun Step 02B first."
        )

    if batch_size < 1:
        raise ValueError(
            "plotting.umap_marker_batch_size must be at least 1."
        )

    marker_umaps_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    saved_files: list[str] = []

    for batch_start in range(
        0,
        len(available_markers),
        batch_size,
    ):
        marker_batch = available_markers[
            batch_start:batch_start + batch_size
        ]

        batch_number = (
            batch_start // batch_size
        ) + 1

        output_path = marker_umaps_dir / (
            f"rna_marker_umaps_{cluster_key_token}"
            f"_batch_{batch_number:02d}.png"
        )

        print(
            f"Creating marker UMAP batch {batch_number}: "
            f"{marker_batch}"
        )

        sc.pl.umap(
            adata,
            color=marker_batch,
            use_raw=use_raw,
            frameon=False,
            ncols=4,
            show=False,
        )

        save_current_figure(
            output_path=output_path,
            dpi=dpi,
            display_figures=display_figures,
        )

        saved_files.append(
            str(output_path)
        )

    return saved_files


# ============================================================
# 5. Initialize active project
# ============================================================

paths = ProjectPaths(
    refactor_root=refactor_dir,
)

paths.ensure_output_dirs()

print("\n============================================================")
print("Step 03A — RNA-informed marker annotation evidence")
print("============================================================")

print("\nRefactored root:")
print(paths.refactor_root)

print("\nActive project id:")
print(paths.project_id)

print("\nActive project root:")
print(paths.project_root)


# ============================================================
# 6. Load Step 03A marker config
# ============================================================

marker_config_path = (
    paths.configs / "rna_marker_validation.yaml"
)

if not marker_config_path.exists():
    raise FileNotFoundError(
        "Step 03A marker config was not found:\n"
        f"{marker_config_path}"
    )

assert_path_inside_project(
    marker_config_path,
    paths.project_root,
    label="rna_marker_validation.yaml",
)

marker_config = load_yaml(
    marker_config_path,
)

input_cfg = require_config_key(
    marker_config,
    "input",
    "rna_marker_validation.yaml",
)

marker_dict = require_config_key(
    marker_config,
    "marker_genes",
    "rna_marker_validation.yaml",
)

outputs_cfg = require_config_key(
    marker_config,
    "outputs",
    "rna_marker_validation.yaml",
)

if not isinstance(marker_dict, dict):
    raise TypeError(
        "marker_genes in rna_marker_validation.yaml must be a dictionary."
    )

if not marker_dict:
    raise ValueError(
        "marker_genes in rna_marker_validation.yaml is empty."
    )

plotting_cfg = marker_config.get(
    "plotting",
    {},
)

use_raw_requested = bool(
    plotting_cfg.get(
        "use_raw",
        True,
    )
)

dotplot_standard_scale = plotting_cfg.get(
    "dotplot_standard_scale",
    "var",
)

if dotplot_standard_scale in {
    "",
    "none",
    "None",
    None,
}:
    dotplot_standard_scale = None

umap_marker_batch_size = int(
    plotting_cfg.get(
        "umap_marker_batch_size",
        12,
    )
)

dpi = int(
    plotting_cfg.get(
        "dpi",
        300,
    )
)

display_figures = bool(
    plotting_cfg.get(
        "display_figures",
        True,
    )
)

print("\nLoaded marker config:")
print(marker_config_path)

print("\nMarker groups requested:")
for marker_group, genes in marker_dict.items():
    print(
        f"- {marker_group}: {len(genes)} markers"
    )


# ============================================================
# 7. Load preprocessing_rna.yaml
# ============================================================

preprocessing_config_path = (
    paths.configs / "preprocessing_rna.yaml"
)

if not preprocessing_config_path.exists():
    raise FileNotFoundError(
        "RNA preprocessing config was not found:\n"
        f"{preprocessing_config_path}"
    )

assert_path_inside_project(
    preprocessing_config_path,
    paths.project_root,
    label="preprocessing_rna.yaml",
)

preprocessing_config = load_yaml(
    preprocessing_config_path,
)

dim_cfg = require_config_key(
    preprocessing_config,
    "dimensionality_reduction",
    "preprocessing_rna.yaml",
)

required_dimensionality_keys = [
    "leiden_key",
    "n_pcs",
    "n_neighbors",
    "leiden_resolution",
]

missing_dimensionality_keys = [
    key
    for key in required_dimensionality_keys
    if key not in dim_cfg
]

if missing_dimensionality_keys:
    raise KeyError(
        "Missing required dimensionality_reduction keys in "
        "preprocessing_rna.yaml:\n"
        f"{missing_dimensionality_keys}"
    )

cluster_key = str(
    dim_cfg["leiden_key"]
)

selected_n_pcs = int(
    dim_cfg["n_pcs"]
)

selected_n_neighbors = int(
    dim_cfg["n_neighbors"]
)

selected_leiden_resolution = float(
    dim_cfg["leiden_resolution"]
)

cluster_key_token = safe_filename_token(
    cluster_key
)

print("\nSelected Step 02B clustering source:")
print(f"cluster_key: {cluster_key}")
print(f"n_pcs: {selected_n_pcs}")
print(f"n_neighbors: {selected_n_neighbors}")
print(
    f"leiden_resolution: "
    f"{selected_leiden_resolution}"
)


# ============================================================
# 8. Resolve input
# ============================================================

rna_clustered_relative = require_config_key(
    input_cfg,
    "rna_clustered_h5ad",
    "rna_marker_validation.yaml::input",
)

rna_clustered_path = paths.prepare_output_path(
    rna_clustered_relative,
    label="rna_clustered_h5ad_input",
)

assert_path_inside_project(
    rna_clustered_path,
    paths.project_root,
    label="rna_clustered_h5ad_input",
)

if not rna_clustered_path.exists():
    raise FileNotFoundError(
        "The final clustered RNA object was not found:\n"
        f"{rna_clustered_path}\n\n"
        "Check input.rna_clustered_h5ad in "
        "configs/rna_marker_validation.yaml and confirm that "
        "Step 02B completed successfully."
    )


# ============================================================
# 9. Resolve outputs
# ============================================================

marker_availability_out = paths.prepare_output_path(
    outputs_cfg.get(
        "marker_availability_table",
        (
            "reports/tables/rna/step03A/"
            "rna_marker_gene_availability.csv"
        ),
    ),
    label="marker_availability_table",
)

mean_expression_out = paths.prepare_output_path(
    outputs_cfg.get(
        "marker_mean_expression_table",
        (
            "reports/tables/rna/step03A/"
            "rna_cluster_marker_mean_expression.csv"
        ),
    ),
    label="marker_mean_expression_table",
)

pct_expressing_out = paths.prepare_output_path(
    outputs_cfg.get(
        "marker_pct_expressing_table",
        (
            "reports/tables/rna/step03A/"
            "rna_cluster_marker_pct_expressing.csv"
        ),
    ),
    label="marker_pct_expressing_table",
)

marker_group_scores_out = paths.prepare_output_path(
    outputs_cfg.get(
        "marker_group_scores_table",
        (
            "reports/tables/rna/step03A/"
            "rna_cluster_marker_group_scores.csv"
        ),
    ),
    label="marker_group_scores_table",
)

marker_group_pct_scores_out = paths.prepare_output_path(
    outputs_cfg.get(
        "marker_group_pct_scores_table",
        (
            "reports/tables/rna/step03A/"
            "rna_cluster_marker_group_pct_scores.csv"
        ),
    ),
    label="marker_group_pct_scores_table",
)

annotation_template_out = paths.prepare_output_path(
    outputs_cfg.get(
        "annotation_template_table",
        (
            "reports/tables/rna/step03A/"
            "rna_informed_annotation_template.csv"
        ),
    ),
    label="annotation_template_table",
)

cluster_sizes_out = paths.prepare_output_path(
    outputs_cfg.get(
        "cluster_sizes_table",
        (
            "reports/tables/rna/step03A/"
            "rna_step03A_cluster_sizes.csv"
        ),
    ),
    label="cluster_sizes_table",
)

marker_dotplot_dir = paths.prepare_output_path(
    outputs_cfg.get(
        "marker_dotplot_dir",
        "reports/figures/rna/step03A/annotation",
    ),
    label="marker_dotplot_dir",
)

marker_umaps_base_dir = paths.prepare_output_path(
    outputs_cfg.get(
        "marker_umaps_dir",
        (
            "reports/figures/rna/step03A/"
            "annotation/marker_umaps"
        ),
    ),
    label="marker_umaps_dir",
)

manifest_out = paths.prepare_output_path(
    outputs_cfg.get(
        "manifest",
        (
            "reports/manifests/"
            "step03A_rna_informed_annotation_manifest.json"
        ),
    ),
    label="step03A_manifest",
)

marker_umaps_dir = (
    marker_umaps_base_dir / cluster_key_token
)

marker_dotplot_out = marker_dotplot_dir / (
    f"rna_marker_dotplot_by_{cluster_key_token}.png"
)

marker_matrixplot_out = marker_dotplot_dir / (
    f"rna_marker_matrixplot_by_{cluster_key_token}.png"
)

cluster_umap_out = marker_dotplot_dir / (
    f"rna_selected_clusters_{cluster_key_token}.png"
)

all_paths = {
    "rna_clustered_input": rna_clustered_path,
    "marker_availability": marker_availability_out,
    "mean_expression": mean_expression_out,
    "pct_expressing": pct_expressing_out,
    "marker_group_scores": marker_group_scores_out,
    "marker_group_pct_scores": marker_group_pct_scores_out,
    "annotation_template": annotation_template_out,
    "cluster_sizes": cluster_sizes_out,
    "marker_dotplot_dir": marker_dotplot_dir,
    "marker_umaps_base_dir": marker_umaps_base_dir,
    "marker_umaps_dir": marker_umaps_dir,
    "marker_dotplot": marker_dotplot_out,
    "marker_matrixplot": marker_matrixplot_out,
    "cluster_umap": cluster_umap_out,
    "manifest": manifest_out,
}

for label, path in all_paths.items():
    assert_path_inside_project(
        path,
        paths.project_root,
        label=label,
    )

for output_file in [
    marker_availability_out,
    mean_expression_out,
    pct_expressing_out,
    marker_group_scores_out,
    marker_group_pct_scores_out,
    annotation_template_out,
    cluster_sizes_out,
    marker_dotplot_out,
    marker_matrixplot_out,
    cluster_umap_out,
    manifest_out,
]:
    ensure_parent_dir(output_file)

marker_dotplot_dir.mkdir(
    parents=True,
    exist_ok=True,
)

marker_umaps_dir.mkdir(
    parents=True,
    exist_ok=True,
)

print("\nResolved Step 03A paths:")
for label, path in all_paths.items():
    print(f"{label}: {path}")


# ============================================================
# 10. Load final clustered RNA object
# ============================================================

print("\nLoading final clustered RNA object:")
print(rna_clustered_path)

rna = sc.read_h5ad(
    rna_clustered_path,
)

print("\nLoaded RNA object:")
print(rna)


# ============================================================
# 11. Validate object and selected clustering
# ============================================================

if rna.n_obs == 0:
    raise ValueError(
        "The clustered RNA object contains zero cells."
    )

if rna.n_vars == 0:
    raise ValueError(
        "The clustered RNA object contains zero genes."
    )

if not rna.obs_names.is_unique:
    duplicated_cells = int(
        rna.obs_names.duplicated().sum()
    )

    raise ValueError(
        "The clustered RNA object has duplicated cell barcodes.\n"
        f"Duplicated barcodes: {duplicated_cells}"
    )

if cluster_key not in rna.obs.columns:
    raise ValueError(
        "The selected Leiden key from preprocessing_rna.yaml "
        "was not found in rna.obs.\n\n"
        f"Selected key: {cluster_key}\n"
        f"Available obs columns:\n"
        f"{list(rna.obs.columns)}\n\n"
        "This usually means preprocessing_rna.yaml changed after "
        "Step 02B or the wrong clustered RNA object was loaded."
    )

if "X_umap" not in rna.obsm:
    raise KeyError(
        "X_umap was not found in rna.obsm.\n"
        "Step 03A must not recompute UMAP. "
        "Rerun or correct Step 02B."
    )

rna.obs[cluster_key] = (
    rna.obs[cluster_key]
    .astype(str)
    .astype("category")
)

cluster_sizes = compute_cluster_sizes(
    adata=rna,
    cluster_key=cluster_key,
)

print("\nSelected cluster counts:")
print(cluster_sizes)

print("\nRNA object validation:")
print(f"n_cells: {rna.n_obs}")
print(f"n_genes_current_matrix: {rna.n_vars}")
print(f"n_clusters: {cluster_sizes.shape[0]}")
print(f"has_raw: {rna.raw is not None}")
print(f"has_X_umap: {'X_umap' in rna.obsm}")


# ============================================================
# 12. Determine marker expression source
# ============================================================

use_raw = bool(
    use_raw_requested
    and rna.raw is not None
)

if use_raw_requested and rna.raw is None:
    print(
        "\nWARNING:"
        "\nplotting.use_raw=true, but rna.raw is missing."
        "\nThe script will fall back to rna.X/rna.var_names."
        "\nIf the current RNA object is HVG-only, some markers may "
        "be missing and annotation evidence may be incomplete."
    )

if use_raw:
    print(
        "\nMarker expression source: rna.raw"
    )
    print(
        f"Genes available in raw matrix: "
        f"{rna.raw.n_vars}"
    )
else:
    print(
        "\nMarker expression source: current rna.X"
    )
    print(
        f"Genes available in current matrix: "
        f"{rna.n_vars}"
    )


# ============================================================
# 13. Marker availability
# ============================================================

marker_genes = flatten_marker_dict(
    marker_dict,
)

available_markers, missing_markers = (
    split_available_missing_markers(
        adata=rna,
        marker_genes=marker_genes,
        use_raw=use_raw,
    )
)

if not available_markers:
    raise ValueError(
        "None of the requested marker genes were found.\n"
        "Check gene symbols, rna.raw, gene-name formatting and "
        "the Step 02B output object."
    )

marker_dict_available = (
    filter_marker_dict_to_available(
        marker_dict=marker_dict,
        available_markers=available_markers,
    )
)

marker_availability = (
    build_marker_availability_table(
        marker_dict=marker_dict,
        available_markers=available_markers,
    )
)

write_csv(
    marker_availability,
    marker_availability_out,
    index=False,
)

marker_availability_fraction = (
    len(available_markers)
    / len(marker_genes)
)

print("\nMarker availability summary:")
print(
    f"Requested unique markers: "
    f"{len(marker_genes)}"
)
print(
    f"Available unique markers: "
    f"{len(available_markers)}"
)
print(
    f"Missing unique markers: "
    f"{len(missing_markers)}"
)
print(
    f"Availability fraction: "
    f"{marker_availability_fraction:.3f}"
)

print("\nAvailable marker groups:")
for marker_group, genes in marker_dict_available.items():
    print(
        f"- {marker_group}: "
        f"{len(genes)} available markers"
    )

print("\nMissing markers:")
print(missing_markers)

print("\nSaved marker availability table:")
print(marker_availability_out)


# ============================================================
# 14. Save cluster-size table
# ============================================================

write_csv(
    cluster_sizes,
    cluster_sizes_out,
    index=False,
)

print("\nSaved cluster-size table:")
print(cluster_sizes_out)


# ============================================================
# 15. Save selected cluster UMAP
# ============================================================

print("\nSaving selected cluster UMAP...")

save_cluster_umap(
    adata=rna,
    cluster_key=cluster_key,
    output_path=cluster_umap_out,
    dpi=dpi,
    display_figures=display_figures,
)

print("\nSaved selected cluster UMAP:")
print(cluster_umap_out)


# ============================================================
# 16. Marker dotplot
# ============================================================

save_marker_dotplot(
    adata=rna,
    marker_dict_available=marker_dict_available,
    cluster_key=cluster_key,
    use_raw=use_raw,
    standard_scale=dotplot_standard_scale,
    output_path=marker_dotplot_out,
    dpi=dpi,
    display_figures=display_figures,
)

print("\nSaved marker dotplot:")
print(marker_dotplot_out)


# ============================================================
# 17. Marker matrixplot
# ============================================================

save_marker_matrixplot(
    adata=rna,
    marker_dict_available=marker_dict_available,
    cluster_key=cluster_key,
    use_raw=use_raw,
    standard_scale=dotplot_standard_scale,
    output_path=marker_matrixplot_out,
    dpi=dpi,
    display_figures=display_figures,
)

print("\nSaved marker matrixplot:")
print(marker_matrixplot_out)


# ============================================================
# 18. Marker UMAP batches
# ============================================================

print("\nCreating marker-expression UMAP batches...")

marker_umap_files = save_marker_umap_batches(
    adata=rna,
    available_markers=available_markers,
    marker_umaps_dir=marker_umaps_dir,
    batch_size=umap_marker_batch_size,
    use_raw=use_raw,
    dpi=dpi,
    cluster_key_token=cluster_key_token,
    display_figures=display_figures,
)

print("\nSaved marker UMAP files:")
for marker_umap_file in marker_umap_files:
    print(marker_umap_file)


# ============================================================
# 19. Marker mean-expression table
# ============================================================

print(
    "\nComputing mean marker expression per cluster..."
)

mean_expression = compute_mean_marker_expression(
    adata=rna,
    marker_genes=available_markers,
    cluster_key=cluster_key,
    use_raw=use_raw,
)

write_csv(
    mean_expression,
    mean_expression_out,
    index=False,
)

print("\nSaved mean-expression table:")
print(mean_expression_out)


# ============================================================
# 20. Marker percent-expressing table
# ============================================================

print(
    "\nComputing percent of cells expressing each marker..."
)

pct_expressing = compute_pct_marker_expressing(
    adata=rna,
    marker_genes=available_markers,
    cluster_key=cluster_key,
    use_raw=use_raw,
)

write_csv(
    pct_expressing,
    pct_expressing_out,
    index=False,
)

print("\nSaved percent-expressing table:")
print(pct_expressing_out)


# ============================================================
# 21. Marker-group score tables
# ============================================================

print("\nComputing marker-group mean-expression scores...")

marker_group_scores = compute_marker_group_scores(
    mean_expression=mean_expression,
    marker_dict_available=marker_dict_available,
)

write_csv(
    marker_group_scores,
    marker_group_scores_out,
    index=False,
)

print("\nSaved marker-group mean-expression scores:")
print(marker_group_scores_out)


print("\nComputing marker-group percent-expressing scores...")

marker_group_pct_scores = compute_marker_group_pct_scores(
    pct_expression=pct_expressing,
    marker_dict_available=marker_dict_available,
)

write_csv(
    marker_group_pct_scores,
    marker_group_pct_scores_out,
    index=False,
)

print("\nSaved marker-group percent-expressing scores:")
print(marker_group_pct_scores_out)


# ============================================================
# 22. Manual annotation worksheet
# ============================================================

print("\nBuilding manual annotation worksheet...")

annotation_template = build_annotation_template(
    cluster_sizes=cluster_sizes,
    marker_group_scores=marker_group_scores,
    marker_group_pct_scores=marker_group_pct_scores,
)

write_csv(
    annotation_template,
    annotation_template_out,
    index=False,
)

print("\nSaved manual annotation worksheet:")
print(annotation_template_out)

print("\nAnnotation worksheet preview:")
print(annotation_template.head(20))


# ============================================================
# 23. Manifest
# ============================================================

manifest = {
    "step": "03A",
    "script": "run_03A_rna_informed_annotation.py",
    "description": (
        "RNA-informed marker-based annotation evidence for "
        "the final selected Step 02B clusters. No final cell-type "
        "labels are assigned in this step."
    ),
    "timestamp": datetime.now().isoformat(
        timespec="seconds",
    ),
    "project": paths.as_dict(),
    "configs": {
        "marker_config": str(marker_config_path),
        "preprocessing_config": str(
            preprocessing_config_path
        ),
    },
    "clustering_source": {
        "source_config": str(
            preprocessing_config_path
        ),
        "source_key": (
            "dimensionality_reduction.leiden_key"
        ),
        "cluster_key": cluster_key,
        "n_pcs": selected_n_pcs,
        "n_neighbors": selected_n_neighbors,
        "leiden_resolution": (
            selected_leiden_resolution
        ),
    },
    "inputs": {
        "rna_clustered_h5ad": str(
            rna_clustered_path
        ),
    },
    "outputs": {
        "marker_availability_table": str(
            marker_availability_out
        ),
        "cluster_sizes_table": str(
            cluster_sizes_out
        ),
        "selected_cluster_umap": str(
            cluster_umap_out
        ),
        "marker_dotplot": str(
            marker_dotplot_out
        ),
        "marker_matrixplot": str(
            marker_matrixplot_out
        ),
        "marker_umaps_dir": str(
            marker_umaps_dir
        ),
        "marker_umap_files": (
            marker_umap_files
        ),
        "marker_mean_expression_table": str(
            mean_expression_out
        ),
        "marker_pct_expressing_table": str(
            pct_expressing_out
        ),
        "marker_group_scores_table": str(
            marker_group_scores_out
        ),
        "marker_group_pct_scores_table": str(
            marker_group_pct_scores_out
        ),
        "annotation_template_table": str(
            annotation_template_out
        ),
        "manifest": str(
            manifest_out
        ),
    },
    "validation": {
        "cluster_key_found_in_rna_obs": bool(
            cluster_key in rna.obs.columns
        ),
        "umap_found_in_rna_obsm": bool(
            "X_umap" in rna.obsm
        ),
        "obs_names_unique": bool(
            rna.obs_names.is_unique
        ),
        "all_paths_project_local": True,
        "requested_use_raw": bool(
            use_raw_requested
        ),
        "actual_use_raw": bool(
            use_raw
        ),
        "has_raw": bool(
            rna.raw is not None
        ),
        "display_figures": bool(
            display_figures
        ),
        "preprocessing_rerun": False,
        "pca_rerun": False,
        "neighbors_rerun": False,
        "umap_rerun": False,
        "leiden_rerun": False,
        "final_annotation_labels_written": False,
        "annotated_h5ad_written": False,
    },
    "summary": {
        "n_cells": int(
            rna.n_obs
        ),
        "n_genes_current_matrix": int(
            rna.n_vars
        ),
        "n_genes_raw_matrix": (
            int(rna.raw.n_vars)
            if rna.raw is not None
            else None
        ),
        "n_clusters": int(
            cluster_sizes.shape[0]
        ),
        "clusters": (
            cluster_sizes["cluster"]
            .astype(str)
            .tolist()
        ),
        "n_marker_groups_requested": int(
            len(marker_dict)
        ),
        "n_marker_groups_available": int(
            len(marker_dict_available)
        ),
        "marker_groups_requested": list(
            marker_dict.keys()
        ),
        "marker_groups_available": list(
            marker_dict_available.keys()
        ),
        "n_marker_genes_requested_unique": int(
            len(marker_genes)
        ),
        "n_marker_genes_available_unique": int(
            len(available_markers)
        ),
        "n_marker_genes_missing_unique": int(
            len(missing_markers)
        ),
        "marker_availability_fraction": float(
            marker_availability_fraction
        ),
        "available_markers": (
            available_markers
        ),
        "missing_markers": (
            missing_markers
        ),
    },
    "interpretation_note": (
        "Cell types must be assigned using combinations of "
        "positive and negative markers, expression prevalence, "
        "cluster size and UMAP localization. The strongest "
        "marker-group score is a review aid and not a final "
        "cell-type prediction. Mixed, low-quality, activated or "
        "small clusters should remain visible for manual review."
    ),
}

write_json(
    manifest,
    manifest_out,
)

print("\nSaved Step 03A manifest:")
print(manifest_out)


# ============================================================
# 24. Final summary
# ============================================================

print("\n============================================================")
print("Step 03A completed successfully")
print("============================================================")

print(f"Active project: {paths.project_id}")
print(f"RNA cells: {rna.n_obs}")
print(f"Selected cluster key: {cluster_key}")
print(f"Number of clusters: {cluster_sizes.shape[0]}")
print(f"Requested markers: {len(marker_genes)}")
print(f"Available markers: {len(available_markers)}")
print(f"Missing markers: {len(missing_markers)}")
print(f"Actual use_raw: {use_raw}")

print("\nNo final labels were assigned.")
print("No annotated h5ad was written.")
print("No PCA, neighbors, UMAP or Leiden analysis was rerun.")

print(
    "\nNext action:"
    "\nInspect the cluster UMAP, dotplot, matrixplot, marker UMAPs, "
    "mean-expression table and percent-expressing table."
    "\nThen review:"
    f"\n{annotation_template_out}"
)

